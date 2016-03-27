/*
    This file is part of darktable,
    copyright (c) 2009--2011 johannes hanika.
    copyright (c) 2012 henrik andersson.
    copyright (c) 2012 aldric renaudin.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bauhaus/bauhaus.h"
#include "common/debug.h"
#include "common/interpolation.h"
#include "common/opencl.h"
#include "control/conf.h"
#include "control/control.h"
#include "develop/develop.h"
#include "develop/imageop.h"
#include "develop/tiling.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "gui/guides.h"
#include "gui/presets.h"
#include "iop/iop_api.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gtk/gtk.h>
#include <inttypes.h>
#include <gdk/gdkkeysyms.h>
#include <assert.h>

DT_MODULE_INTROSPECTION(1, dt_iop_rotate_params_t)

#define CLAMPF(a, mn, mx) ((a) < (mn) ? (mn) : ((a) > (mx) ? (mx) : (a)))

/** flip H/V, rotate an image, then clip the buffer. */
typedef enum dt_iop_rotate_flags_t
{
  FLAG_FLIP_NONE = 0,
  FLAG_FLIP_HORIZONTAL,
  FLAG_FLIP_VERTICAL,
  FLAG_FLIP_BOTH
} dt_iop_rotate_flags_t;

typedef struct dt_iop_rotate_params_t
{
  float angle;
  int crop_auto;
  dt_iop_rotate_flags_t flip;

  // old keystone values for backward compatibility
  float key_h, key_v;
} dt_iop_rotate_params_t;

typedef struct dt_iop_rotate_data_t
{
  float angle;
  dt_image_orientation_t orientation;
  int crop_auto;
  dt_iop_rotate_flags_t flip;

  float m[4]; // rotation matrix

  // old keystone values for backward compatibility
  float ki_h, ki_v;
} dt_iop_rotate_data_t;

typedef struct dt_iop_rotate_gui_data_t
{
  GtkWidget *angle;
  GtkWidget *hvflip;
  GtkWidget *crop_auto;
  GtkWidget *old_keystone;

  float cropw, croph;

  float rotate_ini[2];
} dt_iop_rotate_gui_data_t;

typedef struct dt_iop_rotate_global_data_t
{
  int kernel_rotate_bilinear;
  int kernel_rotate_bicubic;
  int kernel_rotate_lanczos2;
  int kernel_rotate_lanczos3;
  int kernel_flip;
} dt_iop_rotate_global_data_t;

typedef struct dt_iop_rotate_iop_clipping_params_t
{
  float angle, cx, cy, cw, ch, k_h, k_v;
  float kxa, kya, kxb, kyb, kxc, kyc, kxd, kyd;
  int k_type, k_sym;
  int k_apply, crop_auto;
  int ratio_n, ratio_d;
} dt_iop_rotate_iop_clipping_params_t;



const char *name()
{
  return _("rotate");
}

int groups()
{
  return IOP_GROUP_BASIC;
}

int flags()
{
  return IOP_FLAGS_ALLOW_TILING | IOP_FLAGS_TILING_FULL_ROI | IOP_FLAGS_ONE_INSTANCE;
}

int operation_tags()
{
  return IOP_TAG_DISTORT | IOP_TAG_CLIPPING;
}

int operation_tags_filter()
{
  // switch off watermark, it gets confused.
  return IOP_TAG_DECORATION;
}

int accept_extern_params(struct dt_iop_module_t *self, char *iop_name, int params_version)
{
  if(strcmp(iop_name, "clipping") == 0 && params_version == 5) return 1;
  return 0;
}

int handle_extern_params(struct dt_iop_module_t *self, char *iop_name, void *previous_params,
                         void *extern_params, void *new_params)
{
  if(strcmp(iop_name, "clipping") == 0)
  {
    dt_iop_rotate_params_t *p_params = NULL;
    if(previous_params) p_params = (dt_iop_rotate_params_t *)previous_params;
    dt_iop_rotate_params_t *n_params = (dt_iop_rotate_params_t *)new_params;
    dt_iop_rotate_iop_clipping_params_t *clipping_params
        = (dt_iop_rotate_iop_clipping_params_t *)extern_params;

    // new params initialisation
    if(p_params)
    {
      memcpy(n_params, p_params, sizeof(dt_iop_rotate_params_t));
    }
    else
    {
      dt_iop_rotate_params_t tmp = (dt_iop_rotate_params_t){ 0.0f, 1, 0, 0.0f, 0.0f };
      memcpy(n_params, &tmp, sizeof(dt_iop_rotate_params_t));
    }

    // and we include clipping iop params
    n_params->angle += clipping_params->angle;
    n_params->crop_auto = clipping_params->crop_auto;
    n_params->key_h = clipping_params->k_h;
    n_params->key_v = clipping_params->k_v;

    if(clipping_params->cw < 0 && clipping_params->ch < 0)
    {
      if(n_params->flip == FLAG_FLIP_NONE)
        n_params->flip = FLAG_FLIP_BOTH;
      else if(n_params->flip == FLAG_FLIP_HORIZONTAL)
        n_params->flip = FLAG_FLIP_VERTICAL;
      else if(n_params->flip == FLAG_FLIP_VERTICAL)
        n_params->flip = FLAG_FLIP_HORIZONTAL;
      else if(n_params->flip == FLAG_FLIP_BOTH)
        n_params->flip = FLAG_FLIP_NONE;
    }
    else if(clipping_params->cw < 0)
    {
      if(n_params->flip == FLAG_FLIP_NONE)
        n_params->flip = FLAG_FLIP_HORIZONTAL;
      else if(n_params->flip == FLAG_FLIP_HORIZONTAL)
        n_params->flip = FLAG_FLIP_NONE;
      else if(n_params->flip == FLAG_FLIP_VERTICAL)
        n_params->flip = FLAG_FLIP_BOTH;
      else if(n_params->flip == FLAG_FLIP_BOTH)
        n_params->flip = FLAG_FLIP_VERTICAL;
    }
    else if(clipping_params->ch < 0)
    {
      if(n_params->flip == FLAG_FLIP_NONE)
        n_params->flip = FLAG_FLIP_VERTICAL;
      else if(n_params->flip == FLAG_FLIP_HORIZONTAL)
        n_params->flip = FLAG_FLIP_BOTH;
      else if(n_params->flip == FLAG_FLIP_VERTICAL)
        n_params->flip = FLAG_FLIP_NONE;
      else if(n_params->flip == FLAG_FLIP_BOTH)
        n_params->flip = FLAG_FLIP_HORIZONTAL;
    }

    return 1;
  }
  return 0;
}

static void _iop_rotate_add_orientations(dt_image_orientation_t *orientation,
                                         dt_image_orientation_t added_orientation)
{
  dt_image_orientation_t orientation_corrected = *orientation;

  /*
   * if added orientation has ORIENTATION_SWAP_XY set, then we need
   * to swap ORIENTATION_FLIP_Y and ORIENTATION_FLIP_X bits
   * in ini orientation
   */
  if((added_orientation & ORIENTATION_SWAP_XY) == ORIENTATION_SWAP_XY)
  {
    if((*orientation & ORIENTATION_FLIP_Y) == ORIENTATION_FLIP_Y)
      orientation_corrected |= ORIENTATION_FLIP_X;
    else
      orientation_corrected &= ~ORIENTATION_FLIP_X;

    if((*orientation & ORIENTATION_FLIP_X) == ORIENTATION_FLIP_X)
      orientation_corrected |= ORIENTATION_FLIP_Y;
    else
      orientation_corrected &= ~ORIENTATION_FLIP_Y;

    if((*orientation & ORIENTATION_SWAP_XY) == ORIENTATION_SWAP_XY)
      orientation_corrected |= ORIENTATION_SWAP_XY;
  }

  // and now we can automagically compute new new flip
  *orientation = orientation_corrected ^ added_orientation;
}

static void _iop_rotate_mul_mat_vec_2(const float *m, float *p, float *o)
{
  o[0] = p[0] * m[0] + p[1] * m[1];
  o[1] = p[0] * m[2] + p[1] * m[3];
}

static void _iop_rotate_transform(float *points, float *m, float *center_in, float *center_out,
                                  float scale_in, float scale_out, float *key)
{
  float rt[] = { m[0], -m[1], -m[2], m[3] };

  // we want coordinates in buffer centered space
  float pi[2];
  pi[0] = points[0] - center_in[0];
  pi[1] = points[1] - center_in[1];

  // we want coordinates in non-scaled buffer
  pi[0] /= scale_in;
  pi[1] /= scale_in;

  _iop_rotate_mul_mat_vec_2(rt, pi, points);

  // we apply old keystoning
  points[1] *= (1.0f + points[0] * key[0]);
  points[0] *= (1.0f + points[1] * key[1]);

  // we return values in scaled buffer
  points[0] *= scale_out;
  points[1] *= scale_out;

  // and in "classic" (non centered) space
  points[0] += center_out[0];
  points[1] += center_out[1];
}

static void _iop_rotate_backtransform(float *points, float *m, float *center_in, float *center_out,
                                      float scale_in, float scale_out, float *key)
{
  float pi[2];
  pi[0] = points[0] - center_out[0];
  pi[1] = points[1] - center_out[1];
  pi[0] /= scale_out;
  pi[1] /= scale_out;

  pi[1] /= (1.0f + pi[0] * key[0]);
  pi[0] /= (1.0f + pi[1] * key[1]);

  _iop_rotate_mul_mat_vec_2(m, pi, points);

  points[0] *= scale_in;
  points[1] *= scale_in;
  points[0] += center_in[0];
  points[1] += center_in[1];
}



static void _iop_rotate_get_matrice(float angle_rad, float angle_deg, dt_iop_rotate_flags_t flip, float *m)
{
  // we want to be sure to be precise in certain cases
  if(angle_deg == 0.0f)
  {
    m[0] = m[3] = 1.0f;
    m[1] = m[2] = 0.0f;
  }
  else if(angle_deg == 90.0f)
  {
    m[0] = m[3] = 0.0f;
    m[1] = 1.0f;
    m[2] = -1.0f;
  }
  else if(angle_deg == -90.0f)
  {
    m[0] = m[3] = 0.0f;
    m[1] = -1.0f;
    m[2] = 1.0f;
  }
  else if(angle_deg == 180.0f || angle_deg == -180.0f)
  {
    m[0] = m[3] = -1.0f;
    m[1] = m[2] = 0.0f;
  }
  else
  {
    m[0] = cosf(angle_rad);
    m[1] = sinf(angle_rad);
    m[2] = -sinf(angle_rad);
    m[3] = cosf(angle_rad);
  }

  if(flip == FLAG_FLIP_HORIZONTAL || flip == FLAG_FLIP_BOTH)
  {
    m[0] = -m[0];
    m[2] = -m[2];
  }
  if(flip == FLAG_FLIP_VERTICAL || flip == FLAG_FLIP_BOTH)
  {
    m[1] = -m[1];
    m[3] = -m[3];
  }
}

// this calculate the autocrop area, eventually flipped
// corn_out params can be set to NULL
// iwidth and iheight are input buffer sizes unscaled
// returned values :
// cropscale : cropscale * iwidth is the width of the cropped area on out scale.
// => cropx = (out_width_uncropped - iwidth * crospcale) / 2
static void _iop_rotate_autocrop(float angle, float *m, float *key, int iwidth, int iheight, float *cropscale,
                                 int *flipped, float **corn_out_x, float **corn_out_y)
{
  if(iwidth < 1 || iheight < 1) return;

  float corners[4][2] = { { 0.0f, 0.0f }, { iwidth, 0.0f }, { iwidth, iheight }, { 0.0f, iheight } };
  float center_in[2] = { .5f * iwidth, .5f * iheight };
  float center_out[2] = { 0.0f, 0.0f };

  float cropscale0 = 0.0f;

  for(int flip = 0; flip < 2; flip++)
  {
    float newcropscale = 1.0f;

    // we use width and height according to flip or not
    float w = iwidth;
    float h = iheight;
    if(flip)
    {
      h = iwidth;
      w = iheight;
    }

    float p[2];
    for(int c = 0; c < 4; c++)
    {
      p[0] = corners[c][0];
      p[1] = corners[c][1];
      // we want values in buffer-centered space so we set center_out to 0
      _iop_rotate_transform(p, m, center_in, center_out, 1.0f, 1.0f, key);

      if(fabsf(p[0]) > 0.001f)
      {
        newcropscale = fminf(newcropscale, .5f * w / fabsf(p[0]));
      }
      if(fabsf(p[1]) > 0.001f)
      {
        newcropscale = fminf(newcropscale, .5f * h / fabsf(p[1]));
      }
      // we register corn_out values if set and not fliped
      if(flip == 0 && corn_out_x && corn_out_y)
      {
        (*corn_out_x)[c] = p[0] + .5f * iwidth;
        (*corn_out_y)[c] = p[1] + .5f * iheight;
      }
    }
    if(newcropscale >= cropscale0)
    {
      cropscale0 = newcropscale;
      // rotate and clip to max extent
      *cropscale = newcropscale;
      *flipped = flip;
    }
  }
}

int distort_transform(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, float *points, size_t points_count)
{
  dt_iop_rotate_data_t *d = (dt_iop_rotate_data_t *)piece->data;

  // old keystone system
  float key[2] = { 0.0f, 0.0f };
  if(d->ki_h != 0.0f || d->ki_v != 0.0f)
  {
    // correct keystone correction factors by resolution of this buffer
    const float kc = 1.0f / fminf(piece->buf_in.width, piece->buf_in.height);
    key[0] = d->ki_h * kc;
    key[1] = d->ki_v * kc;
  }

  float center_in[2] = { .5f * piece->buf_in.width, .5f * piece->buf_in.height };
  float center_out[2] = { .5f * piece->buf_out.width, .5f * piece->buf_out.height };

  for(size_t i = 0; i < points_count * 2; i += 2)
  {
    _iop_rotate_transform(points + i, d->m, center_in, center_out, 1.0f, 1.0f, key);
  }

  return 1;
}
int distort_backtransform(dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, float *points,
                          size_t points_count)
{
  dt_iop_rotate_data_t *d = (dt_iop_rotate_data_t *)piece->data;

  // old keystone system
  float key[2] = { 0.0f, 0.0f };
  if(d->ki_h != 0.0f || d->ki_v != 0.0f)
  {
    // correct keystone correction factors by resolution of this buffer
    const float kc = 1.0f / fminf(piece->buf_in.width, piece->buf_in.height);
    key[0] = d->ki_h * kc;
    key[1] = d->ki_v * kc;
  }

  float center_in[2] = { .5f * piece->buf_in.width, .5f * piece->buf_in.height };
  float center_out[2] = { .5f * piece->buf_out.width, .5f * piece->buf_out.height };

  for(size_t i = 0; i < points_count * 2; i += 2)
  {
    _iop_rotate_backtransform(points + i, d->m, center_in, center_out, 1.0f, 1.0f, key);
  }

  return 1;
}

// 1st pass: how large would the output be, given this input roi?
// this is always called with the full buffer before processing.
void modify_roi_out(struct dt_iop_module_t *self, struct dt_dev_pixelpipe_iop_t *piece, dt_iop_roi_t *roi_out,
                    const dt_iop_roi_t *roi_in_orig)
{
  dt_iop_roi_t roi_in_d = *roi_in_orig;
  dt_iop_roi_t *roi_in = &roi_in_d;

  dt_iop_rotate_data_t *d = (dt_iop_rotate_data_t *)piece->data;
  dt_iop_rotate_params_t *params = (dt_iop_rotate_params_t *)self->params;
  dt_iop_rotate_gui_data_t *g = (dt_iop_rotate_gui_data_t *)self->gui_data;

  // old keystone system
  float key[2] = { 0.0f, 0.0f };
  if(d->ki_h != 0.0f || d->ki_v != 0.0f)
  {
    // correct keystone correction factors by resolution of this buffer
    const float kc = 1.0f / fminf(roi_in->width, roi_in->height);
    key[0] = d->ki_h * kc;
    key[1] = d->ki_v * kc;
  }

  float cropscale = 1.0f;
  int crop_flip = 0;
  float *corn_out_x = (float *)calloc(4, sizeof(float));
  float *corn_out_y = (float *)calloc(4, sizeof(float));
  _iop_rotate_autocrop(d->angle, d->m, key, piece->buf_in.width, piece->buf_in.height, &cropscale, &crop_flip,
                       &corn_out_x, &corn_out_y);

  *roi_out = *roi_in;

  // we have corn_out values in "full" buffer, but we want them in roi_in one
  // so we have to apply some scale
  for(int c = 0; c < 4; c++)
  {
    // and we set the values
    corn_out_x[c] *= roi_in->scale;
    corn_out_y[c] *= roi_in->scale;
  }

  // we search new width and height of the image (without crop)
  // as it's just rotation, we only check diagonals
  float v1 = fabsf(corn_out_x[0] - corn_out_x[2]);
  float v2 = fabsf(corn_out_x[1] - corn_out_x[3]);
  float new_w, new_h;
  if(v1 > v2)
  {
    new_w = v1;
    new_h = fabsf(corn_out_y[1] - corn_out_y[3]);
  }
  else
  {
    new_w = v2;
    new_h = fabsf(corn_out_y[0] - corn_out_y[2]);
  }

  free(corn_out_x);
  free(corn_out_y);

  // now we crop
  if(d->crop_auto)
  {
    if(!crop_flip)
    {
      roi_out->width = (float)piece->buf_in.width * cropscale;
      roi_out->height = (float)piece->buf_in.height * cropscale;
    }
    else
    {
      roi_out->width = (float)piece->buf_in.height * cropscale;
      roi_out->height = (float)piece->buf_in.width * cropscale;
    }
    roi_out->x = (new_w - (float)roi_out->width) * 0.5f;
    roi_out->y = (new_h - (float)roi_out->height) * 0.5f;
  }
  else
  {
    roi_out->width = new_w;
    roi_out->height = new_h;
  }

  // we eventually set gui values here to display the cropped area
  if(self->dev->gui_attached && g && piece->pipe->type == DT_DEV_PIXELPIPE_PREVIEW)
  {
    if(!d->crop_auto && params->crop_auto)
    {
      if(!crop_flip)
      {
        g->cropw = (float)piece->buf_in.width * cropscale;
        g->croph = (float)piece->buf_in.height * cropscale;
      }
      else
      {
        g->cropw = (float)piece->buf_in.height * cropscale;
        g->croph = (float)piece->buf_in.width * cropscale;
      }
    }
    else
    {
      g->cropw = -1.0f;
      g->croph = -1.0f;
    }
  }
}

// 2nd pass: which roi would this operation need as input to fill the given output region?
void modify_roi_in(struct dt_iop_module_t *self, struct dt_dev_pixelpipe_iop_t *piece,
                   const dt_iop_roi_t *roi_out, dt_iop_roi_t *roi_in)
{
  dt_iop_rotate_data_t *d = (dt_iop_rotate_data_t *)piece->data;
  dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)self->params;
  *roi_in = *roi_out;

  const float in_w = piece->buf_in.width * roi_in->scale;
  const float in_h = piece->buf_in.height * roi_in->scale;
  const float out_w = piece->buf_out.width * roi_out->scale;
  const float out_h = piece->buf_out.height * roi_out->scale;

  // old keystone system
  float key[2] = { 0.0f, 0.0f };
  if(d->ki_h != 0.0f || d->ki_v != 0.0f)
  {
    // correct keystone correction factors by resolution of this buffer
    const float kc = 1.0f / fminf(in_w, in_h);
    key[0] = d->ki_h * kc;
    key[1] = d->ki_v * kc;
  }

  // corner points
  float corners[4][2] = { { roi_out->x, roi_out->y },
                          { roi_out->x + roi_out->width, roi_out->y },
                          { roi_out->x + roi_out->width, roi_out->y + roi_out->height },
                          { roi_out->x, roi_out->y + roi_out->height } };

  float center_in[2] = { .5f * in_w, .5f * in_h };
  float center_out[2] = { .5f * out_w, .5f * out_h };

  for(int c = 0; c < 4; c++)
  {
    _iop_rotate_backtransform(corners[c], d->m, center_in, center_out, roi_in->scale, roi_out->scale, key);
  }

  // we search the min and max of all 4 corners in x and y
  float new_l, new_t, new_r, new_b;
  new_l = fminf(fminf(fminf(corners[0][0], corners[1][0]), corners[2][0]), corners[3][0]);
  new_t = fminf(fminf(fminf(corners[0][1], corners[1][1]), corners[2][1]), corners[3][1]);
  new_r = fmaxf(fmaxf(fmaxf(corners[0][0], corners[1][0]), corners[2][0]), corners[3][0]);
  new_b = fmaxf(fmaxf(fmaxf(corners[0][1], corners[1][1]), corners[2][1]), corners[3][1]);

  // adjust roi_in to minimally needed region
  if(p->angle == 0.0f || p->angle == 90.0f || p->angle == 180.0f || p->angle == -90.0f || p->angle == 180.0f)
  {
    roi_in->x = new_l;
    roi_in->y = new_t;
    roi_in->width = new_r - new_l;
    roi_in->height = new_b - new_t;
  }
  else
  {
    roi_in->x = new_l - 1;
    roi_in->y = new_t - 1;
    roi_in->width = new_r - new_l + 2;
    roi_in->height = new_b - new_t + 2;
  }

  // sanity check.
  roi_in->x = CLAMP(roi_in->x, 0, (int)floorf(in_w) - 2);
  roi_in->y = CLAMP(roi_in->y, 0, (int)floorf(in_h) - 2);
  roi_in->width = CLAMP(roi_in->width, 1, (int)floorf(in_w) - roi_in->x);
  roi_in->height = CLAMP(roi_in->height, 1, (int)floorf(in_h) - roi_in->y);
}

// 3rd (final) pass: you get this input region (may be different from what was requested above),
// do your best to fill the output region!
void process(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, const void *const ivoid,
             void *const ovoid, const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  dt_iop_rotate_data_t *d = (dt_iop_rotate_data_t *)piece->data;

  const int ch = piece->colors;
  const int ch_width = ch * roi_in->width;
  assert(ch == 4);

  // we take care of "special" cases, which can have a fast path
  if(d->orientation != ORIENTATION_NULL)
  {
    const int bpp = sizeof(float) * piece->colors;
    const int stride = bpp * roi_in->width;
    dt_imageio_flip_buffers((char *)ovoid, (const char *)ivoid, bpp, roi_in->width, roi_in->height,
                            roi_in->width, roi_in->height, stride, d->orientation);
    return;
  }

  const struct dt_interpolation *interpolation = dt_interpolation_new(DT_INTERPOLATION_USERPREF);
  const float in_w = piece->buf_in.width * roi_in->scale;
  const float in_h = piece->buf_in.height * roi_in->scale;
  const float out_w = piece->buf_out.width * roi_out->scale;
  const float out_h = piece->buf_out.height * roi_out->scale;
  // old keystone system
  float key[2] = { 0.0f, 0.0f };
  if(d->ki_h != 0.0f || d->ki_v != 0.0f)
  {
    // correct keystone correction factors by resolution of this buffer
    const float kc = 1.0f / fminf(in_w, in_h);
    key[0] = d->ki_h * kc;
    key[1] = d->ki_v * kc;
  }
  float center_in[2] = { .5f * in_w, .5f * in_h };
  float center_out[2] = { .5f * out_w, .5f * out_h };

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(d, interpolation, key, center_in, center_out)
#endif
  // (slow) point-by-point transformation.
  // TODO: optimize with scanlines and linear steps between?
  for(int j = 0; j < roi_out->height; j++)
  {
    float *out = ((float *)ovoid) + (size_t)ch * j * roi_out->width;
    for(int i = 0; i < roi_out->width; i++, out += ch)
    {
      float p[2];

      p[0] = roi_out->x + i + 0.5f;
      p[1] = roi_out->y + j + 0.5f;
      _iop_rotate_backtransform(p, d->m, center_in, center_out, roi_in->scale, roi_out->scale, key);
      p[0] -= roi_in->x + .5f;
      p[1] -= roi_in->y + .5f;

      dt_interpolation_compute_pixel4c(interpolation, (float *)ivoid, out, p[0], p[1], roi_in->width,
                                       roi_in->height, ch_width);
    }
  }
}

#ifdef HAVE_OPENCL
int process_cl(struct dt_iop_module_t *self, dt_dev_pixelpipe_iop_t *piece, cl_mem dev_in, cl_mem dev_out,
               const dt_iop_roi_t *const roi_in, const dt_iop_roi_t *const roi_out)
{
  dt_iop_rotate_data_t *d = (dt_iop_rotate_data_t *)piece->data;
  dt_iop_rotate_global_data_t *gd = (dt_iop_rotate_global_data_t *)self->data;

  cl_int err = -999;
  const int devid = piece->pipe->devid;

  if(d->orientation != ORIENTATION_NULL)
  {
    // fast path for special rotation values (90, -90, 180)
    const int width = roi_in->width;
    const int height = roi_in->height;

    const uint32_t orientation = d->orientation;
    size_t sizes[] = { ROUNDUPWD(width), ROUNDUPWD(height), 1 };

    dt_opencl_set_kernel_arg(devid, gd->kernel_flip, 0, sizeof(cl_mem), (void *)&dev_in);
    dt_opencl_set_kernel_arg(devid, gd->kernel_flip, 1, sizeof(cl_mem), (void *)&dev_out);
    dt_opencl_set_kernel_arg(devid, gd->kernel_flip, 2, sizeof(int), (void *)&width);
    dt_opencl_set_kernel_arg(devid, gd->kernel_flip, 3, sizeof(int), (void *)&height);
    dt_opencl_set_kernel_arg(devid, gd->kernel_flip, 4, sizeof(uint32_t), (void *)&orientation);
    err = dt_opencl_enqueue_kernel_2d(devid, gd->kernel_flip, sizes);
  }
  else
  {
    const int width = roi_out->width;
    const int height = roi_out->height;
    int crkernel = -1;

    const struct dt_interpolation *interpolation = dt_interpolation_new(DT_INTERPOLATION_USERPREF);

    switch(interpolation->id)
    {
      case DT_INTERPOLATION_BILINEAR:
        crkernel = gd->kernel_rotate_bilinear;
        break;
      case DT_INTERPOLATION_BICUBIC:
        crkernel = gd->kernel_rotate_bicubic;
        break;
      case DT_INTERPOLATION_LANCZOS2:
        crkernel = gd->kernel_rotate_lanczos2;
        break;
      case DT_INTERPOLATION_LANCZOS3:
        crkernel = gd->kernel_rotate_lanczos3;
        break;
      default:
        return FALSE;
    }

    int roi[2] = { roi_in->x, roi_in->y };
    float roo[2] = { roi_out->x, roi_out->y };
    float center_in[2]
        = { piece->buf_in.width * roi_in->scale * .5f, piece->buf_in.height * roi_in->scale * .5f };
    float center_out[2]
        = { piece->buf_out.width * roi_in->scale * .5f, piece->buf_out.height * roi_in->scale * .5f };
    float m[4] = { d->m[0], d->m[1], d->m[2], d->m[3] };

    // old keystone system
    float key[2] = { 0.0f, 0.0f };
    if(d->ki_h != 0.0f || d->ki_v != 0.0f)
    {
      // correct keystone correction factors by resolution of this buffer
      const float kc
          = 1.0f / fminf(piece->buf_in.width * roi_in->scale, piece->buf_in.height * roi_in->scale);
      key[0] = d->ki_h * kc;
      key[1] = d->ki_v * kc;
    }

    size_t sizes[3];

    sizes[0] = ROUNDUPWD(width);
    sizes[1] = ROUNDUPHT(height);
    sizes[2] = 1;
    dt_opencl_set_kernel_arg(devid, crkernel, 0, sizeof(cl_mem), &dev_in);
    dt_opencl_set_kernel_arg(devid, crkernel, 1, sizeof(cl_mem), &dev_out);
    dt_opencl_set_kernel_arg(devid, crkernel, 2, sizeof(int), &width);
    dt_opencl_set_kernel_arg(devid, crkernel, 3, sizeof(int), &height);
    dt_opencl_set_kernel_arg(devid, crkernel, 4, sizeof(int), &roi_in->width);
    dt_opencl_set_kernel_arg(devid, crkernel, 5, sizeof(int), &roi_in->height);
    dt_opencl_set_kernel_arg(devid, crkernel, 6, 2 * sizeof(int), &roi);
    dt_opencl_set_kernel_arg(devid, crkernel, 7, 2 * sizeof(float), &roo);
    dt_opencl_set_kernel_arg(devid, crkernel, 8, sizeof(float), &roi_in->scale);
    dt_opencl_set_kernel_arg(devid, crkernel, 9, sizeof(float), &roi_out->scale);
    dt_opencl_set_kernel_arg(devid, crkernel, 10, 2 * sizeof(float), &center_in);
    dt_opencl_set_kernel_arg(devid, crkernel, 11, 2 * sizeof(float), &center_out);
    dt_opencl_set_kernel_arg(devid, crkernel, 12, 2 * sizeof(float), &key);
    dt_opencl_set_kernel_arg(devid, crkernel, 13, 4 * sizeof(float), &m);
    err = dt_opencl_enqueue_kernel_2d(devid, crkernel, sizes);
  }
  if(err != CL_SUCCESS) goto error;

  return TRUE;

error:
  dt_print(DT_DEBUG_OPENCL, "[opencl_clipping] couldn't enqueue kernel! %d\n", err);
  return FALSE;
}
#endif


void init_global(dt_iop_module_so_t *module)
{
  const int program = 18; // rotate.cl from programs.conf
  dt_iop_rotate_global_data_t *gd
      = (dt_iop_rotate_global_data_t *)malloc(sizeof(dt_iop_rotate_global_data_t));
  module->data = gd;
  gd->kernel_rotate_bilinear = dt_opencl_create_kernel(program, "rotate_bilinear");
  gd->kernel_rotate_bicubic = dt_opencl_create_kernel(program, "rotate_bicubic");
  gd->kernel_rotate_lanczos2 = dt_opencl_create_kernel(program, "rotate_lanczos2");
  gd->kernel_rotate_lanczos3 = dt_opencl_create_kernel(program, "rotate_lanczos3");
  const int program2 = 2; // basic.cl from programs.conf
  gd->kernel_flip = dt_opencl_create_kernel(program2, "flip");
}


void cleanup_global(dt_iop_module_so_t *module)
{
  dt_iop_rotate_global_data_t *gd = (dt_iop_rotate_global_data_t *)module->data;
  dt_opencl_free_kernel(gd->kernel_rotate_bilinear);
  dt_opencl_free_kernel(gd->kernel_rotate_bicubic);
  dt_opencl_free_kernel(gd->kernel_rotate_lanczos2);
  dt_opencl_free_kernel(gd->kernel_rotate_lanczos3);
  dt_opencl_free_kernel(gd->kernel_flip);
  free(module->data);
  module->data = NULL;
}


void commit_params(struct dt_iop_module_t *self, dt_iop_params_t *p1, dt_dev_pixelpipe_t *pipe,
                   dt_dev_pixelpipe_iop_t *piece)
{
  dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)p1;
  dt_iop_rotate_data_t *d = (dt_iop_rotate_data_t *)piece->data;

  // copy unchanged values
  d->angle = M_PI / 180.0 * p->angle;
  d->flip = p->flip;

  // orientation value for fast path
  if(p->angle == 0.0f || p->angle == 90.0f || p->angle == 180.0f || p->angle == -90.0f || p->angle == 180.0f)
  {
    dt_image_orientation_t orientation = ORIENTATION_NONE;
    if(p->angle == -90.0f)
      orientation = ORIENTATION_ROTATE_CCW_90_DEG;
    else if(p->angle == 90.0f)
      orientation = ORIENTATION_ROTATE_CW_90_DEG;
    else if(p->angle == 180.0f || p->angle == -180.0f)
      orientation = ORIENTATION_ROTATE_180_DEG;

    if(d->flip == FLAG_FLIP_HORIZONTAL || d->flip == FLAG_FLIP_BOTH)
      _iop_rotate_add_orientations(&orientation, ORIENTATION_FLIP_HORIZONTALLY);
    if(d->flip == FLAG_FLIP_VERTICAL || d->flip == FLAG_FLIP_BOTH)
      _iop_rotate_add_orientations(&orientation, ORIENTATION_FLIP_VERTICALLY);

    d->orientation = orientation;
  }
  else
    d->orientation = ORIENTATION_NULL;

  _iop_rotate_get_matrice(d->angle, p->angle, d->flip, d->m);

  // for crop values, we don't want them if module has focus
  if(self->dev->gui_module == self)
    d->crop_auto = 0;
  else
    d->crop_auto = p->crop_auto;

  // old keystone values
  if(p->key_h >= -1.0 && p->key_h <= 1.0)
    d->ki_h = p->key_h;
  else
    d->ki_h = 0.0f;
  if(p->key_v >= -1.0 && p->key_v <= 1.0)
    d->ki_v = p->key_v;
  else
    d->ki_v = 0.0f;
}

void gui_focus(struct dt_iop_module_t *self, gboolean in)
{
  if(self->enabled) dt_dev_reprocess_all(self->dev);
}


void init_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  piece->data = malloc(sizeof(dt_iop_rotate_data_t));
  self->commit_params(self, self->default_params, pipe, piece);
}

void cleanup_pipe(struct dt_iop_module_t *self, dt_dev_pixelpipe_t *pipe, dt_dev_pixelpipe_iop_t *piece)
{
  free(piece->data);
  piece->data = NULL;
}

void reload_defaults(dt_iop_module_t *self)
{
  dt_iop_rotate_params_t tmp = (dt_iop_rotate_params_t){ 0.0f, 1, 0, 0.0f, 0.0f };
  memcpy(self->params, &tmp, sizeof(dt_iop_rotate_params_t));
  memcpy(self->default_params, &tmp, sizeof(dt_iop_rotate_params_t));
  self->default_enabled = 0;
}

void gui_update(struct dt_iop_module_t *self)
{
  dt_iop_rotate_gui_data_t *g = (dt_iop_rotate_gui_data_t *)self->gui_data;
  dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)self->params;

  /* update ui elements */
  dt_bauhaus_slider_set(g->angle, -p->angle);
  dt_bauhaus_combobox_set(g->hvflip, p->flip);
  dt_bauhaus_combobox_set(g->crop_auto, p->crop_auto);
  if(p->key_h != 0.0f || p->key_v != 0.0f)
  {
    dt_bauhaus_combobox_set(g->old_keystone, 0);
    gtk_widget_show(g->old_keystone);
  }
  else
  {
    gtk_widget_hide(g->old_keystone);
  }
}

void init(dt_iop_module_t *module)
{
  // module->data = malloc(sizeof(dt_iop_rotate_data_t));
  module->params = calloc(1, sizeof(dt_iop_rotate_params_t));
  module->default_params = calloc(1, sizeof(dt_iop_rotate_params_t));
  module->default_enabled = 0;
  module->params_size = sizeof(dt_iop_rotate_params_t);
  module->gui_data = NULL;
  module->priority = 446; // module order created by iop_dependencies.py, do not edit!
}

void cleanup(dt_iop_module_t *module)
{
  free(module->params);
  module->params = NULL;
}

static void _iop_rotate_angle_callback(GtkWidget *slider, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)self->params;
  if(dt_bauhaus_slider_get(slider) == -p->angle) return; // no change
  p->angle = -dt_bauhaus_slider_get(slider);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void _iop_rotate_hvflip_callback(GtkWidget *widget, dt_iop_module_t *self)
{
  if(self->dt->gui->reset) return;
  dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)self->params;
  if(dt_bauhaus_combobox_get(widget) == p->flip) return; // no change
  p->flip = dt_bauhaus_combobox_get(widget);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void _iop_rotate_crop_auto_callback(GtkWidget *combo, dt_iop_module_t *self)
{
  dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)self->params;
  if(dt_bauhaus_combobox_get(combo) == p->crop_auto) return; // no change
  p->crop_auto = dt_bauhaus_combobox_get(combo);

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

static void _iop_rotate_keystone_callback(GtkWidget *combo, dt_iop_module_t *self)
{
  dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)self->params;
  if(dt_bauhaus_combobox_get(combo) == 1)
  {
    p->key_h = p->key_v = 0.0f;
    gtk_widget_hide(combo);
  }

  dt_dev_add_history_item(darktable.develop, self, TRUE);
}

void gui_init(struct dt_iop_module_t *self)
{
  self->gui_data = calloc(1, sizeof(dt_iop_rotate_gui_data_t));
  dt_iop_rotate_gui_data_t *g = (dt_iop_rotate_gui_data_t *)self->gui_data;
  dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)self->params;

  g->cropw = -1.0f;
  g->croph = -1.0f;
  g->rotate_ini[0] = NAN;
  g->rotate_ini[1] = NAN;

  self->widget = gtk_box_new(GTK_ORIENTATION_VERTICAL, DT_BAUHAUS_SPACE);
  g->hvflip = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->hvflip, NULL, _("flip"));
  dt_bauhaus_combobox_add(g->hvflip, _("none"));
  dt_bauhaus_combobox_add(g->hvflip, _("horizontal"));
  dt_bauhaus_combobox_add(g->hvflip, _("vertical"));
  dt_bauhaus_combobox_add(g->hvflip, _("both"));
  g_signal_connect(G_OBJECT(g->hvflip), "value-changed", G_CALLBACK(_iop_rotate_hvflip_callback), self);
  gtk_widget_set_tooltip_text(g->hvflip, _("mirror image horizontally and/or vertically"));
  gtk_box_pack_start(GTK_BOX(self->widget), g->hvflip, TRUE, TRUE, 0);


  g->angle = dt_bauhaus_slider_new_with_range(self, -180.0, 180.0, 0.25, p->angle, 2);
  dt_bauhaus_widget_set_label(g->angle, NULL, _("angle"));
  dt_bauhaus_slider_set_format(g->angle, "%.02f°");
  g_signal_connect(G_OBJECT(g->angle), "value-changed", G_CALLBACK(_iop_rotate_angle_callback), self);
  gtk_widget_set_tooltip_text(g->angle,
                              _("right-click and drag a line on the image to drag a straight line"));
  gtk_box_pack_start(GTK_BOX(self->widget), g->angle, TRUE, TRUE, 0);

  g->crop_auto = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->crop_auto, NULL, _("automatic cropping"));
  dt_bauhaus_combobox_add(g->crop_auto, _("no"));
  dt_bauhaus_combobox_add(g->crop_auto, _("yes"));
  gtk_widget_set_tooltip_text(g->crop_auto, _("automatically crop to avoid black edges"));
  g_signal_connect(G_OBJECT(g->crop_auto), "value-changed", G_CALLBACK(_iop_rotate_crop_auto_callback), self);
  gtk_box_pack_start(GTK_BOX(self->widget), g->crop_auto, TRUE, TRUE, 0);

  // old keystone system (deprecated, user can only set it off)
  g->old_keystone = dt_bauhaus_combobox_new(self);
  dt_bauhaus_widget_set_label(g->old_keystone, NULL, _("keystone (deprecated)"));
  dt_bauhaus_combobox_add(g->old_keystone, _("on"));
  dt_bauhaus_combobox_add(g->old_keystone, _("reset values"));
  g_signal_connect(G_OBJECT(g->old_keystone), "value-changed", G_CALLBACK(_iop_rotate_keystone_callback),
                   self);
  gtk_widget_set_tooltip_text(
      g->old_keystone,
      _("this feature is deprecated and can't be change. use manual perspective module instead."));
  gtk_box_pack_start(GTK_BOX(self->widget), g->old_keystone, TRUE, TRUE, 0);
  if(p->key_h == 0.0f && p->key_v == 0.0f)
  {
    gtk_widget_hide(g->old_keystone);
  }
}

void gui_cleanup(struct dt_iop_module_t *self)
{
}

int button_released(struct dt_iop_module_t *self, double x, double y, int which, uint32_t state)
{
  if(which == 3)
  {
    dt_iop_rotate_gui_data_t *g = (dt_iop_rotate_gui_data_t *)self->gui_data;
    dt_iop_rotate_params_t *p = (dt_iop_rotate_params_t *)self->params;

    float pzx, pzy;
    dt_dev_get_pointer_zoom_pos(self->dev, x, y, &pzx, &pzy);
    pzx += 0.5f;
    pzy += 0.5f;

    float dx = pzx - g->rotate_ini[0];
    float dy = pzy - g->rotate_ini[1];
    if(dx < 0)
    {
      dx = -dx;
      dy = -dy;
    }
    float angle = atan2f(dy, dx);
    if(!(angle >= -M_PI / 2.0 && angle <= M_PI / 2.0)) angle = 0.0f;
    if(angle > M_PI / 4.0)
      angle = M_PI / 2.0 - angle;
    else if(angle < -M_PI / 4.0)
      angle = -M_PI / 2.0 - angle;
    else
      angle = -angle;
    angle = 180.0 / M_PI * angle + p->angle;
    if(angle < -180.0) angle += 360.0;
    if(angle > 180.0) angle -= 360.0;
    dt_bauhaus_slider_set(g->angle, -angle);

    g->rotate_ini[0] = g->rotate_ini[1] = NAN;
    dt_control_change_cursor(GDK_LEFT_PTR);
    return 1;
  }
  return 0;
}

int button_pressed(struct dt_iop_module_t *self, double x, double y, double pressure, int which, int type,
                   uint32_t state)
{
  if(which == 3)
  {
    dt_iop_rotate_gui_data_t *g = (dt_iop_rotate_gui_data_t *)self->gui_data;
    float pzx, pzy;
    dt_dev_get_pointer_zoom_pos(self->dev, x, y, &pzx, &pzy);
    pzx += 0.5f;
    pzy += 0.5f;

    g->rotate_ini[0] = pzx;
    g->rotate_ini[1] = pzy;
    dt_control_change_cursor(GDK_CROSSHAIR);
    return 1;
  }
  return 0;
}

int mouse_moved(struct dt_iop_module_t *self, double x, double y, double pressure, int which)
{
  if(darktable.control->button_down && darktable.control->button_down_which == 3)
  {
    dt_control_queue_redraw_center();
  }
  return 0;
}

void gui_post_expose(struct dt_iop_module_t *self, cairo_t *cr, int32_t width, int32_t height,
                     int32_t pointerx, int32_t pointery)
{
  dt_develop_t *dev = self->dev;
  dt_iop_rotate_gui_data_t *g = (dt_iop_rotate_gui_data_t *)self->gui_data;

  // the usual rescaling stuff
  float wd = dev->preview_pipe->backbuf_width;
  float ht = dev->preview_pipe->backbuf_height;
  if(wd < 1.0 || ht < 1.0) return;
  float zoom_y = dt_control_get_dev_zoom_y();
  float zoom_x = dt_control_get_dev_zoom_x();
  dt_dev_zoom_t zoom = dt_control_get_dev_zoom();
  int closeup = dt_control_get_dev_closeup();
  float zoom_scale = dt_dev_get_zoom_scale(dev, zoom, closeup ? 2 : 1, 1);

  cairo_translate(cr, width / 2.0, height / 2.0f);
  cairo_scale(cr, zoom_scale, zoom_scale);
  cairo_translate(cr, -.5f * wd - zoom_x * wd, -.5f * ht - zoom_y * ht);

  // we draw the cropping area
  if(self->enabled && g->cropw > 0.0f && g->croph > 0.0f && (g->cropw < wd || g->croph < ht))
  {
    cairo_set_source_rgba(cr, .2, .2, .2, .8);
    cairo_set_fill_rule(cr, CAIRO_FILL_RULE_EVEN_ODD);
    cairo_rectangle(cr, 0, 0, wd, ht);
    cairo_rectangle(cr, (wd - g->cropw) * .5f, (ht - g->croph) * .5f, g->cropw, g->croph);
    cairo_fill(cr);
    cairo_set_source_rgb(cr, .7, .7, .7);
    cairo_rectangle(cr, (wd - g->cropw) * .5f, (ht - g->croph) * .5f, g->cropw, g->croph);
    cairo_stroke(cr);
  }

  // we draw the rotation line if needed
  if(!isnan(g->rotate_ini[0]))
  {
    float pzx, pzy;
    dt_dev_get_pointer_zoom_pos(self->dev, pointerx, pointery, &pzx, &pzy);
    pzx += 0.5f;
    pzy += 0.5f;

    PangoRectangle ink;
    PangoLayout *layout;
    PangoFontDescription *desc = pango_font_description_copy_static(darktable.bauhaus->pango_font_desc);
    pango_font_description_set_weight(desc, PANGO_WEIGHT_BOLD);
    pango_font_description_set_absolute_size(desc, DT_PIXEL_APPLY_DPI(16) * PANGO_SCALE / zoom_scale);
    layout = pango_cairo_create_layout(cr);
    pango_layout_set_font_description(layout, desc);
    cairo_arc(cr, g->rotate_ini[0] * wd, g->rotate_ini[1] * ht, DT_PIXEL_APPLY_DPI(3), 0, 2.0 * M_PI);
    cairo_stroke(cr);
    cairo_arc(cr, pzx * wd, pzy * ht, DT_PIXEL_APPLY_DPI(3), 0, 2.0 * M_PI);
    cairo_stroke(cr);
    cairo_move_to(cr, g->rotate_ini[0] * wd, g->rotate_ini[1] * ht);
    cairo_line_to(cr, pzx * wd, pzy * ht);
    cairo_stroke(cr);

    // show rotation angle
    float dx = pzx * wd - g->rotate_ini[0] * wd, dy = pzy * ht - g->rotate_ini[1] * ht;
    if(dx < 0)
    {
      dx = -dx;
      dy = -dy;
    }
    float angle = atan2f(dy, dx);
    angle = angle * 180 / M_PI;
    if(angle > 45.0) angle -= 90;
    if(angle < -45.0) angle += 90;

    char view_angle[16] = { 0 };
    snprintf(view_angle, sizeof(view_angle), "%.2f°", angle);
    pango_layout_set_text(layout, view_angle, -1);
    pango_layout_get_pixel_extents(layout, &ink, NULL);
    cairo_set_source_rgb(cr, .7, .7, .7);
    cairo_move_to(cr, pzx * wd + DT_PIXEL_APPLY_DPI(20) / zoom_scale, pzy * ht - ink.height);
    pango_cairo_show_layout(cr, layout);
    pango_font_description_free(desc);
    g_object_unref(layout);
  }
}

#undef PHI
#undef INVPHI

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-space on;
