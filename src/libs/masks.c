/*
    This file is part of darktable,
    copyright (c) 2011 henrik andersson.

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

#include "common/darktable.h"
#include "develop/develop.h"
#include "control/control.h"
#include "libs/lib.h"
#include "gui/draw.h"

DT_MODULE(1)

#define DT_LIB_MASKS_CURVE_RES 64
#define DT_LIB_MASKS_MAX_CURVE_POINTS 128


typedef struct _curve_t
{
  float anchors[2][DT_LIB_MASKS_MAX_CURVE_POINTS];
  float ctrl1[2][DT_LIB_MASKS_MAX_CURVE_POINTS];
  float ctrl2[2][DT_LIB_MASKS_MAX_CURVE_POINTS];
  int states[DT_LIB_MASKS_MAX_CURVE_POINTS]; // state of the point : 0=not set; 1=auto; 9=manual
  gboolean clockwise; // counter in fact, to avoid intersection with the transition border
  //float feather_unit_vectors[2][DT_LIB_MASKS_MAX_CURVE_POINTS];
  //float feather_weights[DT_LIB_MASKS_MAX_CURVE_POINTS];
  uint32_t anchor_count;
} _curve_t;

typedef struct dt_lib_masks_t
{
  uint32_t viewport_width, viewport_height;

  _curve_t curve;
  gboolean is_curve_edit_mode;
  gboolean is_dragging_anchor;
  int32_t drag_anchor;
  int32_t current_anchor;

  int32_t highlight_anchor;
  gboolean is_feather_highlight;
  int32_t highlight_segment;
  gboolean is_inside;
  
}
dt_lib_masks_t;

const char* name()
{
  return _("masks");
}

uint32_t views()
{
  return DT_VIEW_DARKROOM;
}

uint32_t container()
{
  return DT_UI_CONTAINER_PANEL_LEFT_CENTER;
}

void gui_reset(dt_lib_module_t *self)
{
}

int position()
{
  return 999;
}

static void _curve_resize(_curve_t *curve, float amount)
{
  //get the center of gravity of the form (like if it was a simple polygon)
  float bx = 0.0f;
  float by = 0.0f;
  float surf = 0.0f;
  for(int k = 0; k < curve->anchor_count; k++)
  {
    int k2 = (k+1)%curve->anchor_count;
    surf += curve->anchors[0][k]*curve->anchors[1][k2] - curve->anchors[0][k2]*curve->anchors[1][k];
    
    bx += (curve->anchors[0][k] + curve->anchors[0][k2])*(curve->anchors[0][k]*curve->anchors[1][k2] - curve->anchors[0][k2]*curve->anchors[1][k]);
    by += (curve->anchors[1][k] + curve->anchors[1][k2])*(curve->anchors[0][k]*curve->anchors[1][k2] - curve->anchors[0][k2]*curve->anchors[1][k]);
  }
  bx /= 3.0*surf;
  by /= 3.0*surf;
  
  //first, we have to be sure that the shape is not too small to be resized
  if (amount < 1.0)
  {
    for(int k = 0; k < curve->anchor_count; k++)
    {
      float l = (curve->anchors[0][k]-bx)*(curve->anchors[0][k]-bx) + (curve->anchors[1][k]-by)*(curve->anchors[1][k]-by);
      if ( l < 0.0005f) return;
    }
  }
  //now we move each point
  for(int k = 0; k < curve->anchor_count; k++)
  {
    float x = (curve->anchors[0][k]-bx)*amount;
    float y = (curve->anchors[1][k]-by)*amount;
    
    //we stretch ctrl points
    float ct1x = (curve->ctrl1[0][k]-curve->anchors[0][k])*amount;
    float ct1y = (curve->ctrl1[1][k]-curve->anchors[1][k])*amount;
    float ct2x = (curve->ctrl2[0][k]-curve->anchors[0][k])*amount;
    float ct2y = (curve->ctrl2[1][k]-curve->anchors[1][k])*amount;
    
    //and we set the new points
    curve->anchors[0][k] = bx + x;
    curve->anchors[1][k] = by + y;
    curve->ctrl1[0][k] = curve->anchors[0][k] + ct1x;
    curve->ctrl1[1][k] = curve->anchors[1][k] + ct1y;
    curve->ctrl2[0][k] = curve->anchors[0][k] + ct2x;
    curve->ctrl2[1][k] = curve->anchors[1][k] + ct2y;   
  }  
}

//Get the control points of a segment to match exactly a catmull-rom spline
static void catmull_to_bezier(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4,
                                float* bx1, float* by1, float* bx2, float* by2)
{
  *bx1 = (-x1 + 6*x2 + x3) / 6;
  *by1 = (-y1 + 6*y2 + y3) / 6;
  *bx2 = ( x2 + 6*x3 - x4) / 6;
  *by2 = ( y2 + 6*y3 - y4) / 6;
}

static void _curve_init_ctrl_points (_curve_t *curve, float wd, float ht)
{
  //if we have less that 3 points, what to do ??
  if (curve->anchor_count < 2)
  {
    return;
  }
  
  for(int k = 0; k < curve->anchor_count; k++)
  {
    //if the point as not be set manually, we redfine it
    if (curve->states[k] < 9)
    {
      int k1,k2,k3,k4,k5;
      
      k1 = (k-2)<0?curve->anchor_count+(k-2):k-2;
      k2 = (k-1)<0?curve->anchor_count-1:k-1;
      k3 = k;
      k4 = (k+1)%curve->anchor_count;
      k5 = (k+2)%curve->anchor_count;
      
      float bx1,by1,bx2,by2;
      catmull_to_bezier(curve->anchors[0][k1],curve->anchors[1][k1],
                        curve->anchors[0][k2],curve->anchors[1][k2],
                        curve->anchors[0][k3],curve->anchors[1][k3],
                        curve->anchors[0][k4],curve->anchors[1][k4],
                        &bx1,&by1,&bx2,&by2);
      if (curve->ctrl2[0][k2] == -1.0) curve->ctrl2[0][k2] = bx1;
      if (curve->ctrl2[1][k2] == -1.0) curve->ctrl2[1][k2] = by1;
      curve->ctrl1[0][k] = bx2;
      curve->ctrl1[1][k] = by2;
      catmull_to_bezier(curve->anchors[0][k2],curve->anchors[1][k2],
                        curve->anchors[0][k3],curve->anchors[1][k3],
                        curve->anchors[0][k4],curve->anchors[1][k4],
                        curve->anchors[0][k5],curve->anchors[1][k5],
                        &bx1,&by1,&bx2,&by2);
      if (curve->ctrl1[0][k4] == -1.0) curve->ctrl1[0][k4] = bx2;
      if (curve->ctrl1[1][k4] == -1.0) curve->ctrl1[1][k4] = by2;
      curve->ctrl2[0][k] = bx1;
      curve->ctrl2[1][k] = by1;
    }
  }

  //clockwise ???
  if (curve->anchor_count > 2)
  {
    float sum = 0.0f;
    for(int k = 0; k < curve->anchor_count; k++)
    {      
      int k2 = (k+1)%curve->anchor_count;
      //edge k
      sum += (curve->anchors[0][k2]-curve->anchors[0][k])*(curve->anchors[1][k2]+curve->anchors[1][k]);
    }
    curve->clockwise = (sum < 0);
  }
}

void gui_init(dt_lib_module_t *self)
{
  dt_lib_masks_t *d = (dt_lib_masks_t *)g_malloc(sizeof(dt_lib_masks_t));
  self->data = (void *)d;
  memset(d, 0, sizeof(dt_lib_masks_t));
  
  self->widget = gtk_label_new("The masks...");
}

void gui_cleanup(dt_lib_module_t *self)
{
  g_free(self->data);
  self->data = NULL;
}

#define CURVE_SEGMENT_RES 32

static int _curve_get_segment_near(_curve_t *curve, float ptx, float pty, float mini, gboolean *inside)
{
  int nb = 0;;
  float oldy = -1.0f;
  int rep = -1;
  for (int k=0; k < curve->anchor_count; k++)
  {
    int k0 = (k-1)<0?curve->anchor_count-1:k-1;
    
    if (oldy == -1.0f) oldy = curve->anchors[1][k0] - 0.0005f;

    //we determine the "length" of the curve (upper approx)
    float l = sqrtf((curve->ctrl2[0][k0]-curve->anchors[0][k0])*(curve->ctrl2[0][k0]-curve->anchors[0][k0]) + 
                    (curve->ctrl2[1][k0]-curve->anchors[1][k0])*(curve->ctrl2[1][k0]-curve->anchors[1][k0]));
    l+= sqrtf((curve->ctrl1[0][k]-curve->anchors[0][k])*(curve->ctrl1[0][k]-curve->anchors[0][k]) + 
                    (curve->ctrl1[1][k]-curve->anchors[1][k])*(curve->ctrl1[1][k]-curve->anchors[1][k]));
    l = l*0.8 + sqrtf((curve->anchors[0][k]-curve->anchors[0][k0])*(curve->anchors[0][k]-curve->anchors[0][k0]) + 
                    (curve->anchors[1][k]-curve->anchors[1][k0])*(curve->anchors[1][k]-curve->anchors[1][k0]));
    //so we get the nomber of iterations needed
    int iter = (int) (l*0.5/mini);

    //we go along the segment 
    for (int t=0; t < iter; t++)
    {
      //we can enhance perf here by avoiding to calculate both values each time
      float tt = t/(float)iter;
      float cx = curve->anchors[0][k0]*(1.0-tt)*(1.0-tt)*(1.0-tt) +
                  curve->ctrl2[0][k0]*3*tt*(1-tt)*(1-tt) +
                  curve->ctrl1[0][k]*3*tt*tt*(1-tt) +
                  curve->anchors[0][k]*tt*tt*tt;
      float cy = curve->anchors[1][k0]*(1.0-tt)*(1.0-tt)*(1.0-tt) +
                  curve->ctrl2[1][k0]*3*tt*(1-tt)*(1-tt) +
                  curve->ctrl1[1][k]*3*tt*tt*(1-tt) +
                  curve->anchors[1][k]*tt*tt*tt;
      if (ptx > cx-mini && ptx < cx+mini && pty > cy-mini && pty < cy+mini) rep = k;
      
      //for the "inside" part
      if (((pty > oldy && pty <= cy) || (pty < oldy && pty >= cy)) && ptx < cx) nb++;
      oldy = cy;
    }
  }

  *inside = (nb & 1);
  return rep;  
}

//feather calculating (must be in "real" coordinate, to be sure everything is orthonormal)
static void _curve_ctrl2_to_feather(int ptx,int pty, int ctrlx, int ctrly, int *fx, int *fy, gboolean clockwise)
{
  if (clockwise)
  {
    *fx = ptx + ctrly - pty;
    *fy = pty + ptx - ctrlx;
  }
  else
  {
    *fx = ptx - ctrly + pty;
    *fy = pty - ptx + ctrlx;
  }
}

static void _curve_feather_to_ctrl(int ptx,int pty, int fx, int fy, int *ctrl1x, int *ctrl1y, int *ctrl2x, int *ctrl2y, gboolean clockwise)
{
  if (clockwise)
  {
    *ctrl2x = ptx + pty - fy;
    *ctrl2y = pty + fx - ptx;
    *ctrl1x = ptx - pty + fy;
    *ctrl1y = pty - fx + ptx;
  }
  else
  {
    *ctrl1x = ptx + pty - fy;
    *ctrl1y = pty + fx - ptx;
    *ctrl2x = ptx - pty + fy;
    *ctrl2y = pty - fx + ptx;
  }
}

void gui_post_expose(dt_lib_module_t *self, cairo_t *cri, int32_t width, int32_t height, int32_t pointerx, int32_t pointery)
{
  dt_lib_masks_t *d = (dt_lib_masks_t *)self->data;
  d->viewport_width = width;
  d->viewport_height = height;

  /* get zoomed pointer coords*/
  float pzx,pzy;
  dt_dev_get_pointer_zoom_pos(darktable.develop, pointerx, pointery, &pzx, &pzy);
  pzx += 0.5f;
  pzy += 0.5f;
  
  /* get darkroom zoom and scale expose */
  int32_t zoom, closeup;
  float zoom_x, zoom_y;
  float wd = darktable.develop->preview_pipe->backbuf_width;
  float ht = darktable.develop->preview_pipe->backbuf_height;
  
  /* get zoom scale */
  DT_CTL_GET_GLOBAL(zoom_y, dev_zoom_y);
  DT_CTL_GET_GLOBAL(zoom_x, dev_zoom_x);
  DT_CTL_GET_GLOBAL(closeup, dev_closeup);
  DT_CTL_GET_GLOBAL(zoom, dev_zoom);
  float zoom_scale = dt_dev_get_zoom_scale(darktable.develop, zoom, closeup ? 2 : 1, 1);

  /* scale and translate cairo */
  cairo_translate(cri, width/2.0, height/2.0f);
  cairo_scale(cri, zoom_scale, zoom_scale);
  cairo_translate(cri, -.5f*wd-zoom_x*wd, -.5f*ht-zoom_y*ht);

  double dashed[] = {4.0, 2.0};
  dashed[0] /= zoom_scale;
  dashed[1] /= zoom_scale;
  int len  = sizeof(dashed) / sizeof(dashed[0]);
  
  cairo_set_line_width(cri, 1.25/zoom_scale);
  if (d->curve.anchor_count)
  {
    /* draw each curve segment */
    for (int k=0; k < d->curve.anchor_count; k++)
    {
      int k1,k2;
    
      k1 = (k-1)<0?d->curve.anchor_count-1:k-1;
      k2 = k;
      
      /* highlight segment */
      if (d->highlight_segment == k) cairo_set_source_rgba(cri, 1.0, 1.0, 1.0, 0.9);
      else cairo_set_source_rgba(cri, 1.0, 1.0, 1.0, 0.5);  

      cairo_move_to(cri, d->curve.anchors[0][k1]*wd, d->curve.anchors[1][k1]*ht);
      cairo_set_dash(cri, dashed, len, 0);
      if (d->curve.anchor_count > 2)
        cairo_curve_to(cri, d->curve.ctrl2[0][k1]*wd, d->curve.ctrl2[1][k1]*ht,
                            d->curve.ctrl1[0][k2]*wd, d->curve.ctrl1[1][k2]*ht,
                            d->curve.anchors[0][k2]*wd, d->curve.anchors[1][k2]*ht);
      else cairo_line_to(cri,d->curve.anchors[0][k2]*wd, d->curve.anchors[1][k2]*ht);
      cairo_stroke_preserve(cri);

      /* draw the segment shadow */
      cairo_set_dash(cri, dashed, len, 4);
      if (d->highlight_segment == k) cairo_set_source_rgba(cri, 0.0, 0.0, 0.0, 0.9);
      else cairo_set_source_rgba(cri, 0.0, 0.0, 0.0, 0.5);  
      
      cairo_stroke(cri);
    }
    
    /* draw the anchors */
    cairo_set_dash(cri, dashed, 0, 0);
    float anchor_size = 5.0f / zoom_scale;
    for(int k = 0; k < d->curve.anchor_count; k++)
    {
      if (k == d->highlight_anchor && !d->is_feather_highlight)
      {
        anchor_size = 6.0f / zoom_scale;
        cairo_set_source_rgba(cri, 1.0, 1.0, 1.0, 0.9);
      }
      else cairo_set_source_rgba(cri, 1.0, 1.0, 1.0, 0.5);
      
      cairo_rectangle(cri, 
		      (d->curve.anchors[0][k]*wd) - (anchor_size*0.5), 
		      (d->curve.anchors[1][k]*ht) - (anchor_size*0.5), 
		      anchor_size, anchor_size);
      cairo_fill_preserve(cri);

      if (k == d->highlight_anchor && !d->is_feather_highlight) cairo_set_source_rgba(cri, 0.0, 0.0, 0.0, 0.9);
      else cairo_set_source_rgba(cri, 0.0, 0.0, 0.0, 0.5);
      cairo_stroke(cri);
    }
    
    cairo_set_dash(cri, dashed, 0, 0);
    /* draw the feather weights normals */
    cairo_set_source_rgba(cri, 1.0, 1.0, 1.0, 0.3);
    for(int k = 0; k < d->curve.anchor_count; k++)
    {
      cairo_set_source_rgba(cri, 1.0, 1.0, 1.0, 0.3);
      //uncomment this part if you want to see "real" control points
      /*cairo_move_to(cri, d->curve.anchors[0][k]*wd,d->curve.anchors[1][k]*ht);
      cairo_line_to(cri, d->curve.ctrl1[0][k]*wd,d->curve.ctrl1[1][k]*ht);
      cairo_stroke(cri);
      cairo_move_to(cri, d->curve.anchors[0][k]*wd,d->curve.anchors[1][k]*ht);
      cairo_line_to(cri, d->curve.ctrl2[0][k]*wd,d->curve.ctrl2[1][k]*ht);
      cairo_stroke(cri);*/

      
      int ffx, ffy;
      _curve_ctrl2_to_feather(d->curve.anchors[0][k]*wd, d->curve.anchors[1][k]*ht, 
                              d->curve.ctrl2[0][k]*wd, d->curve.ctrl2[1][k]*ht, &ffx, &ffy, d->curve.clockwise);
      cairo_move_to(cri, d->curve.anchors[0][k]*wd,d->curve.anchors[1][k]*ht);
      cairo_line_to(cri,ffx,ffy);
      cairo_stroke(cri);
      if (d->highlight_anchor == k && d->is_feather_highlight) cairo_set_source_rgba(cri, 1.0, 1.0, 1.0, 0.9);
      cairo_arc (cri, ffx,ffy, 1.5f / zoom_scale, 0, 2.0*M_PI);
      cairo_fill_preserve(cri);

      if (k == d->highlight_anchor && d->is_feather_highlight) cairo_set_source_rgba(cri, 0.0, 0.0, 0.0, 0.9);
      else cairo_set_source_rgba(cri, 0.0, 0.0, 0.0, 0.5);
      cairo_stroke(cri);
    }
  }
}

int mouse_moved(dt_lib_module_t *self, double x, double y, int which)
{
  dt_lib_masks_t *d = (dt_lib_masks_t *)self->data;
  /* get zoomed pointer coords*/
  float pzx,pzy;
  dt_dev_get_pointer_zoom_pos(darktable.develop, x, y, &pzx, &pzy); 
  pzx += 0.5f;
  pzy += 0.5f;

  float wd = darktable.develop->preview_pipe->backbuf_width;
  float ht = darktable.develop->preview_pipe->backbuf_height;
  
  /* if dragging handle update */
  if (d->is_dragging_anchor)
  {
    if (d->is_feather_highlight)
    {
      int ct1x,ct1y,ct2x,ct2y;
      _curve_feather_to_ctrl(d->curve.anchors[0][d->drag_anchor]*wd, d->curve.anchors[1][d->drag_anchor]*ht, pzx*wd, pzy*ht,
                              &ct1x, &ct1y, &ct2x, &ct2y, d->curve.clockwise);
      d->curve.ctrl2[0][d->drag_anchor] = (float) ct2x/wd;
      d->curve.ctrl2[1][d->drag_anchor] = (float) ct2y/ht;
      d->curve.ctrl1[0][d->drag_anchor] = (float) ct1x/wd;
      d->curve.ctrl1[1][d->drag_anchor] = (float) ct1y/ht;
      d->curve.states[d->drag_anchor] = 9;
    }
    else
    {
      d->curve.ctrl1[0][d->drag_anchor] += pzx - d->curve.anchors[0][d->drag_anchor];
      d->curve.ctrl2[0][d->drag_anchor] += pzx - d->curve.anchors[0][d->drag_anchor];
      d->curve.ctrl1[1][d->drag_anchor] += pzy - d->curve.anchors[1][d->drag_anchor];
      d->curve.ctrl2[1][d->drag_anchor] += pzy - d->curve.anchors[1][d->drag_anchor];
      d->curve.anchors[0][d->drag_anchor] = pzx;
      d->curve.anchors[1][d->drag_anchor] = pzy;
      _curve_init_ctrl_points(&d->curve,wd,ht);
    }
    dt_control_queue_redraw_center();
    return 1;
  }
  
  /* reset any higlights */
  int seg = d->highlight_segment;
  int anch = d->highlight_anchor;
  gboolean feath = d->is_feather_highlight;
  d->highlight_segment = -1;
  d->highlight_anchor = -1;
  d->is_feather_highlight = FALSE;
  d->is_inside = FALSE;
  
  /* is mouse over anchor */
  double as = (2.0 / d->viewport_width) * 2.0; 
  for (int k=0;k<d->curve.anchor_count;k++)
  {
    //anchor
    if ( (pzx > d->curve.anchors[0][k] - as && pzx < d->curve.anchors[0][k] + as) && (pzy > d->curve.anchors[1][k] - as && pzy < d->curve.anchors[1][k] + as))
    {
      d->highlight_anchor = k;
      d->is_feather_highlight = FALSE;
    }
    //feather
    int fx, fy;
    _curve_ctrl2_to_feather(d->curve.anchors[0][k]*wd, d->curve.anchors[1][k]*ht, 
                            d->curve.ctrl2[0][k]*wd, d->curve.ctrl2[1][k]*ht, &fx, &fy, d->curve.clockwise);
    float ffx = fx/wd;
    float ffy = fy/ht;
    if ( (pzx > ffx - as && pzx < ffx + as) && (pzy > ffy - as && pzy < ffy + as))
    {
      d->highlight_anchor = k;
      d->is_feather_highlight = TRUE;
    }
  }


  /* is mouse near segment */
  d->highlight_segment = _curve_get_segment_near(&d->curve,pzx,pzy,as*1.5, &d->is_inside);
  
  //if some values change, we redraw
  if (seg != d->highlight_segment || anch != d->highlight_anchor || feath != d->is_feather_highlight)
  {
    dt_control_queue_redraw_center();
    return 1;    
  }
  return 0;
}


int button_released(struct dt_lib_module_t *self, double x, double y, int which, uint32_t state)
{
  dt_lib_masks_t *d = (dt_lib_masks_t *)self->data;

float wd = darktable.develop->preview_pipe->backbuf_width;
  float ht = darktable.develop->preview_pipe->backbuf_height;
  
  if(which == 3)
  {
    if (d->curve.anchor_count < 3)
    {
      /* clear curve */
      d->curve.anchor_count = 0;
      return 1;
    }
    else
    {
      /* remove anchor */
      if (d->current_anchor != -1)
      {
        int k = d->current_anchor;
        memcpy(&d->curve.anchors[0][k],&d->curve.anchors[0][k+1], sizeof(float) * d->curve.anchor_count-k);
        memcpy(&d->curve.anchors[1][k],&d->curve.anchors[1][k+1], sizeof(float) * d->curve.anchor_count-k);
        d->curve.anchor_count--;
        _curve_init_ctrl_points(&d->curve,wd,ht);
        return 1;
      }
    }
  }

  /* end dragging */
  if (d->is_dragging_anchor && !d->is_curve_edit_mode)
  {
    d->is_dragging_anchor = FALSE;
    return 1;
  }

  return 0;
}

int button_pressed (struct dt_lib_module_t *self, double x, double y, int which, int type, uint32_t state)
{
  dt_lib_masks_t *d = (dt_lib_masks_t *)self->data;

  float wd = darktable.develop->preview_pipe->backbuf_width;
  float ht = darktable.develop->preview_pipe->backbuf_height;

  /* get zoomed pointer coords*/
  float pzx,pzy;
  dt_dev_get_pointer_zoom_pos(darktable.develop, x, y, &pzx, &pzy);
  pzx += 0.5f;
  pzy += 0.5f;

  if (d->curve.anchor_count < 3)
    d->is_curve_edit_mode = TRUE;

  if (d->is_curve_edit_mode && which == 1)
  {
    int pos = d->curve.anchor_count==0 ? 0 : d->curve.anchor_count-1;
    /* adding new anchor to curve */
    d->curve.anchors[0][pos] = pzx;
    d->curve.anchors[1][pos] = pzy;
    if (d->curve.anchor_count == 0) d->curve.anchor_count++;
    d->curve.anchor_count++;
    d->curve.anchors[0][d->curve.anchor_count-1] = pzx;
    d->curve.anchors[1][d->curve.anchor_count-1] = pzy;
    d->is_dragging_anchor = TRUE;
    d->current_anchor = d->curve.anchor_count-1;
    d->drag_anchor = d->curve.anchor_count-1;
    d->is_feather_highlight = FALSE;
    _curve_init_ctrl_points(&d->curve,wd,ht);
    dt_control_queue_redraw_center();
    return 1;
  } 
  else
  {

    /* check if we should end editmode */
    if (which == 3 && d->curve.anchor_count >= 3)
    {
      d->is_curve_edit_mode = FALSE;
      d->is_dragging_anchor = FALSE;
      dt_control_queue_redraw_center();
      return 1;
    }

    /* is mouse over anchor */
    d->current_anchor = -1;
    if (d->highlight_anchor >= 0 && which == 1)
    {
      d->current_anchor = d->highlight_anchor;
      d->drag_anchor = d->highlight_anchor;
      d->is_dragging_anchor = TRUE;
      dt_control_queue_redraw_center();
      return 1;
    }

    /* is mouse near segment lets insert anchor and begin drag */
    if (which ==1 && d->highlight_segment != -1)
    {
      int k = d->highlight_segment;
      
      /* insert anchor on segment at px,py */
      d->curve.anchor_count++;
      for (int l = d->curve.anchor_count-1; l >= k; l--)
      {
        d->curve.anchors[0][l] = d->curve.anchors[0][l-1]; 
        d->curve.anchors[1][l] = d->curve.anchors[1][l-1];
        d->curve.ctrl1[0][l] = d->curve.ctrl1[0][l-1]; 
        d->curve.ctrl1[1][l] = d->curve.ctrl1[1][l-1];
        d->curve.ctrl2[0][l] = d->curve.ctrl2[0][l-1]; 
        d->curve.ctrl2[1][l] = d->curve.ctrl2[1][l-1];
        d->curve.states[l] = d->curve.states[l-1]; 
      }
      d->curve.anchors[0][k] = pzx;
      d->curve.anchors[1][k] = pzy;
      d->curve.ctrl1[0][k] = pzx; 
      d->curve.ctrl1[1][k] = pzy;
      d->curve.ctrl2[0][k] = pzx; 
      d->curve.ctrl2[1][k] = pzy;
      d->curve.states[k] = 0;

      d->drag_anchor = k;
      d->is_dragging_anchor = TRUE;
      
      _curve_init_ctrl_points(&d->curve,wd,ht);
      dt_control_queue_redraw_center();
      return 1;
    }
  }

  dt_control_queue_redraw_center();

  return 0;
}

int scrolled(struct dt_lib_module_t *self, double x, double y, int up, uint32_t state)
{
  dt_lib_masks_t *d = (dt_lib_masks_t *)self->data;
  if (!d->is_inside) return 0;
  double as = 1.05;
  if (!up) as = 0.95;
  _curve_resize(&d->curve,as);
  dt_control_queue_redraw_center();
  return 1;
}
