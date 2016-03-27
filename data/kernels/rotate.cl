/*
    This file is part of darktable,
    copyright (c) 2016 aldric renaudin.

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

#include "common.h"

void
mul_mat_vec_2(const float4 m, const float2 *p, float2 *o)
{
  (*o).x = (*p).x*m.x + (*p).y*m.y;
  (*o).y = (*p).x*m.z + (*p).y*m.w;
}

void
backtransform(float2 *p, float2 *o, const float4 m, const float2 k)
{
  (*p).y /= (1.0f + (*p).x*k.x);
  (*p).x /= (1.0f + (*p).y*k.y);
  mul_mat_vec_2(m, p, o);
}


float
interpolation_func_bicubic(float t)
{
  float r;
  t = fabs(t);

  r = (t >= 2.0f) ? 0.0f : ((t > 1.0f) ? (0.5f*(t*(-t*t + 5.0f*t - 8.0f) + 4.0f)) : (0.5f*(t*(3.0f*t*t - 5.0f*t) + 2.0f)));

  return r;
}

#define DT_LANCZOS_EPSILON (1e-9f)

#if 0
float
interpolation_func_lanczos(float width, float t)
{
  float ta = fabs(t);

  float r = (ta > width) ? 0.0f : ((ta < DT_LANCZOS_EPSILON) ? 1.0f : width*native_sin(M_PI_F*t)*native_sin(M_PI_F*t/width)/(M_PI_F*M_PI_F*t*t));

  return r;
}
#else
float
sinf_fast(float t)
{
  const float a = 4.0f/(M_PI_F*M_PI_F);
  const float p = 0.225f;

  t = a*t*(M_PI_F - fabs(t));

  return p*(t*fabs(t) - t) + t;
}

float
interpolation_func_lanczos(float width, float t)
{
  /* Compute a value for sinf(pi.t) in [-pi pi] for which the value will be
   * correct */
  int a = (int)t;
  float r = t - (float)a;

  // Compute the correct sign for sinf(pi.r)
  union { float f; unsigned int i; } sign;
  sign.i = ((a&1)<<31) | 0x3f800000;

  return (DT_LANCZOS_EPSILON + width*sign.f*sinf_fast(M_PI_F*r)*sinf_fast(M_PI_F*t/width))/(DT_LANCZOS_EPSILON + M_PI_F*M_PI_F*t*t);
}
#endif


/* kernel for clip&rotate: bilinear interpolation */
__kernel void
rotate_bilinear(read_only image2d_t in, write_only image2d_t out, const int width, const int height,
            const int in_width, const int in_height,
            const int2 roi_in, const float2 roi_out, const float scale_in, const float scale_out,
            const float2 center_in, const float2 center_out, const float2 k, const float4 mat)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);

  if(x >= width || y >= height) return;

  float2 pi, po;

  pi.x = roi_out.x + x + 0.5f;
  pi.y = roi_out.y + y + 0.5f;

  pi.x -= center_out.x;
  pi.y -= center_out.y;

  pi /= scale_out;
  backtransform(&pi, &po, mat, k);
  po *= scale_in;

  po.x += center_in.x;
  po.y += center_in.y;

  po.x -= roi_in.x + 0.5f;
  po.y -= roi_in.y + 0.5f;

  const int ii = (int)po.x;
  const int jj = (int)po.y;

  float4 o = (ii >=0 && jj >= 0 && ii <= in_width-2 && jj <= in_height-2) ? read_imagef(in, samplerf, po) : (float4)0.0f;

  write_imagef (out, (int2)(x, y), o);
}



/* kernel for clip&rotate: bicubic interpolation */
__kernel void
rotate_bicubic(read_only image2d_t in, write_only image2d_t out, const int width, const int height,
            const int in_width, const int in_height,
            const int2 roi_in, const float2 roi_out, const float scale_in, const float scale_out,
            const float2 center_in, const float2 center_out, const float2 k, const float4 mat)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);

  const int kwidth = 2;

  if(x >= width || y >= height) return;

  float2 pi, po;

  pi.x = roi_out.x + x + 0.5f;
  pi.y = roi_out.y + y + 0.5f;

  pi.x -= center_out.x;
  pi.y -= center_out.y;

  pi /= scale_out;
  backtransform(&pi, &po, mat, k);
  po *= scale_in;

  po.x += center_in.x;
  po.y += center_in.y;

  po.x -= roi_in.x + 0.5f;
  po.y -= roi_in.y + 0.5f;

  int tx = po.x;
  int ty = po.y;

  float4 pixel = (float4)0.0f;
  float weight = 0.0f;

  for(int jj = 1 - kwidth; jj <= kwidth; jj++)
    for(int ii= 1 - kwidth; ii <= kwidth; ii++)
  {
    const int i = tx + ii;
    const int j = ty + jj;

    float wx = interpolation_func_bicubic((float)i - po.x);
    float wy = interpolation_func_bicubic((float)j - po.y);
    float w = wx * wy;

    pixel += read_imagef(in, samplerc, (int2)(i, j)) * w;
    weight += w;
  }

  pixel = weight > 0.0f ? pixel / weight : (float4)0.0f;

  write_imagef (out, (int2)(x, y), pixel);
}


/* kernel for clip&rotate: lanczos2 interpolation */
__kernel void
rotate_lanczos2(read_only image2d_t in, write_only image2d_t out, const int width, const int height,
            const int in_width, const int in_height,
            const int2 roi_in, const float2 roi_out, const float scale_in, const float scale_out,
            const float2 center_in, const float2 center_out, const float2 k, const float4 mat)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);

  const int kwidth = 2;

  if(x >= width || y >= height) return;

  float2 pi, po;

  pi.x = roi_out.x + x + 0.5f;
  pi.y = roi_out.y + y + 0.5f;

  pi.x -= center_out.x;
  pi.y -= center_out.y;

  pi /= scale_out;
  backtransform(&pi, &po, mat, k);
  po *= scale_in;

  po.x += center_in.x;
  po.y += center_in.y;

  po.x -= roi_in.x + 0.5f;
  po.y -= roi_in.y + 0.5f;

  int tx = po.x;
  int ty = po.y;

  float4 pixel = (float4)0.0f;
  float weight = 0.0f;

  for(int jj = 1 - kwidth; jj <= kwidth; jj++)
    for(int ii= 1 - kwidth; ii <= kwidth; ii++)
  {
    const int i = tx + ii;
    const int j = ty + jj;

    float wx = interpolation_func_lanczos(2, (float)i - po.x);
    float wy = interpolation_func_lanczos(2, (float)j - po.y);
    float w = wx * wy;

    pixel += read_imagef(in, samplerc, (int2)(i, j)) * w;
    weight += w;
  }

  pixel = weight > 0.0f ? pixel / weight : (float4)0.0f;

  write_imagef (out, (int2)(x, y), pixel);
}



/* kernel for clip&rotate: lanczos3 interpolation */
__kernel void
rotate_lanczos3(read_only image2d_t in, write_only image2d_t out, const int width, const int height,
            const int in_width, const int in_height,
            const int2 roi_in, const float2 roi_out, const float scale_in, const float scale_out,
            const float2 center_in, const float2 center_out, const float2 k, const float4 mat)
{
  const int x = get_global_id(0);
  const int y = get_global_id(1);

  const int kwidth = 3;

  if(x >= width || y >= height) return;

  float2 pi, po;

  pi.x = roi_out.x + x + 0.5f;
  pi.y = roi_out.y + y + 0.5f;

  pi.x -= center_out.x;
  pi.y -= center_out.y;

  pi /= scale_out;
  backtransform(&pi, &po, mat, k);
  po *= scale_in;

  po.x += center_in.x;
  po.y += center_in.y;

  po.x -= roi_in.x + 0.5f;
  po.y -= roi_in.y + 0.5f;

  int tx = (int)po.x;
  int ty = (int)po.y;

  float4 pixel = (float4)0.0f;
  float weight = 0.0f;

  for(int jj = 1 - kwidth; jj <= kwidth; jj++)
    for(int ii= 1 - kwidth; ii <= kwidth; ii++)
  {
    const int i = tx + ii;
    const int j = ty + jj;

    float wx = interpolation_func_lanczos(3, (float)i - po.x);
    float wy = interpolation_func_lanczos(3, (float)j - po.y);
    float w = wx * wy;

    pixel += read_imagef(in, samplerc, (int2)(i, j)) * w;
    weight += w;
  }

  pixel = weight > 0.0f ? pixel / weight : (float4)0.0f;

  write_imagef (out, (int2)(x, y), pixel);
}
