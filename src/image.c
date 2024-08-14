/*****************************************************************************
 * This file is part of uvg266 VVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "image.h"

#include <limits.h>
#include <stdlib.h>

#include "strategies/strategies-ipol.h"
#include "strategies/strategies-picture.h"
#include "threads.h"

/**
* \brief Allocate a new image with 420.
* This function signature is part of the libkvz API.
* \return image pointer or NULL on failure
*/
uvg_picture * uvg_image_alloc_420(const int32_t width, const int32_t height)
{
  return uvg_image_alloc(UVG_CSP_420, width, height);
}

/**
 * \brief Allocate a new image.
 * \return image pointer or NULL on failure
 */
uvg_picture * uvg_image_alloc(enum uvg_chroma_format chroma_format, const int32_t width, const int32_t height)
{
  //Assert that we have a well defined image
  assert((width % 2) == 0);
  assert((height % 2) == 0);

  const size_t simd_padding_width = 64;

  uvg_picture *im = MALLOC(uvg_picture, 1);
  if (!im) return NULL;

  //Add 4 pixel boundary to each side of luma for ALF
  //This results also 2 pixel boundary for chroma
  unsigned int luma_size = (width + FRAME_PADDING_LUMA) * (height + FRAME_PADDING_LUMA);

  unsigned chroma_sizes[] = { 0, luma_size / 4, luma_size / 2, luma_size };
  unsigned chroma_size = chroma_sizes[chroma_format];

  im->chroma_format = chroma_format;

  //Allocate memory, pad the full data buffer from both ends
  im->fulldata_buf = MALLOC_SIMD_PADDED(uvg_pixel, (luma_size + 2 * chroma_size), simd_padding_width * 2);
  if (!im->fulldata_buf) {
    free(im);
    return NULL;
  }

  //Shift the image to allow ALF filtering
  im->refcount = 1; //We give a reference to caller
  im->width = width;
  im->height = height;
  im->stride = width + FRAME_PADDING_LUMA;
  im->chroma_format = chroma_format;
  const int padding_before_first_pixel_luma = (FRAME_PADDING_LUMA / 2) * (im->stride) + FRAME_PADDING_LUMA / 2;
  const int padding_before_first_pixel_chroma = (FRAME_PADDING_CHROMA / 2) * (im->stride/2) + FRAME_PADDING_CHROMA / 2;
  im->fulldata = &im->fulldata_buf[padding_before_first_pixel_luma] + simd_padding_width / sizeof(uvg_pixel);
  im->base_image = im;

  im->y = im->data[COLOR_Y] = &im->fulldata[0];

  if (chroma_format == UVG_CSP_400) {
    im->u = im->data[COLOR_U] = NULL;
    im->v = im->data[COLOR_V] = NULL;
  } else {
    im->u = im->data[COLOR_U] = &im->fulldata[luma_size - padding_before_first_pixel_luma + padding_before_first_pixel_chroma];
    im->v = im->data[COLOR_V] = &im->fulldata[luma_size - padding_before_first_pixel_luma + chroma_size + padding_before_first_pixel_chroma];
  }

  im->pts = 0;
  im->dts = 0;

  im->interlacing = UVG_INTERLACING_NONE;

  im->roi.roi_array = NULL;
  im->roi.width = 0;
  im->roi.height = 0;

  return im;
}

/**
 * \brief Free an image.
 *
 * Decrement reference count of the image and deallocate associated memory
 * if no references exist any more.
 *
 * \param im image to free
 */
void uvg_image_free(uvg_picture *const im)
{
  if (im == NULL) return;

  int32_t new_refcount = UVG_ATOMIC_DEC(&(im->refcount));
  if (new_refcount > 0) {
    // There are still references so we don't free the data yet.
    return;
  }

  if (im->base_image != im) {
    // Free our reference to the base image.
    uvg_image_free(im->base_image);
  } else {
    free(im->fulldata_buf);
    if (im->roi.roi_array) FREE_POINTER(im->roi.roi_array);
  }

  // Make sure freed data won't be used.
  im->base_image = NULL;
  im->fulldata_buf = NULL;
  im->fulldata = NULL;
  im->y = im->u = im->v = NULL;
  im->data[COLOR_Y] = im->data[COLOR_U] = im->data[COLOR_V] = NULL;
  free(im);
}

/**
 * \brief Get a new pointer to an image.
 *
 * Increment reference count and return the image.
 */
uvg_picture *uvg_image_copy_ref(uvg_picture *im)
{
  int32_t new_refcount = UVG_ATOMIC_INC(&im->refcount);
  // The caller should have had another reference and we added one
  // reference so refcount should be at least 2.
  assert(new_refcount >= 2);
  return im;
}

uvg_picture *uvg_image_make_subimage(uvg_picture *const orig_image,
                             const unsigned x_offset,
                             const unsigned y_offset,
                             const unsigned width,
                             const unsigned height)
{
  // Assert that we have a well defined image
  assert((width % 2) == 0);
  assert((height % 2) == 0);

  assert((x_offset % 2) == 0);
  assert((y_offset % 2) == 0);

  assert(x_offset + width <= orig_image->width);
  assert(y_offset + height <= orig_image->height);

  uvg_picture *im = MALLOC(uvg_picture, 1);
  if (!im) return NULL;

  im->base_image = uvg_image_copy_ref(orig_image->base_image);
  im->refcount = 1; // We give a reference to caller
  im->width = width;
  im->height = height;
  im->stride = orig_image->stride;
  im->chroma_format = orig_image->chroma_format;

  im->y = im->data[COLOR_Y] = &orig_image->y[x_offset + y_offset * orig_image->stride];
  if (orig_image->chroma_format != UVG_CSP_400) {
    im->u = im->data[COLOR_U] = &orig_image->u[x_offset / 2 + y_offset / 2 * orig_image->stride / 2];
    im->v = im->data[COLOR_V] = &orig_image->v[x_offset / 2 + y_offset / 2 * orig_image->stride / 2];
  }

  im->pts = 0;
  im->dts = 0;

  im->roi = orig_image->roi;

  return im;
}

yuv_t * uvg_yuv_t_alloc(int luma_size, int chroma_size)
{
  yuv_t *yuv = (yuv_t *)malloc(sizeof(*yuv));
  yuv->size = luma_size;

  // Get buffers with separate mallocs in order to take advantage of
  // automatic buffer overrun checks.
  yuv->y = (uvg_pixel *)malloc(luma_size * sizeof(*yuv->y));
  if (chroma_size == 0) {
    yuv->u = NULL;
    yuv->v = NULL;
  } else {
    yuv->u = (uvg_pixel *)malloc(chroma_size * sizeof(*yuv->u));
    yuv->v = (uvg_pixel *)malloc(chroma_size * sizeof(*yuv->v));
  }
  
  return yuv;
}

void uvg_yuv_t_free(yuv_t *yuv)
{
  if (yuv) {
    FREE_POINTER(yuv->y);
    FREE_POINTER(yuv->u);
    FREE_POINTER(yuv->v);
  }
  FREE_POINTER(yuv);
}

hi_prec_buf_t * uvg_hi_prec_buf_t_alloc(int luma_size)
{
  // Get buffers with separate mallocs in order to take advantage of
  // automatic buffer overrun checks.
  hi_prec_buf_t *yuv = (hi_prec_buf_t *)malloc(sizeof(*yuv));
  yuv->y = (int16_t *)malloc(luma_size * sizeof(*yuv->y));
  yuv->u = (int16_t *)malloc(luma_size / 2 * sizeof(*yuv->u));
  yuv->v = (int16_t *)malloc(luma_size / 2 * sizeof(*yuv->v));
  yuv->joint_u = (int16_t *)malloc(luma_size / 2 * sizeof(*yuv->u));
  yuv->joint_v = (int16_t *)malloc(luma_size / 2 * sizeof(*yuv->v));
  yuv->size = luma_size;

  return yuv;
}

void uvg_hi_prec_buf_t_free(hi_prec_buf_t * yuv)
{
  free(yuv->y);
  free(yuv->u);
  free(yuv->v);
  free(yuv->joint_v);
  free(yuv->joint_u);
  free(yuv);
}

static INLINE uint32_t reg_sad_maybe_optimized(const uvg_pixel * const data1, const uvg_pixel * const data2,
                                  const int32_t width, const int32_t height, const uint32_t stride1,
                                  const uint32_t stride2, optimized_sad_func_ptr_t optimized_sad)
{
  if (optimized_sad != NULL)
    return optimized_sad(data1, data2, height, stride1, stride2);
  else
    return uvg_reg_sad(data1, data2, width, height, stride1, stride2);
}

/**
 * \brief Diagonally interpolate SAD outside the frame.
 *
 * \param data1   Starting point of the first picture.
 * \param data2   Starting point of the second picture.
 * \param width   Width of the region for which SAD is calculated.
 * \param height  Height of the region for which SAD is calculated.
 * \param width  Width of the pixel array.
 *
 * \returns Sum of Absolute Differences
 */
static unsigned cor_sad(const uvg_pixel *pic_data, const uvg_pixel *ref_data,
                        int block_width, int block_height, unsigned pic_stride)
{
  uvg_pixel ref = *ref_data;
  int x, y;
  unsigned sad = 0;

  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sad += abs(pic_data[y * pic_stride + x] - ref);
    }
  }

  return sad;
}


/**
 * \brief  Handle special cases of comparing blocks that are not completely
 *         inside the frame.
 *
 * \param pic  First frame.
 * \param ref  Second frame.
 * \param pic_x  X coordinate of the first block.
 * \param pic_y  Y coordinate of the first block.
 * \param ref_x  X coordinate of the second block.
 * \param ref_y  Y coordinate of the second block.
 * \param block_width  Width of the blocks.
 * \param block_height  Height of the blocks.
 */
static unsigned image_interpolated_sad(const uvg_picture *pic, const uvg_picture *ref,
                                 int pic_x, int pic_y, int ref_x, int ref_y,
                                 int block_width, int block_height,
                                 optimized_sad_func_ptr_t optimized_sad)
{
  uvg_pixel *pic_data, *ref_data;

  int left, right, top, bottom;
  int result = 0;

  // Change the movement vector to point right next to the frame. This doesn't
  // affect the result but removes some special cases.
  if (ref_x > ref->width)            ref_x = ref->width;
  if (ref_y > ref->height)           ref_y = ref->height;
  if (ref_x + block_width < 0)  ref_x = -block_width;
  if (ref_y + block_height < 0) ref_y = -block_height;

  // These are the number of pixels by how far the movement vector points
  // outside the frame. They are always >= 0. If all of them are 0, the
  // movement vector doesn't point outside the frame.
  left   = (ref_x < 0) ? -ref_x : 0;
  top    = (ref_y < 0) ? -ref_y : 0;
  right  = (ref_x + block_width  > ref->width)  ? ref_x + block_width  - ref->width  : 0;
  bottom = (ref_y + block_height > ref->height) ? ref_y + block_height - ref->height : 0;

  // Center picture to the current block and reference to the point where
  // movement vector is pointing to. That point might be outside the buffer,
  // but that is ok because we project the movement vector to the buffer
  // before dereferencing the pointer.
  pic_data = &pic->y[pic_y * pic->stride + pic_x];
  ref_data = &ref->y[ref_y * ref->stride + ref_x];

  // The handling of movement vectors that point outside the picture is done
  // in the following way.
  // - Correct the index of ref_data so that it points to the top-left
  //   of the area we want to compare against.
  // - Correct the index of pic_data to point inside the current block, so
  //   that we compare the right part of the block to the ref_data.
  // - Reduce block_width and block_height so that the the size of the area
  //   being compared is correct.
  //
  // NOTE: No more correct since hor_sad was modified to be a separate
  // strategy
  if (top && left) {
    result += cor_sad(pic_data,
                      &ref_data[top * ref->stride + left],
                      left, top, pic->stride);
    result += uvg_ver_sad(&pic_data[left],
                      &ref_data[top * ref->stride + left],
                      block_width - left, top, pic->stride);

    result += uvg_hor_sad(pic_data + top * pic->stride,
                          ref_data + top * ref->stride,
                          block_width, block_height - top,
                          pic->stride, ref->stride,
                          left, right);

  } else if (top && right) {
    result += uvg_ver_sad(pic_data,
                      &ref_data[top * ref->stride],
                      block_width - right, top, pic->stride);
    result += cor_sad(&pic_data[block_width - right],
                      &ref_data[top * ref->stride + (block_width - right - 1)],
                      right, top, pic->stride);

    result += uvg_hor_sad(pic_data + top * pic->stride,
                          ref_data + top * ref->stride,
                          block_width, block_height - top,
                          pic->stride, ref->stride,
                          left, right);

  } else if (bottom && left) {
    result += uvg_hor_sad(pic_data, ref_data, block_width, block_height - bottom,
                          pic->stride, ref->stride, left, right);

    result += cor_sad(&pic_data[(block_height - bottom) * pic->stride],
                      &ref_data[(block_height - bottom - 1) * ref->stride + left],
                      left, bottom, pic->stride);
    result += uvg_ver_sad(&pic_data[(block_height - bottom) * pic->stride + left],
                      &ref_data[(block_height - bottom - 1) * ref->stride + left],
                      block_width - left, bottom, pic->stride);
  } else if (bottom && right) {
    result += uvg_hor_sad(pic_data, ref_data, block_width, block_height - bottom,
                          pic->stride, ref->stride, left, right);

    result += uvg_ver_sad(&pic_data[(block_height - bottom) * pic->stride],
                      &ref_data[(block_height - bottom - 1) * ref->stride],
                      block_width - right, bottom, pic->stride);
    result += cor_sad(&pic_data[(block_height - bottom) * pic->stride + block_width - right],
                      &ref_data[(block_height - bottom - 1) * ref->stride + block_width - right - 1],
                      right, bottom, pic->stride);
  } else if (top) {
    result += uvg_ver_sad(pic_data,
                      &ref_data[top * ref->stride],
                      block_width, top, pic->stride);
    result += reg_sad_maybe_optimized(&pic_data[top * pic->stride],
                      &ref_data[top * ref->stride],
                      block_width, block_height - top, pic->stride, ref->stride,
                      optimized_sad);
  } else if (bottom) {
    result += reg_sad_maybe_optimized(pic_data,
                      ref_data,
                      block_width, block_height - bottom, pic->stride, ref->stride,
                      optimized_sad);
    result += uvg_ver_sad(&pic_data[(block_height - bottom) * pic->stride],
                      &ref_data[(block_height - bottom - 1) * ref->stride],
                      block_width, bottom, pic->stride);
  } else if (left | right) {
    result += uvg_hor_sad(pic_data, ref_data,
                          block_width, block_height, pic->stride,
                          ref->stride, left, right);
  } else {
    result += reg_sad_maybe_optimized(pic_data, ref_data,
                                      block_width, block_height,
                                      pic->stride, ref->stride,
                                      optimized_sad);
  }
  return result;
}

/**
* \brief Calculate interpolated SAD between two blocks.
*
* \param pic        Image for the block we are trying to find.
* \param ref        Image where we are trying to find the block.
*
* \returns          Sum of absolute differences
*/
unsigned uvg_image_calc_sad(const uvg_picture *pic,
                            const uvg_picture *ref,
                            int pic_x,
                            int pic_y,
                            int ref_x,
                            int ref_y,
                            int block_width,
                            int block_height,
                            optimized_sad_func_ptr_t optimized_sad)
{
  assert(pic_x >= 0 && pic_x <= pic->width - block_width);
  assert(pic_y >= 0 && pic_y <= pic->height - block_height);

  uint32_t res;

  if (ref_x >= 0 && ref_x <= ref->width  - block_width &&
      ref_y >= 0 && ref_y <= ref->height - block_height)
  {
    // Reference block is completely inside the frame, so just calculate the
    // SAD directly. This is the most common case, which is why it's first.
    const uvg_pixel *pic_data = &pic->y[pic_y * pic->stride + pic_x];
    const uvg_pixel *ref_data = &ref->y[ref_y * ref->stride + ref_x];

    res = reg_sad_maybe_optimized(pic_data,
                                  ref_data,
                                  block_width,
                                  block_height,
                                  pic->stride,
                                  ref->stride,
                                  optimized_sad);
  } else {
    // Call a routine that knows how to interpolate pixels outside the frame.
    res = image_interpolated_sad(pic, ref, pic_x, pic_y, ref_x, ref_y, block_width, block_height, optimized_sad);
  }
  return res >> (UVG_BIT_DEPTH - 8);
}


/**
* \brief Calculate interpolated SATD between two blocks.
*
* \param pic        Image for the block we are trying to find.
* \param ref        Image where we are trying to find the block.
*/
unsigned uvg_image_calc_satd(const uvg_picture *pic,
                             const uvg_picture *ref,
                             int pic_x,
                             int pic_y,
                             int ref_x,
                             int ref_y,
                             int block_width,
                             int block_height,
                             uint8_t ref_wraparound)
{
  assert(pic_x >= 0 && pic_x <= pic->width - block_width);
  assert(pic_y >= 0 && pic_y <= pic->height - block_height);

  if (ref_x >= 0 && ref_x <= ref->width  - block_width &&
      ref_y >= 0 && ref_y <= ref->height - block_height)
  {
    // Reference block is completely inside the frame, so just calculate the
    // SAD directly. This is the most common case, which is why it's first.
    const uvg_pixel *pic_data = &pic->y[pic_y * pic->stride + pic_x];
    const uvg_pixel *ref_data = &ref->y[ref_y * ref->stride + ref_x];
    return uvg_satd_any_size(block_width,
                             block_height,
                             pic_data,
                             pic->stride,
                             ref_data,
                             ref->stride);
  } else {
    // Extrapolate pixels from outside the frame.

    // Space for extrapolated pixels and the part from the picture
    // The extrapolation function will set the pointers and stride.
    uvg_pixel ext_buffer[LCU_LUMA_SIZE];
    uvg_pixel *ext = NULL;
    uvg_pixel *ext_origin = NULL;
    int ext_s = 0;
    uvg_epol_args epol_args = {
      .src = ref->y,
      .src_w = ref->width,
      .src_h = ref->height,
      .src_s = ref->stride,
      .blk_x = ref_x,
      .blk_y = ref_y,
      .blk_w = block_width,
      .blk_h = block_height,
      .pad_l = 0,
      .pad_r = 0,
      .pad_t = 0,
      .pad_b = 0,
      .pad_b_simd = 0,
    };

    // Initialize separately. Gets rid of warning
    // about using nonstandard extension.
    epol_args.buf = ext_buffer;
    epol_args.ext = &ext;
    epol_args.ext_origin = &ext_origin;
    epol_args.ext_s = &ext_s;

    if (ref_wraparound) {
      uvg_get_extended_block_wraparound(&epol_args);
    } else {
      uvg_get_extended_block(&epol_args);
    }

    const uvg_pixel *pic_data = &pic->y[pic_y * pic->stride + pic_x];

    unsigned satd = uvg_satd_any_size(block_width,
      block_height,
      pic_data,
      pic->stride,
      ext_origin,
      ext_s);

    return satd;
  }
}




/**
 * \brief BLock Image Transfer from one buffer to another.
 *
 * It's a stupidly simple loop that copies pixels.
 *
 * \param orig  Start of the originating buffer.
 * \param dst  Start of the destination buffer.
 * \param width  Width of the copied region.
 * \param height  Height of the copied region.
 * \param orig_stride  Width of a row in the originating buffer.
 * \param dst_stride  Width of a row in the destination buffer.
 *
 * This should be inlined, but it's defined here for now to see if Visual
 * Studios LTCG will inline it.
 */
#define BLIT_PIXELS_CASE(n) case n:\
  for (y = 0; y < n; ++y) {\
    memcpy(&dst[y*dst_stride], &orig[y*orig_stride], n * sizeof(uvg_pixel));\
  }\
  break;

void uvg_pixels_blit(const uvg_pixel * const orig, uvg_pixel * const dst,
                         const unsigned width, const unsigned height,
                         const unsigned orig_stride, const unsigned dst_stride)
{
  unsigned y;
  //There is absolutely no reason to have a width greater than the source or the destination stride.
  assert(width <= orig_stride);
  assert(width <= dst_stride);

#ifdef CHECKPOINTS
  char *buffer = malloc((3 * width + 1) * sizeof(char));
  for (y = 0; y < height; ++y) {
    int p;
    for (p = 0; p < width; ++p) {
      sprintf((buffer + 3*p), "%02X ", orig[y*orig_stride]);
    }
    buffer[3*width] = 0;
    CHECKPOINT("uvg_pixels_blit_avx2: %04d: %s", y, buffer);
  }
  FREE_POINTER(buffer);
#endif //CHECKPOINTS

  if (width == orig_stride && width == dst_stride) {
    memcpy(dst, orig, width * height * sizeof(uvg_pixel));
    return;
  }

  int nxn_width = (width == height) ? width : 0;
  switch (nxn_width) {
    BLIT_PIXELS_CASE(4)
    BLIT_PIXELS_CASE(8)
    BLIT_PIXELS_CASE(16)
    BLIT_PIXELS_CASE(32)
    BLIT_PIXELS_CASE(64)
  default:

    if (orig == dst) {
      //If we have the same array, then we should have the same stride
      assert(orig_stride == dst_stride);
      return;
    }
    assert(orig != dst || orig_stride == dst_stride);

    for (y = 0; y < height; ++y) {
      memcpy(&dst[y*dst_stride], &orig[y*orig_stride], width * sizeof(uvg_pixel));
    }
    break;
  }
}
