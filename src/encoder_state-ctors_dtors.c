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

#include "encoder_state-ctors_dtors.h"

#include <stdio.h>
#include <stdlib.h>

#include "bitstream.h"
#include "cabac.h"
#include "cu.h"
#include "debug.h"
#include "encoder.h"
#include "encoder_state-geometry.h"
#include "encoderstate.h"
#include "image.h"
#include "imagelist.h"
#include "uvg266.h"
#include "threadqueue.h"
#include "videoframe.h"
#include "rate_control.h"
#include "alf.h"
#include "reshape.h"


static int encoder_state_config_frame_init(encoder_state_t * const state) {
  state->frame->ref = uvg_image_list_alloc(MAX_REF_PIC_COUNT);
  if(!state->frame->ref) {
    fprintf(stderr, "Failed to allocate the picture list!\n");
    return 0;
  }
  state->frame->ref_list = REF_PIC_LIST_0;
  state->frame->num = 0;
  state->frame->poc = 0;
  state->frame->total_bits_coded = 0;
  state->frame->cur_frame_bits_coded = 0;
  state->frame->cur_gop_bits_coded = 0;
  state->frame->prepared = 0;
  state->frame->done = 1;

  state->frame->rc_alpha = 3.2003;
  state->frame->rc_beta = -1.367;
  state->frame->icost = 0;

  const encoder_control_t * const encoder = state->encoder_control;
  const int num_lcus = encoder->in.width_in_lcu * encoder->in.height_in_lcu;
  state->frame->lcu_stats = calloc(num_lcus, sizeof(lcu_stats_t));
  state->frame->aq_offsets = MALLOC(double, num_lcus);

  for (int y = 0; y < encoder->in.height_in_lcu; y++) {
    for (int x = 0; x < encoder->in.width_in_lcu; x++) {
      int temp = MIN(encoder->cfg.width - x * 64, 64) * MIN(encoder->cfg.height - y * 64, 64);
      state->frame->lcu_stats[x + y * encoder->in.width_in_lcu].pixels = temp;
    }
  }

  state->frame->c_para = malloc(sizeof(double) * num_lcus);
  if(state->frame->c_para == NULL) {
    return 0;
  }
  state->frame->k_para = malloc(sizeof(double) * num_lcus);
  if (state->frame->k_para == NULL) {
    return 0;
  }

  pthread_mutex_init(&state->frame->rc_lock, NULL);

  state->frame->new_ratecontrol = uvg_get_rc_data(NULL);

  return 1;
}

static void encoder_state_config_frame_finalize(encoder_state_t * const state) {
  if (state->frame == NULL) return;

  pthread_mutex_destroy(&state->frame->rc_lock);
  if (state->frame->c_para) FREE_POINTER(state->frame->c_para);
  if (state->frame->k_para) FREE_POINTER(state->frame->k_para);

  uvg_image_list_destroy(state->frame->ref);
  FREE_POINTER(state->frame->lcu_stats);
  FREE_POINTER(state->frame->aq_offsets);

}

static int encoder_state_config_tile_init(encoder_state_t * const state, 
                                          const int lcu_offset_x, const int lcu_offset_y,
                                          const int width, const int height, const int width_in_lcu, const int height_in_lcu) {
  
  const encoder_control_t * const encoder = state->encoder_control;
  state->tile->frame = uvg_videoframe_alloc(width, height, state->encoder_control->chroma_format, encoder->cfg.alf_type, encoder->cfg.cclm);
  
  state->tile->frame->hmvp_lut = malloc(sizeof(cu_info_t) * height_in_lcu * MAX_NUM_HMVP_CANDS);
  state->tile->frame->hmvp_size = calloc(1, sizeof(uint8_t) * height_in_lcu);

  state->tile->frame->rec = NULL;
  
  state->tile->frame->source = NULL;

  if (!state->tile->frame) {
    printf("Error allocating videoframe!\r\n");
    return 0;
  }

  state->tile->lcu_offset_x = lcu_offset_x;
  state->tile->lcu_offset_y = lcu_offset_y;
  state->tile->offset_x     = lcu_offset_x * LCU_WIDTH;
  state->tile->offset_y     = lcu_offset_y * LCU_WIDTH;

  state->tile->lcu_offset_in_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_offset_x + lcu_offset_y * encoder->in.width_in_lcu];
  
  // hor_buf_search and ver_buf_search store single row/col from each LCU row/col.
  // Because these lines are independent, the chroma subsampling only matters in one
  // of the directions, .
  unsigned luma_size = LCU_WIDTH * state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  unsigned chroma_sizes_hor[] = { 0, luma_size / 2, luma_size / 2, luma_size };
  unsigned chroma_sizes_ver[] = { 0, luma_size / 2, luma_size, luma_size };
  unsigned chroma_size_hor = chroma_sizes_hor[state->encoder_control->chroma_format];
  unsigned chroma_size_ver = chroma_sizes_ver[state->encoder_control->chroma_format];

  state->tile->hor_buf_search = uvg_yuv_t_alloc(luma_size, chroma_size_hor);
  state->tile->ver_buf_search = uvg_yuv_t_alloc(luma_size, chroma_size_ver);

  if (encoder->cfg.sao_type) {
    state->tile->hor_buf_before_sao = uvg_yuv_t_alloc(luma_size, chroma_size_hor);
    state->tile->ver_buf_before_sao = uvg_yuv_t_alloc(luma_size, chroma_size_ver);
  } else {
    state->tile->hor_buf_before_sao = NULL;
    state->tile->ver_buf_before_sao = NULL;
  }

  if (encoder->cfg.wpp) {
    int num_jobs = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
    state->tile->wf_jobs = MALLOC(threadqueue_job_t*, num_jobs);
    state->tile->wf_recon_jobs = MALLOC(threadqueue_job_t*, num_jobs);
    for (int i = 0; i < num_jobs; ++i) {
      state->tile->wf_jobs[i] = NULL;
      state->tile->wf_recon_jobs[i] = NULL;
    }
    if (!state->tile->wf_jobs) {
      printf("Error allocating wf_jobs array!\n");
      return 0;
    }
  } else {
    state->tile->wf_jobs = NULL;
    state->tile->wf_recon_jobs = NULL;
  }
  state->tile->id = encoder->tiles_tile_id[state->tile->lcu_offset_in_ts];
  return 1;
}

static void encoder_state_config_tile_finalize(encoder_state_t * const state) {
  if (state->tile == NULL) return;

  uvg_yuv_t_free(state->tile->hor_buf_search);
  uvg_yuv_t_free(state->tile->ver_buf_search);
  uvg_yuv_t_free(state->tile->hor_buf_before_sao);
  uvg_yuv_t_free(state->tile->ver_buf_before_sao);

  if (state->encoder_control->cfg.wpp) {
    int num_jobs = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
    for (int i = 0; i < num_jobs; ++i) {
      uvg_threadqueue_free_job(&state->tile->wf_jobs[i]);
      uvg_threadqueue_free_job(&state->tile->wf_recon_jobs[i]);
    }
  }

  FREE_POINTER(state->tile->frame->hmvp_lut);
  FREE_POINTER(state->tile->frame->hmvp_size);

  uvg_videoframe_free(state->tile->frame);
  state->tile->frame = NULL;
  FREE_POINTER(state->tile->wf_jobs);
  FREE_POINTER(state->tile->wf_recon_jobs);
}

static int encoder_state_config_slice_init(encoder_state_t * const state,
                                           const int start_address_in_ts,
                                           const int end_address_in_ts)
{
  state->slice->id = -1;
  for (int i = 0; i < state->encoder_control->slice_count; ++i) {
    if (state->encoder_control->slice_addresses_in_ts[i] == start_address_in_ts) {
      state->slice->id = i;
      break;
    }
  }
  assert(state->slice->id != -1);

  state->slice->start_in_ts = start_address_in_ts;
  state->slice->end_in_ts = end_address_in_ts;
  
  state->slice->start_in_rs = state->encoder_control->tiles_ctb_addr_ts_to_rs[start_address_in_ts];
  state->slice->end_in_rs = state->encoder_control->tiles_ctb_addr_ts_to_rs[end_address_in_ts];
  
  return 1;
}

static int encoder_state_config_wfrow_init(encoder_state_t * const state, 
                                          const int lcu_offset_y) {
  
  state->wfrow->lcu_offset_y = lcu_offset_y;
  return 1;
}

/**
 * \brief Initializer for main thread related things
   mostly arrays that are only needed one per frame
 * \param state encoder state
 * \returns int
 */
static int encoder_state_main_init(encoder_state_t* const state) {

  uint32_t lcus_in_frame = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  state->tile->frame->lmcs_aps = calloc(1, sizeof(lmcs_aps));
  state->tile->frame->lmcs_avg_processed = calloc(1, lcus_in_frame * sizeof(int8_t));
  state->tile->frame->lmcs_avg = calloc(1, lcus_in_frame * sizeof(int32_t));

  if (state->encoder_control->cfg.alf_type) {
    state->slice->alf = malloc(sizeof(*state->slice->alf));

    state->slice->alf->apss = malloc(sizeof(alf_aps) * ALF_CTB_MAX_NUM_APS);
    state->slice->alf->tile_group_luma_aps_id = malloc(ALF_CTB_MAX_NUM_APS * sizeof(int8_t));
    state->slice->alf->cc_filter_param = malloc(sizeof(*state->slice->alf->cc_filter_param));
    for (int aps_idx = 0; aps_idx < ALF_CTB_MAX_NUM_APS; aps_idx++) {
      state->slice->alf->tile_group_luma_aps_id[aps_idx] = -1;
    }
    state->slice->alf->tile_group_num_aps = -1;
    state->slice->alf->tile_group_chroma_aps_id = -1;
    state->slice->alf->tile_group_cc_alf_cb_enabled_flag = 0;
    state->slice->alf->tile_group_cc_alf_cr_enabled_flag = 0;
    state->slice->alf->tile_group_cc_alf_cb_aps_id = -1;
    state->slice->alf->tile_group_cc_alf_cr_aps_id = -1;
    state->slice->alf->num_of_param_sets = 0;
    memset(state->slice->alf->tile_group_alf_enabled_flag, 0, sizeof(state->slice->alf->tile_group_alf_enabled_flag));
    if (state->encoder_control->cfg.alf_type == UVG_ALF_FULL) {
      uvg_reset_cc_alf_aps_param(state->slice->alf->cc_filter_param);
    }

    state->tile->frame->alf_info = MALLOC(alf_info_t, 1);
    uvg_alf_create(state->tile->frame, state->encoder_control->chroma_format);
    uvg_set_aps_map(state->tile->frame, state->encoder_control->cfg.alf_type);
  }

  return 1;
}

static int encoder_state_main_finalize(encoder_state_t* const state) {

  FREE_POINTER(state->tile->frame->lmcs_aps);
  FREE_POINTER(state->tile->frame->lmcs_avg_processed);
  FREE_POINTER(state->tile->frame->lmcs_avg);

  if (state->encoder_control->cfg.alf_type) {
    if (state->slice->alf->apss != NULL) {
      FREE_POINTER(state->slice->alf->apss);
    }
    if (state->slice->alf->tile_group_luma_aps_id != NULL) {
      FREE_POINTER(state->slice->alf->tile_group_luma_aps_id);
    }
    if (state->slice->alf->cc_filter_param != NULL) {
      FREE_POINTER(state->slice->alf->cc_filter_param);
    }
    FREE_POINTER(state->slice->alf);

    uvg_alf_destroy(state->tile->frame);
    FREE_POINTER(state->tile->frame->alf_info);
    FREE_POINTER(state->tile->frame->alf_param_set_map);
  }
  return 1;
}

int uvg_encoder_state_init(encoder_state_t * const child_state, encoder_state_t * const parent_state) {
  //We require that, if parent_state is NULL:
  //child_state->encoder_control is set
  //
  //If parent_state is not NULL, the following variable should either be set to NULL,
  //in order to inherit from parent, or should point to a valid structure:
  //child_state->frame
  //child_state->tile
  //child_state->slice
  //child_state->wfrow
  
  child_state->parent = parent_state;
  child_state->children = MALLOC(encoder_state_t, 1);
  child_state->children[0].encoder_control = NULL;
  child_state->must_code_qp_delta = false;
  child_state->tqj_bitstream_written = NULL;
  child_state->tqj_recon_done = NULL;
  child_state->tqj_alf_process = NULL;
  
  if (!parent_state) {
    const encoder_control_t * const encoder = child_state->encoder_control;
    child_state->type = ENCODER_STATE_TYPE_MAIN;
    assert(child_state->encoder_control);
    child_state->frame = MALLOC(encoder_state_config_frame_t, 1);
    if (!child_state->frame || !encoder_state_config_frame_init(child_state)) {
      fprintf(stderr, "Could not initialize encoder_state->frame!\n");
      return 0;
    }
    child_state->tile = MALLOC(encoder_state_config_tile_t, 1);
    if (!child_state->tile || !encoder_state_config_tile_init(child_state, 0, 0, encoder->in.width, encoder->in.height, encoder->in.width_in_lcu, encoder->in.height_in_lcu)) {
      fprintf(stderr, "Could not initialize encoder_state->tile!\n");
      return 0;
    }

    child_state->slice = MALLOC(encoder_state_config_slice_t, 1);
    if (!child_state->slice || !encoder_state_config_slice_init(child_state, 0, encoder->in.width_in_lcu * encoder->in.height_in_lcu - 1)) {
      fprintf(stderr, "Could not initialize encoder_state->slice!\n");
      return 0;
    }
    child_state->wfrow = MALLOC(encoder_state_config_wfrow_t, 1);
    if (!child_state->wfrow || !encoder_state_config_wfrow_init(child_state, 0)) {
      fprintf(stderr, "Could not initialize encoder_state->wfrow!\n");
      return 0;
    }
  } else {
    child_state->encoder_control = parent_state->encoder_control;
    if (!child_state->frame) child_state->frame = parent_state->frame;
    if (!child_state->tile) child_state->tile = parent_state->tile;
    if (!child_state->slice) child_state->slice = parent_state->slice;
    if (!child_state->wfrow) child_state->wfrow = parent_state->wfrow;
  }

  // Intialization of the constraint structure
  child_state->constraint = uvg_init_constraint(child_state->constraint, child_state->encoder_control);

  uvg_bitstream_init(&child_state->stream);
  
  // Set CABAC output bitstream
  child_state->cabac.stream = &child_state->stream;
  
  //Create sub-encoders
  {
    const encoder_control_t * const encoder = child_state->encoder_control;
    uint32_t child_count = 0;
    //We first check the type of this element.
    //If it's a MAIN, it can allow both slices or tiles as child
    //If it's a TILE, it can allow slices as child, if its parent is not a slice, or wavefront rows if there is no other children
    //If it's a SLICE, it can allow tiles as child, if its parent is not a tile, or wavefront rows if there is no other children
    //If it's a WAVEFRONT_ROW, it doesn't allow any children
    int children_allow_wavefront_row = 0;
    int children_allow_slice = 0;
    int children_allow_tile = 0;
    int range_start;
    
    // First index of this encoder state in tile scan order.
    int start_in_ts;
    // Index of the first LCU after this state in tile scan order.
    int end_in_ts;
    
    switch(child_state->type) {
      case ENCODER_STATE_TYPE_MAIN:
        children_allow_slice = 1;
        children_allow_tile = 1;
        start_in_ts = 0;
        end_in_ts = child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu;

        encoder_state_main_init(child_state);

        break;
      case ENCODER_STATE_TYPE_SLICE:
        assert(child_state->parent);
        if (child_state->parent->type != ENCODER_STATE_TYPE_TILE) children_allow_tile = 1;
        start_in_ts = child_state->slice->start_in_ts;
        end_in_ts = child_state->slice->end_in_ts + 1;
        int num_wpp_rows = (end_in_ts - start_in_ts) / child_state->tile->frame->width_in_lcu;
        children_allow_wavefront_row = encoder->cfg.wpp && num_wpp_rows > 1;
        break;
      case ENCODER_STATE_TYPE_TILE:
        assert(child_state->parent);
        if (child_state->parent->type != ENCODER_STATE_TYPE_SLICE) children_allow_slice = 1;
        children_allow_wavefront_row =
          encoder->cfg.wpp && child_state->tile->frame->height_in_lcu > 1;
        start_in_ts = child_state->tile->lcu_offset_in_ts;
        end_in_ts = child_state->tile->lcu_offset_in_ts + child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu;
        break;
      case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
        //GCC tries to be too clever...
        start_in_ts = -1;
        end_in_ts = -1;
        break;
      default:
        fprintf(stderr, "Invalid encoder_state->type %d!\n", child_state->type);
        assert(0);
        return 0;
    }
    
    range_start = start_in_ts;
    //printf("%c-%p: start_in_ts=%d, end_in_ts=%d\n",child_state->type, child_state, start_in_ts, end_in_ts);
    while (range_start < end_in_ts && (children_allow_slice || children_allow_tile)) {
      encoder_state_t *new_child = NULL;
      int range_end_slice = range_start; //Will be incremented to get the range of the "thing"
      int range_end_tile = range_start; //Will be incremented to get the range of the "thing"
      
      int tile_allowed = uvg_lcu_at_tile_start(encoder, range_start) && children_allow_tile;
      int slice_allowed = uvg_lcu_at_slice_start(encoder, range_start) && children_allow_slice;
      
      //Find the smallest structure following the cursor
      if (slice_allowed) {
        while(!uvg_lcu_at_slice_end(encoder, range_end_slice)) {
          ++range_end_slice;
        }
      }
      
      if (tile_allowed) {
        while(!uvg_lcu_at_tile_end(encoder, range_end_tile)) {
          ++range_end_tile;
        }
      }
      
      //printf("range_start=%d, range_end_slice=%d, range_end_tile=%d, tile_allowed=%d, slice_allowed=%d end_in_ts=%d\n",range_start,range_end_slice,range_end_tile,tile_allowed,slice_allowed,end_in_ts);
      
      if ((!tile_allowed || (range_end_slice >= range_end_tile)) && !new_child && slice_allowed) {
        //Create a slice
        new_child = &child_state->children[child_count];
        new_child->encoder_control = encoder;
        new_child->type  = ENCODER_STATE_TYPE_SLICE;
        new_child->frame = child_state->frame;
        new_child->tile  = child_state->tile;
        new_child->wfrow = child_state->wfrow;
        new_child->slice = MALLOC(encoder_state_config_slice_t, 1);
        if (!new_child->slice || !encoder_state_config_slice_init(new_child, range_start, range_end_slice)) {
          fprintf(stderr, "Could not initialize encoder_state->slice!\n");
          return 0;
        }
      }
      
      if ((!slice_allowed || (range_end_slice < range_end_tile)) && !new_child && tile_allowed) {
        //Create a tile
        int tile_id = encoder->tiles_tile_id[range_start];
        int tile_x = tile_id % encoder->cfg.tiles_width_count;
        int tile_y = tile_id / encoder->cfg.tiles_width_count;
        
        int lcu_offset_x = encoder->tiles_col_bd[tile_x];
        int lcu_offset_y = encoder->tiles_row_bd[tile_y];
        int width_in_lcu = encoder->tiles_col_bd[tile_x+1]-encoder->tiles_col_bd[tile_x];
        int height_in_lcu = encoder->tiles_row_bd[tile_y+1]-encoder->tiles_row_bd[tile_y];
        int width = MIN(width_in_lcu * LCU_WIDTH, encoder->in.width - lcu_offset_x * LCU_WIDTH);
        int height = MIN(height_in_lcu * LCU_WIDTH, encoder->in.height - lcu_offset_y * LCU_WIDTH);
        
        new_child = &child_state->children[child_count];
        new_child->encoder_control = encoder;
        new_child->type  = ENCODER_STATE_TYPE_TILE;
        new_child->frame = child_state->frame;
        new_child->tile  = MALLOC(encoder_state_config_tile_t, 1);
        new_child->slice = child_state->slice;
        new_child->wfrow = child_state->wfrow;
        
        if (!new_child->tile || !encoder_state_config_tile_init(new_child, lcu_offset_x, lcu_offset_y, width, height, width_in_lcu, height_in_lcu)) {
          fprintf(stderr, "Could not initialize encoder_state->tile!\n");
          return 0;
        }
      }
      
      if (new_child) {
        child_state->children = realloc(child_state->children, sizeof(encoder_state_t) * (2+child_count));
        if (!child_state->children) {
          fprintf(stderr, "Failed to allocate memory for children...\n");
          return 0;
        }

        child_state->children[1 + child_count].encoder_control = NULL;

        //Fix children parent (since we changed the address), except for the last one which is not ready yet
        {
          uint32_t i, j;
          for (i = 0; child_state->children[i].encoder_control && i < child_count; ++i) {
            for (j = 0; child_state->children[i].children[j].encoder_control; ++j) {
              child_state->children[i].children[j].parent = &child_state->children[i];
            }
            for (j = 0; j < child_state->children[i].lcu_order_count; ++j) {
              child_state->children[i].lcu_order[j].encoder_state = &child_state->children[i];
            }
            child_state->children[i].cabac.stream = &child_state->children[i].stream;
          }
        }
        
        if (!uvg_encoder_state_init(&child_state->children[child_count], child_state)) {
          fprintf(stderr, "Unable to init child...\n");
          return 0;
        }
        child_count += 1;
      }
      
      range_start = MAX(range_end_slice, range_end_tile) + 1;
    }
    
    //We create wavefronts only if we have no children
    if (children_allow_wavefront_row && child_count == 0) {
      int first_row = encoder->tiles_ctb_addr_ts_to_rs[start_in_ts] / encoder->in.width_in_lcu;
      int last_row = encoder->tiles_ctb_addr_ts_to_rs[start_in_ts] / encoder->in.width_in_lcu;
      int num_rows;
      int i;
      
      assert(!(children_allow_slice || children_allow_tile));
      assert(child_count == 0);
      
      for (i=start_in_ts; i<end_in_ts; ++i) {
        const int row = encoder->tiles_ctb_addr_ts_to_rs[i] / encoder->in.width_in_lcu;
        if (row < first_row) first_row = row;
        if (row > last_row) last_row = row;
      }
      
      num_rows = last_row - first_row + 1;
      
      //When entropy_coding_sync_enabled_flag is equal to 1 and the first coding tree block in a slice is not the first coding
      //tree block of a row of coding tree blocks in a tile, it is a requirement of bitstream conformance that the last coding tree
      //block in the slice shall belong to the same row of coding tree blocks as the first coding tree block in the slice.
      
      if (encoder->tiles_ctb_addr_ts_to_rs[start_in_ts] % encoder->in.width_in_lcu != child_state->tile->lcu_offset_x) {
        if (num_rows > 1) {
          fprintf(stderr, "Invalid: first CTB in slice %d is not at the tile %d edge, and the slice spans on more than one row.\n", child_state->slice->id, child_state->tile->id);
          return 0;
        }
      }
      
      //FIXME Do the same kind of check if we implement slice segments
    
      child_count = num_rows;
      child_state->children = realloc(child_state->children, sizeof(encoder_state_t) * (num_rows + 1));
      child_state->children[num_rows].encoder_control = NULL;
      
      for (i=0; i < num_rows; ++i) {
        encoder_state_t *new_child = &child_state->children[i];
        
        new_child->encoder_control = encoder;
        new_child->type  = ENCODER_STATE_TYPE_WAVEFRONT_ROW;
        new_child->frame = child_state->frame;
        new_child->tile  = child_state->tile;
        new_child->slice = child_state->slice;
        new_child->wfrow = MALLOC(encoder_state_config_wfrow_t, 1);
        
        if (!new_child->wfrow || !encoder_state_config_wfrow_init(new_child, i)) {
          fprintf(stderr, "Could not initialize encoder_state->wfrow!\n");
          return 0;
        }
        
        if (!uvg_encoder_state_init(new_child, child_state)) {
          fprintf(stderr, "Unable to init child...\n");
          return 0;
        }
      }
    }
    
    child_state->is_leaf = (child_count == 0);
    //This node is a leaf, compute LCU-order
    if (child_state->is_leaf) {
      //All LCU computations are relative to the tile
      //Remark: this could be optimized, but since it's run only once, it's better to do it in a understandable way.
      
      //By default, the full tile
      int lcu_id;
      int lcu_start = 0;
      //End is the element AFTER the end (iterate < lcu_end)
      int lcu_end = child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu;
      
      //Restrict to the current slice if needed
      lcu_start = MAX(lcu_start, child_state->slice->start_in_ts - child_state->tile->lcu_offset_in_ts);
      lcu_end = MIN(lcu_end, child_state->slice->end_in_ts - child_state->tile->lcu_offset_in_ts + 1);
      
      //Restrict to the current wavefront row if needed
      if (child_state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
        lcu_start = MAX(lcu_start, (child_state->wfrow->lcu_offset_y) * child_state->tile->frame->width_in_lcu);
        lcu_end = MIN(lcu_end, (child_state->wfrow->lcu_offset_y + 1) * child_state->tile->frame->width_in_lcu);
      }
      
      child_state->lcu_order_count = lcu_end - lcu_start;
      child_state->lcu_order = MALLOC(lcu_order_element_t, child_state->lcu_order_count);
      assert(child_state->lcu_order);
      
      for (uint32_t i = 0; i < child_state->lcu_order_count; ++i) {
        lcu_id = lcu_start + i;
        child_state->lcu_order[i].encoder_state = child_state;
        child_state->lcu_order[i].id = lcu_id;
        child_state->lcu_order[i].index = i;
        child_state->lcu_order[i].position.x = lcu_id % child_state->tile->frame->width_in_lcu;
        child_state->lcu_order[i].position.y = lcu_id / child_state->tile->frame->width_in_lcu;
        child_state->lcu_order[i].position_px.x = child_state->lcu_order[i].position.x * LCU_WIDTH;
        child_state->lcu_order[i].position_px.y = child_state->lcu_order[i].position.y * LCU_WIDTH;
        child_state->lcu_order[i].size.x = MIN(LCU_WIDTH, encoder->in.width - (child_state->tile->lcu_offset_x * LCU_WIDTH + child_state->lcu_order[i].position_px.x));
        child_state->lcu_order[i].size.y = MIN(LCU_WIDTH, encoder->in.height - (child_state->tile->lcu_offset_y * LCU_WIDTH + child_state->lcu_order[i].position_px.y));
        child_state->lcu_order[i].first_row = uvg_lcu_in_first_row(child_state, child_state->tile->lcu_offset_in_ts + lcu_id);
        child_state->lcu_order[i].last_row = uvg_lcu_in_last_row(child_state, child_state->tile->lcu_offset_in_ts + lcu_id);
        child_state->lcu_order[i].first_column = uvg_lcu_in_first_column(child_state, child_state->tile->lcu_offset_in_ts + lcu_id);
        child_state->lcu_order[i].last_column = uvg_lcu_in_last_column(child_state, child_state->tile->lcu_offset_in_ts + lcu_id);
        
        child_state->lcu_order[i].above = NULL;
        child_state->lcu_order[i].below = NULL;
        child_state->lcu_order[i].left = NULL;
        child_state->lcu_order[i].right = NULL;
        
        if (!child_state->lcu_order[i].first_row) {
          //Find LCU above
          if (child_state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
            uint32_t j;
            //For all previous wavefront rows
            for (j=0; &child_state->parent->children[j] != child_state && child_state->parent->children[j].encoder_control; ++j) {
              if (child_state->parent->children[j].wfrow->lcu_offset_y == child_state->wfrow->lcu_offset_y - 1) {
                uint32_t k;
                for (k=0; k < child_state->parent->children[j].lcu_order_count; ++k) {
                  if (child_state->parent->children[j].lcu_order[k].position.x == child_state->lcu_order[i].position.x) {
                    assert(child_state->parent->children[j].lcu_order[k].position.y == child_state->lcu_order[i].position.y - 1);
                    child_state->lcu_order[i].above = &child_state->parent->children[j].lcu_order[k];
                  }
                }
              }
            }
          } else {
            child_state->lcu_order[i].above = &child_state->lcu_order[i-child_state->tile->frame->width_in_lcu];
          }
          assert(child_state->lcu_order[i].above);
          child_state->lcu_order[i].above->below = &child_state->lcu_order[i];
        }
        if (!child_state->lcu_order[i].first_column) {
          child_state->lcu_order[i].left = &child_state->lcu_order[i-1];
          assert(child_state->lcu_order[i].left->position.x == child_state->lcu_order[i].position.x - 1);
          child_state->lcu_order[i].left->right = &child_state->lcu_order[i];
        }
      }
    } else {
      child_state->lcu_order_count = 0;
      child_state->lcu_order = NULL;
    }
  }
  
  //Validate the structure
  if (child_state->type == ENCODER_STATE_TYPE_TILE) {
    if (child_state->tile->lcu_offset_in_ts < child_state->slice->start_in_ts) {
      fprintf(stderr, "Tile %d starts before slice %d, in which it should be included!\n", child_state->tile->id, child_state->slice->id);
      return 0;
    }
    if (child_state->tile->lcu_offset_in_ts + child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu - 1 > child_state->slice->end_in_ts) {
      fprintf(stderr, "Tile %d ends after slice %d, in which it should be included!\n", child_state->tile->id, child_state->slice->id);
      return 0;
    }
  }
  
  if (child_state->type == ENCODER_STATE_TYPE_SLICE) {
    if (child_state->slice->start_in_ts < child_state->tile->lcu_offset_in_ts) {
      fprintf(stderr, "Slice %d starts before tile %d, in which it should be included!\n", child_state->slice->id, child_state->tile->id);
      return 0;
    }
    if (child_state->slice->end_in_ts > child_state->tile->lcu_offset_in_ts + child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu - 1) {
      fprintf(stderr, "Slice %d ends after tile %d, in which it should be included!\n", child_state->slice->id, child_state->tile->id);
      return 0;
    }
  }
  
#ifdef UVG_DEBUG_PRINT_THREADING_INFO
  if (!parent_state) uvg_dbg_encoder_state_dump_graphviz(child_state);
#endif //UVG_DEBUG_PRINT_THREADING_INFO
  return 1;
}

void uvg_encoder_state_finalize(encoder_state_t * const state) {
  if (state->children) {
    int i=0;
    for (i = 0; state->children[i].encoder_control; ++i) {
      uvg_encoder_state_finalize(&state->children[i]);
    }

    FREE_POINTER(state->children);
  }

  if (state->type == ENCODER_STATE_TYPE_MAIN) {
    encoder_state_main_finalize(state);
  }
  
  FREE_POINTER(state->lcu_order);
  state->lcu_order_count = 0;
  
  if (!state->parent || (state->parent->wfrow != state->wfrow)) {
    FREE_POINTER(state->wfrow);
  }
  
  if (!state->parent || (state->parent->slice != state->slice)) {
    FREE_POINTER(state->slice);
  }
  
  if (!state->parent || (state->parent->tile != state->tile)) {
    encoder_state_config_tile_finalize(state);
    FREE_POINTER(state->tile);
  }
  
  if (!state->parent || (state->parent->frame != state->frame)) {
    encoder_state_config_frame_finalize(state);
    FREE_POINTER(state->frame);
  }
  
  if (state->constraint) {
    // End of the constraint structure
    uvg_constraint_free(state);
  }

  uvg_bitstream_finalize(&state->stream);

  uvg_threadqueue_free_job(&state->tqj_recon_done);
  uvg_threadqueue_free_job(&state->tqj_bitstream_written);
  if (state->encoder_control->cfg.alf_type && state->encoder_control->cfg.wpp) {
    encoder_state_t* parent = state;
    while (parent->parent) parent = parent->parent;
    uvg_threadqueue_free_job(&parent->tqj_alf_process);
  }
}
