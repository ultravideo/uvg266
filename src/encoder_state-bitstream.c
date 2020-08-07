/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "encoder_state-bitstream.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alf.h"
#include "bitstream.h"
#include "cabac.h"
#include "checkpoint.h"
#include "cu.h"
#include "encoder.h"
#include "encoder_state-geometry.h"
#include "encoderstate.h"
#include "imagelist.h"
#include "kvazaar.h"
#include "kvz_math.h"
#include "nal.h"
#include "scalinglist.h"
#include "tables.h"
#include "threadqueue.h"
#include "videoframe.h"

#define JVET_S0266_VUI_length 1
#define LUMA_ADAPTIVE_DEBLOCKING_FILTER_QP_OFFSET 1
#define JVET_S0076_ASPECT1 1
#define JVET_S0193_NO_OUTPUT_PRIOR_PIC 1
#define JVET_S0179_CONDITIONAL_SIGNAL_GCI 1
#define JVET_S0105_GCI_REORDER_IN_CATEGORY 1
#define JVET_S0138_GCI_PTL 1
#define JVET_S_SUB_PROFILE 1

static void encoder_state_write_bitstream_aud(encoder_state_t * const state)
{
  bitstream_t * const stream = &state->stream;
  kvz_nal_write(stream, KVZ_NAL_AUD_NUT, 0, 1);

  WRITE_U(stream, 1, 1, "aud_irap_or_gdr_au_flag");

  uint8_t pic_type = state->frame->slicetype == KVZ_SLICE_I ? 0
                   : state->frame->slicetype == KVZ_SLICE_P ? 1
                   :                                       2;
  WRITE_U(stream, pic_type, 3, "pic_type");

  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encoder_state_write_bitstream_PTL(bitstream_t *stream,
                                              encoder_state_t * const state)
{
  // PTL
  // Profile Tier
  // Main 10 profile == 1
  WRITE_U(stream, 1, 7, "general_profile_idc");
  WRITE_U(stream, state->encoder_control->cfg.high_tier, 1, "general_tier_flag");

#if !JVET_S0179_CONDITIONAL_SIGNAL_GCI
#if JVET_S0179_CONDITIONAL_SIGNAL_GCI
  WRITE_U(stream, 0, 1, "gci_present_flag");
  if (0) //if gci_present_flag
  {
#endif
#if JVET_S0105_GCI_REORDER_IN_CATEGORY
    WRITE_U(stream, 1, 1, "general_progressive_source_flag");
    WRITE_U(stream, state->encoder_control->in.source_scan_type != 0, 1, "general_interlaced_source_flag");

    // Constraint flags
    WRITE_U(stream, 0, 1, "general_non_packed_constraint_flag");
    WRITE_U(stream, 0, 1, "general_frame_only_constraint_flag");
    WRITE_U(stream, 0, 1, "general_non_projected_constraint_flag");

    WRITE_U(stream, 0, 1, "intra_only_constraint_flag");
    WRITE_U(stream, 0, 4, "max_bitdepth_constraint_idc");
    WRITE_U(stream, 0, 2, "max_chroma_format_constraint_idc");
    WRITE_U(stream, 0, 1, "no_res_change_in_clvs_constraint_flag");
    WRITE_U(stream, 0, 1, "one_tile_per_pic_constraint_flag");
    WRITE_U(stream, 0, 1, "one_slice_per_pic_constraint_flag");
    WRITE_U(stream, 0, 1, "one_subpic_per_pic_constraint_flag");
    WRITE_U(stream, 0, 1, "no_qtbtt_dual_tree_intra_constraint_flag");
    WRITE_U(stream, 0, 1, "no_partition_constraints_override_constraint_flag");
    WRITE_U(stream, 0, 1, "no_sao_constraint_flag");
    WRITE_U(stream, 0, 1, "no_alf_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ccalf_constraint_flag");
    WRITE_U(stream, 0, 1, "no_joint_cbcr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ref_wraparound_constraint_flag");
    WRITE_U(stream, 0, 1, "no_temporal_mvp_constraint_flag");
    WRITE_U(stream, 0, 1, "no_sbtmvp_constraint_flag");
    WRITE_U(stream, 0, 1, "no_amvr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_bdof_constraint_flag");
    WRITE_U(stream, 0, 1, "no_dmvr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_cclm_constraint_flag");
    WRITE_U(stream, 0, 1, "no_mts_constraint_flag");
    WRITE_U(stream, 0, 1, "no_sbt_constraint_flag");
    WRITE_U(stream, 0, 1, "no_affine_motion_constraint_flag");
    WRITE_U(stream, 0, 1, "no_bcw_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ibc_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ciip_constraint_flag");
    WRITE_U(stream, 0, 1, "no_fpel_mmvd_constraint_flag");
    WRITE_U(stream, 0, 1, "no_gpm_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ladf_constraint_flag");
    WRITE_U(stream, 0, 1, "no_transform_skip_constraint_flag");
    WRITE_U(stream, 0, 1, "no_bdpcm_constraint_flag");
    WRITE_U(stream, 0, 1, "no_qp_delta_constraint_flag");
    WRITE_U(stream, 0, 1, "no_dep_quant_constraint_flag");
    WRITE_U(stream, 0, 1, "no_sign_data_hiding_constraint_flag");
    WRITE_U(stream, 0, 1, "no_mixed_nalu_types_in_pic_constraint_flag");
    WRITE_U(stream, 0, 1, "no_trail_constraint_flag");
    WRITE_U(stream, 0, 1, "no_stsa_constraint_flag");
    WRITE_U(stream, 0, 1, "no_rasl_constraint_flag");
    WRITE_U(stream, 0, 1, "no_radl_constraint_flag");
    WRITE_U(stream, 0, 1, "no_idr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_cra_constraint_flag");
    WRITE_U(stream, 0, 1, "no_gdr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_aps_constraint_flag");
#endif
  }

  kvz_bitstream_align_zero(stream);
#endif

  // end Profile Tier

  uint8_t level = state->encoder_control->cfg.level;
  // ToDo: level hardcoded to 5.2
  WRITE_U(stream, 86, 8, "general_level_idc");


  WRITE_U(stream, 0, 1, "ptl_frame_only_constraint_flag");
  WRITE_U(stream, 0, 1, "ptl_multilayer_enabled_flag");

#if JVET_S0179_CONDITIONAL_SIGNAL_GCI
#if JVET_S0179_CONDITIONAL_SIGNAL_GCI
  WRITE_U(stream, 0, 1, "gci_present_flag");
  if (0) //if gci_present_flag
  {
#endif
#if JVET_S0105_GCI_REORDER_IN_CATEGORY
    WRITE_U(stream, 1, 1, "general_progressive_source_flag");
    WRITE_U(stream, state->encoder_control->in.source_scan_type != 0, 1, "general_interlaced_source_flag");

    // Constraint flags
    WRITE_U(stream, 0, 1, "general_non_packed_constraint_flag");
    WRITE_U(stream, 0, 1, "general_frame_only_constraint_flag");
    WRITE_U(stream, 0, 1, "general_non_projected_constraint_flag");

    WRITE_U(stream, 0, 1, "intra_only_constraint_flag");
    WRITE_U(stream, 0, 4, "max_bitdepth_constraint_idc");
    WRITE_U(stream, 0, 2, "max_chroma_format_constraint_idc");
    WRITE_U(stream, 0, 1, "no_res_change_in_clvs_constraint_flag");
    WRITE_U(stream, 0, 1, "one_tile_per_pic_constraint_flag");
    WRITE_U(stream, 0, 1, "one_slice_per_pic_constraint_flag");
    WRITE_U(stream, 0, 1, "one_subpic_per_pic_constraint_flag");
    WRITE_U(stream, 0, 1, "no_qtbtt_dual_tree_intra_constraint_flag");
    WRITE_U(stream, 0, 1, "no_partition_constraints_override_constraint_flag");
    WRITE_U(stream, 0, 1, "no_sao_constraint_flag");
    WRITE_U(stream, 0, 1, "no_alf_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ccalf_constraint_flag");
    WRITE_U(stream, 0, 1, "no_joint_cbcr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ref_wraparound_constraint_flag");
    WRITE_U(stream, 0, 1, "no_temporal_mvp_constraint_flag");
    WRITE_U(stream, 0, 1, "no_sbtmvp_constraint_flag");
    WRITE_U(stream, 0, 1, "no_amvr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_bdof_constraint_flag");
    WRITE_U(stream, 0, 1, "no_dmvr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_cclm_constraint_flag");
    WRITE_U(stream, 0, 1, "no_mts_constraint_flag");
    WRITE_U(stream, 0, 1, "no_sbt_constraint_flag");
    WRITE_U(stream, 0, 1, "no_affine_motion_constraint_flag");
    WRITE_U(stream, 0, 1, "no_bcw_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ibc_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ciip_constraint_flag");
    WRITE_U(stream, 0, 1, "no_fpel_mmvd_constraint_flag");
    WRITE_U(stream, 0, 1, "no_gpm_constraint_flag");
    WRITE_U(stream, 0, 1, "no_ladf_constraint_flag");
    WRITE_U(stream, 0, 1, "no_transform_skip_constraint_flag");
    WRITE_U(stream, 0, 1, "no_bdpcm_constraint_flag");
    WRITE_U(stream, 0, 1, "no_qp_delta_constraint_flag");
    WRITE_U(stream, 0, 1, "no_dep_quant_constraint_flag");
    WRITE_U(stream, 0, 1, "no_sign_data_hiding_constraint_flag");
    WRITE_U(stream, 0, 1, "no_mixed_nalu_types_in_pic_constraint_flag");
    WRITE_U(stream, 0, 1, "no_trail_constraint_flag");
    WRITE_U(stream, 0, 1, "no_stsa_constraint_flag");
    WRITE_U(stream, 0, 1, "no_rasl_constraint_flag");
    WRITE_U(stream, 0, 1, "no_radl_constraint_flag");
    WRITE_U(stream, 0, 1, "no_idr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_cra_constraint_flag");
    WRITE_U(stream, 0, 1, "no_gdr_constraint_flag");
    WRITE_U(stream, 0, 1, "no_aps_constraint_flag");
#endif
  }

  kvz_bitstream_align_zero(stream);
#endif

#if !JVET_S_SUB_PROFILE
  WRITE_U(stream, 1, 8, "num_sub_profiles");
  WRITE_U(stream, 0, 32, "general_sub_profile_idc");
#endif

  WRITE_U(stream, 0, 1, "sub_layer_level_present_flag");

  kvz_bitstream_align_zero(stream);


  WRITE_U(stream, 1, 8, "ptl_num_sub_profiles");
  WRITE_U(stream, 0, 32, "general_sub_profile_idc");


  // end PTL
}


static void encoder_state_write_bitstream_vid_parameter_set(bitstream_t* stream,
                                                            encoder_state_t * const state)
{
#ifdef KVZ_DEBUG
  printf("=========== Video Parameter Set ID: 0 ===========\n");
#endif

  WRITE_U(stream, 1, 4, "vps_video_parameter_set_id");
  WRITE_U(stream, 0, 6, "vps_max_layers_minus1" );
  WRITE_U(stream, 1, 3, "vps_max_sub_layers_minus1");

  //for each layer
  for (int i = 0; i < 1; i++) {
    WRITE_U(stream, 0, 6, "vps_layer_id");
  }

  kvz_bitstream_align_zero(stream);
  encoder_state_write_bitstream_PTL(stream, state);
  WRITE_U(stream, 0, 1, "vps_extension_flag")

  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

/*
static void encoder_state_write_bitstream_scaling_list(bitstream_t *stream,
                                                       encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
  uint32_t size_id;
  for (size_id = 0; size_id < SCALING_LIST_SIZE_NUM; size_id++) {
    int32_t list_id;
    for (list_id = 0; list_id < kvz_g_scaling_list_num[size_id]; list_id++) {
      uint8_t scaling_list_pred_mode_flag = 1;
      int32_t pred_list_idx;
      int32_t i;
      uint32_t ref_matrix_id = UINT32_MAX;

      for (pred_list_idx = list_id; pred_list_idx >= 0; pred_list_idx--) {
        const int32_t * const pred_list  = (list_id == pred_list_idx) ?
                                     kvz_scalinglist_get_default(size_id, pred_list_idx) :
                                     encoder->scaling_list.scaling_list_coeff[size_id][pred_list_idx];

        if (!memcmp(encoder->scaling_list.scaling_list_coeff[size_id][list_id], pred_list, sizeof(int32_t) * MIN(8, kvz_g_scaling_list_size[size_id])) &&
            ((size_id < SCALING_LIST_16x16) ||
             (encoder->scaling_list.scaling_list_dc[size_id][list_id] == encoder->scaling_list.scaling_list_dc[size_id][pred_list_idx]))) {
          ref_matrix_id = pred_list_idx;
          scaling_list_pred_mode_flag = 0;
          break;
        }
      }
      WRITE_U(stream, scaling_list_pred_mode_flag, 1, "scaling_list_pred_mode_flag" );

      if (!scaling_list_pred_mode_flag) {
        WRITE_UE(stream, list_id - ref_matrix_id, "scaling_list_pred_matrix_id_delta");
      } else {
        int32_t delta;
        const int32_t coef_num = MIN(MAX_MATRIX_COEF_NUM, kvz_g_scaling_list_size[size_id]);
        const uint32_t * const scan_cg = (size_id == 0) ? g_sig_last_scan_16x16 : g_sig_last_scan_32x32;
        int32_t next_coef = 8;
        const int32_t * const coef_list = encoder->scaling_list.scaling_list_coeff[size_id][list_id];

        if (size_id >= SCALING_LIST_16x16) {
          WRITE_SE(stream, encoder->scaling_list.scaling_list_dc[size_id][list_id] - 8, "scaling_list_dc_coef_minus8");
          next_coef = encoder->scaling_list.scaling_list_dc[size_id][list_id];
        }

        for (i = 0; i < coef_num; i++) {
          delta     = coef_list[scan_cg[i]] - next_coef;
          next_coef = coef_list[scan_cg[i]];
          if (delta > 127)
            delta -= 256;
          if (delta < -128)
            delta += 256;

          WRITE_SE(stream, delta, "scaling_list_delta_coef");
        }
      }
    }
  }
}
*/

static void encoder_state_write_bitstream_VUI(bitstream_t *stream,
                                              encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
#ifdef KVZ_DEBUG
  printf("=========== VUI Set ID: 0 ===========\n");
#endif

  WRITE_U(stream, 0, 1, "vui_progressive_source_flag");
  WRITE_U(stream, 0, 1, "vui_interlaced_source_flag");
#if JVET_S0266_VUI_length
  WRITE_U(stream, 0, 1, "vui_non_packed_constraint_flag");
  WRITE_U(stream, 0, 1, "vui_non_projected_constraint_flag");
#endif

  if (encoder->cfg.vui.sar_width > 0 && encoder->cfg.vui.sar_height > 0) {
    int i;
    static const struct
    {
      uint8_t width;
      uint8_t height;
      uint8_t idc;
    } sar[] = {
      // aspect_ratio_idc = 0 -> unspecified
      {  1,  1, 1 }, { 12, 11, 2 }, { 10, 11, 3 }, { 16, 11, 4 },
      { 40, 33, 5 }, { 24, 11, 6 }, { 20, 11, 7 }, { 32, 11, 8 },
      { 80, 33, 9 }, { 18, 11, 10}, { 15, 11, 11}, { 64, 33, 12},
      {160, 99, 13}, {  4,  3, 14}, {  3,  2, 15}, {  2,  1, 16},
      // aspect_ratio_idc = [17..254] -> reserved
      { 0, 0, 255 }
    };

    for (i = 0; sar[i].idc != 255; i++)
      if (sar[i].width  == encoder->cfg.vui.sar_width &&
          sar[i].height == encoder->cfg.vui.sar_height)
        break;

    WRITE_U(stream, 1, 1, "aspect_ratio_info_present_flag");
    WRITE_U(stream, sar[i].idc, 8, "aspect_ratio_idc");
    if (sar[i].idc == 255) {
      // EXTENDED_SAR
      WRITE_U(stream, encoder->cfg.vui.sar_width, 16, "sar_width");
      WRITE_U(stream, encoder->cfg.vui.sar_height, 16, "sar_height");
    }
  } else
    WRITE_U(stream, 0, 1, "aspect_ratio_info_present_flag");

  if (encoder->cfg.vui.overscan > 0) {
    WRITE_U(stream, 1, 1, "overscan_info_present_flag");
    WRITE_U(stream, encoder->cfg.vui.overscan - 1, 1, "overscan_appropriate_flag");
  }
  else
    WRITE_U(stream, 0, 1, "overscan_info_present_flag");

  WRITE_U(stream, 0, 1, "colour_description_present_flag");

  WRITE_U(stream, 0, 1, "chroma_loc_info_present_flag");

#if JVET_S0266_VUI_length
  if ((stream->cur_bit & 7) != 0) {
    WRITE_U(stream, 1, 0, "vui_payload_bit_equal_to_one");
    while ((stream->cur_bit & 7) != 0)
    {
      WRITE_U(stream, 0, 0, "vui_payload_bit_equal_to_zero");
    }
  }
#endif


  //WRITE_U(stream, 1, 1, "video_signal_type_present_flag");
  //WRITE_U(stream, encoder->cfg.vui.fullrange, 1, "video_full_range_flag");

  //IF overscan info
  //ENDIF

  /*
  if (encoder->cfg.vui.videoformat != 5 ||
      encoder->cfg.vui.fullrange   != 0 ||
      encoder->cfg.vui.colorprim   != 2 ||
      encoder->cfg.vui.transfer    != 2 ||
      encoder->cfg.vui.colormatrix != 2) {
    WRITE_U(stream, 1, 1, "video_signal_type_present_flag");
    WRITE_U(stream, encoder->cfg.vui.videoformat, 3, "chroma_format");
    WRITE_U(stream, encoder->cfg.vui.fullrange, 1, "video_full_range_flag");

    if (encoder->cfg.vui.colorprim   != 2 ||
        encoder->cfg.vui.transfer    != 2 ||
        encoder->cfg.vui.colormatrix != 2) {
      WRITE_U(stream, 1, 1, "colour_description_present_flag");
      WRITE_U(stream, encoder->cfg.vui.colorprim, 8, "colour_primaries");
      WRITE_U(stream, encoder->cfg.vui.transfer, 8, "transfer_characteristics");
      WRITE_U(stream, encoder->cfg.vui.colormatrix, 8, "matrix_coeffs");
    } else
      WRITE_U(stream, 0, 1, "colour_description_present_flag");
  } else
    WRITE_U(stream, 0, 1, "video_signal_type_present_flag");

  //IF video type
  //ENDIF

  if (encoder->cfg.vui.chroma_loc > 0) {
    WRITE_U(stream, 1, 1, "chroma_loc_info_present_flag");
    WRITE_UE(stream, encoder->cfg.vui.chroma_loc, "chroma_sample_loc_type_top_field");
    WRITE_UE(stream, encoder->cfg.vui.chroma_loc, "chroma_sample_loc_type_bottom_field");
  } else
    WRITE_U(stream, 0, 1, "chroma_loc_info_present_flag");

  //IF chroma loc info
  //ENDIF

  WRITE_U(stream, 0, 1, "neutral_chroma_indication_flag");
  WRITE_U(stream, encoder->vui.field_seq_flag, 1, "field_seq_flag"); // 0: frames, 1: fields
  WRITE_U(stream, encoder->vui.frame_field_info_present_flag, 1, "frame_field_info_present_flag");
  WRITE_U(stream, 0, 1, "default_display_window_flag");

  //IF default display window
  //ENDIF

  WRITE_U(stream, encoder->vui.timing_info_present_flag, 1, "vui_timing_info_present_flag");
  if (encoder->vui.timing_info_present_flag) {
    WRITE_U(stream, encoder->vui.num_units_in_tick, 32, "vui_num_units_in_tick");
    WRITE_U(stream, encoder->vui.time_scale, 32, "vui_time_scale");

    WRITE_U(stream, 0, 1, "vui_poc_proportional_to_timing_flag");
    WRITE_U(stream, 0, 1, "vui_hrd_parameters_present_flag");    
  }
  
  WRITE_U(stream, 0, 1, "bitstream_restriction_flag");

  //IF bitstream restriction
  //ENDIF

  */
}


static void encoder_state_write_bitstream_SPS_extension(bitstream_t *stream,
                                                        encoder_state_t * const state)
{
    WRITE_U(stream, 1, 1, "sps_extension_present_flag");

    WRITE_U(stream, 1, 1, "sps_range_extension_flag");
    WRITE_U(stream, 0, 1, "sps_multilayer_extension_flag");
    WRITE_U(stream, 0, 1, "sps_extension_6bits");
    WRITE_U(stream, 0, 1, "sps_extension_6bits");
    WRITE_U(stream, 0, 1, "sps_extension_6bits");
    WRITE_U(stream, 0, 1, "sps_extension_6bits");
    WRITE_U(stream, 0, 1, "sps_extension_6bits");
    WRITE_U(stream, 0, 1, "sps_extension_6bits");

    // Range Extension
    WRITE_U(stream, 0, 1, "transform_skip_rotation_enabled_flag");
    WRITE_U(stream, 0, 1, "transform_skip_context_enabled_flag");

    WRITE_U(stream, 0, 1, "extended_precision_processing_flag");
    WRITE_U(stream, state->encoder_control->cfg.intra_smoothing_disabled, 1, "intra_smoothing_disabled_flag");
    WRITE_U(stream, 0, 1, "high_precision_offsets_enabled_flag");
    WRITE_U(stream, 0, 1, "persistent_rice_adaptation_enabled_flag");
    WRITE_U(stream, 0, 1, "cabac_bypass_alignment_enabled_flag");
 
    
}

static void encoder_state_write_bitstream_seq_parameter_set(bitstream_t* stream,
                                                            encoder_state_t * const state)
{
  const encoder_control_t * encoder = state->encoder_control;

#ifdef KVZ_DEBUG
  printf("=========== Sequence Parameter Set ID: 0 ===========\n");
#endif
  /*


  */

  WRITE_U(stream, 0, 4, "sps_decoding_parameter_set_id");
  WRITE_U(stream, 0, 4, "sps_video_parameter_set_id");

  WRITE_U(stream, 1, 3, "sps_max_sub_layers_minus1");

  WRITE_U(stream, encoder->chroma_format, 2, "chroma_format_idc");
  WRITE_U(stream, kvz_math_floor_log2(LCU_WIDTH) - 5, 2, "sps_log2_ctu_size_minus5");


  WRITE_U(stream, 1, 1, "sps_ptl_dpb_hrd_params_present_flag");

  encoder_state_write_bitstream_PTL(stream, state);

  WRITE_U(stream, 0, 1, "gdr_enabled_flag");

  WRITE_U(stream, 0, 1, "ref_pic_resampling_enabled_flag");


  WRITE_UE(stream, encoder->in.width, "pic_width_max_in_luma_samples");
  WRITE_UE(stream, encoder->in.height, "pic_height_max_in_luma_samples");

  bool use_conformance_window = encoder->in.width != encoder->in.real_width || encoder->in.height != encoder->in.real_height;

  WRITE_U(stream, use_conformance_window, 1, "conformance_window_flag");
  if (use_conformance_window) {
    // The standard does not seem to allow setting conf_win values such that
    // the number of luma samples is not a multiple of 2. Options are to either
    // hide one line or show an extra line of non-video. Neither seems like a
    // very good option, so let's not even try.
    assert(!(encoder->in.width % 2));
    WRITE_UE(stream, 0, "conf_win_left_offset");
    WRITE_UE(stream, (encoder->in.width - encoder->in.real_width) >> 1,
      "conf_win_right_offset");
    WRITE_UE(stream, 0, "conf_win_top_offset");
    WRITE_UE(stream, (encoder->in.height - encoder->in.real_height) >> 1,
      "conf_win_bottom_offset");
  }

  WRITE_U(stream, 0, 1, "subpic_info_present_flag");

  WRITE_UE(stream, encoder->bitdepth-8, "bit_depth_minus8");

  WRITE_U(stream, encoder->cfg.wpp, 1, "sps_entropy_coding_sync_enabled_flag");
  WRITE_U(stream, encoder->tiles_enable || encoder->cfg.wpp, 1, "sps_entry_point_offsets_present_flag");

  WRITE_U(stream, 1, 4, "log2_max_pic_order_cnt_lsb_minus4");
  WRITE_U(stream, 0, 1, "sps_poc_msb_flag");
  WRITE_U(stream, 0, 2, "num_extra_ph_bits_bytes");
  WRITE_U(stream, 0, 2, "num_extra_sh_bits_bytes");

  WRITE_U(stream, 0, 1, "sps_sublayer_dpb_params_flag");

  //for each layer
  if (encoder->cfg.gop_lowdelay) {
    WRITE_UE(stream, encoder->cfg.ref_frames, "sps_max_dec_pic_buffering_minus1");
    WRITE_UE(stream, 0, "sps_max_num_reorder_pics");
  } else {
    WRITE_UE(stream, encoder->cfg.ref_frames + encoder->cfg.gop_len, "sps_max_dec_pic_buffering_minus1");
    WRITE_UE(stream, encoder->cfg.gop_len, "sps_max_num_reorder_pics");
  }
  WRITE_UE(stream, 0, "sps_max_latency_increase_plus1");
  //end for

  WRITE_UE(stream, MIN_SIZE-2, "log2_min_luma_coding_block_size_minus2"); // Min size 2^3 = 8x8
  // if(!no_partition_constraints_override_constraint_flag)
    WRITE_U(stream, 0, 1, "partition_constraints_override_enabled_flag");
  WRITE_UE(stream, 0, "sps_log2_diff_min_qt_min_cb_intra_slice_luma");
  WRITE_UE(stream, 0, "sps_max_mtt_hierarchy_depth_intra_slice_luma");  

  if (encoder->chroma_format != KVZ_CSP_400)
  {
    WRITE_U(stream, 0, 1, "qtbtt_dual_tree_intra_flag");
  }
  
  WRITE_UE(stream, 0, "sps_log2_diff_min_qt_min_cb_inter_slice");
  WRITE_UE(stream, 0, "sps_max_mtt_hierarchy_depth_inter_slice");  


#if 0 // mtt depth intra
  if (max_mtt_depth_intra != 0) {
    WRITE_UE(stream, 0, "sps_log2_diff_max_bt_min_qt_intra_tile_group_luma");
    WRITE_UE(stream, 0, "sps_log2_diff_max_tt_min_qt_intra_tile_group_luma");
  }
#endif
#if 0 // mtt depth inter
  if (max_mtt_depth_inter != 0) {
    WRITE_UE(stream, 0, "sps_log2_diff_max_bt_min_qt_inter_tile_group");
    WRITE_UE(stream, 0, "sps_log2_diff_max_tt_min_qt_inter_tile_group");
  }
#endif
#if 0 // Dual Tree
  if (encoder->cfg.dual_i_tree) {
    WRITE_UE(stream, 0, "sps_log2_diff_min_qt_min_cb_intra_tile_group_chroma");
    WRITE_UE(stream, 0, "sps_max_mtt_hierarchy_depth_intra_tile_group_chroma");

    if (max_mtt_depth_intra != 0) {
      WRITE_UE(stream, 0, "sps_log2_diff_max_bt_min_qt_intra_tile_group_chroma");
      WRITE_UE(stream, 0, "sps_log2_diff_max_tt_min_qt_intra_tile_group_chroma");
    }
  }
#endif

  if (LCU_WIDTH > 32)
    WRITE_U(stream, (TR_MAX_LOG2_SIZE - 5) ? 1 : 0, 1, "sps_max_luma_transform_size_64_flag");


  WRITE_U(stream, 0, 1, "sps_transform_skip_enabled_flag");
  WRITE_U(stream, 0, 1, "sps_mts_enabled_flag");
  WRITE_U(stream, 0, 1, "sps_lfnst_enabled_flag");


  WRITE_U(stream, 0, 1, "sps_joint_cbcr_enabled_flag");

  if (encoder->chroma_format != KVZ_CSP_400) {    
    WRITE_U(stream, 1, 1, "same_qp_table_for_chroma"); //TODO: Enable chroma QP scaling and fix kvz_get_scaled_qp()

    WRITE_SE(stream, 0, "qp_table_starts_minus26");    
    WRITE_UE(stream, 0, "num_points_in_qp_table_minus1");

      WRITE_UE(stream, 0, "delta_qp_in_val_minus1");
      WRITE_UE(stream, 1, "delta_qp_diff_val");

  }

  // if(!no_sao_constraint_flag)
    WRITE_U(stream, encoder->cfg.sao_type ? 1 : 0, 1, "sps_sao_enabled_flag");
  // if(!no_alf_constraint_flag)
    WRITE_U(stream, encoder->cfg.alf_enable, 1, "sps_alf_enable_flag");

    WRITE_U(stream, 0, 1, "sps_lmcs_enable_flag");

    WRITE_U(stream, 0, 1, "sps_weighted_pred_flag");           // Use of Weighting Prediction (P_SLICE)
    WRITE_U(stream, 0, 1, "sps_weighted_bipred_flag");        // Use of Weighting Bi-Prediction (B_SLICE)

    WRITE_U(stream, 0, 1, "long_term_ref_pics_flag");

    WRITE_U(stream, 0, 1, "sps_idr_rpl_present_flag");
    WRITE_U(stream, 0, 1, "rpl1_copy_from_rpl0_flag");

    WRITE_UE(stream, 0, "num_ref_pic_lists_in_sps[0]");
    WRITE_UE(stream, 0, "num_ref_pic_lists_in_sps[0]");

    WRITE_U(stream, 0, 1, "sps_ref_wraparound_enabled_flag");



  // if(!no_temporal_mvp_constraint_flag)
    WRITE_U(stream, state->encoder_control->cfg.tmvp_enable, 1,
        "sps_temporal_mvp_enabled_flag");
  if (state->encoder_control->cfg.tmvp_enable /* && !no_sbtmvp_constraint_flag */) {
    WRITE_U(stream, 0, 1, "sps_sbtmvp_enabled_flag");
  }

  // if(!no_amvr_constraint_flag)
    WRITE_U(stream, 0, 1, "sps_amvr_enabled_flag");
  // if(!no_bdof_constraint_flag)
    WRITE_U(stream, 0, 1, "sps_bdof_enabled_flag");
  // if(!no_dmvr_constraint_flag)
    WRITE_U(stream, 0, 1, "sps_smvd_enabled_flag");

  // if(!no_dmvr_constraint_flag)
    WRITE_U(stream, 0, 1, "sps_dmvr_enabled_flag");
  WRITE_U(stream, 0, 1, "sps_mmvd_enabled_flag");  


  WRITE_UE(stream, MRG_MAX_NUM_CANDS - 6, "six_minus_max_num_merge_cand");
  WRITE_U(stream, 0, 1, "sps_sbt_enabled_flag");
  WRITE_U(stream, 0, 1, "sps_affine_enabled_flag");

  WRITE_U(stream, 0, 1, "sps_bcw_enabled_flag");

  WRITE_U(stream, 0, 1, "sps_ciip_enabled_flag");
  if (6 /*MAX_NUM_MERGE_CAND*/ >= 2)
  {
    WRITE_U(stream, 0, 1, "sps_gpm_enabled_flag");
  }

  WRITE_UE(stream, 0, "log2_parallel_merge_level_minus2");

  WRITE_U(stream, 0, 1, "sps_isp_enabled_flag");
  WRITE_U(stream, 0, 1, "sps_mrl_enabled_flag");
  WRITE_U(stream, 0, 1, "sps_mip_enabled_flag");
  // if(!no_cclm_constraint_flag)
    WRITE_U(stream, 0, 1, "sps_cclm_enabled_flag");

  WRITE_U(stream, 0, 1, "sps_chroma_horizontal_collocated_flag");
  WRITE_U(stream, 0, 1, "sps_chroma_vertical_collocated_flag");

  WRITE_U(stream, 0, 1, "sps_palette_enabled_flag");
  WRITE_U(stream, 0, 1, "sps_ibc_enabled_flag");

#if LUMA_ADAPTIVE_DEBLOCKING_FILTER_QP_OFFSET
  // if(!no_ladf_constraint_flag)
  WRITE_U(stream, 0, 1, "sps_ladf_enabled_flag");
#endif

  WRITE_U(stream, 0, 1, "scaling_list_enabled_flag");

  WRITE_U(stream, 0, 1, "pic_dep_quant_enabled_flag");

  WRITE_U(stream, encoder->cfg.signhide_enable, 1, "pic_sign_data_hiding_enabled_flag");

  WRITE_U(stream, 0, 1, "sps_virtual_boundaries_enabled_flag");

  WRITE_U(stream, 0, 1, "general_hrd_parameters_present_flag");
  /*
  WRITE_U(stream, encoder->vui.timing_info_present_flag, 1, "general_hrd_parameters_present_flag");
  if (encoder->vui.timing_info_present_flag) {
    WRITE_U(stream, encoder->vui.num_units_in_tick, 32, "num_units_in_tick");
    WRITE_U(stream, encoder->vui.time_scale, 32, "time_scale");

    WRITE_U(stream, 0, 1, "general_nal_hrd_parameters_present_flag");
    WRITE_U(stream, 0, 1, "general_vcl_hrd_parameters_present_flag");
    WRITE_U(stream, 0, 1, "general_same_pic_timing_in_all_ols_flag");
    WRITE_U(stream, 0, 1, "general_decoding_unit_hrd_params_present_flag");
    
    WRITE_U(stream, 0, 4, "bit_rate_scale");
    WRITE_U(stream, 0, 4, "cpb_size_scale");

    WRITE_UE(stream, 0, "hrd_cpb_cnt_minus1");
  }
  */
  WRITE_U(stream, 0, 1, "field_seq_flag");
  WRITE_U(stream, 0, 1, "vui_parameters_present_flag");
  
  // ToDo: Check and enable
  //encoder_state_write_bitstream_VUI(stream, state);

  encoder_state_write_bitstream_SPS_extension(stream, state);

  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encoder_state_write_bitstream_pic_parameter_set(bitstream_t* stream,
                                                            encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
#ifdef KVZ_DEBUG
  printf("=========== Picture Parameter Set ID: 0 ===========\n");
#endif
  WRITE_U(stream, 0, 6, "pps_pic_parameter_set_id");

  WRITE_U(stream, 0, 4, "pps_seq_parameter_set_id");

  WRITE_U(stream, 0, 1, "mixed_nalu_types_in_pic_flag");

  WRITE_UE(stream, encoder->in.width, "pic_width_in_luma_samples");
  WRITE_UE(stream, encoder->in.height, "pic_height_in_luma_samples");

  bool use_conformance_window = encoder->in.width != encoder->in.real_width || encoder->in.height != encoder->in.real_height;

  WRITE_U(stream, use_conformance_window, 1, "conformance_window_flag");
  if (use_conformance_window) {
    // The standard does not seem to allow setting conf_win values such that
    // the number of luma samples is not a multiple of 2. Options are to either
    // hide one line or show an extra line of non-video. Neither seems like a
    // very good option, so let's not even try.
    assert(!(encoder->in.width % 2));
    WRITE_UE(stream, 0, "conf_win_left_offset");
    WRITE_UE(stream, (encoder->in.width - encoder->in.real_width) >> 1,
      "conf_win_right_offset");
    WRITE_UE(stream, 0, "conf_win_top_offset");
    WRITE_UE(stream, (encoder->in.height - encoder->in.real_height) >> 1,
      "conf_win_bottom_offset");
  }
  WRITE_U(stream, 0, 1, "scaling_window_flag");
  WRITE_U(stream, 0, 1, "output_flag_present_flag");
  WRITE_U(stream, 1, 1, "pps_no_pic_partition_flag");
  WRITE_U(stream, 0, 1, "subpic_id_mapping_in_pps_flag");

  /*
  WRITE_U(stream, encoder->tiles_enable ? 0 : 1, 1, "single_tile_in_pic_flag");

  if (encoder->tiles_enable) {

    WRITE_U(stream, encoder->tiles_uniform_spacing_flag, 1, "uniform_spacing_flag");

    if (!encoder->tiles_uniform_spacing_flag) {
      int i;
      for (i = 0; i < encoder->cfg.tiles_width_count - 1; ++i) {
        WRITE_UE(stream, encoder->tiles_col_width[i] - 1, "column_width_minus1[...]");
      }
      for (i = 0; i < encoder->cfg.tiles_height_count - 1; ++i) {
        WRITE_UE(stream, encoder->tiles_row_height[i] - 1, "row_height_minus1[...]");
      }
    }
    else {
      WRITE_UE(stream, encoder->cfg.tiles_width_count - 1, "num_tile_columns_minus1");
      WRITE_UE(stream, encoder->cfg.tiles_height_count - 1, "num_tile_rows_minus1");
      // ToDo: Signal the tiles properly
    }
    WRITE_U(stream, 0, 1, "brick_splitting_present_flag");
    WRITE_U(stream, 1, 1, "single_brick_per_slice_flag");

    WRITE_U(stream, 0, 1, "loop_filter_across_bricks_enabled_flag");
    // if loop_filter_across_bricks_enabled_flag
    //WRITE_U(stream, 0, 1, "loop_filter_across_tiles_enabled_flag");
  }
  */
  
  WRITE_U(stream, 0, 1, "cabac_init_present_flag");

  WRITE_UE(stream, 0, "num_ref_idx_l0_default_active_minus1");
  WRITE_UE(stream, 0, "num_ref_idx_l1_default_active_minus1");
  WRITE_U(stream, 0, 1, "rpl1_idx_present_flag");

  WRITE_U(stream, 0, 1, "weighted_pred_flag");   // Use of Weighting Prediction (P_SLICE)
  WRITE_U(stream, 0, 1, "weighted_bipred_flag");  // Use of Weighting Bi-Prediction (B_SLICE)
  WRITE_U(stream, 0, 1, "pps_ref_wraparound_enabled_flag");

  WRITE_SE(stream, ((int8_t)encoder->cfg.qp) - 26, "init_qp_minus26");
  WRITE_U(stream, encoder->max_qp_delta_depth >= 0 ? 1:0, 1, "cu_qp_delta_enabled_flag");
  if (encoder->max_qp_delta_depth >= 0) {
    // Use separate QP for each LCU when rate control is enabled.    
    WRITE_UE(stream, encoder->max_qp_delta_depth, "diff_cu_qp_delta_depth");
  }

  WRITE_U(stream, 0,1, "pps_chroma_tool_offsets_present_flag");
  /* // If chroma_tool_offsets_present
  //TODO: add QP offsets
  WRITE_SE(stream, 0, "pps_cb_qp_offset");
  WRITE_SE(stream, 0, "pps_cr_qp_offset");
  WRITE_U(stream, 0, 1, "pps_joint_cbcr_qp_offset_present_flag");
    // If pps_joint_cbcr_qp_offset_present_flag
    //WRITE_SE(stream, 0, "pps_joint_cbcr_qp_offset");
  WRITE_U(stream, 0, 1, "pps_slice_chroma_qp_offsets_present_flag");

  WRITE_U(stream, 0, 1, "cu_chroma_qp_offset_enabled_flag");
  */

  WRITE_U(stream, 1, 1, "deblocking_filter_control_present_flag");

  //IF deblocking_filter
  WRITE_U(stream, 0, 1, "deblocking_filter_override_enabled_flag");
  WRITE_U(stream, encoder->cfg.deblock_enable ? 0 : 1, 1,
          "pps_disable_deblocking_filter_flag");

  //IF !disabled
  if (encoder->cfg.deblock_enable) {
      WRITE_SE(stream, encoder->cfg.deblock_beta, "pps_beta_offset_div2");
      WRITE_SE(stream, encoder->cfg.deblock_tc, "pps_tc_offset_div2");
  }

  WRITE_U(stream, 0, 1, "picture_header_extension_present_flag");
  WRITE_U(stream, 0, 1, "slice_header_extension_present_flag");
  WRITE_U(stream, 0, 1, "pps_extension_present_flag");

  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encoder_state_write_bitstream_prefix_sei_version(encoder_state_t * const state)
{
#define STR_BUF_LEN 1000
  bitstream_t * const stream = &state->stream;
  int i, length;
  char buf[STR_BUF_LEN] = { 0 };
  char *s = buf + 16;
  const kvz_config * const cfg = &state->encoder_control->cfg;

  // random uuid_iso_iec_11578 generated with www.famkruithof.net/uuid/uuidgen
  static const uint8_t uuid[16] = {
    0x32, 0xfe, 0x46, 0x6c, 0x98, 0x41, 0x42, 0x69,
    0xae, 0x35, 0x6a, 0x91, 0x54, 0x9e, 0xf3, 0xf1
  };
  memcpy(buf, uuid, 16);

  // user_data_payload_byte
  s += sprintf(s, "uvg266 VVC Encoder v. " VERSION_STRING " - "
                  "Copyleft 2020- - http://ultravideo.fi/ - options:");
  s += sprintf(s, " %dx%d", cfg->width, cfg->height);
  s += sprintf(s, " deblock=%d:%d:%d", cfg->deblock_enable,
               cfg->deblock_beta, cfg->deblock_tc);
  s += sprintf(s, " sao=%d", cfg->sao_type);
  s += sprintf(s, " intra_period=%d", cfg->intra_period);
  s += sprintf(s, " qp=%d", cfg->qp);
  s += sprintf(s, " ref=%d", cfg->ref_frames);

  length = (int)(s - buf + 1);  // length, +1 for \0

  // Assert this so that in the future if the message gets longer, we remember
  // to increase the buf len. Divide by 2 for margin.
  assert(length < STR_BUF_LEN / 2);

  // payloadType = 5 -> user_data_unregistered
  WRITE_U(stream, 5, 8, "last_payload_type_byte");

  // payloadSize
  for (i = 0; i <= length - 255; i += 255)
    WRITE_U(stream, 255, 8, "ff_byte");
  WRITE_U(stream, length - i, 8, "last_payload_size_byte");

  for (i = 0; i < length; i++)
    WRITE_U(stream, ((uint8_t *)buf)[i], 8, "sei_payload");

  // The bitstream is already aligned, but align it anyway.
  kvz_bitstream_align(stream);

#undef STR_BUF_LEN
}

/*
static void encoder_state_write_active_parameter_sets_sei_message(encoder_state_t * const state) {

  const encoder_control_t * const encoder = state->encoder_control;
  bitstream_t * const stream = &state->stream;

  int i = 0;

  int active_vps_id = 0;
  int self_contained_cvs_flag = 0;
  int no_parameter_set_update_flag = 0;
  int num_sps_ids_minus1 = 0;
  int layer_sps_idx = 0;
  int active_seq_parameter_set_id = 0;
  int vps_base_layer_internal_flag = 0;

  int max_layers_minus1 = 0;

  WRITE_U(stream, 129, 8, "last_payload_type_byte"); //active_parameter_sets
  WRITE_U(stream, 2, 8, "last_payload_size_byte");
  WRITE_U(stream, active_vps_id, 4, "active_video_parameter_set_id");
  WRITE_U(stream, self_contained_cvs_flag, 1, "self_contained_cvs_flag");
  WRITE_U(stream, no_parameter_set_update_flag, 1, "no_parameter_set_update_flag");
  WRITE_UE(stream, num_sps_ids_minus1, "num_sps_ids_minus1");
  //for (i = 0; i <= num_sps_ids_minus1; ++i) {
  WRITE_UE(stream, active_seq_parameter_set_id, "active_seq_parameter_set_id");
  //}
  // for (i = vps_base_layer_internal_flag; i <= max_layers_minus1; ++i){
  WRITE_UE(stream, layer_sps_idx, "layer_sps_idx");
  //}

  kvz_bitstream_rbsp_trailing_bits(stream); //rbsp_trailing_bits
}
*/

static void encoder_state_write_picture_timing_sei_message(encoder_state_t * const state) {

  bitstream_t * const stream = &state->stream;

  if (state->encoder_control->vui.frame_field_info_present_flag){

    int8_t odd_picture = state->frame->num % 2;
    int8_t pic_struct = 0; //0: progressive picture, 1: top field, 2: bottom field, 3...
    int8_t source_scan_type = 1; //0: interlaced, 1: progressive

    switch (state->tile->frame->source->interlacing){
    case 0: //Progressive frame
      pic_struct = 0;
      source_scan_type = 1;
      break;
    case 1: //Top field first
      pic_struct = odd_picture ? 2 : 1;
      source_scan_type = 0;
      break;
    case 2: //Bottom field first
      pic_struct = odd_picture ? 1 : 2;
      source_scan_type = 0;
      break;
    default:
      assert(0); //Should never execute
      break;
    }

    WRITE_U(stream, 1, 8, "last_payload_type_byte"); //pic_timing
    WRITE_U(stream, 1, 8, "last_payload_size_byte");
    WRITE_U(stream, pic_struct, 4, "pic_struct");
    WRITE_U(stream, source_scan_type, 2, "source_scan_type");
    WRITE_U(stream, 0, 1, "duplicate_flag");

    kvz_bitstream_align(stream);
  }
}


// ToDo: Enable tiles/wpp
static void encoder_state_entry_points_explore(const encoder_state_t * const state, int * const r_count, int * const r_max_length) {
  int i;
  for (i = 0; state->children[i].encoder_control; ++i) {
    if (state->children[i].is_leaf) {
      const int my_length = kvz_bitstream_tell(&state->children[i].stream)/8;
      ++(*r_count);
      if (my_length > *r_max_length) {
        *r_max_length = my_length;
      }
    } else {
      encoder_state_entry_points_explore(&state->children[i], r_count, r_max_length);
    }
  }
}

static void encoder_state_write_bitstream_entry_points_write(bitstream_t * const stream, const encoder_state_t * const state, const int num_entry_points, const int write_length, int * const r_count) {
  int i;
  for (i = 0; state->children[i].encoder_control; ++i) {
    if (state->children[i].is_leaf) {
      const int my_length = kvz_bitstream_tell(&state->children[i].stream)/8;
      ++(*r_count);
      //Don't write the last one
      if (*r_count < num_entry_points) {
        WRITE_U(stream, my_length - 1, write_length, "entry_point_offset-minus1")
      }
    } else {
      encoder_state_write_bitstream_entry_points_write(stream, &state->children[i], num_entry_points, write_length, r_count);
    }
  }
}

static void kvz_encoder_state_write_bitstream_picture_header(
    struct bitstream_t * const stream,
    struct encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;

#ifdef KVZ_DEBUG
  printf("=========== Picture Header ===========\n");
#endif

  if (state->frame->pictype == KVZ_NAL_IDR_W_RADL
    || state->frame->pictype == KVZ_NAL_IDR_N_LP) {
    WRITE_U(stream, 1, 1, "ph_gdr_or_irap_pic_flag");
#if JVET_S0076_ASPECT1
    WRITE_U(stream, 0, 1, "ph_non_ref_pic_flag");
#endif
    WRITE_U(stream, 0, 1, "ph_gdr_pic_flag");
    WRITE_U(stream, 0, 1, "ph_inter_slice_allowed_flag");
  }
  else {
    WRITE_U(stream, 0, 1, "ph_gdr_or_irap_pic_flag");
#if JVET_S0076_ASPECT1
    WRITE_U(stream, 0, 1, "ph_non_ref_pic_flag");
#endif
    WRITE_U(stream, 1, 1, "ph_inter_slice_allowed_flag");
    WRITE_U(stream, 1, 1, "ph_intra_slice_allowed_flag");
  }

#if !JVET_S0076_ASPECT1
  WRITE_U(stream, 0, 1, "non_reference_picture_flag");
#endif
  WRITE_UE(stream, 0, "ph_pic_parameter_set_id");

  WRITE_U(stream, state->frame->poc & 0x1f, 5, "ph_pic_order_cnt_lsb");

  if (state->frame->pictype == KVZ_NAL_IDR_W_RADL
    || state->frame->pictype == KVZ_NAL_IDR_N_LP) {
  }
  else {
    // ToDo: ALF flag
    //WRITE_U(stream, state->encoder_control->cfg.tmvp_enable, 1, "ph_pic_temporal_mvp_enabled_flag");
    WRITE_U(stream, 0, 1, "ph_mvd_l1_zero_flag");
  }

  

  // alf enable flags and aps IDs
  if (encoder->cfg.alf_enable)
  {
    if (encoder->cfg.alf_info_in_ph_flag)
    {
     /* WRITE_FLAG(picHeader->getAlfEnabledFlag(COMPONENT_Y), "ph_alf_enabled_flag");
      if (picHeader->getAlfEnabledFlag(COMPONENT_Y))
      {
        WRITE_CODE(picHeader->getNumAlfAps(), 3, "ph_num_alf_aps_ids_luma");
        const std::vector<int>&   apsId = picHeader->getAlfAPSs();
        for (int i = 0; i < picHeader->getNumAlfAps(); i++)
        {
          WRITE_CODE(apsId[i], 3, "ph_alf_aps_id_luma");
        }

        const int alfChromaIdc = picHeader->getAlfEnabledFlag(COMPONENT_Cb) + picHeader->getAlfEnabledFlag(COMPONENT_Cr) * 2;
        if (sps->getChromaFormatIdc() != CHROMA_400)
        {
          WRITE_CODE(alfChromaIdc, 2, "ph_alf_chroma_idc");
        }
        if (alfChromaIdc)
        {
          WRITE_CODE(picHeader->getAlfApsIdChroma(), 3, "ph_alf_aps_id_chroma");
        }
        if (sps->getCCALFEnabledFlag())
        {
          WRITE_FLAG(picHeader->getCcAlfEnabledFlag(COMPONENT_Cb), "ph_cc_alf_cb_enabled_flag");
          if (picHeader->getCcAlfEnabledFlag(COMPONENT_Cb))
          {
            WRITE_CODE(picHeader->getCcAlfCbApsId(), 3, "ph_cc_alf_cb_aps_id");
          }
          WRITE_FLAG(picHeader->getCcAlfEnabledFlag(COMPONENT_Cr), "ph_cc_alf_cr_enabled_flag");
          if (picHeader->getCcAlfEnabledFlag(COMPONENT_Cr))
          {
            WRITE_CODE(picHeader->getCcAlfCrApsId(), 3, "ph_cc_alf_cr_aps_id");
          }
        }
      }*/
    }
    else
    {
      state->tile->frame->alf_info->g_ctu_enable_flag[COMPONENT_Y] = true;
      state->tile->frame->alf_info->g_ctu_enable_flag[COMPONENT_Cb] = true;
      state->tile->frame->alf_info->g_ctu_enable_flag[COMPONENT_Cr] = true;
      state->tile->frame->alf_info->g_alf_cc_enable_flag[COMPONENT_Cb] = encoder->cfg.alf_cc_enabled_flag;
      state->tile->frame->alf_info->g_alf_cc_enable_flag[COMPONENT_Cr] = encoder->cfg.alf_cc_enabled_flag;
    }
  }
  else
  {
    state->tile->frame->alf_info->g_ctu_enable_flag[COMPONENT_Y] = false;
    state->tile->frame->alf_info->g_ctu_enable_flag[COMPONENT_Cb] = false;
    state->tile->frame->alf_info->g_ctu_enable_flag[COMPONENT_Cr] = false;
    state->tile->frame->alf_info->g_alf_cc_enable_flag[COMPONENT_Cb] = false;
    state->tile->frame->alf_info->g_alf_cc_enable_flag[COMPONENT_Cr] = false;
  }

  WRITE_U(stream, state->encoder_control->cfg.tmvp_enable, 1, "pic_temporal_mvp_enabled_flag");
  WRITE_U(stream, 0, 1, "pic_mvd_l1_zero_flag");

  if (encoder->cfg.sao_type) {
    WRITE_U(stream, 1, 1, "slice_sao_luma_flag");
    if (encoder->chroma_format != KVZ_CSP_400) {
      WRITE_U(stream, 1, 1, "slice_sao_chroma_flag");
    }
  }

  // getDeblockingFilterControlPresentFlag

  // END PICTURE HEADER

}

static void kvz_encoder_state_write_bitstream_ref_pic_list(
  struct bitstream_t* const stream,
  struct encoder_state_t* const state)
{
  int j;
  int ref_negative = 0;
  int ref_positive = 0;
  const encoder_control_t* const encoder = state->encoder_control;
  if (encoder->cfg.gop_len) {
    for (j = 0; j < state->frame->ref->used_size; j++) {
      if (state->frame->ref->pocs[j] < state->frame->poc) {
        ref_negative++;
      }
      else {
        ref_positive++;
      }
    }
  }
  else ref_negative = state->frame->ref->used_size;

  int last_poc = 0;
  int poc_shift = 0;


  WRITE_UE(stream, ref_negative, "num_ref_entries[0]");  
  for (j = 0; j < ref_negative; j++) {
    int8_t delta_poc = 0;

    if (encoder->cfg.gop_len) {
      int8_t found = 0;
      do {
        delta_poc = encoder->cfg.gop[state->frame->gop_offset].ref_neg[j + poc_shift];
        for (int i = 0; i < state->frame->ref->used_size; i++) {
          if (state->frame->ref->pocs[i] == state->frame->poc - delta_poc) {
            found = 1;
            break;
          }
        }
        if (!found) poc_shift++;
        if (j + poc_shift == ref_negative) {
          fprintf(stderr, "Failure, reference not found!");
          exit(EXIT_FAILURE);
        }
      } while (!found);
    }
    /*
    WRITE_U(stream, j, 1, "inter_layer_ref_pic_flag");
    if (j) {
      WRITE_UE(stream, j, "ilrp_idx");
    }
    */
    WRITE_UE(stream, delta_poc, "abs_delta_poc_st");
    if (delta_poc+1) WRITE_U(stream, 1, 1, "strp_entry_sign_flag");
    last_poc = delta_poc;
    
  }
  last_poc = 0;
  poc_shift = 0;
  WRITE_UE(stream, ref_positive, "num_ref_entries[1]");
  for (j = 0; j < ref_positive; j++) {
    int8_t delta_poc = 0;

    if (encoder->cfg.gop_len) {
      int8_t found = 0;
      do {
        delta_poc = encoder->cfg.gop[state->frame->gop_offset].ref_pos[j + poc_shift];
        for (int i = 0; i < state->frame->ref->used_size; i++) {
          if (state->frame->ref->pocs[i] == state->frame->poc + delta_poc) {
            found = 1;
            break;
          }
        }
        if (!found) poc_shift++;
        if (j + poc_shift == ref_positive) {
          fprintf(stderr, "Failure, reference not found!");
          exit(EXIT_FAILURE);
        }
      } while (!found);
    }
    /*
    WRITE_U(stream, j, 1, "inter_layer_ref_pic_flag");
    if (j) {
      WRITE_UE(stream, j, "ilrp_idx");
    }
    */

    WRITE_UE(stream, delta_poc, "abs_delta_poc_st");
    if (delta_poc+1) WRITE_U(stream, 0, 1, "strp_entry_sign_flag");
    last_poc = delta_poc;
  }

}

void kvz_encoder_state_write_bitstream_slice_header(
    struct bitstream_t * const stream,
    struct encoder_state_t * const state,
    bool independent)
{
  const encoder_control_t * const encoder = state->encoder_control;

#ifdef KVZ_DEBUG
  printf("=========== Slice ===========\n");
#endif

  bool first_slice_segment_in_pic = (state->slice->start_in_rs == 0);
  if ((state->encoder_control->cfg.slices & KVZ_SLICES_WPP)
      && state->wfrow->lcu_offset_y > 0)
  {
    first_slice_segment_in_pic = false;
  }

  //WRITE_U(stream, first_slice_segment_in_pic, 1, "first_slice_segment_in_pic_flag");

  WRITE_U(stream, 1, 1, "picture_header_in_slice_header_flag");

  kvz_encoder_state_write_bitstream_picture_header(stream, state);

  if (state->frame->pictype != KVZ_NAL_IDR_W_RADL
    && state->frame->pictype != KVZ_NAL_IDR_N_LP) {
    WRITE_UE(stream, state->frame->slicetype, "sh_slice_type");
  }

  if (state->frame->pictype == KVZ_NAL_CRA_NUT || state->frame->pictype == KVZ_NAL_IDR_N_LP || state->frame->pictype == KVZ_NAL_IDR_W_RADL || state->frame->pictype == KVZ_NAL_GDR_NUT)
  {
    WRITE_U(stream, 0, 1, "sh_no_output_of_prior_pics_flag");
  }

    if (alf_enabled)
    {
      WRITE_U(stream, 0/*state->slice->tile_group_num_aps - 1*/, 3, "slice_num_alf_aps_ids_luma");
      const int* aps_ids = state->slice->tile_group_luma_aps_id;
      for (int i = 0; i < 0/*state->slice->tile_group_num_aps - 1*/; i++)
      {
        WRITE_U(stream, aps_ids[i], 3, "slice_alf_aps_id_luma");
      }
      const int alf_chroma_idc = state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] + state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] * 2;
      if (encoder->chroma_format != KVZ_CSP_400)
      {
        WRITE_U(stream, alf_chroma_idc, 2, "slice_alf_chroma_idc");
      }
      if (alf_chroma_idc)
      {
        WRITE_U(stream, state->slice->tile_group_chroma_aps_id, 3, "slice_alf_aps_id_chroma");
      }

      if (encoder->cfg.alf_cc_enabled_flag)
      {
        /*CcAlfFilterParam &filterParam = pcSlice->m_ccAlfFilterParam;
        WRITE_FLAG(filterParam.ccAlfFilterEnabled[COMPONENT_Cb - 1] ? 1 : 0, "slice_cc_alf_cb_enabled_flag");
        if (filterParam.ccAlfFilterEnabled[COMPONENT_Cb - 1])
        {
          // write CC ALF Cb APS ID
          WRITE_CODE(pcSlice->getTileGroupCcAlfCbApsId(), 3, "slice_cc_alf_cb_aps_id");
        }
        // Cr
        WRITE_FLAG(filterParam.ccAlfFilterEnabled[COMPONENT_Cr - 1] ? 1 : 0, "slice_cc_alf_cr_enabled_flag");
        if (filterParam.ccAlfFilterEnabled[COMPONENT_Cr - 1])
        {
          // write CC ALF Cr APS ID
          WRITE_CODE(pcSlice->getTileGroupCcAlfCrApsId(), 3, "slice_cc_alf_cr_aps_id");
        }*/
      }
    }
  }

  if (state->encoder_control->cfg.tmvp_enable) {
    //WRITE_U(stream, ref_negative ? 1 : 0, 1, "slice_temporal_mvp_enabled_flag");
    WRITE_U(stream, 0, 1, "sh_collocated_from_l0_flag");
  }

  int slice_qp_delta = state->frame->QP - encoder->cfg.qp;
  WRITE_SE(stream, slice_qp_delta, "sh_qp_delta");


  if (encoder->cfg.sao_type) {
    WRITE_U(stream, 1, 1, "sh_sao_luma_flag");
    if (encoder->chroma_format != KVZ_CSP_400) {
      WRITE_U(stream, 1, 1, "sh_sao_chroma_flag");
    }
  }

  if (state->frame->slicetype != KVZ_SLICE_I) {

    // BT Size set only with non-I-frames, in I-frames the size is 32x32
    // but in other frames it is CTU size >> <this value>
    //WRITE_UE(stream, 0, "max_binary_tree_unit_size"); // Max BT size == CTU size

  }

  if (encoder->tiles_enable || encoder->cfg.wpp) {
    int num_entry_points = 0;
    int max_length_seen = 0;

    if (state->is_leaf) {
      num_entry_points = 1;
    }
    else {
      encoder_state_entry_points_explore(state, &num_entry_points, &max_length_seen);
    }

    int num_offsets = num_entry_points - 1;

    //WRITE_UE(stream, num_offsets, "num_entry_point_offsets");
    if (num_offsets > 0) {
      int entry_points_written = 0;
      int offset_len = kvz_math_floor_log2(max_length_seen) + 1;
      WRITE_UE(stream, offset_len - 1, "offset_len_minus1");
      encoder_state_write_bitstream_entry_points_write(stream, state, num_entry_points, offset_len, &entry_points_written);
    }
  }

  //WRITE_U(stream, 0, 1, "slice_ts_residual_coding_disabled_flag");

  //kvz_bitstream_align(stream);
}


/**
 * \brief Add a checksum SEI message to the bitstream.
 * \param encoder The encoder.
 * \returns Void
 */
static void add_checksum(encoder_state_t * const state)
{
  bitstream_t * const stream = &state->stream;
  const videoframe_t * const frame = state->tile->frame;
  unsigned char checksum[3][SEI_HASH_MAX_LENGTH];

  kvz_nal_write(stream, KVZ_NAL_SUFFIX_SEI_NUT, 0, 0);

  WRITE_U(stream, 132, 8, "sei_type");

  int num_colors = (state->encoder_control->chroma_format == KVZ_CSP_400 ? 1 : 3);

  switch (state->encoder_control->cfg.hash)
  {
  case KVZ_HASH_CHECKSUM:
    kvz_image_checksum(frame->rec, checksum, state->encoder_control->bitdepth);

    WRITE_U(stream, 2 + num_colors * 4, 8, "size");
    WRITE_U(stream, 2, 8, "hash_type");  // 2 = checksum
    WRITE_U(stream, 0, 1, "dph_sei_single_component_flag");
    WRITE_U(stream, 0, 7, "dph_sei_reserved_zero_7bits");

    for (int i = 0; i < num_colors; ++i) {
      uint32_t checksum_val = (
        (checksum[i][0] << 24) + (checksum[i][1] << 16) +
        (checksum[i][2] << 8) + (checksum[i][3]));
      WRITE_U(stream, checksum_val, 32, "picture_checksum");
      CHECKPOINT("checksum[%d] = %u", i, checksum_val);
    }

    break;

  case KVZ_HASH_MD5:
    kvz_image_md5(frame->rec, checksum, state->encoder_control->bitdepth);

    WRITE_U(stream, 2 + num_colors * 16, 8, "size");
    WRITE_U(stream, 0, 8, "hash_type");  // 0 = md5
    WRITE_U(stream, 0, 1, "dph_sei_single_component_flag");
    WRITE_U(stream, 0, 7, "dph_sei_reserved_zero_7bits");

    for (int i = 0; i < num_colors; ++i) {
      for (int b = 0; b < 16; ++b) {
        WRITE_U(stream, checksum[i][b], 8, "picture_md5");
      }
    }

    break;

  case KVZ_HASH_NONE:
    // Means we shouldn't be writing this SEI.
    assert(0);
  }

  kvz_bitstream_align(stream);

  // spec:sei_rbsp() rbsp_trailing_bits
  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encoder_state_write_slice_header(
  bitstream_t * stream,
  encoder_state_t * state,
  bool independent)
{
  kvz_nal_write(stream, state->frame->pictype, (state->frame->pictype==KVZ_NAL_STSA)?1:0, state->frame->first_nal);
  state->frame->first_nal = false;

  kvz_encoder_state_write_bitstream_slice_header(stream, state, independent);
  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

/**
 * \brief Move child state bitstreams to the parent stream.
 */
static void encoder_state_write_bitstream_children(encoder_state_t * const state)
{
  // Write Slice headers to the parent stream instead of the child stream
  // in case the child stream is a leaf with something in it already.
  for (int i = 0; state->children[i].encoder_control; ++i) {
    if (state->children[i].type == ENCODER_STATE_TYPE_SLICE) {
      encoder_state_write_slice_header(&state->stream, &state->children[i], true);
    } else if (state->children[i].type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
      if ((state->encoder_control->cfg.slices & KVZ_SLICES_WPP) && i != 0) {
        // Add header for dependent WPP row slice.
        encoder_state_write_slice_header(&state->stream, &state->children[i], false);
      }
    }
    kvz_encoder_state_write_bitstream(&state->children[i]);
    kvz_bitstream_move(&state->stream, &state->children[i].stream);
  }
  
}

static void encoder_state_write_bitstream_main(encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
  bitstream_t * const stream = &state->stream;
  uint64_t curpos = kvz_bitstream_tell(stream);

  // The first NAL unit of the access unit must use a long start code.
  state->frame->first_nal = true;

  // Access Unit Delimiter (AUD)
  if (encoder->cfg.aud_enable) {
    state->frame->first_nal = false;
    encoder_state_write_bitstream_aud(state);
  }

  if (encoder_state_must_write_vps(state)) {
    state->frame->first_nal = false;
    kvz_encoder_state_write_parameter_sets(&state->stream, state);
  }

  // Send Kvazaar version information only in the first frame.
  if (state->frame->num == 0 && encoder->cfg.add_encoder_info) {
    kvz_nal_write(stream, KVZ_NAL_PREFIX_SEI_NUT, 0, state->frame->first_nal);
    state->frame->first_nal = false;
    encoder_state_write_bitstream_prefix_sei_version(state);

    // spec:sei_rbsp() rbsp_trailing_bits
    kvz_bitstream_add_rbsp_trailing_bits(stream);
  }

  //SEI messages for interlacing
  if (encoder->vui.frame_field_info_present_flag) {
    // These should be optional, needed for earlier versions
    // of HM decoder to accept bitstream
    //kvz_nal_write(stream, KVZ_NAL_PREFIX_SEI_NUT, 0, 0);
    //encoder_state_write_active_parameter_sets_sei_message(state);
    //kvz_bitstream_rbsp_trailing_bits(stream);

    kvz_nal_write(stream, KVZ_NAL_PREFIX_SEI_NUT, 0, state->frame->first_nal);
    state->frame->first_nal = false;
    encoder_state_write_picture_timing_sei_message(state);

    // spec:sei_rbsp() rbsp_trailing_bits
    kvz_bitstream_add_rbsp_trailing_bits(stream);
  }

  // Adaptation parameter set (APS)
  //send LMCS APS when LMCSModel is updated. It can be updated even current slice does not enable reshaper.
  //For example, in RA, update is on intra slice, but intra slice may not use reshaper
  if (0 /*pcSlice->getSPS()->getUseLmcs()*/)
  {
    /*//only 1 LMCS data for 1 picture
    int apsId = picHeader->getLmcsAPSId();
    ParameterSetMap<APS> *apsMap = m_pcEncLib->getApsMap();
    APS* aps = apsMap->getPS((apsId << NUM_APS_TYPE_LEN) + LMCS_APS);
    bool writeAPS = aps && apsMap->getChangedFlag((apsId << NUM_APS_TYPE_LEN) + LMCS_APS);
    if (writeAPS)
    {
    actualTotalBits += xWriteAPS(accessUnit, aps, m_pcEncLib->getLayerId(), true);
    apsMap->clearChangedFlag((apsId << NUM_APS_TYPE_LEN) + LMCS_APS);
    CHECK(aps != picHeader->getLmcsAPS(), "Wrong LMCS APS pointer in compressGOP");
    }*/
  }

  // Adaptation parameter set (APS)
  // only 1 SCALING LIST data for 1 picture
  if (0 /*pcSlice->getSPS()->getScalingListFlag() && (m_pcCfg->getUseScalingListId() == SCALING_LIST_FILE_READ)*/)
  {
    /*int apsId = picHeader->getScalingListAPSId();
    ParameterSetMap<APS> *apsMap = m_pcEncLib->getApsMap();
    APS* aps = apsMap->getPS((apsId << NUM_APS_TYPE_LEN) + SCALING_LIST_APS);
    bool writeAPS = aps && apsMap->getChangedFlag((apsId << NUM_APS_TYPE_LEN) + SCALING_LIST_APS);
    if (writeAPS)
    {
      actualTotalBits += xWriteAPS(accessUnit, aps, m_pcEncLib->getLayerId(), true);
      apsMap->clearChangedFlag((apsId << NUM_APS_TYPE_LEN) + SCALING_LIST_APS);
      CHECK(aps != picHeader->getScalingListAPS(), "Wrong SCALING LIST APS pointer in compressGOP");
    }*/
  }

  // Adaptation parameter set (APS)
  for (int i = 0; state->children[i].encoder_control; ++i) {
    if (state->children[i].type == ENCODER_STATE_TYPE_SLICE) {
      if (encoder->cfg.alf_enable && (state->children[i].slice->tile_group_alf_enabled_flag[COMPONENT_Y] || state->children[i].slice->tile_group_cc_alf_cb_enabled_flag || state->children[i].slice->tile_group_cc_alf_cr_enabled_flag))
      {
        for (int aps_id = ALF_CTB_MAX_NUM_APS - 0/*state->children[i].slice->tile_group_num_aps + 1*/; aps_id < ALF_CTB_MAX_NUM_APS; aps_id++) //HD: shouldn't this be looping over slice_alf_aps_id_luma[ i ]? By looping over MAX_NUM_APS, it is possible unused ALF APS is written. Please check!
        {
          //ParameterSetMap<APS> *apsMap = m_pcEncLib->getApsMap();
          param_set_map *aps_map = state->children[i].slice->param_set_map;

          //APS* aps = apsMap->getPS((apsId << NUM_APS_TYPE_LEN) + ALF_APS);
          alf_aps aps = aps_map[aps_id + T_ALF_APS].parameter_set;

          //bool writeAPS = aps && apsMap->getChangedFlag((apsId << NUM_APS_TYPE_LEN) + ALF_APS);
          bool writeAPS = aps_map[aps_id + T_ALF_APS].b_changed;
          if (!writeAPS && state->children[i].slice->apss && state->children[i].
            slice->apss[aps_id].aps_id >= 0 && state->children[i].slice->apss[aps_id].aps_id < 8)
          {
            writeAPS = true;
            //aps = pcSlice->getAlfAPSs()[apsId]; // use asp from slice header
            aps = state->children[i].slice->apss[aps_id];

            //*apsMap->allocatePS(apsId) = *aps; //allocate and cpy
            assert(aps_id < ALF_CTB_MAX_NUM_APS);
            bool found = false;
            for (int ID = 0; ID < ALF_CTB_MAX_NUM_APS; ID++)
            {
              if (aps_map[ID + T_ALF_APS].parameter_set.aps_id == ID) {
                found = true;
              }
            }
            if (!found)
            {
              aps_map[aps_id + T_ALF_APS].b_changed = true;
              //apsMap[apsId].p_nalu_data = 0;
            }
            aps_map[aps_id + T_ALF_APS].parameter_set.aps_id = aps_id;
            copy_alf_param(&aps_map[aps_id + T_ALF_APS].parameter_set, &aps);
            //

            //m_pcALF->setApsIdStart(apsId);
            g_aps_id_start = aps_id;
          }
          else if (state->children[i].slice->tile_group_cc_alf_cb_enabled_flag && !writeAPS && aps_id == state->children[i].slice->tile_group_cc_alf_cb_aps_id)
          {
            writeAPS = true;
            aps = aps_map[(state->children[i].slice->tile_group_cc_alf_cb_aps_id << NUM_APS_TYPE_LEN) + T_ALF_APS].parameter_set;
          }
          else if (state->children[i].slice->tile_group_cc_alf_cr_enabled_flag && !writeAPS && aps_id == state->children[i].slice->tile_group_cc_alf_cr_aps_id)
          {
            writeAPS = true;
            aps = aps_map[(state->children[i].slice->tile_group_cc_alf_cr_aps_id << NUM_APS_TYPE_LEN) + T_ALF_APS].parameter_set;
          }

          if (writeAPS)
          {
            //actualTotalBits += xWriteAPS(accessUnit, aps);
            kvz_nal_write(stream, NAL_UNIT_PREFIX_APS, 0, state->frame->first_nal);
            state->frame->first_nal = false;
            encoder_state_write_adaptation_parameter_set(state, &aps);

            //apsMap->clearChangedFlag((apsId << NUM_APS_TYPE_LEN) + ALF_APS);
            aps_map[aps_id + T_ALF_APS].b_changed = false;

            //CHECK(aps != pcSlice->getAlfAPSs()[apsId] && apsId != pcSlice->getTileGroupCcAlfCbApsId() && apsId != pcSlice->getTileGroupCcAlfCrApsId(), "Wrong APS pointer in compressGOP");
            assert(aps.aps_id == state->children[i].slice->apss[aps_id].aps_id
                  && aps_id != state->children[i].slice->tile_group_cc_alf_cr_aps_id
                  && aps_id != state->children[i].slice->tile_group_cc_alf_cr_aps_id); //"Wrong APS id");
          }
        }
      }
    }
  }

  encoder_state_write_bitstream_children(state);

  if (state->encoder_control->cfg.hash != KVZ_HASH_NONE) {
    // Calculate checksum
    add_checksum(state);
  }

  //kvz_nal_write(stream, KVZ_NAL_EOS_NUT, 0, 1);
  //kvz_nal_write(stream, KVZ_NAL_EOS_NUT, 0, 1);
  

  //Get bitstream length for stats
  uint64_t newpos = kvz_bitstream_tell(stream);
  state->stats_bitstream_length = (newpos >> 3) - (curpos >> 3);

  if (state->frame->num > 0) {
    state->frame->total_bits_coded = state->previous_encoder_state->frame->total_bits_coded;
  }
  state->frame->total_bits_coded += newpos - curpos;

    state->frame->cur_gop_bits_coded = state->previous_encoder_state->frame->cur_gop_bits_coded;
  state->frame->cur_gop_bits_coded += newpos - curpos;
}

void kvz_encoder_state_write_bitstream(encoder_state_t * const state)
{
  if (!state->is_leaf) {
    switch (state->type) {
      case ENCODER_STATE_TYPE_MAIN:
        encoder_state_write_bitstream_main(state);
        break;
      case ENCODER_STATE_TYPE_TILE:
      case ENCODER_STATE_TYPE_SLICE:
        encoder_state_write_bitstream_children(state);
        break;
      default:
        fprintf(stderr, "Unsupported node type %c!\n", state->type);
        assert(0);
    }
  }
}

void kvz_encoder_state_worker_write_bitstream(void * opaque)
{
  kvz_encoder_state_write_bitstream((encoder_state_t *) opaque);
}

void kvz_encoder_state_write_parameter_sets(bitstream_t *stream,
                                            encoder_state_t * const state)
{
  // Video Parameter Set (VPS)  
  //kvz_nal_write(stream, KVZ_NAL_VPS_NUT, 0, 1);
  //encoder_state_write_bitstream_vid_parameter_set(stream, state);
  

  // Sequence Parameter Set (SPS)
  kvz_nal_write(stream, KVZ_NAL_SPS_NUT, 0, 1);
  encoder_state_write_bitstream_seq_parameter_set(stream, state);

  // Picture Parameter Set (PPS)
  kvz_nal_write(stream, KVZ_NAL_PPS_NUT, 0, 1);
  encoder_state_write_bitstream_pic_parameter_set(stream, state);
}


