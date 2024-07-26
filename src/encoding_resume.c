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

#include "encoding_resume.h"
#if defined(UVG_ENCODING_RESUME)

#define PROCESS_VAR(type, num, var, file, read_mode) {\
  if (read_mode){\
    fread(var, sizeof(type), num, file);\
  } else {\
    fwrite(var, sizeof(type), num, file);\
  }\
}\

#define PROCESS_FIELD(type, s_field, file, read_mode) {\
  type tmp_var;\
  if (read_mode){\
    fread(&tmp_var, sizeof(type), 1, file);\
    s_field = tmp_var;\
  } else {\
    tmp_var = s_field;\
    fwrite(&tmp_var, sizeof(type), 1, file);\
  }\
}\

static const char* const resume_fname_format = "%s/%s_poc%d_(%d,%d)%s";
//#define RESUME_MAX_NUM_PREC "9999"
//#define RESUME_MAX_FTYPE_LEN 6

#define GET_RESUME_ID(chroma_only) (chroma_only) ? "enc_resume_c" : "enc_resume"
#define GET_REC_FILENAME(fname, id, poc, posx, posy) snprintf(fname, FILENAME_MAX, resume_fname_format, RESUME_DIRNAME, id, poc, posx, posy, ".rec")
#define GET_COEFF_FILENAME(fname, id, poc, posx, posy) snprintf(fname, FILENAME_MAX, resume_fname_format, RESUME_DIRNAME, id, poc, posx, posy, ".coeff")
#define GET_CU_FILENAME(fname, id, poc, posx, posy) snprintf(fname, FILENAME_MAX, resume_fname_format, RESUME_DIRNAME, id, poc, posx, posy, ".cu")

static char tmp_fname[FILENAME_MAX];

#define OPEN_RESUME_FILE(file, fname_macro, id, poc, posx, posy, mode){\
  fname_macro(tmp_fname, id, poc, posx, posy);\
  file = fopen(tmp_fname, mode);\
}\


static bool use_resume(const encoder_state_t * const state, const int x, const int y, const bool chroma_only) {
#ifdef RESUME_SLICE_COND
  if (state->frame->slicetype != RESUME_SLICE_COND) return false;
#endif // RESUME_SLICE_COND

#ifdef RESUME_POC_LTE_COND
  if (state->frame->poc > RESUME_POC_LTE_COND) return false;
#endif // RESUME_POC_LTE_COND

#ifdef RESUME_X_LTE_COND
  if (x > RESUME_X_LTE_COND) return false;
#endif // RESUME_X_LTE_COND

#ifdef RESUME_Y_LTE_COND
  if (y > RESUME_Y_LTE_COND) return false;
#endif // RESUME_Y_LTE_COND

#ifdef RESUME_NO_CHROMA_COND
  if (chroma_only) return false;
#endif // RESUME_CHROMA_COND

  return true;
}

//Check if resume data is available
bool uvg_can_resume_encoding(const encoder_state_t * const state, const int x, const int y, const bool chroma_only) {
  //char* tmp_fname = malloc(strlen(resume_fname_format) + strlen(RESUME_DIRNAME) + strlen(id) + strlen(RESUME_MAX_NUM_PREC) * 2 + RESUME_MAX_FTYPE_LEN + 1);
  
  if (!use_resume(state, x, y, chroma_only)) return false;

  const char* const id = GET_RESUME_ID(chroma_only);
  const int32_t poc = state->frame->poc;
  const int posx = x >> LOG2_LCU_WIDTH;
  const int posy = y >> LOG2_LCU_WIDTH;//state->wfrow->lcu_offset_y;
  FILE* tmp_file;

    //Check if data files exist
  OPEN_RESUME_FILE(tmp_file, GET_REC_FILENAME, id, poc, posx, posy, "r");
  if (!tmp_file) return false;
  fclose(tmp_file);

  OPEN_RESUME_FILE(tmp_file, GET_COEFF_FILENAME, id, poc, posx, posy, "r");
  if (!tmp_file) return false;
  fclose(tmp_file);

  OPEN_RESUME_FILE(tmp_file, GET_CU_FILENAME, id, poc, posx, posy, "r");
  if (!tmp_file) return false;
  fclose(tmp_file);

  //free(tmp_fname);
  return true;
}

static void process_rec(FILE* const file, lcu_yuv_t* const rec, bool read_mode) {
  PROCESS_VAR(uvg_pixel, LCU_LUMA_SIZE, rec->y, file, read_mode);
  PROCESS_VAR(uvg_pixel, LCU_CHROMA_SIZE, rec->u, file, read_mode);
  PROCESS_VAR(uvg_pixel, LCU_CHROMA_SIZE, rec->v, file, read_mode);
  PROCESS_VAR(uvg_pixel, LCU_CHROMA_SIZE, rec->joint_u, file, read_mode);
  PROCESS_VAR(uvg_pixel, LCU_CHROMA_SIZE, rec->joint_v, file, read_mode);
}

static void process_coeff(FILE* const file, lcu_coeff_t* const coeff, bool read_mode) {
  PROCESS_VAR(coeff_t, LCU_LUMA_SIZE, coeff->y, file, read_mode);
  PROCESS_VAR(coeff_t, LCU_CHROMA_SIZE, coeff->u, file, read_mode);
  PROCESS_VAR(coeff_t, LCU_CHROMA_SIZE, coeff->v, file, read_mode);
  PROCESS_VAR(coeff_t, LCU_CHROMA_SIZE, coeff->joint_uv, file, read_mode);
}

static void process_cu(FILE* const file, lcu_t* const lcu, bool read_mode) {
  for (unsigned y = 0; y < LCU_WIDTH; y += SCU_WIDTH) {
    for (unsigned x = 0; x < LCU_WIDTH; x += SCU_WIDTH) {
      cu_info_t* cu = LCU_GET_CU_AT_PX(lcu, x, y);

      PROCESS_FIELD(uint8_t, cu->type, file, read_mode);
      if (cu->type == CU_NOTSET) continue;

      PROCESS_FIELD(uint8_t, cu->skipped, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->merged, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->merge_idx, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->tr_skip, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->tr_idx, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->joint_cb_cr, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->log2_width, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->log2_height, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->log2_chroma_width, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->log2_chroma_height, file, read_mode);
      PROCESS_VAR(uint16_t, 1, &(cu->cbf), file, read_mode);
      PROCESS_VAR(uint8_t, 1, &(cu->root_cbf), file, read_mode);
      PROCESS_FIELD(uint32_t, cu->split_tree, file, read_mode);
      PROCESS_VAR(uint8_t, 1, &(cu->qp), file, read_mode);
      PROCESS_FIELD(uint8_t, cu->bdpcmMode, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->violates_mts_coeff_constraint, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->mts_last_scan_pos, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->violates_lfnst_constrained_luma, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->violates_lfnst_constrained_chroma, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->lfnst_last_scan_pos, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->lfnst_idx, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->cr_lfnst_idx, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->luma_deblocking, file, read_mode);
      PROCESS_FIELD(uint8_t, cu->chroma_deblocking, file, read_mode);

      if (cu->type == CU_INTER) {
        PROCESS_VAR(mv_t   , 2*2, cu->inter.mv, file, read_mode);  
        PROCESS_VAR(uint8_t, 2, cu->inter.mv_ref, file, read_mode);
        PROCESS_FIELD(uint8_t, cu->inter.mv_cand0, file, read_mode);
        PROCESS_FIELD(uint8_t, cu->inter.mv_cand1, file, read_mode);
        PROCESS_FIELD(uint8_t, cu->inter.mv_dir, file, read_mode);
        PROCESS_FIELD(uint8_t, cu->inter.imv, file, read_mode);
      } else {
        PROCESS_VAR(int8_t, 1, &(cu->intra.mode), file, read_mode);
        PROCESS_VAR(int8_t, 1, &(cu->intra.mode_chroma), file, read_mode);
        PROCESS_VAR(uint8_t, 1, &(cu->intra.multi_ref_idx), file, read_mode);
        PROCESS_VAR(int8_t, 1, &(cu->intra.mip_flag), file, read_mode);
        PROCESS_VAR(int8_t, 1, &(cu->intra.mip_is_transposed), file, read_mode);
        PROCESS_VAR(int8_t, 1, &(cu->intra.isp_mode), file, read_mode);
        PROCESS_FIELD(uint8_t, cu->intra.isp_cbfs, file, read_mode);
        PROCESS_FIELD(uint8_t, cu->intra.isp_index, file, read_mode);
      }
    }
  }
}


void uvg_process_resume_encoding(const encoder_state_t * const state, const int x, const int y, const bool chroma_only, lcu_t * const lcu, bool read_mode)
{
  if (!use_resume(state, x, y, chroma_only)) return;

  FILE* file;
  const char* const id = GET_RESUME_ID(chroma_only);
  const int32_t poc = state->frame->poc;
  const int posx = x >> LOG2_LCU_WIDTH;
  const int posy = y >> LOG2_LCU_WIDTH;//state->wfrow->lcu_offset_y;
  const char* const mode = read_mode ? "rb" : "wb";
  
  OPEN_RESUME_FILE(file, GET_REC_FILENAME, id, poc, posx, posy, mode);
  process_rec(file, &(lcu->rec), read_mode);
  assert(!ferror(file) && "Processing failed");
  //lcu->rec.chroma_format = state->encoder_control->chroma_format;
  fclose(file);
  
  OPEN_RESUME_FILE(file, GET_COEFF_FILENAME, id, poc, posx, posy, mode);
  process_coeff(file, &(lcu->coeff), read_mode);
  assert(!ferror(file) && "Processing failed");
  fclose(file);

  OPEN_RESUME_FILE(file, GET_CU_FILENAME, id, poc, posx, posy, mode);
  process_cu(file, lcu, read_mode);
  assert(!ferror(file) && "Processing failed");
  fclose(file);

}

#endif