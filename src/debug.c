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

 /**
  * \file
  * Functions for debugging purposes.
  */

#include "encoderstate.h"
#include "global.h"

#ifdef KVZ_DEBUG_PRINT_YUVIEW_CSV

FILE* yuview_output;

void kvz_dbg_yuview_init(const encoder_control_t* const encoder, char* filename, char* sequence) {

  char buf[100];
  snprintf(buf, 100, "%s.csv", filename);

  yuview_output = fopen(buf, "wb");

  fprintf(yuview_output, "%%;syntax-version;v1.22\r\n");
  fprintf(yuview_output, "%%;%%width and height of Chroma TUs are scaled to luma resolution\r\n");
  fprintf(yuview_output, "%%;seq-specs;;%s;layer2;;%d;%d;%d\r\n", sequence, encoder->in.width, encoder->in.height, (int)encoder->cfg.framerate);
  fprintf(yuview_output, "%%;type;0;PCM;map\r\n");
  fprintf(yuview_output, "%%;mapColor;0;0;0;0;120\r\n");
  fprintf(yuview_output, "%%;type;1;IntraDirLuma;range\r\n");
  fprintf(yuview_output, "%%;defaultRange;0;35;Autumn\r\n");
  fprintf(yuview_output, "%%;type;2;IntraDirChroma;range\r\n");
  fprintf(yuview_output, "%%;defaultRange;0;36;Autumn\r\n");
  fprintf(yuview_output, "%%;type;3;RefIdxSkipL0;range\r\n");
  fprintf(yuview_output, "%%;range;0;4;0;65;75;154;160;255;255;255\r\n");
  fprintf(yuview_output, "%%;type;4;RefIdxSkipL1;range\r\n");
  fprintf(yuview_output, "%%;range;0;4;53;128;0;66;160;255;255;255\r\n");
  fprintf(yuview_output, "%%;type;5;RefIdxMergeL0;range\r\n");
  fprintf(yuview_output, "%%;range;0;4;75;154;160;255;0;65;255;255\r\n");
  fprintf(yuview_output, "%%;type;6;RefIdxMergeL1;range\r\n");
  fprintf(yuview_output, "%%;range;0;4;0;65;160;255;56;132;255;255\r\n");
  fprintf(yuview_output, "%%;type;7;RefIdxInterL0;range\r\n");
  fprintf(yuview_output, "%%;range;0;4;160;255;113;199;0;65;255;255\r\n");
  fprintf(yuview_output, "%%;type;8;RefIdxInterL1;range\r\n");
  fprintf(yuview_output, "%%;range;0;4;146;239;160;255;0;65;255;255\r\n");
  fprintf(yuview_output, "%%;type;9;MVSkipL0;vector\r\n");
  fprintf(yuview_output, "%%;vectorColor;0;0;0;255\r\n");
  fprintf(yuview_output, "%%;scaleFactor;4\r\n");
  fprintf(yuview_output, "%%;type;10;MVSkipL1;vector\r\n");
  fprintf(yuview_output, "%%;vectorColor;255;255;255;255\r\n");
  fprintf(yuview_output, "%%;scaleFactor;4\r\n");
  fprintf(yuview_output, "%%;type;11;MVMergeL0;vector\r\n");
  fprintf(yuview_output, "%%;vectorColor;0;0;0;255\r\n");
  fprintf(yuview_output, "%%;scaleFactor;4\r\n");
  fprintf(yuview_output, "%%;type;12;MVMergeL1;vector\r\n");
  fprintf(yuview_output, "%%;vectorColor;255;255;255;255\r\n");
  fprintf(yuview_output, "%%;scaleFactor;4\r\n");
  fprintf(yuview_output, "%%;type;13;MVInterL0;vector\r\n");
  fprintf(yuview_output, "%%;vectorColor;0;0;0;255\r\n");
  fprintf(yuview_output, "%%;scaleFactor;4\r\n");
  fprintf(yuview_output, "%%;type;14;MVInterL1;vector\r\n");
  fprintf(yuview_output, "%%;vectorColor;255;255;255;255\r\n");
  fprintf(yuview_output, "%%;scaleFactor;4\r\n");
  fprintf(yuview_output, "%%;type;15;NumSigCoeffY;range\r\n");
  fprintf(yuview_output, "%%;range;0;1024;0;255;0;255;0;255;255;255\r\n");
  fprintf(yuview_output, "%%;scaleToBlockSize;1\r\n");
  fprintf(yuview_output, "%%;type;16;NumSigCoeffU;range\r\n");
  fprintf(yuview_output, "%%;range;0;1024;0;255;0;255;0;255;255;255\r\n");
  fprintf(yuview_output, "%%;scaleToBlockSize;1\r\n");
  fprintf(yuview_output, "%%;type;17;NumSigCoeffV;range\r\n");
  fprintf(yuview_output, "%%;range;0;1024;0;255;0;255;0;255;255;255\r\n");
  fprintf(yuview_output, "%%;scaleToBlockSize;1\r\n");
  
}

void kvz_dbg_yuview_add_vector(const encoder_state_t* const state, int x, int y, int width, int height, int type, int x_vec, int y_vec) {
  fprintf(yuview_output, "%d;%d;%d;%d;%d;%d;%d;%d\r\n", state->tile->frame->poc, x, y, width, height, type, x_vec, y_vec);
}

void kvz_dbg_yuview_add(const encoder_state_t* const state, int x, int y, int width, int height, int type, int val) {
  fprintf(yuview_output, "%d;%d;%d;%d;%d;%d;%d\r\n", state->tile->frame->poc, x, y, width, height, type, val);
}

void kvz_dbg_yuview_cleanup() {
  fclose(yuview_output);
  yuview_output = NULL;
}

#endif //KVZ_DEBUG_PRINT_YUVIEW_CSV

#ifdef KVZ_DEBUG_PRINT_THREADING_INFO
void kvz_dbg_encoder_state_dump_graphviz(const encoder_state_t* const state) {
  int i;

  if (!state->parent) {
    const encoder_control_t* const encoder = state->encoder_control;
    int y, x;
    //Empty lines (easier to copy-paste)
    printf("\n\n\n\n\n");
    //Some styling...
    printf("digraph EncoderStates {\n");
    printf(" fontname = \"Bitstream Vera Sans\"\n");
    printf(" fontsize = 8\n\n");
    printf(" node [\n");
    printf("  fontname = \"Bitstream Vera Sans\"\n");
    printf("  fontsize = 8\n");
    printf("  shape = \"record\"\n");
    printf(" ]\n\n");
    printf(" edge [\n");
    printf("  arrowtail = \"empty\"\n");
    printf(" ]\n\n");

    printf(" \"Map\" [\n");
    printf("  shape=plaintext\n");
    printf("  label = <<table cellborder=\"1\" cellspacing=\"0\" border=\"0\">");
    printf("<tr><td colspan=\"%d\" height=\"20\" valign=\"bottom\"><b>RS Map</b></td></tr>", encoder->in.width_in_lcu);
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      printf("<tr>");
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;

        printf("<td>%d</td>", lcu_id_rs);
      }
      printf("</tr>");
    }
    printf("<tr><td colspan=\"%d\" height=\"20\" valign=\"bottom\"><b>TS Map</b></td></tr>", encoder->in.width_in_lcu);
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      printf("<tr>");
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        const int lcu_id_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_id_rs];

        printf("<td>%d</td>", lcu_id_ts);
      }
      printf("</tr>");
    }
    printf("<tr><td colspan=\"%d\" height=\"20\" valign=\"bottom\"><b>Tile map</b></td></tr>", encoder->in.width_in_lcu);
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      printf("<tr>");
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        const int lcu_id_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_id_rs];

        printf("<td>%d</td>", encoder->tiles_tile_id[lcu_id_ts]);
      }
      printf("</tr>");
    }
    printf("<tr><td colspan=\"%d\" height=\"20\" valign=\"bottom\"><b>Slice map</b></td></tr>", encoder->in.width_in_lcu);
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      printf("<tr>");
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        const int lcu_id_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_id_rs];
        int slice_id = 0;

        //Not efficient, but who cares
        for (i = 0; i < encoder->slice_count; ++i) {
          if (encoder->slice_addresses_in_ts[i] <= lcu_id_ts) {
            slice_id = i;
          }
        }

        printf("<td>%d</td>", slice_id);
      }
      printf("</tr>");
    }
    printf("</table>>\n ]\n");
  }

  printf(" \"%p\" [\n", state);
  printf("  label = \"{encoder_state|");
  printf("+ type=%c\\l", state->type);
  if (!state->parent || state->frame != state->parent->global) {
    printf("|+ global\\l");
  }
  if (!state->parent || state->tile != state->parent->tile) {
    printf("|+ tile\\l");
    printf(" - id = %d\\l", state->tile->id);
    printf(" - lcu_offset_x = %d\\l", state->tile->lcu_offset_x);
    printf(" - lcu_offset_y = %d\\l", state->tile->lcu_offset_y);
    printf(" - lcu_offset_in_ts = %d\\l", state->tile->lcu_offset_in_ts);
  }
  if (!state->parent || state->slice != state->parent->slice) {
    printf("|+ slice\\l");
    printf(" - id = %d\\l", state->slice->id);
    printf(" - start_in_ts = %d\\l", state->slice->start_in_ts);
    printf(" - end_in_ts = %d\\l", state->slice->end_in_ts);
    printf(" - start_in_rs = %d\\l", state->slice->start_in_rs);
    printf(" - end_in_rs = %d\\l", state->slice->end_in_rs);
  }
  if (!state->parent || state->wfrow != state->parent->wfrow) {
    printf("|+ wfrow\\l");
    printf(" - lcu_offset_y = %d\\l", state->wfrow->lcu_offset_y);
  }
  printf("}\"\n");
  printf(" ]\n");

  if (state->parent) {
    printf(" \"%p\" -> \"%p\"\n", state->parent, state);
  }

  for (i = 0; state->children[i].encoder_control; ++i) {
    encoder_state_dump_graphviz(&state->children[i]);
  }

  if (!state->parent) {
    printf("}\n");
    //Empty lines (easier to copy-paste)
    printf("\n\n\n\n\n");
  }
}
#endif //KVZ_DEBUG_PRINT_THREADING_INFO