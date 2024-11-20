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

#include "cabac.h"

#include "encoder.h"
#include "encoderstate.h"
#include "uvg266.h"

#ifdef UVG_DEBUG_PRINT_CABAC
uint32_t uvg_cabac_bins_count = 0;
bool uvg_cabac_bins_verbose = true;
#endif

const uint8_t uvg_g_auc_renorm_table[32] =
{
  6, 5, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

static const uint8_t uvg_tb_max[257] = { 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8 };


/**
 * \brief Initialize struct cabac_data.
 */
void uvg_cabac_start(cabac_data_t * const data)
{
  data->low = 0;
  data->range = 510;
  data->bits_left = 23;
  data->num_buffered_bytes = 0;
  data->buffered_byte = 0xff;
  data->only_count = 0; // By default, write bits out
  data->update = 0; 
}

/**
 * \brief
 */
void uvg_cabac_encode_bin(cabac_data_t * const data, const uint32_t bin_value)
{
  uint32_t lps = CTX_LPS(data->cur_ctx, data->range);

  data->range -= lps;

  // Not the Most Probable Symbol?
  if ((bin_value ? 1 : 0) != CTX_MPS(data->cur_ctx)) {
    int num_bits = uvg_g_auc_renorm_table[lps >> 3];
    data->low = (data->low + data->range) << num_bits;
    data->range = lps << num_bits;
    
    data->bits_left -= num_bits;
    if (data->bits_left < 12) {
      uvg_cabac_write(data);
    }
  } else {    
    if (data->range < 256) {
      data->low <<= 1;
      data->range <<= 1;
      data->bits_left--;

      if (data->bits_left < 12) {
        uvg_cabac_write(data);
      }
    }
  }
  CTX_UPDATE(data->cur_ctx, bin_value);
}

/**
 * \brief
 */
void uvg_cabac_write(cabac_data_t * const data)
{
  uint32_t lead_byte = data->low >> (24 - data->bits_left);
  data->bits_left += 8;
  data->low &= 0xffffffffu >> data->bits_left;

  // Binary counter mode
  if(data->only_count) {
    data->num_buffered_bytes++;
    return;
  }

  if (lead_byte == 0xff) {
    data->num_buffered_bytes++;
  } else {
    if (data->num_buffered_bytes > 0) {
      uint32_t carry = lead_byte >> 8;
      uint32_t byte = data->buffered_byte + carry;
      data->buffered_byte = lead_byte & 0xff;
      uvg_bitstream_put_byte(data->stream, byte);

      byte = (0xff + carry) & 0xff;
      while (data->num_buffered_bytes > 1) {
        uvg_bitstream_put_byte(data->stream, byte);
        data->num_buffered_bytes--;
      }
    } else {
      data->num_buffered_bytes = 1;
      data->buffered_byte = lead_byte;
    }
  }
}

/**
 * \brief
 */
void uvg_cabac_finish(cabac_data_t * const data)
{
  assert(data->bits_left <= 32);

  if (data->low >> (32 - data->bits_left)) {
    uvg_bitstream_put_byte(data->stream, data->buffered_byte + 1);
    while (data->num_buffered_bytes > 1) {
      uvg_bitstream_put_byte(data->stream, 0);
      data->num_buffered_bytes--;
    }
    data->low -= 1 << (32 - data->bits_left);
  } else {
    if (data->num_buffered_bytes > 0) {
      uvg_bitstream_put_byte(data->stream, data->buffered_byte);
    }
    while (data->num_buffered_bytes > 1) {
      uvg_bitstream_put_byte(data->stream, 0xff);
      data->num_buffered_bytes--;
    }
  }

  {
    uint8_t bits = (uint8_t)(24 - data->bits_left);
    uvg_bitstream_put(data->stream, data->low >> 8, bits);
  }
}

/*!
  \brief Encode terminating bin
  \param binValue bin value
*/
void uvg_cabac_encode_bin_trm(cabac_data_t * const data, const uint8_t bin_value)
{
  data->range -= 2;
  if(bin_value) {
    data->low += data->range;
    data->low <<= 7;
    data->range = 2 << 7;
    data->bits_left -= 7;
  } else if (data->range >= 256) {
    return;
  } else {
    data->low <<= 1;
    data->range <<= 1;
    data->bits_left--;
  }

  if (data->bits_left < 12) {
    uvg_cabac_write(data);
  }
}

/**
 * \brief encode truncated binary code
 */
void uvg_cabac_encode_trunc_bin(cabac_data_t * const data, const uint32_t bin_value, const uint32_t max_value, double* bits_out) {
  int thresh;
  int symbol = bin_value;
  if (max_value > 256) {
    uint32_t threshVal = 1 << 8;
    thresh = 8;
    while (threshVal <= max_value) {
      thresh++;
      threshVal <<= 1;
    }
    thresh--;
  } else {
    thresh = uvg_tb_max[max_value];
  }

  int val = 1 << thresh;

  int b = max_value - val;
  if (symbol < val - b) {
    CABAC_BINS_EP(data, symbol, thresh, "TruncSymbols");
    if (bits_out) *bits_out += thresh;
  } else {
    symbol += val - b;
    CABAC_BINS_EP(data, symbol, thresh + 1, "TruncSymbols");
    if (bits_out) *bits_out += thresh + 1;
  }
}


/**
 * \brief
 */
void uvg_cabac_encode_bin_ep(cabac_data_t * const data, const uint32_t bin_value)
{
  data->low <<= 1;
  if (bin_value) {
    data->low += data->range;
  }
  data->bits_left--;

  if (data->bits_left < 12) {
    uvg_cabac_write(data);
  }
}

// Import from VTM 4.0
void uvg_cabac_encode_aligned_bins_ep(cabac_data_t * const data, const uint32_t bin_values, int num_bins) 
{
  uint32_t rem_bins = num_bins;
  while (rem_bins > 0) {
    //The process of encoding an EP bin is the same as that of coding a normal
    //bin where the symbol ranges for 1 and 0 are both half the range:
    //
    //  low = (low + range/2) << 1       (to encode a 1)
    //  low =  low            << 1       (to encode a 0)
    //
    //  i.e.
    //  low = (low + (bin * range/2)) << 1
    //
    //  which is equivalent to:
    //
    //  low = (low << 1) + (bin * range)
    //
    //  this can be generalised for multiple bins, producing the following expression:
    //
    unsigned bins_to_code = MIN(rem_bins, 8); //code bytes if able to take advantage of the system's byte-write function
    unsigned bin_mask = (1 << bins_to_code) - 1;
    unsigned new_bins = (bin_values >> (rem_bins - bins_to_code)) & bin_mask;
    data->low = (data->low << bins_to_code) + (new_bins << 8); //range is known to be 256
    rem_bins -= bins_to_code;
    data->bits_left -= bins_to_code;
    if (data->bits_left < 12) {
      uvg_cabac_write(data);
    }
  }
}

/**
 * \brief
 */
void uvg_cabac_encode_bins_ep(cabac_data_t * const data, uint32_t bin_values, int num_bins)
{
  uint32_t pattern;
  
  if (data->range == 256) {
    uvg_cabac_encode_aligned_bins_ep(data, bin_values, num_bins);
    return;
  }
  while (num_bins > 8) {
    num_bins -= 8;
    pattern = bin_values >> num_bins;
    data->low <<= 8;
    data->low += data->range * pattern;
    bin_values -= pattern << num_bins;
    data->bits_left -= 8;

    if(data->bits_left < 12) {
      uvg_cabac_write(data);
    }
  }

  data->low <<= num_bins;
  data->low += data->range * bin_values;
  data->bits_left -= num_bins;

  if (data->bits_left < 12) {
    uvg_cabac_write(data);
  }
}

/**
 * \brief Coding of remainder abs coeff value.
 * \param remainder Value of remaining abs coeff
 * \param rice_param Reference to Rice parameter.
 */
int uvg_cabac_write_coeff_remain(cabac_data_t * const cabac, const uint32_t remainder, const uint32_t rice_param, const unsigned int cutoff)
{
  const unsigned threshold = cutoff << rice_param;
  uint32_t bins = remainder;
  uint32_t bits = 0;
  if (bins < threshold) {
    uint32_t length = (bins >> rice_param) + 1;
    CABAC_BINS_EP(cabac, ((1 << (length)) - 2) , length, "coeff_abs_level_remaining");
    CABAC_BINS_EP(cabac, bins & ((1 << rice_param) - 1), rice_param, "coeff_abs_level_remaining");
    bits += length;
    bits += rice_param;
  } else {
    const unsigned max_prefix_length = 32 - cutoff - 15/*max_log2_tr_dynamic_range*/;
    unsigned prefix_length = 0;
    unsigned code_value = (bins >> rice_param) - cutoff;
    unsigned suffix_length;
    if ((int32_t)code_value >= ((1 << max_prefix_length) - 1)) {
      prefix_length = max_prefix_length;
      suffix_length = 15 /*max_log2_tr_dynamic_range*/;
    } else {
      while ((int32_t)code_value > ((2 << prefix_length) - 2)) {
        prefix_length++;
      }
      suffix_length = prefix_length + rice_param + 1;
    }
    const unsigned total_prefix_length = prefix_length + cutoff;
    const unsigned bit_mask = (1 << rice_param) - 1;
    const unsigned prefix = (1 << total_prefix_length) - 1;
    const unsigned suffix = ((code_value - ((1 << prefix_length) - 1)) << rice_param) | (bins & bit_mask);
    CABAC_BINS_EP(cabac, prefix, total_prefix_length, "coeff_abs_level_remaining");
    CABAC_BINS_EP(cabac, suffix, suffix_length, "coeff_abs_level_remaining");
    bits += total_prefix_length + suffix_length;
  }
  return bits;
}


/**
 * \brief
 */
void uvg_cabac_write_unary_max_symbol(cabac_data_t * const data, 
  cabac_ctx_t * const ctx, 
  uint32_t symbol,
  const int32_t offset,
  const uint32_t max_symbol, 
  double* bits_out)
{
  int8_t code_last = max_symbol > symbol;

  assert(symbol <= max_symbol);

  if (!max_symbol) return;
  
  CABAC_FBITS_UPDATE(data, ctx, symbol, *bits_out, "ums");

  if (!symbol) return;

  data->cur_ctx = &ctx[offset];

  while (--symbol) {
    CABAC_FBITS_UPDATE(data, &ctx[offset], 1, *bits_out, "ums");
  }
  if (code_last) {
    CABAC_FBITS_UPDATE(data, &ctx[offset], 0,*bits_out, "ums");
  }
}

/**
 * This can be used for Truncated Rice binarization with cRiceParam=0.
 */
void uvg_cabac_write_unary_max_symbol_ep(cabac_data_t * const data, unsigned int symbol, const unsigned int max_symbol)
{
  /*if (symbol == 0) {
    CABAC_BIN_EP(data, 0, "ums_ep");
  } else {
    // Make a bit-string of (symbol) times 1 and a single 0, except when
    // symbol == max_symbol.
    unsigned bins = ((1 << symbol) - 1) << (symbol < max_symbol);
    CABAC_BINS_EP(data, bins, symbol + (symbol < max_symbol), "ums_ep");
  }*/

  int8_t code_last = max_symbol > symbol;

  assert(symbol <= max_symbol);

  CABAC_BIN_EP(data, symbol ? 1 : 0, "ums_ep");

  if (!symbol) return;

  while (--symbol) {
    CABAC_BIN_EP(data, 1, "ums_ep");
  }
  if (code_last) {
    CABAC_BIN_EP(data, 0, "ums_ep");
  }
}

/**
 * \brief
 */
uint32_t uvg_cabac_write_ep_ex_golomb(encoder_state_t * const state,
                                  cabac_data_t * const data,
                                  uint32_t symbol,
                                  uint32_t count)
{
  uint32_t bins = 0;
  int32_t num_bins = 0;

  while (symbol >= (uint32_t)(1 << count)) {
    bins = 2 * bins + 1;
    ++num_bins;
    symbol -= 1 << count;
    ++count;
  }
  bins = 2 * bins;
  ++num_bins;

  bins      = (bins << count) | symbol;
  num_bins += count;

  CABAC_BINS_EP(data, bins, num_bins, "ep_ex_golomb");
  return num_bins;
}
