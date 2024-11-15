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

#include "strategies/strategies-picture.h"

#include "strategies/altivec/picture-altivec.h"
#if defined(__AVX512F__)
#include "strategies/avx2/picture-avx2.h"
#endif
#include "strategies/generic/picture-generic.h"
#include "strategies/sse2/picture-sse2.h"
#include "strategies/sse41/picture-sse41.h"
#include "strategies/sse42/picture-sse42.h"
#include "strategyselector.h"


// Define function pointers.
crc32c_4x4_func * uvg_crc32c_4x4 = 0;
crc32c_8x8_func * uvg_crc32c_8x8 = 0;
reg_sad_func * uvg_reg_sad = 0;

cost_pixel_nxn_func * uvg_sad_4x4 = 0;
cost_pixel_nxn_func * uvg_sad_8x8 = 0;
cost_pixel_nxn_func * uvg_sad_16x16 = 0;
cost_pixel_nxn_func * uvg_sad_32x32 = 0;
cost_pixel_nxn_func * uvg_sad_64x64 = 0;

cost_pixel_nxn_func * uvg_satd_4x4 = 0;
cost_pixel_nxn_func * uvg_satd_8x8 = 0;
cost_pixel_nxn_func * uvg_satd_16x16 = 0;
cost_pixel_nxn_func * uvg_satd_32x32 = 0;
cost_pixel_nxn_func * uvg_satd_64x64 = 0;

cost_pixel_nxn_multi_func * uvg_sad_4x4_dual = 0;
cost_pixel_nxn_multi_func * uvg_sad_8x8_dual = 0;
cost_pixel_nxn_multi_func * uvg_sad_16x16_dual = 0;
cost_pixel_nxn_multi_func * uvg_sad_32x32_dual = 0;
cost_pixel_nxn_multi_func * uvg_sad_64x64_dual = 0;

cost_pixel_nxn_multi_func * uvg_satd_4x4_dual = 0;
cost_pixel_nxn_multi_func * uvg_satd_8x8_dual = 0;
cost_pixel_nxn_multi_func * uvg_satd_16x16_dual = 0;
cost_pixel_nxn_multi_func * uvg_satd_32x32_dual = 0;
cost_pixel_nxn_multi_func * uvg_satd_64x64_dual = 0;

cost_pixel_any_size_func * uvg_satd_any_size = 0;
cost_pixel_any_size_func * uvg_satd_any_size_vtm = 0;
cost_pixel_any_size_multi_func * uvg_satd_any_size_quad = 0;

pixels_calc_ssd_func * uvg_pixels_calc_ssd = 0;

inter_recon_bipred_func * uvg_bipred_average = 0;

get_optimized_sad_func *uvg_get_optimized_sad = 0;
ver_sad_func *uvg_ver_sad = 0;
hor_sad_func *uvg_hor_sad = 0;

pixel_var_func *uvg_pixel_var = 0;

generate_residual_func *uvg_generate_residual = 0;




int uvg_strategy_register_picture(void* opaque, uint8_t bitdepth) {
  bool success = true;

  success &= uvg_strategy_register_picture_generic(opaque, bitdepth);

  if (uvg_g_hardware_flags.intel_flags.sse2) {
    success &= uvg_strategy_register_picture_sse2(opaque, bitdepth);
  }
  if (uvg_g_hardware_flags.intel_flags.sse41) {
    success &= uvg_strategy_register_picture_sse41(opaque, bitdepth);
  }
  if (uvg_g_hardware_flags.intel_flags.sse42) {
    success &= uvg_strategy_register_picture_sse42(opaque, bitdepth);
  }
#if defined(__AVX512F__)
  if (uvg_g_hardware_flags.intel_flags.avx2) {
    success &= uvg_strategy_register_picture_avx2(opaque, bitdepth);
  }
#endif
  if (uvg_g_hardware_flags.powerpc_flags.altivec) {
    success &= uvg_strategy_register_picture_altivec(opaque, bitdepth);
  }

  return success;
}


/**
* \brief  Get a function that calculates SATD for NxN block.
*
* \param width  Width of the region for which SATD is calculated.
*
* \returns  Pointer to cost_16bit_nxn_func.
*/
cost_pixel_nxn_func * uvg_pixels_get_satd_func(unsigned width, unsigned height)
{
  if(width == height) {
    switch (width) {
      case 4:
        return uvg_satd_4x4;
      case 8:
        return uvg_satd_8x8;
      case 16:
        return uvg_satd_16x16;
      case 32:
        return uvg_satd_32x32;
      case 64:
        return uvg_satd_64x64;
      default:
        return NULL;
    }
  }
  return NULL;
}


/**
* \brief  Get a function that calculates SAD for NxN block.
*
* \param width  Width of the region for which SAD is calculated.
*
* \returns  Pointer to cost_16bit_nxn_func.
*/
cost_pixel_nxn_func * uvg_pixels_get_sad_func(unsigned width, unsigned height)
{
  if(width == height) {
    switch (width) {
      case 4:
        return uvg_sad_4x4;
      case 8:
        return uvg_sad_8x8;
      case 16:
        return uvg_sad_16x16;
      case 32:
        return uvg_sad_32x32;
      case 64:
        return uvg_sad_64x64;
      default:
        return NULL;
    }
  }
  return NULL;
}

/**
* \brief  Get a function that calculates SATDs for 2 NxN blocks.
*
* \param width  Width of the region for which SATD is calculated.
* \param height  Height of the region for which SATD is calculated.
*
* \returns  Pointer to cost_pixel_nxn_multi_func.
*/
cost_pixel_nxn_multi_func * uvg_pixels_get_satd_dual_func(unsigned width, unsigned height)
{
  if(width == height) {
    switch (width) {
      case 4:
        return uvg_satd_4x4_dual;
      case 8:
        return uvg_satd_8x8_dual;
      case 16:
        return uvg_satd_16x16_dual;
      case 32:
        return uvg_satd_32x32_dual;
      case 64:
        return uvg_satd_64x64_dual;
      default:
        return NULL;
    }
  }
  return NULL;
}


/**
* \brief  Get a function that calculates SADs for 2 NxN blocks.
*
* \param width  Width of the region for which SAD is calculated.
*
* \returns  Pointer to cost_pixel_nxn_multi_func.
*/
cost_pixel_nxn_multi_func * uvg_pixels_get_sad_dual_func(unsigned width, unsigned height)
{
  if(width == height) {
    switch (width) {
      case 4:
        return uvg_sad_4x4_dual;
      case 8:
        return uvg_sad_8x8_dual;
      case 16:
        return uvg_sad_16x16_dual;
      case 32:
        return uvg_sad_32x32_dual;
      case 64:
        return uvg_sad_64x64_dual;
      default:
        return NULL;
    }
  }
  return NULL;
}

// Precomputed CRC32C lookup table for polynomial 0x04C11DB7
const uint32_t uvg_crc_table[256] = {
  0x00000000, 0xf26b8303, 0xe13b70f7, 0x1350f3f4, 0xc79a971f, 0x35f1141c,
  0x26a1e7e8, 0xd4ca64eb, 0x8ad958cf, 0x78b2dbcc, 0x6be22838, 0x9989ab3b,
  0x4d43cfd0, 0xbf284cd3, 0xac78bf27, 0x5e133c24, 0x105ec76f, 0xe235446c,
  0xf165b798, 0x030e349b, 0xd7c45070, 0x25afd373, 0x36ff2087, 0xc494a384,
  0x9a879fa0, 0x68ec1ca3, 0x7bbcef57, 0x89d76c54, 0x5d1d08bf, 0xaf768bbc,
  0xbc267848, 0x4e4dfb4b, 0x20bd8ede, 0xd2d60ddd, 0xc186fe29, 0x33ed7d2a,
  0xe72719c1, 0x154c9ac2, 0x061c6936, 0xf477ea35, 0xaa64d611, 0x580f5512,
  0x4b5fa6e6, 0xb93425e5, 0x6dfe410e, 0x9f95c20d, 0x8cc531f9, 0x7eaeb2fa,
  0x30e349b1, 0xc288cab2, 0xd1d83946, 0x23b3ba45, 0xf779deae, 0x05125dad,
  0x1642ae59, 0xe4292d5a, 0xba3a117e, 0x4851927d, 0x5b016189, 0xa96ae28a,
  0x7da08661, 0x8fcb0562, 0x9c9bf696, 0x6ef07595, 0x417b1dbc, 0xb3109ebf,
  0xa0406d4b, 0x522bee48, 0x86e18aa3, 0x748a09a0, 0x67dafa54, 0x95b17957,
  0xcba24573, 0x39c9c670, 0x2a993584, 0xd8f2b687, 0x0c38d26c, 0xfe53516f,
  0xed03a29b, 0x1f682198, 0x5125dad3, 0xa34e59d0, 0xb01eaa24, 0x42752927,
  0x96bf4dcc, 0x64d4cecf, 0x77843d3b, 0x85efbe38, 0xdbfc821c, 0x2997011f,
  0x3ac7f2eb, 0xc8ac71e8, 0x1c661503, 0xee0d9600, 0xfd5d65f4, 0x0f36e6f7,
  0x61c69362, 0x93ad1061, 0x80fde395, 0x72966096, 0xa65c047d, 0x5437877e,
  0x4767748a, 0xb50cf789, 0xeb1fcbad, 0x197448ae, 0x0a24bb5a, 0xf84f3859,
  0x2c855cb2, 0xdeeedfb1, 0xcdbe2c45, 0x3fd5af46, 0x7198540d, 0x83f3d70e,
  0x90a324fa, 0x62c8a7f9, 0xb602c312, 0x44694011, 0x5739b3e5, 0xa55230e6,
  0xfb410cc2, 0x092a8fc1, 0x1a7a7c35, 0xe811ff36, 0x3cdb9bdd, 0xceb018de,
  0xdde0eb2a, 0x2f8b6829, 0x82f63b78, 0x709db87b, 0x63cd4b8f, 0x91a6c88c,
  0x456cac67, 0xb7072f64, 0xa457dc90, 0x563c5f93, 0x082f63b7, 0xfa44e0b4,
  0xe9141340, 0x1b7f9043, 0xcfb5f4a8, 0x3dde77ab, 0x2e8e845f, 0xdce5075c,
  0x92a8fc17, 0x60c37f14, 0x73938ce0, 0x81f80fe3, 0x55326b08, 0xa759e80b,
  0xb4091bff, 0x466298fc, 0x1871a4d8, 0xea1a27db, 0xf94ad42f, 0x0b21572c,
  0xdfeb33c7, 0x2d80b0c4, 0x3ed04330, 0xccbbc033, 0xa24bb5a6, 0x502036a5,
  0x4370c551, 0xb11b4652, 0x65d122b9, 0x97baa1ba, 0x84ea524e, 0x7681d14d,
  0x2892ed69, 0xdaf96e6a, 0xc9a99d9e, 0x3bc21e9d, 0xef087a76, 0x1d63f975,
  0x0e330a81, 0xfc588982, 0xb21572c9, 0x407ef1ca, 0x532e023e, 0xa145813d,
  0x758fe5d6, 0x87e466d5, 0x94b49521, 0x66df1622, 0x38cc2a06, 0xcaa7a905,
  0xd9f75af1, 0x2b9cd9f2, 0xff56bd19, 0x0d3d3e1a, 0x1e6dcdee, 0xec064eed,
  0xc38d26c4, 0x31e6a5c7, 0x22b65633, 0xd0ddd530, 0x0417b1db, 0xf67c32d8,
  0xe52cc12c, 0x1747422f, 0x49547e0b, 0xbb3ffd08, 0xa86f0efc, 0x5a048dff,
  0x8ecee914, 0x7ca56a17, 0x6ff599e3, 0x9d9e1ae0, 0xd3d3e1ab, 0x21b862a8,
  0x32e8915c, 0xc083125f, 0x144976b4, 0xe622f5b7, 0xf5720643, 0x07198540,
  0x590ab964, 0xab613a67, 0xb831c993, 0x4a5a4a90, 0x9e902e7b, 0x6cfbad78,
  0x7fab5e8c, 0x8dc0dd8f, 0xe330a81a, 0x115b2b19, 0x020bd8ed, 0xf0605bee,
  0x24aa3f05, 0xd6c1bc06, 0xc5914ff2, 0x37faccf1, 0x69e9f0d5, 0x9b8273d6,
  0x88d28022, 0x7ab90321, 0xae7367ca, 0x5c18e4c9, 0x4f48173d, 0xbd23943e,
  0xf36e6f75, 0x0105ec76, 0x12551f82, 0xe03e9c81, 0x34f4f86a, 0xc69f7b69,
  0xd5cf889d, 0x27a40b9e, 0x79b737ba, 0x8bdcb4b9, 0x988c474d, 0x6ae7c44e,
  0xbe2da0a5, 0x4c4623a6, 0x5f16d052, 0xad7d5351,
};
