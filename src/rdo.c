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

#include "rdo.h"

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "cabac.h"
#include "context.h"
#include "encode_coding_tree.h"
#include "encoder.h"
#include "imagelist.h"
#include "inter.h"
#include "uvg_math.h"
#include "scalinglist.h"
#include "strategyselector.h"
#include "tables.h"
#include "transform.h"

#include "strategies/strategies-quant.h"


#define SCAN_SET_SIZE        16
#define LOG2_SCAN_SET_SIZE    4
#define SBH_THRESHOLD         4

#define RD_SAMPLING_MAX_LAST_QP     50

static FILE *fastrd_learning_outfile[RD_SAMPLING_MAX_LAST_QP + 1] = {NULL};
static pthread_mutex_t outfile_mutex[RD_SAMPLING_MAX_LAST_QP + 1];

const uint32_t uvg_g_go_rice_range[5] = { 7, 14, 26, 46, 78 };
const uint32_t uvg_g_go_rice_prefix_len[5] = { 8, 7, 6, 5, 4 };
static const uint32_t g_auiGoRiceParsCoeff[32] =
{
  0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3
};

/**
 * Entropy bits to estimate coded bits in RDO / RDOQ (From VTM 13.0)
 */
const uint32_t uvg_entropy_bits[2*256] = {
  0x0005c, 0x48000 , 0x00116, 0x3b520, 0x001d0, 0x356cb, 0x0028b, 0x318a9,
  0x00346, 0x2ea40 , 0x00403, 0x2c531, 0x004c0, 0x2a658, 0x0057e, 0x28beb,
  0x0063c, 0x274ce , 0x006fc, 0x26044, 0x007bc, 0x24dc9, 0x0087d, 0x23cfc,
  0x0093f, 0x22d96 , 0x00a01, 0x21f60, 0x00ac4, 0x2122e, 0x00b89, 0x205dd,
  0x00c4e, 0x1fa51 , 0x00d13, 0x1ef74, 0x00dda, 0x1e531, 0x00ea2, 0x1db78,
  0x00f6a, 0x1d23c , 0x01033, 0x1c970, 0x010fd, 0x1c10b, 0x011c8, 0x1b903,
  0x01294, 0x1b151 , 0x01360, 0x1a9ee, 0x0142e, 0x1a2d4, 0x014fc, 0x19bfc,
  0x015cc, 0x19564 , 0x0169c, 0x18f06, 0x0176d, 0x188de, 0x0183f, 0x182e8,
  0x01912, 0x17d23 , 0x019e6, 0x1778a, 0x01abb, 0x1721c, 0x01b91, 0x16cd5,
  0x01c68, 0x167b4 , 0x01d40, 0x162b6, 0x01e19, 0x15dda, 0x01ef3, 0x1591e,
  0x01fcd, 0x15480 , 0x020a9, 0x14fff, 0x02186, 0x14b99, 0x02264, 0x1474e,
  0x02343, 0x1431b , 0x02423, 0x13f01, 0x02504, 0x13afd, 0x025e6, 0x1370f,
  0x026ca, 0x13336 , 0x027ae, 0x12f71, 0x02894, 0x12bc0, 0x0297a, 0x12821,
  0x02a62, 0x12494 , 0x02b4b, 0x12118, 0x02c35, 0x11dac, 0x02d20, 0x11a51,
  0x02e0c, 0x11704 , 0x02efa, 0x113c7, 0x02fe9, 0x11098, 0x030d9, 0x10d77,
  0x031ca, 0x10a63 , 0x032bc, 0x1075c, 0x033b0, 0x10461, 0x034a5, 0x10173,
  0x0359b, 0x0fe90 , 0x03693, 0x0fbb9, 0x0378c, 0x0f8ed, 0x03886, 0x0f62b,
  0x03981, 0x0f374 , 0x03a7e, 0x0f0c7, 0x03b7c, 0x0ee23, 0x03c7c, 0x0eb89,
  0x03d7d, 0x0e8f9 , 0x03e7f, 0x0e671, 0x03f83, 0x0e3f2, 0x04088, 0x0e17c,
  0x0418e, 0x0df0e , 0x04297, 0x0dca8, 0x043a0, 0x0da4a, 0x044ab, 0x0d7f3,
  0x045b8, 0x0d5a5 , 0x046c6, 0x0d35d, 0x047d6, 0x0d11c, 0x048e7, 0x0cee3,
  0x049fa, 0x0ccb0 , 0x04b0e, 0x0ca84, 0x04c24, 0x0c85e, 0x04d3c, 0x0c63f,
  0x04e55, 0x0c426 , 0x04f71, 0x0c212, 0x0508d, 0x0c005, 0x051ac, 0x0bdfe,
  0x052cc, 0x0bbfc , 0x053ee, 0x0b9ff, 0x05512, 0x0b808, 0x05638, 0x0b617,
  0x0575f, 0x0b42a , 0x05888, 0x0b243, 0x059b4, 0x0b061, 0x05ae1, 0x0ae83,
  0x05c10, 0x0acaa , 0x05d41, 0x0aad6, 0x05e74, 0x0a907, 0x05fa9, 0x0a73c,
  0x060e0, 0x0a575 , 0x06219, 0x0a3b3, 0x06354, 0x0a1f5, 0x06491, 0x0a03b,
  0x065d1, 0x09e85 , 0x06712, 0x09cd4, 0x06856, 0x09b26, 0x0699c, 0x0997c,
  0x06ae4, 0x097d6 , 0x06c2f, 0x09634, 0x06d7c, 0x09495, 0x06ecb, 0x092fa,
  0x0701d, 0x09162 , 0x07171, 0x08fce, 0x072c7, 0x08e3e, 0x07421, 0x08cb0,
  0x0757c, 0x08b26 , 0x076da, 0x089a0, 0x0783b, 0x0881c, 0x0799f, 0x0869c,
  0x07b05, 0x0851f , 0x07c6e, 0x083a4, 0x07dd9, 0x0822d, 0x07f48, 0x080b9,
  0x080b9, 0x07f48 , 0x0822d, 0x07dd9, 0x083a4, 0x07c6e, 0x0851f, 0x07b05,
  0x0869c, 0x0799f , 0x0881c, 0x0783b, 0x089a0, 0x076da, 0x08b26, 0x0757c,
  0x08cb0, 0x07421 , 0x08e3e, 0x072c7, 0x08fce, 0x07171, 0x09162, 0x0701d,
  0x092fa, 0x06ecb , 0x09495, 0x06d7c, 0x09634, 0x06c2f, 0x097d6, 0x06ae4,
  0x0997c, 0x0699c , 0x09b26, 0x06856, 0x09cd4, 0x06712, 0x09e85, 0x065d1,
  0x0a03b, 0x06491 , 0x0a1f5, 0x06354, 0x0a3b3, 0x06219, 0x0a575, 0x060e0,
  0x0a73c, 0x05fa9 , 0x0a907, 0x05e74, 0x0aad6, 0x05d41, 0x0acaa, 0x05c10,
  0x0ae83, 0x05ae1 , 0x0b061, 0x059b4, 0x0b243, 0x05888, 0x0b42a, 0x0575f,
  0x0b617, 0x05638 , 0x0b808, 0x05512, 0x0b9ff, 0x053ee, 0x0bbfc, 0x052cc,
  0x0bdfe, 0x051ac , 0x0c005, 0x0508d, 0x0c212, 0x04f71, 0x0c426, 0x04e55,
  0x0c63f, 0x04d3c , 0x0c85e, 0x04c24, 0x0ca84, 0x04b0e, 0x0ccb0, 0x049fa,
  0x0cee3, 0x048e7 , 0x0d11c, 0x047d6, 0x0d35d, 0x046c6, 0x0d5a5, 0x045b8,
  0x0d7f3, 0x044ab , 0x0da4a, 0x043a0, 0x0dca8, 0x04297, 0x0df0e, 0x0418e,
  0x0e17c, 0x04088 , 0x0e3f2, 0x03f83, 0x0e671, 0x03e7f, 0x0e8f9, 0x03d7d,
  0x0eb89, 0x03c7c , 0x0ee23, 0x03b7c, 0x0f0c7, 0x03a7e, 0x0f374, 0x03981,
  0x0f62b, 0x03886 , 0x0f8ed, 0x0378c, 0x0fbb9, 0x03693, 0x0fe90, 0x0359b,
  0x10173, 0x034a5 , 0x10461, 0x033b0, 0x1075c, 0x032bc, 0x10a63, 0x031ca,
  0x10d77, 0x030d9 , 0x11098, 0x02fe9, 0x113c7, 0x02efa, 0x11704, 0x02e0c,
  0x11a51, 0x02d20 , 0x11dac, 0x02c35, 0x12118, 0x02b4b, 0x12494, 0x02a62,
  0x12821, 0x0297a , 0x12bc0, 0x02894, 0x12f71, 0x027ae, 0x13336, 0x026ca,
  0x1370f, 0x025e6 , 0x13afd, 0x02504, 0x13f01, 0x02423, 0x1431b, 0x02343,
  0x1474e, 0x02264 , 0x14b99, 0x02186, 0x14fff, 0x020a9, 0x15480, 0x01fcd,
  0x1591e, 0x01ef3 , 0x15dda, 0x01e19, 0x162b6, 0x01d40, 0x167b4, 0x01c68,
  0x16cd5, 0x01b91 , 0x1721c, 0x01abb, 0x1778a, 0x019e6, 0x17d23, 0x01912,
  0x182e8, 0x0183f , 0x188de, 0x0176d, 0x18f06, 0x0169c, 0x19564, 0x015cc,
  0x19bfc, 0x014fc , 0x1a2d4, 0x0142e, 0x1a9ee, 0x01360, 0x1b151, 0x01294,
  0x1b903, 0x011c8 , 0x1c10b, 0x010fd, 0x1c970, 0x01033, 0x1d23c, 0x00f6a,
  0x1db78, 0x00ea2 , 0x1e531, 0x00dda, 0x1ef74, 0x00d13, 0x1fa51, 0x00c4e,
  0x205dd, 0x00b89 , 0x2122e, 0x00ac4, 0x21f60, 0x00a01, 0x22d96, 0x0093f,
  0x23cfc, 0x0087d , 0x24dc9, 0x007bc, 0x26044, 0x006fc, 0x274ce, 0x0063c,
  0x28beb, 0x0057e , 0x2a658, 0x004c0, 0x2c531, 0x00403, 0x2ea40, 0x00346,
  0x318a9, 0x0028b , 0x356cb, 0x001d0, 0x3b520, 0x00116, 0x48000, 0x0005c,
};

// Entropy bits scaled so that 50% probability yields 1 bit.
const float uvg_f_entropy_bits[256*2] =
{
  0.002807617187500, 9.000000000000000, 0.008483886718750, 7.415039062500000, 0.014160156250000, 6.678070068359375, 0.019866943359375, 6.192657470703125,
  0.025573730468750, 5.830078125000000, 0.031341552734375, 5.540557861328125, 0.037109375000000, 5.299560546875000, 0.042907714843750, 5.093109130859375,
  0.048706054687500, 4.912536621093750, 0.054565429687500, 4.752075195312500, 0.060424804687500, 4.607696533203125, 0.066314697265625, 4.476440429687500,
  0.072235107421875, 4.356140136718750, 0.078155517578125, 4.245117187500000, 0.084106445312500, 4.142028808593750, 0.090118408203125, 4.045806884765625,
  0.096130371093750, 3.955596923828125, 0.102142333984375, 3.870727539062500, 0.108215332031250, 3.790557861328125, 0.114318847656250, 3.714599609375000,
  0.120422363281250, 3.642456054687500, 0.126556396484375, 3.573730468750000, 0.132720947265625, 3.508148193359375, 0.138916015625000, 3.445404052734375,
  0.145141601562500, 3.385284423828125, 0.151367187500000, 3.327575683593750, 0.157653808593750, 3.272094726562500, 0.163940429687500, 3.218627929687500,
  0.170288085937500, 3.167114257812500, 0.176635742187500, 3.117370605468750, 0.183013916015625, 3.069274902343750, 0.189422607421875, 3.022705078125000,
  0.195861816406250, 2.977630615234375, 0.202331542968750, 2.933898925781250, 0.208831787109375, 2.891479492187500, 0.215362548828125, 2.850250244140625,
  0.221923828125000, 2.810180664062500, 0.228515625000000, 2.771179199218750, 0.235137939453125, 2.733215332031250, 0.241790771484375, 2.696228027343750,
  0.248443603515625, 2.660156250000000, 0.255157470703125, 2.624969482421875, 0.261901855468750, 2.590606689453125, 0.268676757812500, 2.557067871093750,
  0.275482177734375, 2.524261474609375, 0.282318115234375, 2.492218017578125, 0.289184570312500, 2.460845947265625, 0.296081542968750, 2.430145263671875,
  0.303039550781250, 2.400085449218750, 0.309997558593750, 2.370635986328125, 0.317016601562500, 2.341796875000000, 0.324035644531250, 2.313507080078125,
  0.331115722656250, 2.285766601562500, 0.338226318359375, 2.258544921875000, 0.345367431640625, 2.231811523437500, 0.352539062500000, 2.205596923828125,
  0.359741210937500, 2.179809570312500, 0.367004394531250, 2.154510498046875, 0.374298095703125, 2.129638671875000, 0.381622314453125, 2.105194091796875,
  0.388977050781250, 2.081146240234375, 0.396362304687500, 2.057495117187500, 0.403808593750000, 2.034210205078125, 0.411285400390625, 2.011322021484375,
  0.418792724609375, 1.988769531250000, 0.426361083984375, 1.966583251953125, 0.433959960937500, 1.944732666015625, 0.441589355468750, 1.923187255859375,
  0.449249267578125, 1.901977539062500, 0.456970214843750, 1.881072998046875, 0.464721679687500, 1.860443115234375, 0.472534179687500, 1.840118408203125,
  0.480377197265625, 1.820098876953125, 0.488250732421875, 1.800323486328125, 0.496185302734375, 1.780822753906250, 0.504150390625000, 1.761596679687500,
  0.512145996093750, 1.742614746093750, 0.520233154296875, 1.723876953125000, 0.528320312500000, 1.705383300781250, 0.536468505859375, 1.687103271484375,
  0.544677734375000, 1.669097900390625, 0.552917480468750, 1.651275634765625, 0.561218261718750, 1.633666992187500, 0.569549560546875, 1.616302490234375,
  0.577941894531250, 1.599121093750000, 0.586364746093750, 1.582153320312500, 0.594848632812500, 1.565368652343750, 0.603393554687500, 1.548797607421875,
  0.611968994140625, 1.532409667968750, 0.620635986328125, 1.516174316406250, 0.629302978515625, 1.500152587890625, 0.638061523437500, 1.484313964843750,
  0.646850585937500, 1.468627929687500, 0.655700683593750, 1.453094482421875, 0.664611816406250, 1.437744140625000, 0.673583984375000, 1.422576904296875,
  0.682586669921875, 1.407531738281250, 0.691650390625000, 1.392669677734375, 0.700805664062500, 1.377960205078125, 0.709991455078125, 1.363372802734375,
  0.719238281250000, 1.348937988281250, 0.728546142578125, 1.334655761718750, 0.737915039062500, 1.320526123046875, 0.747344970703125, 1.306518554687500,
  0.756835937500000, 1.292633056640625, 0.766387939453125, 1.278900146484375, 0.776000976562500, 1.265289306640625, 0.785675048828125, 1.251800537109375,
  0.795440673828125, 1.238433837890625, 0.805236816406250, 1.225219726562500, 0.815124511718750, 1.212097167968750, 0.825073242187500, 1.199096679687500,
  0.835083007812500, 1.186218261718750, 0.845184326171875, 1.173461914062500, 0.855346679687500, 1.160797119140625, 0.865570068359375, 1.148254394531250,
  0.875885009765625, 1.135803222656250, 0.886260986328125, 1.123474121093750, 0.896697998046875, 1.111267089843750, 0.907257080078125, 1.099121093750000,
  0.917846679687500, 1.087097167968750, 0.928527832031250, 1.075195312500000, 0.939300537109375, 1.063354492187500, 0.950164794921875, 1.051635742187500,
  0.961090087890625, 1.040008544921875, 0.972106933593750, 1.028442382812500, 0.983184814453125, 1.016998291015625, 0.994384765625000, 1.005645751953125,
  1.005645751953125, 0.994384765625000, 1.016998291015625, 0.983184814453125, 1.028442382812500, 0.972106933593750, 1.040008544921875, 0.961090087890625,
  1.051635742187500, 0.950164794921875, 1.063354492187500, 0.939300537109375, 1.075195312500000, 0.928527832031250, 1.087097167968750, 0.917846679687500,
  1.099121093750000, 0.907257080078125, 1.111267089843750, 0.896697998046875, 1.123474121093750, 0.886260986328125, 1.135803222656250, 0.875885009765625,
  1.148254394531250, 0.865570068359375, 1.160797119140625, 0.855346679687500, 1.173461914062500, 0.845184326171875, 1.186218261718750, 0.835083007812500,
  1.199096679687500, 0.825073242187500, 1.212097167968750, 0.815124511718750, 1.225219726562500, 0.805236816406250, 1.238433837890625, 0.795440673828125,
  1.251800537109375, 0.785675048828125, 1.265289306640625, 0.776000976562500, 1.278900146484375, 0.766387939453125, 1.292633056640625, 0.756835937500000,
  1.306518554687500, 0.747344970703125, 1.320526123046875, 0.737915039062500, 1.334655761718750, 0.728546142578125, 1.348937988281250, 0.719238281250000,
  1.363372802734375, 0.709991455078125, 1.377960205078125, 0.700805664062500, 1.392669677734375, 0.691650390625000, 1.407531738281250, 0.682586669921875,
  1.422576904296875, 0.673583984375000, 1.437744140625000, 0.664611816406250, 1.453094482421875, 0.655700683593750, 1.468627929687500, 0.646850585937500,
  1.484313964843750, 0.638061523437500, 1.500152587890625, 0.629302978515625, 1.516174316406250, 0.620635986328125, 1.532409667968750, 0.611968994140625,
  1.548797607421875, 0.603393554687500, 1.565368652343750, 0.594848632812500, 1.582153320312500, 0.586364746093750, 1.599121093750000, 0.577941894531250,
  1.616302490234375, 0.569549560546875, 1.633666992187500, 0.561218261718750, 1.651275634765625, 0.552917480468750, 1.669097900390625, 0.544677734375000,
  1.687103271484375, 0.536468505859375, 1.705383300781250, 0.528320312500000, 1.723876953125000, 0.520233154296875, 1.742614746093750, 0.512145996093750,
  1.761596679687500, 0.504150390625000, 1.780822753906250, 0.496185302734375, 1.800323486328125, 0.488250732421875, 1.820098876953125, 0.480377197265625,
  1.840118408203125, 0.472534179687500, 1.860443115234375, 0.464721679687500, 1.881072998046875, 0.456970214843750, 1.901977539062500, 0.449249267578125,
  1.923187255859375, 0.441589355468750, 1.944732666015625, 0.433959960937500, 1.966583251953125, 0.426361083984375, 1.988769531250000, 0.418792724609375,
  2.011322021484375, 0.411285400390625, 2.034210205078125, 0.403808593750000, 2.057495117187500, 0.396362304687500, 2.081146240234375, 0.388977050781250,
  2.105194091796875, 0.381622314453125, 2.129638671875000, 0.374298095703125, 2.154510498046875, 0.367004394531250, 2.179809570312500, 0.359741210937500,
  2.205596923828125, 0.352539062500000, 2.231811523437500, 0.345367431640625, 2.258544921875000, 0.338226318359375, 2.285766601562500, 0.331115722656250,
  2.313507080078125, 0.324035644531250, 2.341796875000000, 0.317016601562500, 2.370635986328125, 0.309997558593750, 2.400085449218750, 0.303039550781250,
  2.430145263671875, 0.296081542968750, 2.460845947265625, 0.289184570312500, 2.492218017578125, 0.282318115234375, 2.524261474609375, 0.275482177734375,
  2.557067871093750, 0.268676757812500, 2.590606689453125, 0.261901855468750, 2.624969482421875, 0.255157470703125, 2.660156250000000, 0.248443603515625,
  2.696228027343750, 0.241790771484375, 2.733215332031250, 0.235137939453125, 2.771179199218750, 0.228515625000000, 2.810180664062500, 0.221923828125000,
  2.850250244140625, 0.215362548828125, 2.891479492187500, 0.208831787109375, 2.933898925781250, 0.202331542968750, 2.977630615234375, 0.195861816406250,
  3.022705078125000, 0.189422607421875, 3.069274902343750, 0.183013916015625, 3.117370605468750, 0.176635742187500, 3.167114257812500, 0.170288085937500,
  3.218627929687500, 0.163940429687500, 3.272094726562500, 0.157653808593750, 3.327575683593750, 0.151367187500000, 3.385284423828125, 0.145141601562500,
  3.445404052734375, 0.138916015625000, 3.508148193359375, 0.132720947265625, 3.573730468750000, 0.126556396484375, 3.642456054687500, 0.120422363281250,
  3.714599609375000, 0.114318847656250, 3.790557861328125, 0.108215332031250, 3.870727539062500, 0.102142333984375, 3.955596923828125, 0.096130371093750,
  4.045806884765625, 0.090118408203125, 4.142028808593750, 0.084106445312500, 4.245117187500000, 0.078155517578125, 4.356140136718750, 0.072235107421875,
  4.476440429687500, 0.066314697265625, 4.607696533203125, 0.060424804687500, 4.752075195312500, 0.054565429687500, 4.912536621093750, 0.048706054687500,
  5.093109130859375, 0.042907714843750, 5.299560546875000, 0.037109375000000, 5.540557861328125, 0.031341552734375, 5.830078125000000, 0.025573730468750,
  6.192657470703125, 0.019866943359375, 6.678070068359375, 0.014160156250000, 7.415039062500000, 0.008483886718750, 9.000000000000000, 0.002807617187500,

};


// This struct is for passing data to uvg_rdoq_sign_hiding
struct sh_rates_t {
  // Bit cost of increasing rate by one.
  int32_t inc[32 * 32];
  // Bit cost of decreasing rate by one.
  int32_t dec[32 * 32];
  // Bit cost of going from zero to one.
  int32_t sig_coeff_inc[32 * 32];
  // Coeff minus quantized coeff.
  int32_t quant_delta[32 * 32];
};

int uvg_init_rdcost_outfiles(const char *dir_path)
{
#define RD_SAMPLING_MAX_FN_LENGTH 4095
  static const char *basename_tmpl = "/%02i.txt";
  char fn_template[RD_SAMPLING_MAX_FN_LENGTH + 1];
  char fn[RD_SAMPLING_MAX_FN_LENGTH + 1];
  int rv = 0, qp;

  // As long as QP is a two-digit number, template and produced string should
  // be equal in length ("%i" -> "22")
  assert(RD_SAMPLING_MAX_LAST_QP <= 99);

  strncpy(fn_template, dir_path, RD_SAMPLING_MAX_FN_LENGTH);
  strncat(fn_template, basename_tmpl, RD_SAMPLING_MAX_FN_LENGTH - strlen(dir_path));
  assert(strlen(fn_template) <= RD_SAMPLING_MAX_FN_LENGTH);

  for (qp = 0; qp <= RD_SAMPLING_MAX_LAST_QP; qp++) {
    pthread_mutex_t *curr = outfile_mutex + qp;

    if (pthread_mutex_init(curr, NULL) != 0) {
      fprintf(stderr, "Failed to create mutex\n");
      rv = -1;
      qp--;
      goto out_destroy_mutexes;
    }
  }

  for (qp = 0; qp <= RD_SAMPLING_MAX_LAST_QP; qp++) {
    FILE *curr;

    snprintf(fn, RD_SAMPLING_MAX_FN_LENGTH, fn_template, qp);
    fn[RD_SAMPLING_MAX_FN_LENGTH] = 0;
    curr = fopen(fn, "w");
    if (curr == NULL) {
      fprintf(stderr, "Failed to open %s: %s\n", fn, strerror(errno));
      rv = -1;
      qp--;
      goto out_close_files;
    }
    fastrd_learning_outfile[qp] = curr;
  }
  goto out;

out_close_files:
  for (; qp >= 0; qp--) {
    fclose(fastrd_learning_outfile[qp]);
    fastrd_learning_outfile[qp] = NULL;
  }
  goto out;

out_destroy_mutexes:
  for (; qp >= 0; qp--) {
    pthread_mutex_destroy(outfile_mutex + qp);
  }
  goto out;

out:
  return rv;
#undef RD_SAMPLING_MAX_FN_LENGTH
}


/**
 * \brief Calculate actual (or really close to actual) bitcost for coding
 * coefficients.
 *
 * \param coeff coefficient array
 * \param width coeff block width
 * \param color data type (0 == luma)
 *
 * \returns bits needed to code input coefficients
 */
static INLINE double get_coeff_cabac_cost(
  const encoder_state_t * const state,
  const coeff_t *coeff,
  const cu_loc_t* const cu_loc,
  color_t color,
  int8_t scan_mode,
  int8_t tr_skip,
  cu_info_t* cur_tu)
{
  const int width  = cu_loc->width;
  const int height = cu_loc->height;
  const int sub_coeff_w = color == COLOR_Y ? cu_loc->width  : cu_loc->chroma_width;
  const int sub_coeff_h = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;

  // Make sure there are coeffs present
  bool found = false;
  for (int i = 0; i < sub_coeff_w * sub_coeff_h; i++) {
    if (coeff[i] != 0) {
      found = 1;
      break;
    }
  }
  if (!found) return 0;

  // Take a copy of the CABAC so that we don't overwrite the contexts when
  // counting the bits.
  cabac_data_t cabac_copy;
  memcpy(&cabac_copy, &state->search_cabac, sizeof(cabac_copy));

  // Clear bytes and bits and set mode to "count"
  cabac_copy.only_count = 1;
  cabac_copy.update = 1;
  double bits = 0;

  // Execute the coding function.
  // It is safe to drop the const modifier since state won't be modified
  // when cabac.only_count is set.
  if(!tr_skip) {
    uvg_encode_coeff_nxn((encoder_state_t*) state,
                         &cabac_copy,
                         coeff,
                         cu_loc,
                         color,
                         scan_mode,
                         cur_tu,                   
                         &bits);
  }
  else {
    uvg_encode_ts_residual((encoder_state_t* const)state,
      &cabac_copy,
      coeff,
      width,
      height,
      color,
      scan_mode,
      &bits);
  }
  if(state->search_cabac.update) {
    memcpy((cabac_data_t *)&state->search_cabac, &cabac_copy, sizeof(cabac_copy));
  }
  return bits;
}

static INLINE void save_ccc(int qp, const coeff_t *coeff, int32_t size, double ccc)
{
  pthread_mutex_t *mtx = outfile_mutex + qp;

  assert(sizeof(coeff_t) == sizeof(int16_t));
  assert(qp <= RD_SAMPLING_MAX_LAST_QP);

  pthread_mutex_lock(mtx);

  fwrite(&size,  sizeof(size),     1,    fastrd_learning_outfile[qp]);
  fwrite(&ccc,   sizeof(ccc),      1,    fastrd_learning_outfile[qp]);
  fwrite( coeff, sizeof(coeff_t),  size, fastrd_learning_outfile[qp]);

  pthread_mutex_unlock(mtx);
}

static INLINE void save_accuracy(int qp, double ccc, uint32_t fast_cost)
{
  pthread_mutex_t *mtx = outfile_mutex + qp;

  assert(qp <= RD_SAMPLING_MAX_LAST_QP);

  pthread_mutex_lock(mtx);
  fprintf(fastrd_learning_outfile[qp], "%u %f\n", fast_cost, ccc);
  pthread_mutex_unlock(mtx);
}

/**
 * \brief Estimate bitcost for coding coefficients.
 *
 * \param coeff   coefficient array
 * \param width   coeff block width
 * \param color    data type (0 == luma)
 *
 * \returns       number of bits needed to code coefficients
 */
double uvg_get_coeff_cost(
  const encoder_state_t * const state,
  const coeff_t *coeff,
  cu_info_t* cur_tu,
  const cu_loc_t* const cu_loc,
  color_t color,
  int8_t scan_mode,
  int8_t tr_skip,
  int coeff_order)
{
  uint8_t save_cccs = state->encoder_control->cfg.fastrd_sampling_on;
  uint8_t check_accuracy = state->encoder_control->cfg.fastrd_accuracy_check_on;

  const int width  = color == COLOR_Y ? cu_loc->width  : cu_loc->chroma_width;
  const int height = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;
  int x_local = cu_loc->x % LCU_WIDTH;
  int y_local = cu_loc->y % LCU_WIDTH;
  const int sub_coeff_w = color == COLOR_Y ? cu_loc->width : cu_loc->chroma_width;
  const int sub_coeff_h = color == COLOR_Y ? cu_loc->height : cu_loc->chroma_height;
  const int lcu_width = color == COLOR_Y ? LCU_WIDTH : LCU_WIDTH_C;


  const coeff_t* coeff_ptr = NULL;
  coeff_t sub_coeff[TR_MAX_WIDTH * TR_MAX_WIDTH];

  if (coeff_order == COEFF_ORDER_LINEAR) {
    coeff_ptr = coeff;
  }
  else {
    // Coeff order CU
    uvg_get_sub_coeff(sub_coeff, coeff, x_local, y_local, sub_coeff_w, sub_coeff_h, lcu_width);
    coeff_ptr = sub_coeff;
  }

  if (state->qp < state->encoder_control->cfg.fast_residual_cost_limit &&
      state->qp < MAX_FAST_COEFF_COST_QP && !tr_skip) {
    // TODO: do we need to assert(0) out of the fast-estimation branch if we
    // are to save block costs, or should we just warn about it somewhere
    // earlier (configuration validation I guess)?
    if (save_cccs) {
      assert(0 && "Fast RD sampling does not work with fast-residual-cost");
      return UINT32_MAX; // Hush little compiler don't you cry, not really gonna return anything after assert(0)
    } else {
      uint64_t weights = uvg_fast_coeff_get_weights(state);
      uint32_t fast_cost = uvg_fast_coeff_cost(coeff_ptr, width, height, weights);
      if (check_accuracy) {
        double ccc = get_coeff_cabac_cost(state, coeff_ptr, cu_loc, color, scan_mode, tr_skip, cur_tu);
        save_accuracy(state->qp, ccc, fast_cost);
      }
      return fast_cost;
    }
  } else {
    double ccc = get_coeff_cabac_cost(state, coeff_ptr, cu_loc, color, scan_mode, tr_skip, cur_tu);
    if (save_cccs) {
      save_ccc(state->qp, coeff, width * height, ccc);
    }
    return ccc;
  }
}

#define COEF_REMAIN_BIN_REDUCTION 5
/** Calculates the cost for specific absolute transform level
 * \param abs_level scaled quantized level
 * \param ctx_num_one current ctxInc for coeff_abs_level_greater1 (1st bin of coeff_abs_level_minus1 in AVC)
 * \param ctx_num_abs current ctxInc for coeff_abs_level_greater2 (remaining bins of coeff_abs_level_minus1 in AVC)
 * \param abs_go_rice Rice parameter for coeff_abs_level_minus3
 * \returns cost of given absolute transform level
 * From VTM 13.0
 */
INLINE int32_t uvg_get_ic_rate(encoder_state_t * const state,
                    uint32_t abs_level,
                    uint16_t ctx_num_gt1,
                    uint16_t ctx_num_gt2,
                    uint16_t ctx_num_par,
                    uint16_t abs_go_rice,
                    uint32_t reg_bins,
                    int8_t type,
                    int use_limited_prefix_length)
{
  cabac_data_t * const cabac = &state->cabac;
  int32_t rate = 1 << CTX_FRAC_BITS; // cost of sign bit
  uint32_t base_level  =  4;
  cabac_ctx_t *base_par_ctx = (type == 0) ? &(cabac->ctx.cu_parity_flag_model_luma[0]) : &(cabac->ctx.cu_parity_flag_model_chroma[0]);
  cabac_ctx_t *base_gt1_ctx = (type == 0) ? &(cabac->ctx.cu_gtx_flag_model_luma[1][0]) : &(cabac->ctx.cu_gtx_flag_model_chroma[1][0]);
  cabac_ctx_t* base_gt2_ctx = (type == 0) ? &(cabac->ctx.cu_gtx_flag_model_luma[0][0]) : &(cabac->ctx.cu_gtx_flag_model_chroma[0][0]);
  uint16_t go_rice_zero = 1 << abs_go_rice;
  int maxLog2TrDynamicRange = 15;

  if (reg_bins < 4)
  {
    uint32_t  symbol = (abs_level == 0 ? go_rice_zero : abs_level <= go_rice_zero ? abs_level - 1 : abs_level);
    uint32_t  length;
    const uint32_t threshold = COEF_REMAIN_BIN_REDUCTION;
    if (symbol < (threshold << abs_go_rice))
    {
      length = symbol >> abs_go_rice;
      rate += (length + 1 + abs_go_rice) << CTX_FRAC_BITS;
    } else if(use_limited_prefix_length) {
      const uint32_t maximumPrefixLength = (32 - (COEF_REMAIN_BIN_REDUCTION + maxLog2TrDynamicRange));

      uint32_t prefixLength = 0;
      uint32_t suffix = (symbol >> abs_go_rice) - COEF_REMAIN_BIN_REDUCTION;

      while ((prefixLength < maximumPrefixLength) && ((int32_t)suffix > ((2 << prefixLength) - 2)))
      {
        prefixLength++;
      }

      const uint32_t suffixLength = (prefixLength == maximumPrefixLength) ? (maxLog2TrDynamicRange - abs_go_rice) : (prefixLength + 1/*separator*/);

      rate += (COEF_REMAIN_BIN_REDUCTION + prefixLength + suffixLength + abs_go_rice) << CTX_FRAC_BITS;
    }
    else {
      length = abs_go_rice;
      symbol = symbol - (threshold << abs_go_rice);
      while ((int32_t)symbol >= (1 << length))
      {
        symbol -= (1 << (length++));
      }
      rate += (threshold + length + 1 - abs_go_rice + length) << CTX_FRAC_BITS;
    }
    return rate;
  }

  if ( abs_level >= base_level ) {
    int32_t symbol     = abs_level - base_level;
    int32_t length;
    if (symbol < (COEF_REMAIN_BIN_REDUCTION << abs_go_rice)) {
      length = symbol>>abs_go_rice;
      rate += (length + 1 + abs_go_rice) << CTX_FRAC_BITS;
    }
    else if (use_limited_prefix_length) {
      const uint32_t maximumPrefixLength = (32 - (COEF_REMAIN_BIN_REDUCTION + maxLog2TrDynamicRange));

      uint32_t prefixLength = 0;
      uint32_t suffix = (symbol >> abs_go_rice) - COEF_REMAIN_BIN_REDUCTION;

      while ((prefixLength < maximumPrefixLength) && ((int32_t)suffix > ((2 << prefixLength) - 2)))
      {
        prefixLength++;
      }

      const uint32_t suffixLength = (prefixLength == maximumPrefixLength) ? (maxLog2TrDynamicRange - abs_go_rice) : (prefixLength + 1/*separator*/);

      rate += (COEF_REMAIN_BIN_REDUCTION + prefixLength + suffixLength + abs_go_rice) << CTX_FRAC_BITS;
    }
    else {
      length = abs_go_rice;
      symbol  = symbol - ( COEF_REMAIN_BIN_REDUCTION << abs_go_rice);
      while (symbol >= (1<<length)) {
        symbol -=  (1<<(length++));
      }
      rate += (COEF_REMAIN_BIN_REDUCTION+length+1-abs_go_rice+length) << CTX_FRAC_BITS;
    }

    rate += CTX_ENTROPY_BITS(&base_par_ctx[ctx_num_par], (abs_level - 2) & 1);
    rate += CTX_ENTROPY_BITS(&base_gt1_ctx[ctx_num_gt1], 1);
    rate += CTX_ENTROPY_BITS(&base_gt2_ctx[ctx_num_gt2], 1);

  }
  else if (abs_level == 1)
  {    
    rate += CTX_ENTROPY_BITS(&base_gt1_ctx[ctx_num_gt1], 0);
  }
  else if (abs_level == 2)
  {
    rate += CTX_ENTROPY_BITS(&base_par_ctx[ctx_num_par], 0);
    rate += CTX_ENTROPY_BITS(&base_gt1_ctx[ctx_num_gt1], 1);
    rate += CTX_ENTROPY_BITS(&base_gt2_ctx[ctx_num_gt2], 0);
  }
  else if (abs_level == 3)
  {
    rate += CTX_ENTROPY_BITS(&base_par_ctx[ctx_num_par], 1);
    rate += CTX_ENTROPY_BITS(&base_gt1_ctx[ctx_num_gt1], 1);
    rate += CTX_ENTROPY_BITS(&base_gt2_ctx[ctx_num_gt2], 0);
  }
  else
  {
    rate = 0;
  }

  return rate;
}

/** Get the best level in RD sense
 * \param coded_cost reference to coded cost
 * \param coded_cost0 reference to cost when coefficient is 0
 * \param coded_cost_sig reference to cost of significant coefficient
 * \param level_double reference to unscaled quantized level
 * \param max_abs_level scaled quantized level
 * \param ctx_num_sig current ctxInc for coeff_abs_significant_flag
 * \param ctx_num_one current ctxInc for coeff_abs_level_greater1 (1st bin of coeff_abs_level_minus1 in AVC)
 * \param ctx_num_abs current ctxInc for coeff_abs_level_greater2 (remaining bins of coeff_abs_level_minus1 in AVC)
 * \param abs_go_rice current Rice parameter for coeff_abs_level_minus3
 * \param q_bits quantization step size
 * \param temp correction factor
 * \param last indicates if the coefficient is the last significant
 * \returns best quantized transform level for given scan position
 * This method calculates the best quantized transform level for a given scan position.
 * From VTM 13.0
 */
INLINE uint32_t uvg_get_coded_level( encoder_state_t * const state, double *coded_cost, double *coded_cost0, double *coded_cost_sig,
                           int32_t level_double, uint32_t max_abs_level,
                           uint16_t ctx_num_sig, uint16_t ctx_num_gt1, uint16_t ctx_num_gt2, uint16_t ctx_num_par,
                           uint16_t abs_go_rice,
                           uint32_t reg_bins,
                           int32_t q_bits,double error_scale, int8_t last, int8_t type)
{
  cabac_data_t * const cabac = &state->cabac;
  double cur_cost_sig   = 0;
  uint32_t best_abs_level = 0;
  int32_t abs_level;
  int32_t min_abs_level;
  cabac_ctx_t* base_sig_model = type?(cabac->ctx.cu_sig_model_chroma[0]):(cabac->ctx.cu_sig_model_luma[0]);
  const double lambda = type ? state->c_lambda : state->lambda;

  if( !last && max_abs_level < 3 ) {
    *coded_cost_sig = lambda * CTX_ENTROPY_BITS(&base_sig_model[ctx_num_sig], 0);
    *coded_cost     = *coded_cost0 + *coded_cost_sig;
    if (max_abs_level == 0) return best_abs_level;
  } else {
    *coded_cost = MAX_DOUBLE;
  }

  if( !last ) {
    cur_cost_sig = lambda * CTX_ENTROPY_BITS(&base_sig_model[ctx_num_sig], 1);
  }

  min_abs_level    = ( max_abs_level > 1 ? max_abs_level - 1 : 1 );
  for (abs_level = max_abs_level; abs_level >= min_abs_level ; abs_level-- ) {
    double err       = (double)(level_double - ( abs_level * (1 << q_bits) ) );
    double cur_cost  = err * err * error_scale + lambda *
                       uvg_get_ic_rate( state, abs_level, ctx_num_gt1, ctx_num_gt2, ctx_num_par,
                                    abs_go_rice, reg_bins, type, true);
    cur_cost        += cur_cost_sig;

    if( cur_cost < *coded_cost ) {
      best_abs_level  = abs_level;
      *coded_cost     = cur_cost;
      *coded_cost_sig = cur_cost_sig;
    }
  }

  return best_abs_level;
}


/** Calculates the cost of signaling the last significant coefficient in the block
 * \param pos_x X coordinate of the last significant coefficient
 * \param pos_y Y coordinate of the last significant coefficient
 * \returns cost of last significant coefficient
 * \param uiWidth width of the transform unit (TU)
 *
 * From VTM 13.0
*/
static double get_rate_last(double lambda,
                            const uint32_t  pos_x, const uint32_t pos_y,
                            int32_t* last_x_bits, int32_t* last_y_bits)
{
  uint32_t ctx_x   = g_group_idx[pos_x];
  uint32_t ctx_y   = g_group_idx[pos_y];
  double uiCost = last_x_bits[ ctx_x ] + last_y_bits[ ctx_y ];
  if( ctx_x > 3 ) {
    uiCost += CTX_FRAC_ONE_BIT * ((ctx_x - 2) >> 1);
  }
  if( ctx_y > 3 ) {
    uiCost += CTX_FRAC_ONE_BIT * ((ctx_y - 2) >> 1);
  }
  return lambda * uiCost;
}

static void calc_last_bits(encoder_state_t * const state, int32_t width, int32_t height, int8_t type,
                           int32_t* last_x_bits, int32_t* last_y_bits)
{
  cabac_data_t * const cabac = &state->cabac;
  int32_t bits_x = 0, bits_y = 0;
  int32_t blk_size_offset_x, blk_size_offset_y, shiftX, shiftY;
  int32_t ctx;

  cabac_ctx_t *base_ctx_x = (type ? cabac->ctx.cu_ctx_last_x_chroma : cabac->ctx.cu_ctx_last_x_luma);
  cabac_ctx_t *base_ctx_y = (type ? cabac->ctx.cu_ctx_last_y_chroma : cabac->ctx.cu_ctx_last_y_luma);

  static const int prefix_ctx[8] = { 0, 0, 0, 3, 6, 10, 15, 21 };
  blk_size_offset_x = type ? 0: prefix_ctx[uvg_math_floor_log2(width)];
  blk_size_offset_y = type ? 0: prefix_ctx[uvg_math_floor_log2(height)];
  shiftX = type ? CLIP(0, 2, width>>3) :((uvg_math_floor_log2(width) +1)>>2);
  shiftY = type ? CLIP(0, 2, height>>3) :((uvg_math_floor_log2(height) +1)>>2);


  for (ctx = 0; ctx < g_group_idx[ width - 1 ]; ctx++) {
    int32_t ctx_offset = blk_size_offset_x + (ctx >>shiftX);
    last_x_bits[ ctx ] = bits_x + CTX_ENTROPY_BITS(&base_ctx_x[ ctx_offset ],0);
    bits_x += CTX_ENTROPY_BITS(&base_ctx_x[ ctx_offset ],1);
  }
  last_x_bits[ctx] = bits_x;
  for (ctx = 0; ctx < g_group_idx[ height - 1 ]; ctx++) {
    int32_t ctx_offset = blk_size_offset_y + (ctx >>shiftY);
    last_y_bits[ ctx ] = bits_y + CTX_ENTROPY_BITS(&base_ctx_y[ ctx_offset ],0);
    bits_y +=  CTX_ENTROPY_BITS(&base_ctx_y[ ctx_offset ],1);
  }
  last_y_bits[ctx] = bits_y;
}

/**
 * \brief Select which coefficient to change for sign hiding, and change it.
 *
 * When sign hiding is enabled, the last sign bit of the last coefficient is
 * calculated from the parity of the other coefficients. If the parity is not
 * correct, one coefficient has to be changed by one. This function uses
 * tables generated during RDOQ to select the best coefficient to change.
 */
void uvg_rdoq_sign_hiding(
  const encoder_state_t *const state,
  const int32_t qp_scaled,
  const uint32_t *const scan2raster,
  const struct sh_rates_t *const sh_rates,
  const int32_t last_pos,
  const coeff_t *const coeffs,
  coeff_t *const quant_coeffs,
  const int8_t color,
  const bool need_sqrt_adjust)
{
  const encoder_control_t * const ctrl = state->encoder_control;
  const double lambda = color ? state->c_lambda : state->lambda;

  int inv_quant = uvg_g_inv_quant_scales[need_sqrt_adjust][qp_scaled % 6];
  // This somehow scales quant_delta into fractional bits. Instead of the bits
  // being multiplied by lambda, the residual is divided by it, or something
  // like that.
  const int64_t rd_factor = (const int64_t)(inv_quant * inv_quant * (1 << (2 * (qp_scaled / 6)))
                      / lambda / 16 / (1 << (2 * (ctrl->bitdepth - 8))) + 0.5);
  const int last_cg = (last_pos - 1) >> LOG2_SCAN_SET_SIZE;

  for (int32_t cg_scan = last_cg; cg_scan >= 0; cg_scan--) {
    const int32_t cg_coeff_scan = cg_scan << LOG2_SCAN_SET_SIZE;
    
    // Find positions of first and last non-zero coefficients in the CG.
    int32_t last_nz_scan = -1;
    for (int32_t coeff_i = SCAN_SET_SIZE - 1; coeff_i >= 0; --coeff_i) {
      if (quant_coeffs[scan2raster[coeff_i + cg_coeff_scan]]) {
        last_nz_scan = coeff_i;
        break;
      }
    }
    int32_t first_nz_scan = SCAN_SET_SIZE;
    for (int32_t coeff_i = 0; coeff_i <= last_nz_scan; coeff_i++) {
      if (quant_coeffs[scan2raster[coeff_i + cg_coeff_scan]]) {
        first_nz_scan = coeff_i;
        break;
      }
    }

    if (last_nz_scan - first_nz_scan < SBH_THRESHOLD) {
      continue;
    }

    const int32_t signbit = quant_coeffs[scan2raster[cg_coeff_scan + first_nz_scan]] <= 0;
    unsigned abs_coeff_sum = 0;
    for (int32_t coeff_scan = first_nz_scan; coeff_scan <= last_nz_scan; coeff_scan++) {
      abs_coeff_sum += quant_coeffs[scan2raster[coeff_scan + cg_coeff_scan]];
    }
    if (signbit == (abs_coeff_sum & 0x1)) {
      // Sign already matches with the parity, no need to modify coefficients.
      continue;
    }

    // Otherwise, search for the best coeff to change by one and change it.

    struct {
      int64_t cost;
      int pos;
      int change;
    } current, best = { MAX_INT64, 0, 0 };

    const int last_coeff_scan = (cg_scan == last_cg ? last_nz_scan : SCAN_SET_SIZE - 1);
    for (int coeff_scan = last_coeff_scan; coeff_scan >= 0; --coeff_scan) {
      current.pos = scan2raster[coeff_scan + cg_coeff_scan];
      // Shift the calculation back into original precision to avoid
      // changing the bitstream.
#     define PRECISION_INC (15 - CTX_FRAC_BITS)
      int64_t quant_cost_in_bits = rd_factor * sh_rates->quant_delta[current.pos];

      coeff_t abs_coeff = abs(quant_coeffs[current.pos]);

      if (abs_coeff != 0) {
        // Choose between incrementing and decrementing a non-zero coeff.

        int64_t inc_bits = sh_rates->inc[current.pos];
        int64_t dec_bits = sh_rates->dec[current.pos];
        if (abs_coeff == 1) {
          // We save sign bit and sig_coeff goes to zero.
          dec_bits -= sh_rates->sig_coeff_inc[current.pos];
        }
        if (cg_scan == last_cg && last_nz_scan == coeff_scan && abs_coeff == 1) {
          // Changing the last non-zero bit in the last cg to zero.
          // This might save a lot of bits if the next bits are already
          // zeros, or just a coupple fractional bits if they are not.
          // TODO: Check if calculating the real savings makes sense.
          dec_bits -= 4 * CTX_FRAC_ONE_BIT;
        }

        inc_bits = -quant_cost_in_bits + inc_bits * (1 << PRECISION_INC);
        dec_bits = quant_cost_in_bits + dec_bits * (1 << PRECISION_INC);

        if (inc_bits < dec_bits) {
          current.change = 1;
          current.cost = inc_bits;
        } else {
          current.change = -1;
          current.cost = dec_bits;

          if (coeff_scan == first_nz_scan && abs_coeff == 1) {
            // Don't turn first non-zero coeff into zero.
            // Seems kind of arbitrary. It's probably because it could lead to
            // breaking SBH_THRESHOLD.
            current.cost = MAX_INT64;
          }
        }
      } else {
        // Try incrementing a zero coeff.

        // Add sign bit, other bits and sig_coeff goes to one.
        int bits = CTX_FRAC_ONE_BIT + sh_rates->inc[current.pos] + sh_rates->sig_coeff_inc[current.pos];
        current.cost = -llabs(quant_cost_in_bits) + bits;
        current.change = 1;

        if (coeff_scan < first_nz_scan) {
          if (((coeffs[current.pos] >= 0) ? 0 : 1) != signbit) {
            current.cost = MAX_INT64;
          }
        }
      }

      if (current.cost < best.cost) {
        best = current;
      }
    }

    if (quant_coeffs[best.pos] == 32767 || quant_coeffs[best.pos] == -32768) {
      best.change = -1;
    }

    if (coeffs[best.pos] >= 0) {
      quant_coeffs[best.pos] += best.change;
    } else {
      quant_coeffs[best.pos] -= best.change;
    }
  }
}

static unsigned templateAbsSum(const coeff_t* coeff, int baseLevel, uint32_t  posX, uint32_t  posY, uint32_t width, uint32_t height, uint8_t mts_index)
{
  const coeff_t* pData = coeff + posX + posY * width;
  coeff_t          sum = 0;
  if (posX < width - 1)
  {
    sum += mts_index && posX + 1 >= 16 ? 0 : abs(pData[1]);
    if (posX < width - 2)
    {
      sum += mts_index && posX + 2 >= 16 ? 0 : abs(pData[2]);
    }
    if (posY < height - 1)
    {
      sum += mts_index && (posY + 1 >= 16 || posX + 1 >= 16) ? 0 : abs(pData[width + 1]);
    }
  }
  if (posY < height - 1)
  {
    sum += mts_index && posY + 1 >= 16 ? 0 : abs(pData[width]);
    if (posY < height - 2)
    {
      sum += mts_index && posY + 2 >= 16 ? 0 : abs(pData[width << 1]);
    }
  }
  return MAX(MIN(sum - 5 * baseLevel, 31), 0);
}

static INLINE int x_get_ic_rate_ts(const uint32_t            abs_level,
  const cabac_ctx_t* frac_bits_par,
  const cabac_ctx_t* frac_bits_sign,
  const cabac_ctx_t* frac_bits_gt1,
  const cabac_ctx_t* frac_bits_gtx_ctx,
  int* num_ctx_bins,
  const uint8_t             sign,
  const uint16_t            rice_par,
  const bool                use_limited_prefix_length,
  const int                 max_log2_tr_dynamic_range,
  int rem_reg_bins)
{

  if (rem_reg_bins < 4) // Full by-pass coding
  {
    int rate = abs_level ? (CTX_FRAC_ONE_BIT) : 0; // 1 bit to signal sign of non-zero

    uint32_t symbol = abs_level;

    uint32_t length;
    const int threshold = COEF_REMAIN_BIN_REDUCTION;
    if ((int32_t)symbol < (threshold << rice_par))
    {
      length = symbol >> rice_par;
      rate += (length + 1 + rice_par) << CTX_FRAC_BITS;
    }
    else if (use_limited_prefix_length)
    {
      const uint32_t maximumPrefixLength = (32 - (COEF_REMAIN_BIN_REDUCTION + max_log2_tr_dynamic_range));

      uint32_t prefixLength = 0;
      uint32_t suffix = (symbol >> rice_par) - COEF_REMAIN_BIN_REDUCTION;

      while ((prefixLength < maximumPrefixLength) && ((int32_t)suffix > ((2 << prefixLength) - 2)))
      {
        prefixLength++;
      }

      const uint32_t suffixLength = (prefixLength == maximumPrefixLength) ? (max_log2_tr_dynamic_range - rice_par) : (prefixLength + 1/*separator*/);

      rate += (COEF_REMAIN_BIN_REDUCTION + prefixLength + suffixLength + rice_par) << CTX_FRAC_BITS;
    }
    else
    {
      length = rice_par;
      symbol = symbol - (threshold << rice_par);
      while ((int32_t)symbol >= (1 << length))
      {
        symbol -= (1 << (length++));
      }
      rate += (threshold + length + 1 - rice_par + length) << CTX_FRAC_BITS;
    }

    return rate;
  }

  else if (rem_reg_bins >= 4 && rem_reg_bins < 8) // First pass context coding and all by-pass coding ( Sign flag is not counted here)
  {
    int rate = CTX_ENTROPY_BITS(frac_bits_sign, sign); // frac_bits_sign.intBits[sign]; // sign bits
    if (abs_level)
      (*num_ctx_bins)++;

    if (abs_level > 1)
    {
      rate += CTX_ENTROPY_BITS(frac_bits_gt1, 1); // frac_bits_gt1.intBits[1];
      rate += CTX_ENTROPY_BITS(frac_bits_par, (abs_level - 2) & 1); // frac_bits_par.intBits[(abs_level - 2) & 1];

      (*num_ctx_bins) += 2;

      uint32_t cutoffVal = 2;

      if (abs_level >= cutoffVal)
      {
        uint32_t symbol = (abs_level - cutoffVal) >> 1;
        uint32_t length;
        const int threshold = COEF_REMAIN_BIN_REDUCTION;
        if ((int32_t)symbol < (threshold << rice_par))
        {
          length = symbol >> rice_par;
          rate += (length + 1 + rice_par) << CTX_FRAC_BITS;
        }
        else if (use_limited_prefix_length)
        {
          const uint32_t maximumPrefixLength = (32 - (COEF_REMAIN_BIN_REDUCTION + max_log2_tr_dynamic_range));

          uint32_t prefixLength = 0;
          uint32_t suffix = (symbol >> rice_par) - COEF_REMAIN_BIN_REDUCTION;

          while ((prefixLength < maximumPrefixLength) && ((int32_t)suffix > ((2 << prefixLength) - 2)))
          {
            prefixLength++;
          }

          const uint32_t suffixLength = (prefixLength == maximumPrefixLength) ? (max_log2_tr_dynamic_range - rice_par) : (prefixLength + 1/*separator*/);

          rate += (COEF_REMAIN_BIN_REDUCTION + prefixLength + suffixLength + rice_par) << CTX_FRAC_BITS;
        }
        else
        {
          length = rice_par;
          symbol = symbol - (threshold << rice_par);
          while ((int32_t)symbol >= (1 << length))
          {
            symbol -= (1 << (length++));
          }
          rate += (threshold + length + 1 - rice_par + length) << CTX_FRAC_BITS;
        }
      }
    }
    else if (abs_level == 1)
    {
      rate += CTX_ENTROPY_BITS(frac_bits_gt1, 0);  // frac_bits_gt1.intBits[0];
      num_ctx_bins++;
    }
    else
    {
      rate = 0;
    }
    return rate;

  }
  int rate = CTX_ENTROPY_BITS(frac_bits_sign, sign);

  if (abs_level)
    num_ctx_bins++;

  if (abs_level > 1)
  {
    rate += CTX_ENTROPY_BITS(frac_bits_gt1, 1); // frac_bits_gt1.intBits[1];
    rate += CTX_ENTROPY_BITS(frac_bits_sign, (abs_level - 2) & 1); // frac_bits_par.intBits[(abs_level - 2) & 1];
    num_ctx_bins += 2;

    uint32_t cutoffVal = 2;
    const int numGtBins = 4;
    for (int i = 0; i < numGtBins; i++)
    {
      if (abs_level >= cutoffVal)
      {
        const uint16_t ctxGtX = cutoffVal >> 1;
        // const BinFracBits* fracBitsGtX = fracBitsAccess.getFracBitsArray(ctxGtX);
        unsigned gtX = ((int32_t)abs_level >= (cutoffVal + 2));
        rate += CTX_ENTROPY_BITS(&frac_bits_gtx_ctx[ctxGtX], gtX);// fracBitsGtX.intBits[gtX];
        num_ctx_bins++;
      }
      cutoffVal += 2;
    }

    if (abs_level >= cutoffVal)
    {
      uint32_t symbol = (abs_level - cutoffVal) >> 1;
      uint32_t length;
      const uint32_t threshold = COEF_REMAIN_BIN_REDUCTION;
      if (symbol < (threshold << rice_par))
      {
        length = symbol >> rice_par;
        rate += (length + 1 + rice_par) << CTX_FRAC_BITS;
      }
      else if (use_limited_prefix_length)
      {
        const uint32_t maximumPrefixLength = (32 - (COEF_REMAIN_BIN_REDUCTION + max_log2_tr_dynamic_range));

        uint32_t prefixLength = 0;
        uint32_t suffix = (symbol >> rice_par) - COEF_REMAIN_BIN_REDUCTION;

        while ((prefixLength < maximumPrefixLength) && ((int32_t)suffix > ((2 << prefixLength) - 2)))
        {
          prefixLength++;
        }

        const uint32_t suffixLength = (prefixLength == maximumPrefixLength) ? (max_log2_tr_dynamic_range - rice_par) : (prefixLength + 1/*separator*/);

        rate += (COEF_REMAIN_BIN_REDUCTION + prefixLength + suffixLength + rice_par) << CTX_FRAC_BITS;
      }
      else
      {
        length = rice_par;
        symbol = symbol - (threshold << rice_par);
        while ((int32_t)symbol >= (1 << length))
        {
          symbol -= (1 << (length++));
        }
        rate += (threshold + length + 1 - rice_par + length) << CTX_FRAC_BITS;
      }
    }
  }
  else if (abs_level == 1)
  {
    rate += CTX_ENTROPY_BITS(frac_bits_gt1, 0); // frac_bits_gt1.intBits[0];
    num_ctx_bins++;
  }
  else
  {
    rate = 0;
  }
  return rate;

}

static inline uint32_t get_coded_level_ts_pred(double* coded_cost,
  double* coded_cost0,
  double* coded_cost_sig,
  int    level_double,
  int                 q_bits,
  double              error_scale,
  uint32_t* coeff_levels,
  double* coeff_level_error,
  const cabac_ctx_t* frac_bits_sig,
  const cabac_ctx_t* frac_bits_par,
  const cabac_ctx_t* frac_bits_sign,
  const cabac_ctx_t* frac_bits_gt1,
  const cabac_ctx_t* frac_bits_gtx_ctx,
  const uint8_t      sign,
  int                right_pixel,
  int                below_pixel,
  uint16_t           rice_par,
  bool               is_last,
  bool               use_limited_prefix_length,
  const int          max_log2_tr_dynamic_range,
  int* num_used_ctx_bins,
  int rem_reg_bins,
  int tested_levels,
  double lambda
)
{
  double curr_cost_sig = 0;
  uint32_t   best_abs_level = 0;
  *num_used_ctx_bins = 0;
  int num_best_ctx_bin = 0;
  int bdpcm = 0;
  if (!is_last && coeff_levels[0] < 3)
  {
    if (rem_reg_bins >= 4)
      *coded_cost_sig = lambda * CTX_ENTROPY_BITS(frac_bits_sig, 0);
    else
      *coded_cost_sig = lambda * (1 << CTX_FRAC_BITS);
    *coded_cost = *coded_cost0 + *coded_cost_sig;
    if (rem_reg_bins >= 4)
      (*num_used_ctx_bins)++;
    if (coeff_levels[0] == 0)
    {
      return best_abs_level;
    }
  }
  else
  {
    *coded_cost = MAX_DOUBLE;
  }

  if (!is_last)
  {
    if (rem_reg_bins >= 4)
      curr_cost_sig = lambda * CTX_ENTROPY_BITS(frac_bits_sig, 1);
    else
      curr_cost_sig = lambda * (1 << CTX_FRAC_BITS);
    if (coeff_levels[0] >= 3 && rem_reg_bins >= 4)
      (*num_used_ctx_bins)++;
  }

  for (int errorInd = 1; errorInd <= tested_levels; errorInd++)
  {
    int absLevel = coeff_levels[errorInd - 1];
    double dErr = 0.0;
    dErr = (double)(level_double - ((absLevel) << q_bits));
    coeff_level_error[errorInd] = dErr * dErr * error_scale;
    int modAbsLevel = absLevel;
    if (rem_reg_bins >= 4)
    {
      modAbsLevel = uvg_derive_mod_coeff(right_pixel, below_pixel, absLevel, bdpcm);
    }
    int numCtxBins = 0;
    double dCurrCost = coeff_level_error[errorInd] + lambda * 
      x_get_ic_rate_ts(modAbsLevel, frac_bits_par,  frac_bits_sign, frac_bits_gt1, frac_bits_gtx_ctx,
        &numCtxBins, sign, rice_par, use_limited_prefix_length, max_log2_tr_dynamic_range, rem_reg_bins);

    if (rem_reg_bins >= 4)
      dCurrCost += curr_cost_sig; // if cctx.numCtxBins < 4, xGetICRateTS return rate including sign cost. dont need to add any more

    if (dCurrCost < *coded_cost)
    {
      best_abs_level = absLevel;
      *coded_cost = dCurrCost;
      *coded_cost_sig = curr_cost_sig;
      num_best_ctx_bin = numCtxBins;
    }
  }
  *num_used_ctx_bins += num_best_ctx_bin;
  return best_abs_level;
}


int uvg_ts_rdoq(encoder_state_t* const state, coeff_t* src_coeff, coeff_t* dest_coeff, int32_t width,
                int32_t height, int8_t type, int8_t scan_mode) {
  const encoder_control_t* const encoder = state->encoder_control;
  const cabac_data_t* cabac = &state->cabac;
  

  const bool extended_precision = false;
  const int  max_log2_tr_dynamic_range = 15;
  uint32_t log2_tr_width = uvg_math_floor_log2(width);
  uint32_t log2_tr_height = uvg_math_floor_log2(height);
  const uint32_t log2_block_width  = uvg_g_convert_to_log2[width];
  const uint32_t log2_block_height = uvg_g_convert_to_log2[height];
  const uint32_t log2_cg_width = g_log2_sbb_size[log2_tr_width][log2_tr_height][0];
  const uint32_t log2_cg_height = g_log2_sbb_size[log2_tr_width][log2_tr_height][1];

  const uint32_t log2_cg_size = log2_cg_width + log2_cg_height;
  //TODO: Scaling list

  double   block_uncoded_cost = 0;
  uint32_t cg_num = width * height >> log2_cg_size;

  int32_t qp_scaled = uvg_get_scaled_qp(type, state->qp, (encoder->bitdepth - 8) * 6, encoder->qp_map[0]);
  qp_scaled = MAX(qp_scaled, 4 + 6 * MIN_QP_PRIME_TS);
  int32_t max_num_coeff = width * height;

  // TODO: Scaling list
  
  double cost_coeff[32 * 32];
  double cost_sig[32 * 32];
  double cost_coeff0[32 * 32];

  double   cost_coeffgroup_sig[64];
  uint32_t sig_coeffgroup_flag[64];

  switch (cg_num) {
  case  1: FILL_ARRAY(sig_coeffgroup_flag, 0, 1); FILL_ARRAY(cost_coeffgroup_sig, 0, 1); break;
  case  2: FILL_ARRAY(sig_coeffgroup_flag, 0, 2); FILL_ARRAY(cost_coeffgroup_sig, 0, 2); break;
  case  4: FILL_ARRAY(sig_coeffgroup_flag, 0, 4); FILL_ARRAY(cost_coeffgroup_sig, 0, 4);  break;
  case  8: FILL_ARRAY(sig_coeffgroup_flag, 0, 8); FILL_ARRAY(cost_coeffgroup_sig, 0, 8);  break;
  case 16: FILL_ARRAY(sig_coeffgroup_flag, 0, 16); FILL_ARRAY(cost_coeffgroup_sig, 0, 16);  break;
  case 32: FILL_ARRAY(sig_coeffgroup_flag, 0, 32); FILL_ARRAY(cost_coeffgroup_sig, 0, 32);  break;
  case 64: FILL_ARRAY(sig_coeffgroup_flag, 0, 64); FILL_ARRAY(cost_coeffgroup_sig, 0, 64); break;
  default: assert(0 && "There should be 1, 4, 16 or 64 coefficient groups");
  }

  const bool   needs_sqrt2_scale = false; // from VTM: should always be false - transform-skipped blocks don't require sqrt(2) compensation.
  const int    q_bits = QUANT_SHIFT + qp_scaled / 6  + (needs_sqrt2_scale ? -1 : 0);  // Right shift of non-RDOQ quantizer;  level = (coeff*uiQ + offset)>>q_bits
  const int32_t quant_coeff = uvg_g_quant_scales[needs_sqrt2_scale][qp_scaled % 6];
 
  const double error_scale = (double)(1 << CTX_FRAC_BITS) / quant_coeff / quant_coeff;

  double lambda = type == 0 ? state->lambda : state->c_lambda;

  const coeff_t entropy_coding_maximum = (1 << max_log2_tr_dynamic_range) - 1;

  const uint32_t* const scan = uvg_get_scan_order_table(SCAN_GROUP_4X4, scan_mode, log2_block_width, log2_block_height, 0);
  const uint32_t* const scan_cg = uvg_get_scan_order_table(SCAN_GROUP_UNGROUPED, scan_mode, log2_block_width, log2_block_height, 0);

  uint32_t coeff_levels[3];
  double   coeff_level_error[4];

  const int sbSizeM1 = (1 << log2_cg_size) - 1;
  double    base_cost = 0;
  uint32_t  go_rice_par = 0;
    
  int scan_pos;
  struct {
    double coded_level_and_dist;
    double uncoded_dist;
    double sig_cost;
    double sig_cost_0;
    int32_t nnz_before_pos0;
    int32_t num_sbb_ctx_bins;
  } rd_stats;

  bool any_sig_cg = false;

  int rem_reg_bins = (width * height * 7) >> 2;

  for (uint32_t sbId = 0; sbId < cg_num; sbId++)
  {
    uint32_t cg_blkpos = scan_cg[sbId];

    int no_coeff_coded = 0;
    base_cost = 0.0;
    FILL(rd_stats, 0);

    rd_stats.num_sbb_ctx_bins = 0;

    for (int scan_pos_in_sb = 0; scan_pos_in_sb <= sbSizeM1; scan_pos_in_sb++)
    {
      scan_pos = (sbId << log2_cg_size) + scan_pos_in_sb;
      int last_pos_coded = sbSizeM1;
      uint32_t blkpos = scan[scan_pos];
      uint32_t  pos_y = blkpos >> log2_block_width;
      uint32_t  pos_x = blkpos - (pos_y << log2_block_width); 
      //===== quantization =====

      // set coeff
      const int64_t          tmp_level = (int64_t)(abs(src_coeff[blkpos])) * quant_coeff;
      const int level_double = (const int)MIN(tmp_level, MAX_INT64 - (1ll << ((long long)q_bits - 1ll)));

      uint32_t roundAbsLevel = MIN((uint32_t)(entropy_coding_maximum), (uint32_t)((level_double + ((1) << (q_bits - 1))) >> q_bits));
      uint32_t min_abs_level = (roundAbsLevel > 1 ? roundAbsLevel - 1 : 1);

      uint32_t down_abs_level = MIN((uint32_t)(entropy_coding_maximum), (uint32_t)(level_double >> q_bits));
      uint32_t up_abs_level = MIN((uint32_t)(entropy_coding_maximum), down_abs_level + 1);

      int tested_levels = 0;
      coeff_levels[tested_levels++] = roundAbsLevel;

      if (min_abs_level != roundAbsLevel)
        coeff_levels[tested_levels++] = min_abs_level;

      int right_pixel, below_pixel, pred_pixel;

      right_pixel = pos_x > 0 ? src_coeff[pos_x + pos_y * width - 1] : 0;
      below_pixel = pos_y > 0 ? src_coeff[pos_x + (pos_y - 1) * width] : 0;

      pred_pixel = uvg_derive_mod_coeff(right_pixel, below_pixel, up_abs_level, 0);

      if (up_abs_level != roundAbsLevel && up_abs_level != min_abs_level && pred_pixel == 1)
        coeff_levels[tested_levels++] = up_abs_level;

      double err = (double)(level_double);
      coeff_level_error[0] = err * err * error_scale;

      cost_coeff0[scan_pos] = coeff_level_error[0];
      block_uncoded_cost += cost_coeff0[scan_pos];
      dest_coeff[blkpos] = coeff_levels[0];

      //===== coefficient level estimation =====

      unsigned    ctx_id_sig = uvg_context_get_sig_ctx_idx_abs_ts(dest_coeff, pos_x, pos_y, width);
      uint32_t    c_level;
      const cabac_ctx_t* frac_bits_par = &cabac->ctx.transform_skip_par;

      go_rice_par = 1;
      unsigned ctx_id_sign = uvg_sign_ctx_id_abs_ts(dest_coeff, pos_x, pos_y, width, 0);
      const cabac_ctx_t* frac_bits_sign = &cabac->ctx.transform_skip_res_sign[ctx_id_sign];
      const uint8_t     sign = src_coeff[blkpos] < 0 ? 1 : 0;
      

      unsigned gt1_ctx_id = uvg_lrg1_ctx_id_abs_ts(dest_coeff, pos_x, pos_y, width, 0);
      const cabac_ctx_t* frac_bits_gt1 = &cabac->ctx.transform_skip_gt1[gt1_ctx_id];

      const cabac_ctx_t* frac_bits_sig = &cabac->ctx.transform_skip_sig[ctx_id_sig]; 
      bool is_last = false; //
      if (scan_pos_in_sb == last_pos_coded && no_coeff_coded == 0)
      {
        is_last = true;
      }
      int num_used_ctx_bins = 0;
      c_level = get_coded_level_ts_pred(&cost_coeff[scan_pos], &cost_coeff0[scan_pos], &cost_sig[scan_pos], level_double,
        q_bits, error_scale, coeff_levels, coeff_level_error,
        frac_bits_sig, frac_bits_par, frac_bits_sign, frac_bits_gt1, cabac->ctx.transform_skip_gt2,
        sign, right_pixel, below_pixel, go_rice_par, is_last, extended_precision,
        max_log2_tr_dynamic_range, &num_used_ctx_bins, rem_reg_bins, tested_levels, lambda);

      rem_reg_bins -= num_used_ctx_bins;
      rd_stats.num_sbb_ctx_bins += num_used_ctx_bins;


      if (c_level > 0)
      {
        no_coeff_coded++;
      }

      coeff_t level = c_level;
      dest_coeff[blkpos] = (level != 0 && src_coeff[blkpos] < 0) ? -level : level;
      base_cost += cost_coeff[scan_pos];
      rd_stats.sig_cost += cost_sig[scan_pos];

      if (dest_coeff[blkpos])
      {
        sig_coeffgroup_flag[cg_blkpos] = 1;
        rd_stats.coded_level_and_dist += cost_coeff[scan_pos] - cost_sig[scan_pos];
        rd_stats.uncoded_dist += cost_coeff0[scan_pos];
      }
    } //end for (iScanPosinCG)

    const cabac_ctx_t* fracBitsSigGroup = &cabac->ctx.sig_coeff_group_model[(type == 0 ? 0 : 1) * 2 + 1];
    if (sig_coeffgroup_flag[cg_blkpos])
    {
      base_cost += lambda*CTX_ENTROPY_BITS(fracBitsSigGroup, 0) - rd_stats.sig_cost;
      cost_coeffgroup_sig[sbId] = lambda * CTX_ENTROPY_BITS(fracBitsSigGroup, 0);

      rem_reg_bins += rd_stats.num_sbb_ctx_bins; // skip sub-block
    }
    else if (sbId != cg_num - 1 || any_sig_cg)
    {
      // rd-cost if SigCoeffGroupFlag = 0, initialization
      double cost_zero_sb = base_cost;
      
      base_cost += lambda * CTX_ENTROPY_BITS(fracBitsSigGroup, 1);
      cost_zero_sb += lambda * CTX_ENTROPY_BITS(fracBitsSigGroup, 0);
      cost_coeffgroup_sig[sbId] = lambda * CTX_ENTROPY_BITS(fracBitsSigGroup, 1);

      cost_zero_sb += rd_stats.uncoded_dist;         // distortion for resetting non-zero levels to zero levels
      cost_zero_sb -= rd_stats.coded_level_and_dist;   // distortion and level cost for keeping all non-zero levels
      cost_zero_sb -= rd_stats.sig_cost;             // sig cost for all coeffs, including zero levels and non-zerl levels

      if (cost_zero_sb < base_cost)
      {
        base_cost = cost_zero_sb;
        cost_coeffgroup_sig[sbId] = lambda * CTX_ENTROPY_BITS(fracBitsSigGroup, 0);
        rem_reg_bins += rd_stats.num_sbb_ctx_bins; // skip sub-block
        for (int scanPosInSB = 0; scanPosInSB <= sbSizeM1; scanPosInSB++)
        {
          scan_pos = (sbId << log2_cg_size) + scanPosInSB;
          uint32_t blkPos = scan[scan_pos];

          if (dest_coeff[blkPos])
          {
            dest_coeff[blkPos] = 0;
            cost_coeff[scan_pos] = cost_coeff0[scan_pos];
            cost_sig[scan_pos] = 0;
          }
        }
      }
      else
      {
        any_sig_cg = true;
      }
    }
  }

  int abs_sum = 0;
  //===== estimate last position =====
  for (int scanPos = 0; scanPos < max_num_coeff; scanPos++)
  {
    int blkPos = scan[scanPos];
    coeff_t level = dest_coeff[blkPos];
    abs_sum += abs(level);
  }
  return abs_sum;
}


static uint32_t context_get_sig_ctx_idx_abs(const coeff_t* coeff, uint32_t pos_x, uint32_t pos_y,
                                            uint32_t width, uint32_t height, int8_t color,
                                            int32_t* temp_diag, int32_t* temp_sum, int8_t mts)
{
  const coeff_t* data = coeff + pos_x + pos_y * width;
  const int     diag = pos_x + pos_y;
  int           num_pos = 0;
  int           sum_abs = 0;
#define UPDATE(x) {int a=abs(x);sum_abs+=MIN(4+(a&1),a);num_pos+=(a?1:0);}
  if (pos_x < width - 1)
  {
    UPDATE(mts && pos_x + 1 >= 16 ? 0 : data[1]);
    if (pos_x < width - 2)
    {
      UPDATE(mts && pos_x + 2 >= 16 ? 0 : data[2]);
    }
    if (pos_y < height - 1)
    {
      UPDATE(mts && (pos_y + 1 >= 16 || pos_x + 1 >= 16) ? 0 : data[width + 1]);
    }
  }
  if (pos_y < height - 1)
  {
    UPDATE(mts && pos_x + 1 >= 16 ? 0 : data[width]);
    if (pos_y < height - 2)
    {
      UPDATE(mts && pos_x + 2 >= 16 ? 0 : data[width << 1]);
    }
  }
#undef UPDATE
  int ctx_ofs = MIN((sum_abs + 1) >> 1, 3) + (diag < 2 ? 4 : 0);
  if (color == COLOR_Y)
  {
    ctx_ofs += diag < 5 ? 4 : 0;
  }

  *temp_diag = diag;
  *temp_sum = sum_abs - num_pos;
  return ctx_ofs;
}

/** RDOQ with CABAC
 * \returns void
 * Rate distortion optimized quantization for entropy
 * coding engines using probability models like CABAC
 * From VTM 13.0
 */
void uvg_rdoq(
  encoder_state_t * const state,
  coeff_t *coef,
  coeff_t *dest_coeff,
  int32_t width,
  int32_t height,
  int8_t color,
  int8_t scan_mode,
  int8_t block_type,
  uint16_t cbf,
  uint8_t lfnst_idx, uint8_t mts_idx)
{
  const encoder_control_t * const encoder = state->encoder_control;
  cabac_data_t * const cabac = &state->cabac;
  const uint32_t log2_block_width = uvg_g_convert_to_log2[width];
  const uint32_t log2_block_height = uvg_g_convert_to_log2[height];
  bool needs_block_size_trafo_scale = !false && ((log2_block_width + log2_block_height) % 2 == 1);
  needs_block_size_trafo_scale |= 0; // Non log2 block size

  int32_t  transform_shift   = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - ((log2_block_width + log2_block_height) >> 1);  // Represents scaling through forward transform
  uint16_t go_rice_param     = 0;
  uint32_t reg_bins = (width * height * 28) >> 4;
  
  int32_t  scalinglist_type= (block_type == CU_INTRA ? 0 : 3) + color;

  int32_t qp_scaled = uvg_get_scaled_qp(color, state->qp, (encoder->bitdepth - 8) * 6, encoder->qp_map[0]);
  
  int32_t q_bits = QUANT_SHIFT + qp_scaled/6 + transform_shift - needs_block_size_trafo_scale;

  const double lambda = color ? state->c_lambda : state->lambda;
  const int32_t default_quant_coeff = uvg_g_quant_scales[needs_block_size_trafo_scale][qp_scaled % 6];
  const bool use_scaling_list = state->encoder_control->cfg.scaling_list != UVG_SCALING_LIST_OFF;

  const int32_t *quant_coeff  = encoder->scaling_list.quant_coeff[log2_block_width][log2_block_height][scalinglist_type][qp_scaled%6];
  const double *err_scale     = encoder->scaling_list.error_scale[log2_block_width][log2_block_height][scalinglist_type][qp_scaled%6];

  double block_uncoded_cost = 0;
  
  double cost_coeff [ 32 * 32 ];
  double cost_sig   [ 32 * 32 ];
  double cost_coeff0[ 32 * 32 ];

  struct sh_rates_t sh_rates;

  FILL(sh_rates, 0);

  memset(dest_coeff, 0, sizeof(coeff_t) * width * height);

  const uint32_t log2_cg_size = uvg_g_log2_sbb_size[log2_block_width][log2_block_height][0] + uvg_g_log2_sbb_size[log2_block_width][log2_block_height][1];
  const uint32_t log2_cg_width  = uvg_g_log2_sbb_size[log2_block_width][log2_block_height][0];
  const uint32_t log2_cg_height = uvg_g_log2_sbb_size[log2_block_width][log2_block_height][1];

  const uint32_t cg_width  = (MIN((uint8_t)TR_MAX_WIDTH, width) >> log2_cg_width);
  const uint32_t cg_height = (MIN((uint8_t)TR_MAX_WIDTH, height) >> log2_cg_height);

  const uint32_t * const scan = uvg_get_scan_order_table(SCAN_GROUP_4X4, scan_mode, log2_block_width, log2_block_height, 0);
  const uint32_t * const scan_cg = uvg_get_scan_order_table(SCAN_GROUP_UNGROUPED, scan_mode, log2_block_width, log2_block_height, 0);

  const uint32_t cg_size = 16;
  const int32_t  shift = 4 >> 1;
  const uint32_t num_blk_side = MAX(width >> shift, 1);
  double   cost_coeffgroup_sig[ 64 ];
  uint32_t sig_coeffgroup_flag[ 64 ];

  uint16_t    ctx_set    = 0;
  double      base_cost  = 0;
  int32_t temp_diag = -1;
  int32_t temp_sum = -1;

  int32_t cg_last_scanpos = -1;
  int32_t last_scanpos = -1;

  uint32_t       cg_num          = lfnst_idx > 0 ? 1 : width * height >> 4;

  double         dTransShift = (double)transform_shift + (needs_block_size_trafo_scale ? -0.5 : 0.0);
  // Compensate for scaling of bitcount in Lagrange cost function
  double scale       = CTX_FRAC_ONE_BIT;
  // Compensate for scaling through forward transform
  scale              = scale * pow(2.0, -2.0 * dTransShift);
  const double  default_error_scale = scale / default_quant_coeff / default_quant_coeff;

  // Explicitly tell the only possible numbers of elements to be zeroed.
  // Hope the compiler is able to utilize this information.
  switch (cg_num) {
    case  1: FILL_ARRAY(sig_coeffgroup_flag, 0,  1); break;
    case  2: FILL_ARRAY(sig_coeffgroup_flag, 0,  2); break;
    case  4: FILL_ARRAY(sig_coeffgroup_flag, 0,  4); break;
    case  8: FILL_ARRAY(sig_coeffgroup_flag, 0,  8); break;
    case 16: FILL_ARRAY(sig_coeffgroup_flag, 0, 16); break;
    case 32: FILL_ARRAY(sig_coeffgroup_flag, 0, 32); break;
    case 64: FILL_ARRAY(sig_coeffgroup_flag, 0, 64); break;
    default: assert(0 && "There should be 1, 2, 4, 8, 16, 32 or 64 coefficient groups");
  }

  cabac_ctx_t *base_coeff_group_ctx = &(cabac->ctx.sig_coeff_group_model[color ? 2 : 0]);
  cabac_ctx_t *baseCtx              = (color == 0) ? &(cabac->ctx.cu_sig_model_luma[0][0]) : &(cabac->ctx.cu_sig_model_chroma[0][0]);
  cabac_ctx_t* base_gt1_ctx = (color == 0) ? &(cabac->ctx.cu_gtx_flag_model_luma[1][0]) : &(cabac->ctx.cu_gtx_flag_model_chroma[1][0]);

  struct {
    double coded_level_and_dist;
    double uncoded_dist;
    double sig_cost;
    double sig_cost_0;
    int32_t nnz_before_pos0;
  } rd_stats;

  //Find last cg and last scanpos
  const int max_lfnst_pos = ((height == 4 && width == 4) || (height == 8 && width == 8)) ? 7 : 15;
  int32_t   cg_scanpos;
  uint32_t  max_scan_group_size = lfnst_idx > 0 ? max_lfnst_pos : cg_size - 1;
  for (cg_scanpos = (cg_num - 1); cg_scanpos >= 0; cg_scanpos--)
  {
    uint32_t cg_blkpos = scan_cg[cg_scanpos];
    uint32_t cg_pos_y = cg_blkpos / num_blk_side;
    uint32_t cg_pos_x = cg_blkpos - (cg_pos_y * num_blk_side);
    if (mts_idx != 0 && (cg_pos_y >= 4 || cg_pos_x >= 4)) continue;
    for (int32_t scanpos_in_cg = max_scan_group_size; scanpos_in_cg >= 0; scanpos_in_cg--)
    {
      int32_t  scanpos        = cg_scanpos*cg_size + scanpos_in_cg;
      
      uint32_t blkpos         = scan[scanpos];
      int32_t q               = use_scaling_list ? quant_coeff[blkpos] : default_quant_coeff;
      int32_t level_double    = coef[blkpos];
      level_double            = MIN(abs(level_double) * q, MAX_INT - (1 << (q_bits - 1)));
      uint32_t max_abs_level  = (level_double + (1 << (q_bits - 1))) >> q_bits;

      double err = (double)level_double;

      cost_coeff0[scanpos] = err * err * (use_scaling_list ? err_scale[blkpos] : default_error_scale);      
      
      dest_coeff[blkpos] = max_abs_level;
      if (max_abs_level > 0) {
        last_scanpos    = scanpos;        
        cg_last_scanpos = cg_scanpos;
        sh_rates.sig_coeff_inc[blkpos] = 0;
        break;
      }
      block_uncoded_cost += cost_coeff0[scanpos];
      base_cost += cost_coeff0[scanpos];
    }
    if (last_scanpos != -1) break;
  }

  if (last_scanpos == -1) {
    return;
  }


  for (; cg_scanpos >= 0; cg_scanpos--) cost_coeffgroup_sig[cg_scanpos] = 0;

  int32_t last_x_bits[32], last_y_bits[32];

  for (int32_t cg_scanpos = cg_last_scanpos; cg_scanpos >= 0; cg_scanpos--) {
    uint32_t cg_blkpos  = scan_cg[cg_scanpos];
    uint32_t cg_pos_y   = cg_blkpos / num_blk_side;
    uint32_t cg_pos_x   = cg_blkpos - (cg_pos_y * num_blk_side);

    FILL(rd_stats, 0);
    if (mts_idx != 0 && (cg_pos_y >= 4 || cg_pos_x >= 4)) continue;
    for (int32_t scanpos_in_cg = max_scan_group_size; scanpos_in_cg >= 0; scanpos_in_cg--)  {
      int32_t  scanpos = cg_scanpos*cg_size + scanpos_in_cg;
      if (scanpos > last_scanpos) {
        continue;
      }
      uint32_t blkpos         = scan[scanpos];
      int32_t q               = use_scaling_list ? quant_coeff[blkpos] : default_quant_coeff;
      double temp             = (use_scaling_list ? err_scale[blkpos] : default_error_scale);
      int32_t level_double    = coef[blkpos];
      level_double            = MIN(abs(level_double) * q , MAX_INT - (1 << (q_bits - 1)));
      uint32_t max_abs_level  = (level_double + (1 << (q_bits - 1))) >> q_bits;
      dest_coeff[blkpos] = max_abs_level;
      double err = (double)level_double;

      cost_coeff0[scanpos] = err * err * (use_scaling_list ? err_scale[blkpos] : default_error_scale);

      block_uncoded_cost      += cost_coeff0[ scanpos ];

      if (last_scanpos >= 0) {

        uint32_t  pos_y = blkpos >> log2_block_width;
        uint32_t  pos_x = blkpos - (pos_y << log2_block_width);
        //===== coefficient level estimation =====
        int32_t  level;
        
        uint16_t ctx_sig = 0;
        if (scanpos != last_scanpos) {
          // VVC document 9.3.4.2.8, context for sig_coeff_flag calculated here
          ctx_sig = context_get_sig_ctx_idx_abs(dest_coeff, pos_x, pos_y, width, height, color, &temp_diag, &temp_sum, mts_idx);
        }
        
        if (temp_diag != -1) {
          ctx_set = (MIN(temp_sum, 4) + 1) + (!temp_diag ? ((color == 0) ? 15 : 5) : (color == 0) ? temp_diag < 3 ? 10 : (temp_diag < 10 ? 5 : 0) : 0);
        }
        else ctx_set = 0;

        if (reg_bins < 4) {
          int  sumAll = templateAbsSum(dest_coeff, 0, pos_x, pos_y, width, height,mts_idx);
          go_rice_param = g_auiGoRiceParsCoeff[sumAll];
        }

        uint16_t  gt1_ctx = ctx_set;
        uint16_t  gt2_ctx = ctx_set;
        uint16_t  par_ctx = ctx_set;

        if (scanpos == last_scanpos) {
          level = uvg_get_coded_level(state, &cost_coeff[scanpos], &cost_coeff0[scanpos], &cost_sig[scanpos],
            level_double, max_abs_level, 0, gt1_ctx, gt2_ctx, par_ctx, go_rice_param,
            reg_bins, q_bits, temp, 1, color);          
        }
        else {
          level = uvg_get_coded_level(state, &cost_coeff[scanpos], &cost_coeff0[scanpos], &cost_sig[scanpos],
            level_double, max_abs_level, ctx_sig, gt1_ctx, gt2_ctx, par_ctx, go_rice_param,
            reg_bins, q_bits, temp, 0, color);
          if (encoder->cfg.signhide_enable) {
            int greater_than_zero = CTX_ENTROPY_BITS(&baseCtx[ctx_sig], 1);
            int zero = CTX_ENTROPY_BITS(&baseCtx[ctx_sig], 0);
            sh_rates.sig_coeff_inc[blkpos] = (reg_bins < 4 ? 0 : greater_than_zero - zero);
          }
        }



        if (encoder->cfg.signhide_enable) {
          sh_rates.quant_delta[blkpos] = (level_double - level * (1 << q_bits)) >> (q_bits - 8);
          if (level > 0) {
            int32_t rate_now = uvg_get_ic_rate(state, level, gt1_ctx, gt2_ctx, par_ctx, go_rice_param, reg_bins, color, false);
            sh_rates.inc[blkpos] = uvg_get_ic_rate(state, level + 1, gt1_ctx, gt2_ctx, par_ctx, go_rice_param, reg_bins, color, false) - rate_now;
            sh_rates.dec[blkpos] = uvg_get_ic_rate(state, level - 1, gt1_ctx, gt2_ctx, par_ctx, go_rice_param, reg_bins, color, false) - rate_now;
          }
          else { // level == 0
            if (reg_bins < 4) {
              int32_t rate_now = uvg_get_ic_rate(state, level, gt1_ctx, gt2_ctx, par_ctx, go_rice_param, reg_bins, color, false);
              sh_rates.inc[blkpos] = uvg_get_ic_rate(state, level + 1, gt1_ctx, gt2_ctx, par_ctx, go_rice_param, reg_bins, color, false) - rate_now;
            }
            else {
              sh_rates.inc[blkpos] = CTX_ENTROPY_BITS(&base_gt1_ctx[gt1_ctx], 0);
            }
          }
        }
        dest_coeff[blkpos] = (coeff_t)level;
        base_cost += cost_coeff[scanpos];

        //===== context set update =====
        if ((scanpos % SCAN_SET_SIZE == 0) && scanpos > 0) {
          go_rice_param = 0;
        }
        else if (reg_bins >= 4) {
          reg_bins -= (level < 2 ? level : 3) + (scanpos != last_scanpos);
          int  sumAll = templateAbsSum(coef, 4, pos_x, pos_y, width, height, mts_idx);
          go_rice_param = g_auiGoRiceParsCoeff[sumAll];
        }
      }
      else {
        base_cost += cost_coeff0[scanpos];
      }

      rd_stats.sig_cost += cost_sig[scanpos];
      if ( scanpos_in_cg == 0 ) {
        rd_stats.sig_cost_0 = cost_sig[scanpos];
      }
      if ( dest_coeff[blkpos] )  {
        sig_coeffgroup_flag[cg_blkpos] = 1;
        rd_stats.coded_level_and_dist   += cost_coeff[scanpos] - cost_sig[scanpos];
        rd_stats.uncoded_dist           += cost_coeff0[scanpos];
        if ( scanpos_in_cg != 0 ) {
          rd_stats.nnz_before_pos0++;
        }
      }
    } //end for (scanpos_in_cg)

    if( cg_scanpos ) {
      if (sig_coeffgroup_flag[cg_blkpos] == 0) {
        uint32_t ctx_sig  = uvg_context_get_sig_coeff_group(sig_coeffgroup_flag, cg_pos_x,
                                                        cg_pos_y, cg_width, cg_height);
        cost_coeffgroup_sig[cg_scanpos] = lambda *CTX_ENTROPY_BITS(&base_coeff_group_ctx[ctx_sig],0);
        base_cost += cost_coeffgroup_sig[cg_scanpos]  - rd_stats.sig_cost;
      } else {
        if (cg_scanpos < cg_last_scanpos){
          double cost_zero_cg;
          uint32_t ctx_sig;
          if (rd_stats.nnz_before_pos0 == 0) {
            base_cost -= rd_stats.sig_cost_0;
            rd_stats.sig_cost -= rd_stats.sig_cost_0;
          }
          // rd-cost if SigCoeffGroupFlag = 0, initialization
          cost_zero_cg = base_cost;

          // add SigCoeffGroupFlag cost to total cost
          ctx_sig = uvg_context_get_sig_coeff_group(sig_coeffgroup_flag, cg_pos_x,
            cg_pos_y, cg_width, cg_height);

          cost_coeffgroup_sig[cg_scanpos] = lambda * CTX_ENTROPY_BITS(&base_coeff_group_ctx[ctx_sig], 1);
          base_cost += cost_coeffgroup_sig[cg_scanpos];
          cost_zero_cg += lambda * CTX_ENTROPY_BITS(&base_coeff_group_ctx[ctx_sig], 0);

          // try to convert the current coeff group from non-zero to all-zero
          cost_zero_cg += rd_stats.uncoded_dist;          // distortion for resetting non-zero levels to zero levels
          cost_zero_cg -= rd_stats.coded_level_and_dist;  // distortion and level cost for keeping all non-zero levels
          cost_zero_cg -= rd_stats.sig_cost;              // sig cost for all coeffs, including zero levels and non-zerl levels

          // if we can save cost, change this block to all-zero block
          if (cost_zero_cg < base_cost) {

            sig_coeffgroup_flag[cg_blkpos] = 0;
            base_cost = cost_zero_cg;

            cost_coeffgroup_sig[cg_scanpos] = lambda * CTX_ENTROPY_BITS(&base_coeff_group_ctx[ctx_sig], 0);

            // reset coeffs to 0 in this block
            for (int32_t scanpos_in_cg = max_scan_group_size; scanpos_in_cg >= 0; scanpos_in_cg--) {
              int32_t  scanpos = cg_scanpos*cg_size + scanpos_in_cg;
              uint32_t blkpos = scan[scanpos];
              if (dest_coeff[blkpos]){
                dest_coeff[blkpos] = 0;
                cost_coeff[scanpos] = cost_coeff0[scanpos];
                cost_sig[scanpos] = 0;
              }
            }
          } // end if ( cost_all_zeros < base_cost )
        }
      } // end if if (sig_coeffgroup_flag[ cg_blkpos ] == 0)
    } else {
      sig_coeffgroup_flag[cg_blkpos] = 1;
    }
  } //end for (cg_scanpos)

  //===== estimate last position =====
  double  best_cost        = 0;
  int32_t ctx_cbf          = 0;
  int8_t found_last        = 0;
  int32_t best_last_idx_p1 = 0;

  if( block_type != CU_INTRA && !color ) {
    best_cost  = block_uncoded_cost +  lambda * CTX_ENTROPY_BITS(&(cabac->ctx.cu_qt_root_cbf_model),0);
    base_cost +=   lambda * CTX_ENTROPY_BITS(&(cabac->ctx.cu_qt_root_cbf_model),1);
  } else {
    cabac_ctx_t* base_cbf_model = NULL;
    switch (color) {
      case COLOR_Y:
        base_cbf_model = cabac->ctx.qt_cbf_model_luma;
        break;
      case COLOR_U:
        base_cbf_model = cabac->ctx.qt_cbf_model_cb;
        break;
      case COLOR_V:
        base_cbf_model = cabac->ctx.qt_cbf_model_cr;
        break;
      default:
        assert(0);
    }
    // This cbf should work even with non-square blocks
    ctx_cbf    = ( color != COLOR_V ? 0 : cbf_is_set(cbf, COLOR_U));
    best_cost  = block_uncoded_cost +  lambda * CTX_ENTROPY_BITS(&base_cbf_model[ctx_cbf],0);
    base_cost +=   lambda * CTX_ENTROPY_BITS(&base_cbf_model[ctx_cbf],1);
  }

  calc_last_bits(state, width, height, color, last_x_bits, last_y_bits);
  for ( int32_t cg_scanpos = cg_last_scanpos; cg_scanpos >= 0; cg_scanpos--) {
    uint32_t cg_blkpos = scan_cg[cg_scanpos];
    base_cost -= cost_coeffgroup_sig[cg_scanpos];

    if (sig_coeffgroup_flag[ cg_blkpos ]) {
      for ( int32_t scanpos_in_cg = max_scan_group_size; scanpos_in_cg >= 0; scanpos_in_cg--) {
        int32_t  scanpos = cg_scanpos*cg_size + scanpos_in_cg;
        if (scanpos > last_scanpos) continue;
        uint32_t blkpos  = scan[scanpos];

        if( dest_coeff[ blkpos ] ) {
          uint32_t   pos_y = blkpos >> log2_block_width;
          uint32_t   pos_x = blkpos - ( pos_y << log2_block_width );

          double cost_last = get_rate_last(lambda, pos_x, pos_y, last_x_bits,last_y_bits );
          double totalCost = base_cost + cost_last - cost_sig[ scanpos ];

          if( totalCost < best_cost ) {
            best_last_idx_p1 = scanpos + 1;
            best_cost        = totalCost;
          }
          if( dest_coeff[ blkpos ] > 1 ) {
            found_last = 1;
            break;
          }
          base_cost -= cost_coeff[scanpos];
          base_cost += cost_coeff0[scanpos];
        } else {
          base_cost -= cost_sig[scanpos];
        }
      } //end for
      if (found_last) break;
    } // end if (sig_coeffgroup_flag[ cg_blkpos ])
  } // end for

  uint32_t abs_sum = 0;
  if(!mts_idx || (width < 32 && height < 32)) {
    for ( int32_t scanpos = 0; scanpos < best_last_idx_p1; scanpos++) {
      int32_t blkPos     = scan[scanpos];
      int32_t level      = dest_coeff[blkPos];
      abs_sum            += level;
      dest_coeff[blkPos] = (coeff_t)(( coef[blkPos] < 0 ) ? -level : level);
    }
  }
  else {
    for ( int32_t scanpos = 0; scanpos < best_last_idx_p1; scanpos++) {
      int32_t blkPos     = scan[scanpos];
      int32_t blk_x = blkPos & (width - 1);
      int32_t blk_y = blkPos >> log2_block_width;
      int32_t level      = blk_x >= 16 || blk_y >= 16 ? 0 : dest_coeff[blkPos];
      abs_sum            += level;
      dest_coeff[blkPos] = (coeff_t)((level != 0 && coef[blkPos] < 0) ? -level : level);
    }
  }
  //===== clean uncoded coefficients =====
  for ( int32_t scanpos = best_last_idx_p1; scanpos <= last_scanpos; scanpos++) {
    dest_coeff[scan[scanpos]] = 0;
  }

  if (encoder->cfg.signhide_enable && abs_sum >= 2) {
    uvg_rdoq_sign_hiding(state, qp_scaled, scan, &sh_rates, best_last_idx_p1, coef, dest_coeff, color, needs_block_size_trafo_scale);
  }
}


/**
 * Calculate cost of actual motion vectors using CABAC coding
 */
double uvg_get_mvd_coding_cost_cabac(const encoder_state_t* state,
                                     const cabac_data_t* cabac,
                                     const int32_t mvd_hor,
                                     const int32_t mvd_ver)
{
  cabac_data_t cabac_copy = *cabac;
  cabac_copy.only_count = 1;
  double bits = 0;
  // It is safe to drop const here because cabac->only_count is set.
  uvg_encode_mvd((encoder_state_t*) state, &cabac_copy, mvd_hor, mvd_ver, &bits);

  return bits;
}


/** MVD cost calculation with CABAC
* \returns int
* Calculates Motion Vector cost and related costs using CABAC coding
*/
double uvg_calc_ibc_mvd_cost_cabac(const encoder_state_t * state,
                               int x,
                               int y,
                               int mv_shift,
                               mv_t mv_cand[2][2],
                               inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                               int16_t num_cand,
                               int32_t ref_idx,
                               double* bitcost)
{
  cabac_data_t state_cabac_copy;
  cabac_data_t* cabac;
  uint32_t merge_idx;
  vector2d_t mvd = { 0, 0 };
  int8_t merged = 0;
  int8_t cur_mv_cand = 0;

  x *= 1 << mv_shift;
  y *= 1 << mv_shift;

  // Check every candidate to find a match
  for (merge_idx = 0; merge_idx < (uint32_t)num_cand; merge_idx++) {
    if (merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][0] == x &&
      merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][1] == y)
    {
      merged = 1;
      break;
    }
  }

  // Store cabac state and contexts
  memcpy(&state_cabac_copy, &state->search_cabac, sizeof(cabac_data_t));

  // Clear bytes and bits and set mode to "count"
  state_cabac_copy.only_count = 1;

  cabac = &state_cabac_copy;
  double bits = 0;

  if (!merged) {
    vector2d_t mvd1 = {
      x - mv_cand[0][0],
      y - mv_cand[0][1],
    };
    vector2d_t mvd2 = {
      x - mv_cand[1][0],
      y - mv_cand[1][1],
    };

    uvg_change_precision_vector2d(INTERNAL_MV_PREC, 2, &mvd1);
    uvg_change_precision_vector2d(INTERNAL_MV_PREC, 2, &mvd2);

    double cand1_cost = uvg_get_mvd_coding_cost_cabac(state, cabac, mvd1.x, mvd1.y);
    double cand2_cost = uvg_get_mvd_coding_cost_cabac(state, cabac, mvd2.x, mvd2.y);

    // Select candidate 1 if it has lower cost
    if (cand2_cost < cand1_cost) {
      cur_mv_cand = 1;
      mvd = mvd2;
    } else {
      mvd = mvd1;
    }
  }

  cabac->cur_ctx = &(cabac->ctx.cu_merge_flag_ext_model);

  CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_merge_flag_ext_model), merged, bits, "MergeFlag");
  num_cand = state->encoder_control->cfg.max_merge;
  if (merged) {
    if (num_cand > 1) {
      int32_t ui;
      for (ui = 0; ui < num_cand - 1; ui++) {
        int32_t symbol = (ui != merge_idx);
        if (ui == 0) {
          CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_merge_idx_ext_model), symbol, bits, "MergeIndex");
        } else {
          CABAC_BIN_EP(cabac, symbol, "MergeIndex");
          bits += 1;
        }
        if (symbol == 0) break;
      }
    }
  } else {

    // It is safe to drop const here because cabac->only_count is set.
    uvg_encode_mvd((encoder_state_t*) state, cabac, mvd.x, mvd.y, &bits);

    // Signal which candidate MV to use
    cabac->cur_ctx = &(cabac->ctx.mvp_idx_model);
    CABAC_BIN(cabac, cur_mv_cand, "mvp_flag");
  }

  *bitcost = bits;

  // Store bitcost before restoring cabac
  return *bitcost * state->lambda_sqrt;
}

/** MVD cost calculation with CABAC
* \returns int
* Calculates Motion Vector cost and related costs using CABAC coding
*/
double uvg_calc_mvd_cost_cabac(const encoder_state_t * state,
                               int x,
                               int y,
                               int mv_shift,
                               mv_t mv_cand[2][2],
                               inter_merge_cand_t merge_cand[MRG_MAX_NUM_CANDS],
                               int16_t num_cand,
                               int32_t ref_idx,
                               double* bitcost)
{
  cabac_data_t state_cabac_copy;
  cabac_data_t* cabac;
  uint32_t merge_idx;
  vector2d_t mvd = { 0, 0 };
  int8_t merged = 0;
  int8_t cur_mv_cand = 0;

  x *= 1 << mv_shift;
  y *= 1 << mv_shift;

  // Check every candidate to find a match
  for (merge_idx = 0; merge_idx < (uint32_t)num_cand; merge_idx++) {
    if (merge_cand[merge_idx].dir == 3) continue;
    if (merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][0] == x &&
      merge_cand[merge_idx].mv[merge_cand[merge_idx].dir - 1][1] == y &&
      state->frame->ref_LX[merge_cand[merge_idx].dir - 1][
        merge_cand[merge_idx].ref[merge_cand[merge_idx].dir - 1]
      ] == ref_idx)
    {
      merged = 1;
      break;
    }
  }

  // Store cabac state and contexts
  memcpy(&state_cabac_copy, &state->search_cabac, sizeof(cabac_data_t));

  // Clear bytes and bits and set mode to "count"
  state_cabac_copy.only_count = 1;

  cabac = &state_cabac_copy;
  double bits = 0;

  if (!merged) {
    vector2d_t mvd1 = {
      x - mv_cand[0][0],
      y - mv_cand[0][1],
    };
    vector2d_t mvd2 = {
      x - mv_cand[1][0],
      y - mv_cand[1][1],
    };

    uvg_change_precision_vector2d(INTERNAL_MV_PREC, 2, &mvd1);
    uvg_change_precision_vector2d(INTERNAL_MV_PREC, 2, &mvd2);

    double cand1_cost = uvg_get_mvd_coding_cost_cabac(state, cabac, mvd1.x, mvd1.y);
    double cand2_cost = uvg_get_mvd_coding_cost_cabac(state, cabac, mvd2.x, mvd2.y);

    // Select candidate 1 if it has lower cost
    if (cand2_cost < cand1_cost) {
      cur_mv_cand = 1;
      mvd = mvd2;
    } else {
      mvd = mvd1;
    }
  }

  cabac->cur_ctx = &(cabac->ctx.cu_merge_flag_ext_model);

  CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_merge_flag_ext_model), merged, bits, "MergeFlag");
  num_cand = state->encoder_control->cfg.max_merge;
  if (merged) {
    if (num_cand > 1) {
      int32_t ui;
      for (ui = 0; ui < num_cand - 1; ui++) {
        int32_t symbol = (ui != merge_idx);
        if (ui == 0) {
          CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_merge_idx_ext_model), symbol, bits, "MergeIndex");
        } else {
          CABAC_BIN_EP(cabac, symbol, "MergeIndex");
          bits += 1;
        }
        if (symbol == 0) break;
      }
    }
  } else {
    uint32_t ref_list_idx;
    uint32_t j;
    int ref_list[2] = { 0, 0 };
    for (j = 0; j < state->frame->ref->used_size; j++) {
      if (state->frame->ref->pocs[j] < state->frame->poc) {
        ref_list[0]++;
      } else {
        ref_list[1]++;
      }
    }

    //ToDo: bidir mv support
    for (ref_list_idx = 0; ref_list_idx < 2; ref_list_idx++) {
      if (/*cur_cu->inter.mv_dir*/ 1 & (1 << ref_list_idx)) {
        if (ref_list[ref_list_idx] > 1) {
          // parseRefFrmIdx
          int32_t ref_frame = ref_idx;
          
          CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_ref_pic_model[0]), (ref_frame != 0), bits, "ref_idx_lX");

          if (ref_frame > 0) {
            int32_t i;
            int32_t ref_num = ref_list[ref_list_idx] - 2;
            
            ref_frame--;

            for (i = 0; i < ref_num; ++i) {
              const uint32_t symbol = (i == ref_frame) ? 0 : 1;

              if (i == 0) {
                CABAC_FBITS_UPDATE(cabac, &(cabac->ctx.cu_ref_pic_model[1]), symbol, bits, "ref_idx_lX");
              } else {
                CABAC_BIN_EP(cabac, symbol, "ref_idx_lX");
                bits += 1;
              }
              if (symbol == 0) break;
            }
          }
        }

        // ToDo: Bidir vector support
        if (!(state->frame->ref_list == REF_PIC_LIST_1 && /*cur_cu->inter.mv_dir == 3*/ 0)) {
          // It is safe to drop const here because cabac->only_count is set.
          uvg_encode_mvd((encoder_state_t*) state, cabac, mvd.x, mvd.y, &bits);
        }

        // Signal which candidate MV to use
        cabac->cur_ctx = &(cabac->ctx.mvp_idx_model);
        CABAC_BIN(cabac, cur_mv_cand, "mvp_flag");
      }
    }
  }

  *bitcost = bits;

  // Store bitcost before restoring cabac
  return *bitcost * state->lambda_sqrt;
}

void uvg_close_rdcost_outfiles(void)
{
  int i;

  for (i = 0; i < RD_SAMPLING_MAX_LAST_QP; i++) {
    FILE *curr = fastrd_learning_outfile[i];
    pthread_mutex_t *curr_mtx = outfile_mutex + i;
    if (curr != NULL) {
      fclose(curr);
    }
    if (curr_mtx != NULL) {
      pthread_mutex_destroy(curr_mtx);
    }
  }
}
