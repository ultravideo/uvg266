#ifndef UVG266_H_
#define UVG266_H_
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
 * \ingroup Control
 * \file
 * This file defines the public API of uvg266 when used as a library.
 */

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(UVG_DLL_EXPORTS)
  #if !defined(PIC)
    // Building static uvg266 library.
    #define UVG_PUBLIC
  #elif defined(_WIN32) || defined(__CYGWIN__)
    // Building uvg266 DLL on Windows.
    #define UVG_PUBLIC __declspec(dllexport)
  #elif defined(__GNUC__)
    // Building uvg266 shared library with GCC.
    #define UVG_PUBLIC __attribute__ ((visibility ("default")))
  #else
    #define UVG_PUBLIC
  #endif
#else
  #if defined(UVG_STATIC_LIB)
    // Using static uvg266 library.
    #define UVG_PUBLIC
  #elif defined(_WIN32) || defined(__CYGWIN__)
    // Using uvg266 DLL on Windows.
    #define UVG_PUBLIC __declspec(dllimport)
  #else
    // Using uvg266 shared library and not on Windows.
    #define UVG_PUBLIC
  #endif
#endif

/**
 * Maximum length of a GoP structure.
 */
#define UVG_MAX_GOP_LENGTH 32

 /**
 * Maximum amount of GoP layers.
 */
#define UVG_MAX_GOP_LAYERS 6

/**
 * Size of data chunks.
 */
#define UVG_DATA_CHUNK_SIZE 4096

#ifndef UVG_BIT_DEPTH
#define UVG_BIT_DEPTH 8
#endif

#if UVG_BIT_DEPTH == 8
typedef uint8_t uvg_pixel;
#else
typedef uint16_t uvg_pixel;
#endif

typedef int16_t uvg_pixel_im;  // For intermediate precision (interpolation/bipred).

/**
 * \brief Opaque data structure representing one instance of the encoder.
 */
typedef struct uvg_encoder uvg_encoder;

/**
 * \brief Integer motion estimation algorithms.
 */
enum uvg_ime_algorithm {
  UVG_IME_HEXBS = 0,
  UVG_IME_TZ = 1,
  UVG_IME_FULL = 2,
  UVG_IME_FULL8 = 3, //! \since 3.6.0
  UVG_IME_FULL16 = 4, //! \since 3.6.0
  UVG_IME_FULL32 = 5, //! \since 3.6.0
  UVG_IME_FULL64 = 6, //! \since 3.6.0
  UVG_IME_DIA = 7, // Experimental. TODO: change into a proper doc comment
};

/**
 * \brief Interlacing methods.
 * \since 3.2.0
 */
enum uvg_interlacing
{
  UVG_INTERLACING_NONE = 0,
  UVG_INTERLACING_TFF = 1, // top field first
  UVG_INTERLACING_BFF = 2, // bottom field first
};

/**
* \brief Constrain movement vectors.
* \since 3.3.0
*/
enum uvg_mv_constraint
{
  UVG_MV_CONSTRAIN_NONE = 0,
  UVG_MV_CONSTRAIN_FRAME = 1,  // Don't refer outside the frame.
  UVG_MV_CONSTRAIN_TILE = 2,  // Don't refer to other tiles.
  UVG_MV_CONSTRAIN_FRAME_AND_TILE = 3,  // Don't refer outside the tile.
  UVG_MV_CONSTRAIN_FRAME_AND_TILE_MARGIN = 4,  // Keep enough margin for fractional pixel margins not to refer outside the tile.
};

/**
* \brief Constrain movement vectors.
* \since 3.5.0
*/
enum uvg_hash
{
  UVG_HASH_NONE = 0,
  UVG_HASH_CHECKSUM = 1,
  UVG_HASH_MD5 = 2,
};

/**
* \brief cu split termination mode
* \since since 3.8.0
*/
enum uvg_cu_split_termination
{
  UVG_CU_SPLIT_TERMINATION_ZERO = 0,
  UVG_CU_SPLIT_TERMINATION_OFF = 1
};

/**
* \brief me early termination mode
* \since since 3.8.0
*/
enum uvg_me_early_termination
{
  UVG_ME_EARLY_TERMINATION_OFF = 0,
  UVG_ME_EARLY_TERMINATION_ON = 1,
  UVG_ME_EARLY_TERMINATION_SENSITIVE = 2
};


/**
 * \brief Format the pixels are read in.
 * This is separate from chroma subsampling, because we might want to read
 * interleaved formats in the future.
 * \since 3.12.0
 */
enum uvg_input_format {
  UVG_FORMAT_P400 = 0,
  UVG_FORMAT_P420 = 1,
  UVG_FORMAT_P422 = 2,
  UVG_FORMAT_P444 = 3,
};

/**
* \brief Chroma subsampling format used for encoding.
* \since 3.12.0
*/
enum uvg_chroma_format {
  UVG_CSP_400 = 0,
  UVG_CSP_420 = 1,
  UVG_CSP_422 = 2,
  UVG_CSP_444 = 3,
};

/**
 * \brief Chroma subsampling format used for encoding.
 * \since 3.15.0
 */
enum uvg_slices {
  UVG_SLICES_NONE,
  UVG_SLICES_TILES = (1 << 0), /*!< \brief Put each tile in a slice. */
  UVG_SLICES_WPP   = (1 << 1), /*!< \brief Put each row in a slice. */
};

enum uvg_sao {
  UVG_SAO_OFF = 0,
  UVG_SAO_EDGE = 1,
  UVG_SAO_BAND = 2,
  UVG_SAO_FULL = 3
};

enum uvg_alf {
  UVG_ALF_OFF = 0,
  UVG_ALF_NO_CC_ALF = 1,
  UVG_ALF_FULL = 2
};

enum uvg_mts {
  UVG_MTS_OFF = 0,
  UVG_MTS_INTRA = 1,
  UVG_MTS_INTER = 2,
  UVG_MTS_BOTH = 3,
  UVG_MTS_IMPLICIT = 4,
};


//MTS transform tags
typedef enum tr_type_t {
  DCT2 = 0,
  DCT8 = 1,
  DST7 = 2,
  NUM_TRANS_TYPE = 3,
  DCT2_MTS = 4
} tr_type_t;

enum uvg_scalinglist {
  UVG_SCALING_LIST_OFF = 0,
  UVG_SCALING_LIST_CUSTOM = 1,
  UVG_SCALING_LIST_DEFAULT = 2,  
};

enum uvg_rc_algorithm
{
  UVG_NO_RC = 0,
  UVG_LAMBDA = 1,
  UVG_OBA = 2,
};

enum uvg_file_format
{
  UVG_FORMAT_AUTO = 0,
  UVG_FORMAT_Y4M = 1,
  UVG_FORMAT_YUV = 2
};

enum uvg_amvr_resolution
{
  UVG_IMV_OFF     = 0,
  UVG_IMV_FPEL    = 1,
  UVG_IMV_4PEL    = 2,
  UVG_IMV_HPEL    = 3
};

enum uvg_roi_format
{
  UVG_ROI_TXT = 0,
  UVG_ROI_BIN = 1
};

// Map from input format to chroma format.
#define UVG_FORMAT2CSP(format) ((enum uvg_chroma_format)format)

/**
 * \brief GoP picture configuration.
 */
typedef struct uvg_gop_config {
  double qp_factor;
  int8_t qp_offset;    /*!< \brief QP offset */
  int8_t poc_offset;   /*!< \brief POC offset */
  int8_t layer;        /*!< \brief Current layer */
  int8_t is_ref;       /*!< \brief Flag if this picture is used as a reference */
  int8_t ref_pos_count;/*!< \brief Reference picture count */
  int8_t ref_pos[16];  /*!< \brief reference picture offset list */
  int8_t ref_neg_count;/*!< \brief Reference picture count */
  int8_t ref_neg[16];  /*!< \brief reference picture offset list */
  double qp_model_offset;
  double qp_model_scale;
} uvg_gop_config;

/**
 * \brief Struct which contains all configuration data
 *
 * Functions config_alloc, config_init and config_destroy must be used to
 * maintain ABI compatibility. Do not copy this struct, as the size might
 * change.
 */
typedef struct uvg_config
{
  int32_t qp;        /*!< \brief Quantization parameter */
  int32_t intra_period; /*!< \brief the period of intra frames in stream */

  /** \brief How often the VPS, SPS and PPS are re-sent
   *
   * -1: never
   *  0: first frame only
   *  1: every intra frame
   *  2: every other intra frame
   *  3: every third intra frame
   *  and so on
   */
  int32_t vps_period;

  int32_t width;   /*!< \brief frame width, must be a multiple of 8 */
  int32_t height;  /*!< \brief frame height, must be a multiple of 8 */
  int32_t framerate_num; /*!< \brief Framerate numerator */
  int32_t framerate_denom; /*!< \brief Framerate denominator */
  int32_t lmcs_enable;   /*!< \brief Flag to enable luma mapping with chroma scaling - filter */
  int32_t deblock_enable; /*!< \brief Flag to enable deblocking filter */
  enum uvg_sao sao_type;     /*!< \brief Flag to enable sample adaptive offset filter */
  enum uvg_alf alf_type;     /*!< \brief Flag to enable adaptive loop filter */
  int32_t alf_info_in_ph_flag; /*!< \brief Flag to enable if ALF is applied to all slices in picture */
  int32_t alf_slice_enable_flag[3/*MAX_NUM_COMPONENT*/];
  int32_t alf_non_linear_luma;    /*!< \brief Flag to enable non linear alf for luma */
  int32_t alf_non_linear_chroma;    /*!< \brief Flag to enable non linear alf for chroma */
  int32_t alf_allow_predefined_filters;
  int32_t rdoq_enable;    /*!< \brief Flag to enable RD optimized quantization. */
  int32_t signhide_enable;   /*!< \brief Flag to enable sign hiding. */
  int32_t rdo;            /*!< \brief RD-calculation level (0..2) */
  int32_t full_intra_search; /*!< \brief If true, don't skip modes in intra search. */
  int32_t trskip_enable;    /*!< \brief Flag to enable transform skip. */
  int32_t chroma_trskip_enable; /*!< \brief Flag to enable transform skip for chroma blocks. */
  int32_t trskip_max_size;    /*!< \brief Transform skip max block size. */
  enum uvg_mts mts;        /*< \brief flag to enable multiple transform selection*/
  int32_t mts_implicit;        /*< \brief flag to enable implicit multiple transform selection*/
  enum uvg_ime_algorithm ime_algorithm;  /*!< \brief Integer motion estimation algorithm. */
  int32_t fme_level;      /*!< \brief Fractional pixel motion estimation level (0: disabled, 1: enabled). */
  int8_t source_scan_type; /*!< \brief Source scan type (0: progressive, 1: top field first, 2: bottom field first).*/
  int32_t bipred;         /*!< \brief Bi-prediction (0: disabled, 1: enabled). */
  int32_t deblock_beta;   /*!< \brief (deblocking) beta offset (div 2), range -6...6 */
  int32_t deblock_tc;     /*!< \brief (deblocking) tc offset (div 2), range -6...6 */
  struct
  {
    int32_t sar_width;   /*!< \brief the horizontal size of the sample aspect ratio (in arbitrary units) */
    int32_t sar_height;  /*!< \brief the vertical size of the sample aspect ratio (in the same arbitrary units as sar_width). */
    int8_t overscan;     /*!< \brief Crop overscan setting */
    int8_t videoformat;  /*!< \brief Video format */
    int8_t fullrange;    /*!< \brief Flag to indicate full-range */
    int8_t colorprim;    /*!< \brief Color primaries */
    int8_t transfer;     /*!< \brief Transfer characteristics */
    int8_t colormatrix;  /*!< \brief Color matrix coefficients */
    int32_t chroma_loc;   /*!< \brief Chroma sample location */
  } vui;
  int32_t aud_enable;     /*!< \brief Flag to use access unit delimiters */
  int32_t ref_frames;     /*!< \brief number of reference frames to use */
  char * cqmfile;        /*!< \brief Pointer to custom quantization matrices filename */

  int32_t tiles_width_count;      /*!< \brief number of tiles separation in x direction */
  int32_t tiles_height_count;      /*!< \brief number of tiles separation in y direction */
  int32_t* tiles_width_split;      /*!< \brief tiles split x coordinates (dimension: tiles_width_count) */
  int32_t* tiles_height_split;      /*!< \brief tiles split y coordinates (dimension: tiles_height_count) */

  int wpp;
  int owf;

  int32_t slice_count;
  int32_t* slice_addresses_in_ts;

  int32_t threads;
  int32_t cpuid;

  struct {
    int32_t min[UVG_MAX_GOP_LAYERS];
    int32_t max[UVG_MAX_GOP_LAYERS];
  } pu_depth_inter, pu_depth_intra;

  int32_t add_encoder_info;
  int8_t gop_len;            /*!< \brief length of GOP for the video sequence */
  int8_t gop_lowdelay;       /*!< \brief specifies that the GOP does not use future pictures */
  uvg_gop_config gop[UVG_MAX_GOP_LENGTH];  /*!< \brief Array of GOP settings */

  int32_t target_bitrate;

  int8_t mv_rdo;            /*!< \brief MV RDO calculation in search (0: estimation, 1: RDO). */
  int8_t calc_psnr;         /*!< \since 3.1.0 \brief Print PSNR in CLI. */

  enum uvg_mv_constraint mv_constraint;  /*!< \since 3.3.0 \brief Constrain movement vectors. */
  enum uvg_hash hash;  /*!< \since 3.5.0 \brief What hash algorithm to use. */

  enum uvg_cu_split_termination cu_split_termination; /*!< \since 3.8.0 \brief Mode of cu split termination. */

  enum uvg_me_early_termination me_early_termination; /*!< \since 3.8.0 \brief Mode of me early termination. */
  int32_t intra_rdo_et; /*!< \since 4.1.0 \brief Use early termination in intra rdo. */

  int32_t lossless; /*!< \brief Use lossless coding. */

  int32_t tmvp_enable; /*!> \brief Use Temporal Motion Vector Predictors. */

  int32_t rdoq_skip; /*!< \brief Mode of rdoq skip */

  enum uvg_input_format input_format; /*!< \brief Use Temporal Motion Vector Predictors. */
  int32_t input_bitdepth; /*!< \brief Use Temporal Motion Vector Predictors. */

  struct {
    unsigned d;  // depth
    unsigned t;  // temporal
  } gop_lp_definition;

  int32_t implicit_rdpcm; /*!< \brief Enable implicit residual DPCM. */

  struct {
    char *file_path;
    enum uvg_roi_format format;
  } roi; /*!< \brief Specify delta QPs for region of interest coding. */

  unsigned slices; /*!< \since 3.15.0 \brief How to map slices to frame. */

  /**
   * \brief Use adaptive QP for 360 video with equirectangular projection.
   */
  int32_t erp_aqp;

  /** \brief The HEVC level */
  uint8_t level;
  /** \brief Whether we ignore and just warn from all of the errors about the output not conforming to the level's requirements. */
  uint8_t force_level;
  /** \brief Whether we use the high tier bitrates. Requires the level to be 4 or higher. */
  uint8_t high_tier;
  /** \brief The maximum allowed bitrate for this level and tier. */
  uint32_t max_bitrate;

  /** \brief Maximum steps that hexagonal and diagonal motion estimation can use. -1 to disable */
  uint32_t me_max_steps;

  /** \brief Offset to add to QP for intra frames */
  int8_t intra_qp_offset;
  /** \brief Select intra QP Offset based on GOP length */
  uint8_t intra_qp_offset_auto;

  /** \brief Minimum QP that uses CABAC for residual cost instead of a fast estimate. */
  int8_t fast_residual_cost_limit;

  /** \brief Set QP at CU level keeping pic_init_qp_minus26 in PPS zero */
  int8_t set_qp_in_cu;

  /** \brief Flag to enable/disable open GOP configuration */
  int8_t open_gop;

  int32_t vaq; /** \brief Enable variance adaptive quantization*/

  /** \brief Type of scaling lists to use */
  int8_t scaling_list;

  /** \brief Maximum number of merge cadidates */
  uint8_t max_merge;

  /** \brief Enable Early Skip Mode Decision */
  uint8_t early_skip;

  /** \brief Disable intra smoothing when true */
  uint8_t intra_smoothing_disabled;  /** \brief Enable Machine learning CU depth prediction for Intra encoding. */
  uint8_t ml_pu_depth_intra;  
  
  /** \brief Used for partial frame encoding*/
  struct {
    uint8_t startCTU_x;
    uint8_t startCTU_y;
    uint16_t fullWidth;
    uint16_t fullHeight;
  } partial_coding;

  /** \brief Always consider CU without any quantized residual */
  uint8_t zero_coeff_rdo;

  /** \brief Currently unused parameter for OBA rc */
  int8_t frame_allocation;

  /** \brief used rc scheme, 0 for QP */
  int8_t rc_algorithm;

  /** \brief whether to use hadamard based bit allocation for intra frames or not */
  uint8_t intra_bit_allocation;

  uint8_t clip_neighbour;

  enum uvg_file_format file_format;

  char *stats_file_prefix;

  uint8_t log2_parallel_merge_level;

  char *fast_coeff_table_fn;   /*!< \brief Pointer to fast coeff table filename */

  /** \brief whether we're sampling TBs and their costs for fast cost
   *         estimation training */
  uint8_t rdo_cost_sampling_mode_on;

  /** \brief whether we're running in normal mode, sampling TBs and their cost
   *         for fast estimation training, or comparing estimator accuracy to
   *         CABAC */
  uint8_t fastrd_sampling_on;
  uint8_t fastrd_accuracy_check_on;

  char *fastrd_learning_outdir_fn;

  int8_t num_used_table;
  int8_t qp_table_start_minus26[3];
  int8_t qp_table_length_minus1[3];
  int8_t delta_qp_in_val_minus1[3][16];
  int8_t delta_qp_out_val[3][16];

  int8_t chroma_scale_in[3][17];
  int8_t chroma_scale_out[3][17];

  /** \brief enable use of multiple reference lines in intra prediction */
  int8_t mrl;

  /** \brief enable matrix weighted intra prediction */
  int8_t mip;
  /** \brief enable low frequency non-separable transform */
  int8_t lfnst;

  /** \brief enable intra sub partitions*/
  int8_t isp;

  int8_t jccr;

  int8_t cclm;

  int8_t amvr; /* \brief Adaptive motion vector resolution parameter */

  /** \brief whether to try combining intra cus at the lower depth when search
   *         is not performed at said depth*/
  uint8_t combine_intra_cus;

  uint8_t force_inter;
  char* cabac_debug_file_name;

  uint8_t dual_tree;

  uint8_t min_qt_size[3]; /* intra, inter, dual tree chroma*/
  uint8_t max_bt_size[3]; /* intra, inter, dual tree chroma*/
  uint8_t max_tt_size[3]; /* intra, inter, dual tree chroma*/

  uint8_t max_btt_depth[3]; /* intra, inter, dual tree chroma*/

  uint8_t intra_rough_search_levels;

  uint8_t ibc; /* \brief Intra Block Copy parameter */
  uint8_t dep_quant;

  uint8_t ref_wraparound; /* \brief MV reference wraparound */

} uvg_config;

/**
 * \brief Struct which contains all picture data
 *
 * Function picture_alloc in uvg_api must be used for allocation.
 */
typedef struct uvg_picture {
  uvg_pixel *fulldata_buf;     //!< \brief Allocated buffer with padding (only used in the base_image)
  uvg_pixel *fulldata;         //!< \brief Allocated buffer portion that's actually used

  uvg_pixel *y;                //!< \brief Pointer to luma pixel array.
  uvg_pixel *u;                //!< \brief Pointer to chroma U pixel array.
  uvg_pixel *v;                //!< \brief Pointer to chroma V pixel array.
  uvg_pixel *data[3]; //!< \brief Alternate access method to same data.

  int32_t width;           //!< \brief Luma pixel array width.
  int32_t height;          //!< \brief Luma pixel array height.

  int32_t stride;          //!< \brief Luma pixel array width for the full picture (should be used as stride)

  struct uvg_picture *base_image; //!< \brief Pointer to the picture which owns the pixels
  int32_t refcount;        //!< \brief Number of references to the picture

  int64_t pts;             //!< \brief Presentation timestamp. Should be set for input frames.
  int64_t dts;             //!< \brief Decompression timestamp.

  enum uvg_interlacing interlacing; //!< \since 3.2.0 \brief Field order for interlaced pictures.
  enum uvg_chroma_format chroma_format;

  int32_t ref_pocs[16];

  struct
  {
    int width;
    int height;
    int8_t *roi_array;
  } roi;

} uvg_picture;

/**
 * \brief NAL unit type codes.
 *
 * These are the nal_unit_type codes from Table 7-1 ITU-T H.265 v1.0.
 */
enum uvg_nal_unit_type {

  // Coded slices

  UVG_NAL_TRAIL = 0,
  UVG_NAL_STSA = 1,
  UVG_NAL_RADL = 2,
  UVG_NAL_RASL = 3,

  // Intra random access point pictures
  UVG_NAL_IDR_W_RADL = 7,
  UVG_NAL_IDR_N_LP = 8,
  UVG_NAL_CRA_NUT = 9,
  UVG_NAL_GDR_NUT = 10,


  // non-VCL  
  UVG_NAL_VPS_NUT = 14,
  UVG_NAL_SPS_NUT = 15,
  UVG_NAL_PPS_NUT = 16,
  NAL_UNIT_PREFIX_APS = 17,
  NAL_UNIT_SUFFIX_APS = 18,

  UVG_NAL_AUD_NUT = 20,

  UVG_NAL_EOS_NUT = 21,
  UVG_NAL_EOB_NUT = 22,
  UVG_NAL_PREFIX_SEI_NUT = 23,
  UVG_NAL_SUFFIX_SEI_NUT = 24,
  

};

enum uvg_slice_type {
  UVG_SLICE_B = 0,
  UVG_SLICE_P = 1,
  UVG_SLICE_I = 2,
};

enum uvg_split_mode {
  UVG_NO_SPLIT = 0,
  UVG_QUAD_SPLIT = 1,
  UVG_HORZ_SPLIT = 2,
  UVG_VERT_SPLIT = 3,
};

/**
 * \brief Other information about an encoded frame
 */
typedef struct uvg_frame_info {

  /**
   * \brief Picture order count
   */
  int32_t poc;

  /**
   * \brief Quantization parameter
   */
  int8_t qp;

  /**
   * \brief Type of the NAL VCL unit
   */
  enum uvg_nal_unit_type nal_unit_type;

  /**
   * \brief Type of the slice
   */
  enum uvg_slice_type slice_type;

  /**
   * \brief Reference picture lists
   *
   * The first list contains the reference picture POCs that are less than the
   * POC of this frame and the second one contains those that are greater.
   */
  int ref_list[2][16];

  /**
   * \brief Lengths of the reference picture lists
   */
  int ref_list_len[2];

} uvg_frame_info;

/**
 * \brief A linked list of chunks of data.
 *
 * Used for returning the encoded data.
 */
typedef struct uvg_data_chunk {
  /// \brief Buffer for the data.
  uint8_t data[UVG_DATA_CHUNK_SIZE];

  /// \brief Number of bytes filled in this chunk.
  uint32_t len;

  /// \brief Next chunk in the list.
  struct uvg_data_chunk *next;
} uvg_data_chunk;

typedef struct uvg_api {

  /**
   * \brief Allocate a uvg_config structure.
   *
   * The returned structure should be deallocated by calling config_destroy.
   *
   * \return allocated config, or NULL if allocation failed.
   */
  uvg_config *  (*config_alloc)(void);

  /**
   * \brief Deallocate a uvg_config structure.
   *
   * If cfg is NULL, do nothing. Otherwise, the given structure must have been
   * returned from config_alloc.
   *
   * \param cfg   configuration
   * \return      1 on success, 0 on failure
   */
  int           (*config_destroy)(uvg_config *cfg);

  /**
   * \brief Initialize a config structure
   *
   * Set all fields in the given config to default values.
   *
   * \param cfg   configuration
   * \return      1 on success, 0 on failure
   */
  int           (*config_init)(uvg_config *cfg);

  /**
   * \brief Set an option.
   *
   * \param cfg   configuration
   * \param name  name of the option to set
   * \param value value to set
   * \return      1 on success, 0 on failure
   */
  int           (*config_parse)(uvg_config *cfg, const char *name, const char *value);

  /**
   * \brief Allocate a uvg_picture.
   *
   * The returned uvg_picture should be deallocated by calling picture_free.
   *
   * \param width   width of luma pixel array to allocate
   * \param height  height of luma pixel array to allocate
   * \return        allocated picture, or NULL if allocation failed.
   */
  uvg_picture * (*picture_alloc)(int32_t width, int32_t height);

  /**
   * \brief Deallocate a uvg_picture.
   *
   * If pic is NULL, do nothing. Otherwise, the picture must have been returned
   * from picture_alloc.
   */
  void          (*picture_free)(uvg_picture *pic);

  /**
   * \brief Deallocate a list of data chunks.
   *
   * Deallocates the given chunk and all chunks that follow it in the linked
   * list.
   */
  void          (*chunk_free)(uvg_data_chunk *chunk);

  /**
   * \brief Create an encoder.
   *
   * The returned encoder should be closed by calling encoder_close.
   *
   * Only one encoder may be open at a time.
   *
   * \param cfg   encoder configuration
   * \return      g_created encoder, or NULL if creation failed.
   */
  uvg_encoder * (*encoder_open)(const uvg_config *cfg);

  /**
   * \brief Deallocate an encoder.
   *
   * If encoder is NULL, do nothing. Otherwise, the encoder must have been
   * returned from encoder_open.
   */
  void          (*encoder_close)(uvg_encoder *encoder);

  /**
   * \brief Get parameter sets.
   *
   * Encode the VPS, SPS and PPS.
   *
   * If data_out is set to non-NULL values, the caller is responsible for
   * calling chunk_free on it.
   *
   * A null pointer may be passed in place of the parameter data_out or len_out
   * to skip returning the corresponding value.
   *
   * \param encoder   encoder
   * \param data_out  Returns the encoded parameter sets.
   * \param len_out   Returns number of bytes in the encoded data.
   * \return          1 on success, 0 on error.
   */
  int           (*encoder_headers)(uvg_encoder *encoder,
                                   uvg_data_chunk **data_out,
                                   uint32_t *len_out);

  /**
   * \brief Encode one frame.
   *
   * Add pic_in to the encoding pipeline. If an encoded frame is ready, return
   * the bitstream, length of the bitstream, the reconstructed frame, the
   * original frame and frame info in data_out, len_out, pic_out, src_out and
   * info_out, respectively. Otherwise, set the output parameters to NULL.
   * 
   * Region of interest (ROI) / delta QP map can be specified in the input
   * picture's ROI field but only when a ROI file is not used.
   *
   * After passing all of the input frames, the caller should keep calling this
   * function with pic_in set to NULL, until no more data is returned in the
   * output parameters.
   *
   * The caller must not modify pic_in after passing it to this function.
   *
   * If data_out, pic_out and src_out are set to non-NULL values, the caller is
   * responsible for calling chunk_free and picture_free on them.
   *
   * A null pointer may be passed in place of any of the parameters data_out,
   * len_out, pic_out, src_out or info_out to skip returning the corresponding
   * value.
   *
   * \param encoder   encoder
   * \param pic_in    input frame or NULL
   * \param data_out  Returns the encoded data.
   * \param len_out   Returns number of bytes in the encoded data.
   * \param pic_out   Returns the reconstructed picture.
   * \param src_out   Returns the original picture.
   * \param info_out  Returns information about the encoded picture.
   * \return          1 on success, 0 on error.
   */
  int           (*encoder_encode)(uvg_encoder *encoder,
                                  uvg_picture *pic_in,
                                  uvg_data_chunk **data_out,
                                  uint32_t *len_out,
                                  uvg_picture **pic_out,
                                  uvg_picture **src_out,
                                  uvg_frame_info *info_out);

  /**
   * \brief Allocate a uvg_picture.
   *
   * The returned uvg_picture should be deallocated by calling picture_free.
   *
   * \since 3.12.0
   * \param chroma_fomat  Chroma subsampling to use.
   * \param width   width of luma pixel array to allocate
   * \param height  height of luma pixel array to allocate
   * \return        allocated picture, or NULL if allocation failed.
   */
  uvg_picture * (*picture_alloc_csp)(enum uvg_chroma_format chroma_fomat, int32_t width, int32_t height);
} uvg_api;


UVG_PUBLIC const uvg_api * uvg_api_get(int bit_depth);

#ifdef __cplusplus
}
#endif

#endif // UVG266_H_
