

#include "alf.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cabac.h"
#include "rdo.h"
#include "strategies/strategies-sao.h"


void kvz_alf_init(encoder_state_t *const state,
  encoder_state_config_slice_t *slice,
  alf_info_t *alf)
{
  //Slice alf init. Johonkin muualle?

  /*
  for (int i = 0; i < MAX_NUM_APS; i++)
  {
    slice->param_set_map[i].parameter_set = malloc(sizeof(alf_aps));
  }*/

  //slice->param_set_map->p_nalu_data = malloc();

  //Default clp_rng for a slice
  g_clp_rngs.comp[COMPONENT_Y].min = g_clp_rngs.comp[COMPONENT_Cb].min  = g_clp_rngs.comp[COMPONENT_Cr].min = 0;
  g_clp_rngs.comp[COMPONENT_Y].max = (1<< ALF_NUM_BITS)-1;
  g_clp_rngs.comp[COMPONENT_Y].bd  = ALF_NUM_BITS;
  g_clp_rngs.comp[COMPONENT_Y].n   = 0;
  g_clp_rngs.comp[COMPONENT_Cb].max = g_clp_rngs.comp[COMPONENT_Cr].max = (1<< ALF_NUM_BITS)-1;
  g_clp_rngs.comp[COMPONENT_Cb].bd  = g_clp_rngs.comp[COMPONENT_Cr].bd  = ALF_NUM_BITS;
  g_clp_rngs.comp[COMPONENT_Cb].n   = g_clp_rngs.comp[COMPONENT_Cr].n   = 0;
  g_clp_rngs.used = g_clp_rngs.chroma = false;

  /*
  for (int i = 0; i < MAX_NUM_CHANNEL_TYPE; i++)
  {
    g_alf_covariance_frame[i] = malloc(sizeof(g_alf_covariance_frame[i]));
  }*/
 
  /*alf->classifier = NULL;
  for (int i = 0; i < NUM_DIRECTIONS; i++)
  {
    alf->laplacian[i] = NULL;
  }*/
  alf->created = false;
}

//-------------------------help functions---------------------------

int alf_clip_pixel(const int a, const clp_rng clp_rng)
{
  return MIN(MAX(clp_rng.min, a), clp_rng.max);
}

int alf_clip3(const int minVal, const int maxVal, const int a)
{
  return MIN(MAX(minVal, a), maxVal);
}

void gns_backsubstitution(double r[MAX_NUM_ALF_COEFF][MAX_NUM_ALF_COEFF], double* z, int size, double* a)
{
  size--;
  a[size] = z[size] / r[size][size];

  for (int i = size - 1; i >= 0; i--)
  {
    double sum = 0;

    for (int j = i + 1; j <= size; j++)
    {
      sum += r[i][j] * a[j];
    }

    a[i] = (z[i] - sum) / r[i][i];
  }
}

void gns_transpose_backsubstitution(double u[MAX_NUM_ALF_COEFF][MAX_NUM_ALF_COEFF], double* rhs, double* x, int order)
{
  /* Backsubstitution starts */
  x[0] = rhs[0] / u[0][0];               /* First row of U' */
  for (int i = 1; i < order; i++)
  {         /* For the rows 1..order-1 */

    double sum = 0; //Holds backsubstitution from already handled rows

    for (int j = 0; j < i; j++) /* Backsubst already solved unknowns */
    {
      sum += x[j] * u[j][i];
    }

    x[i] = (rhs[i] - sum) / u[i][i];       /* i'th component of solution vect. */
  }
}

int gns_cholesky_dec(double **inp_matr, double out_matr[MAX_NUM_ALF_COEFF][MAX_NUM_ALF_COEFF], int num_eq)
{
  static double inv_diag[MAX_NUM_ALF_COEFF];  /* Vector of the inverse of diagonal entries of outMatr */

  for (int i = 0; i < num_eq; i++)
  {
    for (int j = i; j < num_eq; j++)
    {
      /* Compute the scaling factor */
      double scale = inp_matr[i][j];
      if (i > 0)
      {
        for (int k = i - 1; k >= 0; k--)
        {
          scale -= out_matr[k][j] * out_matr[k][i];
        }
      }

      /* Compute i'th row of outMatr */
      if (i == j)
      {
        if (scale <= REG_SQR) // if(scale <= 0 )  /* If inpMatr is singular */
        {
          return 0;
        }
        else              /* Normal operation */
          inv_diag[i] = 1.0 / (out_matr[i][i] = sqrt(scale));
      }
      else
      {
        out_matr[i][j] = scale * inv_diag[i]; /* Upper triangular part          */
        out_matr[j][i] = 0.0;              /* Lower triangular part set to 0 */
      }
    }
  }
  return 1; /* Signal that Cholesky factorization is successfully performed */
}

int gns_solve_by_chol(double **lhs, double *rhs, double *x, int num_eq)
{
  static double aux[MAX_NUM_ALF_COEFF];     /* Auxiliary vector */
  static double u[MAX_NUM_ALF_COEFF][MAX_NUM_ALF_COEFF];    /* Upper triangular Cholesky factor of LHS */
  int res = 1;  // Signal that Cholesky factorization is successfully performed

                /* The equation to be solved is LHSx = rhs */

                /* Compute upper triangular U such that U'*U = LHS */
  if (gns_cholesky_dec(lhs, u, num_eq)) /* If Cholesky decomposition has been successful */
  {
    /* Now, the equation is  U'*U*x = rhs, where U is upper triangular
    * Solve U'*aux = rhs for aux
    */
    gns_transpose_backsubstitution(u, rhs, aux, num_eq);

    /* The equation is now U*x = aux, solve it for x (new motion coefficients) */
    gns_backsubstitution(u, aux, num_eq, x);

  }
  else /* LHS was singular */
  {
    res = 0;

    /* Regularize LHS */
    for (int i = 0; i < num_eq; i++)
    {
      lhs[i][i] += REG;
    }

    /* Compute upper triangular U such that U'*U = regularized LHS */
    res = gns_cholesky_dec(lhs, u, num_eq);

    if (!res)
    {
      memset(x, 0, sizeof(double)*num_eq);
      return 0;
    }

    /* Solve  U'*aux = rhs for aux */
    gns_transpose_backsubstitution(u, rhs, aux, num_eq);

    /* Solve U*x = aux for x */
    gns_backsubstitution(u, aux, num_eq, x);
  }
  return res;
}

double calc_error_for_coeffs(double **ee, double *y, const int *coeff, const int num_coeff, const int bit_depth)
{
  double factor = 1 << (bit_depth - 1);
  double error = 0;

  for (int i = 0; i < num_coeff; i++)   //diagonal
  {
    double sum = 0;
    for (int j = i + 1; j < num_coeff; j++)
    {
      // E[j][i] = E[i][j], sum will be multiplied by 2 later
      sum += ee[i][j] * coeff[j];
    }
    error += ((ee[i][i] * coeff[i] + sum * 2) / factor - 2 * y[i]) * coeff[i];
  }

  return error / factor;
}

int get_golomb_k_min(channel_type channel, const int num_filters, int k_min_tab[MAX_NUM_ALF_LUMA_COEFF], int bits_coeff_scan[11/*m_MAX_SCAN_VAL*/][16/*m_MAX_EXP_GOLOMB*/])
{
  int k_start;
  const int max_golomb_idx = channel == CHANNEL_TYPE_LUMA ? 3 : 2;

  int min_bits_k_start = MAX_INT;
  int min_k_start = -1;

  for (int k = 1; k < 8; k++)
  {
    int bits_k_start = 0; k_start = k;
    for (int scan_pos = 0; scan_pos < max_golomb_idx; scan_pos++)
    {
      int k_min = k_start;
      int min_bits = bits_coeff_scan[scan_pos][k_min];

      if (bits_coeff_scan[scan_pos][k_start + 1] < min_bits)
      {
        k_min = k_start + 1;
        min_bits = bits_coeff_scan[scan_pos][k_min];
      }
      k_start = k_min;
      bits_k_start += min_bits;
    }
    if (bits_k_start < min_bits_k_start)
    {
      min_bits_k_start = bits_k_start;
      min_k_start = k;
    }
  }

  k_start = min_k_start;
  for (int scan_pos = 0; scan_pos < max_golomb_idx; scan_pos++)
  {
    int k_min = k_start;
    int min_bits = bits_coeff_scan[scan_pos][k_min];

    if (bits_coeff_scan[scan_pos][k_start + 1] < min_bits)
    {
      k_min = k_start + 1;
      min_bits = bits_coeff_scan[scan_pos][k_min];
    }

    k_min_tab[scan_pos] = k_min;
    k_start = k_min;
  }

  return min_k_start;
}


int length_golomb(int coeff_val, int k)
{
  int m = 2 << (k - 1);
  int q = coeff_val / m;
  if (coeff_val != 0)
  {
    return q + 2 + k;
  }
  else
  {
    return q + 1 + k;
  }
}

double get_dist_force_0(channel_type channel, const int num_filters, double error_tab_force_0_coeff[MAX_NUM_ALF_CLASSES][2], bool* coded_var_bins)
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  const int *golomb_idx = channel == CHANNEL_TYPE_LUMA ? alf_golomb_idx_7 : alf_golomb_idx_5;

  static int bits_var_bin[MAX_NUM_ALF_CLASSES];

  memset(g_bits_coeff_scan, 0, sizeof(g_bits_coeff_scan));
  for (int ind = 0; ind < num_filters; ++ind)
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      int coeff_val = abs(g_filter_coeff_set[ind][i]);
      for (int k = 1; k < 15; k++)
      {
        g_bits_coeff_scan[golomb_idx[i]][k] += length_golomb(coeff_val, k);
      }
    }
  }

  get_golomb_k_min(channel, num_filters, g_k_min_tab, g_bits_coeff_scan);

  for (int ind = 0; ind < num_filters; ++ind)
  {
    bits_var_bin[ind] = 0;
    for (int i = 0; i < num_coeff - 1; i++)
    {
      bits_var_bin[ind] += length_golomb(abs(g_filter_coeff_set[ind][i]), g_k_min_tab[golomb_idx[i]]);
    }
  }

  double dist_force_0 = 0;
  memset(coded_var_bins, 0, sizeof(*coded_var_bins) * MAX_NUM_ALF_CLASSES);

  for (int filtIdx = 0; filtIdx < num_filters; filtIdx++)
  {
    double costDiff = error_tab_force_0_coeff[filtIdx][0] - (error_tab_force_0_coeff[filtIdx][1] + 981.62883931057581/*m_lambda[COMPONENT_Y]*/ * bits_var_bin[filtIdx]);
    coded_var_bins[filtIdx] = costDiff > 0 ? true : false;
    dist_force_0 += error_tab_force_0_coeff[filtIdx][coded_var_bins[filtIdx] ? 1 : 0];
  }

  return dist_force_0;
}

int get_cost_filter_coeff_force_0(channel_type channel, int **p_diff_q_filter_coeff_int_pp, const int num_filters, bool* coded_var_bins)
{
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  const int max_golomb_idx = channel == CHANNEL_TYPE_LUMA ? 3 : 2;
  const int *golomb_idx = channel == CHANNEL_TYPE_LUMA ? alf_golomb_idx_7 : alf_golomb_idx_5;
  memset(g_bits_coeff_scan, 0, sizeof(g_bits_coeff_scan));

  for (int ind = 0; ind < num_filters; ++ind)
  {
    if (!coded_var_bins[ind])
    {
      continue;
    }
    for (int i = 0; i < num_coeff - 1; i++)
    {
      int coeff_val = abs(p_diff_q_filter_coeff_int_pp[ind][i]);
      for (int k = 1; k < 15; k++)
      {
        g_bits_coeff_scan[golomb_idx[i]][k] += length_golomb(coeff_val, k);
      }
    }
  }

  int k_min = get_golomb_k_min(channel, num_filters, g_k_min_tab, g_bits_coeff_scan);

  // Coding parameters
  int len = k_min           //min_golomb_order
    + max_golomb_idx   //golomb_order_increase_flag
    + num_filters;    //filter_coefficient_flag[i]

                     // Filter coefficients
  for (int ind = 0; ind < num_filters; ++ind)
  {
    if (coded_var_bins[ind])
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        len += length_golomb(abs(p_diff_q_filter_coeff_int_pp[ind][i]), g_k_min_tab[golomb_idx[i]]); // alf_coeff_luma_delta[i][j]
      }
    }
  }
  return len;
}

int get_cost_filter_coeff(channel_type channel, int **p_diff_q_filter_coeff_int_pp, const int num_filters)
{
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  const int max_golomb_idx = channel == CHANNEL_TYPE_LUMA ? 3 : 2;
  const int *golomb_idx = channel == CHANNEL_TYPE_LUMA ? alf_golomb_idx_7 : alf_golomb_idx_5;

  memset(g_bits_coeff_scan, 0, sizeof(g_bits_coeff_scan));

  for (int ind = 0; ind < num_filters; ++ind)
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      int coeff_val = abs(p_diff_q_filter_coeff_int_pp[ind][i]);
      for (int k = 1; k < 15; k++)
      {
        g_bits_coeff_scan[golomb_idx[i]][k] += length_golomb(coeff_val, k);
      }
    }
  }

  int k_min = get_golomb_k_min(channel, num_filters, g_k_min_tab, g_bits_coeff_scan);

  // Coding parameters
  int len = k_min           //min_golomb_order
    + max_golomb_idx;  //golomb_order_increase_flag

                     // Filter coefficients

  //len += lengthFilterCoeffs( alfShape, numFilters, pDiffQFilterCoeffIntPP, m_kMinTab );  // alf_coeff_luma_delta[i][j]
  for (int ind = 0; ind < num_filters; ++ind)
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      len += length_golomb(abs(p_diff_q_filter_coeff_int_pp[ind][i]), g_k_min_tab[golomb_idx[i]]);
    }
  }

  return len;
}

int get_tb_length(int ui_symbol, const int ui_max_symbol)
{
  int ui_thresh;
  if (ui_max_symbol > 256)
  {
    int ui_thresh_val = 1 << 8;
    ui_thresh = 8;
    while (ui_thresh_val <= ui_max_symbol)
    {
      ui_thresh++;
      ui_thresh_val <<= 1;
    }
    ui_thresh--;
  }
  else
  {
    ui_thresh = kvz_tb_max[ui_max_symbol];
  }

  int ui_val = 1 << ui_thresh;
  assert(ui_val <= ui_max_symbol);
  assert((ui_val << 1) > ui_max_symbol);
  assert(ui_symbol < ui_max_symbol);
  int b = ui_max_symbol - ui_val;
  assert(b < ui_val);
  if (ui_symbol < ui_val - b)
  {
    return ui_thresh;
  }
  else
  {
    return ui_thresh + 1;
  }
}

int get_non_filter_coeff_rate(alf_aps *aps)
{
  short* filter_coeff_delta_idx = aps->filter_coeff_delta_idx;
  int *num_luma_filters = &aps->num_luma_filters;
  int *fixed_filter_pattern = &aps->fixed_filter_pattern;
  int *fixed_filter_set_index = &aps->fixed_filter_set_index;

  int len = 1   // alf_coefficients_delta_flag
    + 1 //lengthTruncatedUnary(0, 3)    // chroma_idc = 0, it is signalled when ALF is enabled for luma
    + get_tb_length(*num_luma_filters - 1, MAX_NUM_ALF_CLASSES);   //numLumaFilters

  if (*num_luma_filters > 1)
  {
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
    {
      len += get_tb_length((int)filter_coeff_delta_idx[i], *num_luma_filters);  //filter_coeff_delta[i]
    }
  }
  len++; //fixed filter set flag
  if (*fixed_filter_set_index > 0)
  {
    len += get_tb_length(*fixed_filter_set_index - 1, ALF_NUM_FIXED_FILTER_SETS);
    len += 1; //fixed filter flag pattern
    if (*fixed_filter_pattern > 0)
      len += MAX_NUM_ALF_CLASSES;  //"fixed_filter_flag" for each class
  }

  return len;
}


double calculate_error(alf_covariance cov)
{
  static double c[MAX_NUM_ALF_COEFF];

  gns_solve_by_chol(cov.ee, cov.y, c, cov.num_coeff);

  double sum = 0;
  for (int i = 0; i < cov.num_coeff; i++)
  {
    sum += c[i] * cov.y[i];
  }

  return cov.pix_acc - sum;
}

int get_coeff_rate(alf_aps *aps, bool is_chroma)
{
  short *luma_coeff = aps->luma_coeff;
  short *chroma_coeff = aps->chroma_coeff;
  bool *alf_luma_coeff_flag = aps->alf_luma_coeff_flag;
  int *num_luma_filters = &aps->num_luma_filters;
  bool *alf_luma_coeff_delta_flag = &aps->alf_luma_coeff_delta_flag;

  channel_type channel = !is_chroma ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;
  const int *golomb_idx = !is_chroma ? alf_golomb_idx_7 : alf_golomb_idx_5;
  const int num_coeff = !is_chroma ? 13 : 7;

  int i_bits = 0;
  if (!is_chroma)
  {
    i_bits++;                                               // alf_coefficients_delta_flag
    if (!*alf_luma_coeff_delta_flag)
    {
      if (*num_luma_filters > 1)
      {
        i_bits++;                                           // coeff_delta_pred_mode_flag
      }
    }
  }

  memset(g_bits_coeff_scan, 0, sizeof(g_bits_coeff_scan));

  const int max_golomb_idx = is_chroma ? 2 : 3;

  const short* coeff = is_chroma ? chroma_coeff : luma_coeff;
  const int num_filters = is_chroma ? 1 : *num_luma_filters;

  // vlc for all
  for (int ind = 0; ind < num_filters; ++ind)
  {
    if (is_chroma || !*alf_luma_coeff_delta_flag || alf_luma_coeff_flag[ind])
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        int coeff_val = abs(coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i]);

        for (int k = 1; k < 15; k++)
        {
          g_bits_coeff_scan[golomb_idx[i]][k] += length_golomb(coeff_val, k);
        }
      }
    }
  }

  int k_min = get_golomb_k_min(channel, num_filters, g_k_min_tab, g_bits_coeff_scan);

  // Golomb parameters
  //i_bits += lengthUvlc(k_min - 1);  // "min_golomb_order"
  int ui_code = k_min - 1;
  int ui_length = 1;
  int ui_temp = ++ui_code;

  //CHECK(!ui_temp, "Integer overflow");

  while (1 != ui_temp)
  {
    ui_temp >>= 1;
    ui_length += 2;
  }
  // Take care of cases where uiLength > 32
  i_bits += (ui_length >> 1) + ((ui_length + 1) >> 1);

  int golomb_order_increase_flag = 0;

  for (int idx = 0; idx < max_golomb_idx; idx++)
  {
    golomb_order_increase_flag = (g_k_min_tab[idx] != k_min) ? 1 : 0;
    //CHECK(!(g_k_min_tab[idx] <= k_min + 1), "ALF Golomb parameter not consistent");
    i_bits += golomb_order_increase_flag;                           //golomb_order_increase_flag
    k_min = g_k_min_tab[idx];
  }

  if (!is_chroma)
  {
    if (*alf_luma_coeff_delta_flag)
    {
      i_bits += num_filters;             //filter_coefficient_flag[i]
    }
  }

  // Filter coefficients
  for (int ind = 0; ind < num_filters; ++ind)
  {
    if (!is_chroma && !alf_luma_coeff_flag[ind] && *alf_luma_coeff_delta_flag)
    {
      continue;
    }

    for (int i = 0; i < num_coeff - 1; i++)
    {
      i_bits += length_golomb(coeff[ind* MAX_NUM_ALF_LUMA_COEFF + i], g_k_min_tab[golomb_idx[i]]);  // alf_coeff_chroma[i], alf_coeff_luma_delta[i][j]
    }
  }

  return i_bits;
}

int length_truncated_unary(int symbol, int max_symbol)
{
  if (max_symbol == 0)
  {
    return 0;
  }

  bool code_last = (max_symbol > symbol);
  int num_bins = 0;
  while (symbol--)
  {
    num_bins++;
  }
  if (code_last)
  {
    num_bins++;
  }

  return num_bins;
}

double get_unfiltered_distortion_cov_channel(alf_covariance* cov, channel_type channel)
{
  double dist = 0;
  if (channel == CHANNEL_TYPE_LUMA)
  {
    dist = get_unfiltered_distortion_cov_classes(cov, MAX_NUM_ALF_CLASSES);
  }
  else
  {
    dist = get_unfiltered_distortion_cov_classes(cov, 1) + length_truncated_unary(0, 3) * 981.62883931057581/*m_lambda[COMPONENT_Cb]*/;
  }
  return dist;
}

double get_unfiltered_distortion_cov_classes(alf_covariance* cov, const int num_classes)
{
  double dist = 0;
  for (int class_idx = 0; class_idx < num_classes; class_idx++)
  {
    dist += cov[class_idx].pix_acc;
  }
  return dist;
}

void get_frame_stats(channel_type channel, int i_shape_idx, int m_ctuEnableFlag)
{
  int numClasses = channel == CHANNEL_TYPE_LUMA ? MAX_NUM_ALF_CLASSES : 1;
  int numCoeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  for (int i = 0; i < numClasses; i++)
  {
    //g_alf_covariance_frame[channel][iShapeIdx][i].reset();
    g_alf_covariance_frame[channel][i_shape_idx][i].pix_acc = 0;
    memset(g_alf_covariance_frame[channel][i_shape_idx][i].y, 0, sizeof(*g_alf_covariance_frame[channel][i_shape_idx][i].y) * numCoeff);
    for (int ii = 0; ii < numCoeff; ii++)
    {
      memset(g_alf_covariance_frame[channel][i_shape_idx][i].ee[ii], 0, sizeof(*g_alf_covariance_frame[channel][i_shape_idx][i].ee[ii]) * numCoeff);
    }
  }
  if (channel == CHANNEL_TYPE_LUMA)
  {
    get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_LUMA][i_shape_idx], g_alf_covariance[COMPONENT_Y][i_shape_idx], m_ctuEnableFlag, numClasses);
  }
  else
  {
    get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA][i_shape_idx], g_alf_covariance[COMPONENT_Cb][i_shape_idx], m_ctuEnableFlag, numClasses);
    get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA][i_shape_idx], g_alf_covariance[COMPONENT_Cr][i_shape_idx], m_ctuEnableFlag, numClasses);
  }
}

void get_frame_stat(alf_covariance* frame_cov, alf_covariance** ctb_cov, int ctb_enable_flags, const int num_classes)
{
  for (int i = 0; i < num_ctus_in_pic; i++)
  {
    if (ctb_enable_flags)
    {
      for (int j = 0; j < num_classes; j++)
      {
        //frameCov[j] += ctbCov[i][j];
        int numCoeff = frame_cov[j].num_coeff;
        for (int jj = 0; jj < numCoeff; jj++)
        {
          for (int ii = 0; ii < numCoeff; ii++)
          {
            frame_cov[j].ee[jj][ii] += ctb_cov[i][j].ee[jj][ii];
          }
          frame_cov[j].y[jj] += ctb_cov[i][j].y[jj];
        }
        frame_cov[j].pix_acc += ctb_cov[i][j].pix_acc;
      }
    }
  }
}

//-------------------------------------------------------------------

//-------------------------encoding functions------------------------

void kvz_alf_enc_process(encoder_state_t *const state,
  const lcu_order_element_t *const lcu)
{
  if (false/*cs.slice->getPendingRasInit()*/ || true/*cs.slice->isIDRorBLA()*/)
  {
    //memset(cs.slice->getAPSs(), 0, sizeof(*cs.slice->getAPSs())*MAX_NUM_APS);
    memset(state->slice->aps, 0, sizeof(alf_aps) * MAX_NUM_APS);
    
    //state->slice->aps = malloc(sizeof(alf_aps*) * MAX_NUM_APS);
    
    //m_apsMap->clearMap();
    //free(state->slice->param_set_map);
  }

  alf_aps* new_aps = &state->slice->aps[31];
  alf_aps aps;// = malloc(sizeof(alf_aps));
  int width = state->tile->frame->width;
  int height = state->tile->frame->height;
  //state->slice->aps = aps;

  //alfSliceParam.reset();
  /*
  memset(state->encoder_control->cfg.alf_slice_enable_flag, false, sizeof(state->encoder_control->cfg.alf_slice_enable_flag));
  memset(g_luma_coeff, 0, sizeof(g_luma_coeff));
  memset(g_chroma_coeff, 0, sizeof(g_chroma_coeff));
  memset(g_filter_coeff_delta_idx, 0, sizeof(g_filter_coeff_delta_idx));
  memset(g_alf_luma_coeff_flag, true, sizeof(g_alf_luma_coeff_flag));
  g_num_luma_filters = 1;
  g_alf_luma_coeff_delta_flag = false;
  g_alf_luma_coeff_delta_prediction_flag = false;
  g_fixed_filter_pattern = 0;
  memset(g_fixed_filter_idx, 0, sizeof(g_fixed_filter_idx));
  g_fixed_filter_set_index = 0;
  */
  memset(aps.enabled_flag, false, sizeof(aps.enabled_flag));
  memset(aps.luma_coeff, 0, sizeof(aps.luma_coeff));
  memset(aps.chroma_coeff, 0, sizeof(aps.chroma_coeff));
  memset(aps.filter_coeff_delta_idx, 0, sizeof(aps.filter_coeff_delta_idx));
  memset(aps.alf_luma_coeff_flag, true, sizeof(aps.alf_luma_coeff_flag));
  aps.num_luma_filters = 1;
  aps.alf_luma_coeff_delta_flag = false;
  aps.alf_luma_coeff_delta_prediction_flag = false;
  aps.t_layer = 0;
  memset(aps.new_filter_flag, 0, sizeof(aps.new_filter_flag));
  aps.fixed_filter_pattern = 0;
  memset(aps.fixed_filter_idx, 0, sizeof(aps.fixed_filter_idx));
  aps.fixed_filter_set_index = 0;

  //const TempCtx  ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));

  for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++)
  {
    //m_ctuEnableFlag[compIdx] = cs.picture->getAlfCtuEnableFlag(compIdx);
    memset(g_ctu_enable_flag[compIdx], 0, sizeof(uint8_t) * num_ctus_in_pic);    
  }

  
  int shiftLuma = 2 * 0 /*DISTORTION_PRECISION_ADJUSTMENT(m_inputBitDepth[CHANNEL_TYPE_LUMA])*/;
  int shiftChroma = 2 * 0  /*DISTORTION_PRECISION_ADJUSTMENT(m_inputBitDepth[CHANNEL_TYPE_CHROMA])*/;
  g_lambda[COMPONENT_Y] = 981.62883931057581/*lambdas[COMPONENT_Y]*/ * (1 << shiftLuma); // 981.62883931057581
  g_lambda[COMPONENT_Cb] = 981.62883931057581/*lambdas[COMPONENT_Cb]*/ * (1 << shiftChroma); // 981.62883931057581
  g_lambda[COMPONENT_Cr] = 981.62883931057581/*lambdas[COMPONENT_Cr]*/ * (1 << shiftChroma); // 981.62883931057581

  /*
  PelUnitBuf orgYuv = cs.getOrgBuf();
  m_tempBuf.copyFrom(cs.getRecoBuf());
  PelUnitBuf recYuv = m_tempBuf.getBuf(cs.area);
  recYuv.extendBorderPel(MAX_ALF_FILTER_LENGTH >> 1);
  */

  // derive classification
  //kutsutaan suurimmalla blokilla, mahdollisesti siis koko kuvalla
  kvz_alf_derive_classification(state, lcu, width, height, 0, 0);

  kvz_alf_reset_pcm_blk_class_info(state, lcu, width, height, 0, 0);

  // get CTB stats for filtering
  kvz_alf_derive_stats_for_filtering(state, lcu);
  // derive filter (luma)
  kvz_alf_encoder(state, lcu, &aps, CHANNEL_TYPE_LUMA);
  // derive filter (chroma)
  kvz_alf_encoder(state, lcu, &aps, CHANNEL_TYPE_CHROMA);
  //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
  kvz_alf_encoder_ctb(state, lcu, &aps);

  kvz_alf_reconstructor(state, lcu);
}

void kvz_alf_enc_create(encoder_state_t const *state,
  const lcu_order_element_t *lcu)
{
  kvz_alf_create(state, lcu);

  for (int channel_idx = 0; channel_idx < MAX_NUM_CHANNEL_TYPE; channel_idx++)
  {
    channel_type ch_type = (channel_type)channel_idx;
    int num_classes = channel_idx ? 1 : MAX_NUM_ALF_CLASSES;
    int num_coeffs = channel_idx ? 7 : 13;
    //m_alfCovarianceFrame[chType] = new AlfCovariance*[m_filterShapes[chType].size()];
    g_alf_covariance_frame[ch_type] = malloc(sizeof(**g_alf_covariance_frame[ch_type]));
    for (int i = 0; i != 1/*m_filterShapes[chType].size()*/; i++)
    {
      g_alf_covariance_frame[ch_type][i] = malloc(num_classes * sizeof(alf_covariance));
      for (int k = 0; k < num_classes; k++)
      {
        //m_alfCovarianceFrame[chType][i][k].create(m_filterShapes[chType][i].numCoeff);

        g_alf_covariance_frame[ch_type][i][k].num_coeff = num_coeffs;
        g_alf_covariance_frame[ch_type][i][k].y = malloc(num_coeffs * sizeof(*g_alf_covariance_frame[ch_type][i][k].y));
        g_alf_covariance_frame[ch_type][i][k].ee = malloc(num_coeffs * sizeof(**g_alf_covariance_frame[ch_type][i][k].ee));

        for (int ii = 0; ii < num_coeffs; ii++)
        {
          g_alf_covariance_frame[ch_type][i][k].ee[ii] = malloc(num_coeffs * sizeof(*g_alf_covariance_frame[ch_type][i][k].ee[ii]));
        }
      }
    }
  }

  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    g_ctu_enable_flag[comp_idx] = malloc(num_ctus_in_pic * sizeof(*g_ctu_enable_flag[comp_idx]));
    g_ctu_enable_flag_tmp[comp_idx] = malloc(num_ctus_in_pic * sizeof(*g_ctu_enable_flag_tmp[comp_idx]));

    int num_classes = comp_idx ? 1 : MAX_NUM_ALF_CLASSES;
    int num_coeffs = comp_idx ? 7 : 13;

    g_alf_covariance[comp_idx] = malloc(sizeof(***g_alf_covariance[comp_idx]));

    for (int i = 0; i != 1/*m_filterShapes[chType].size()*/; i++)
    {
      g_alf_covariance[comp_idx][i] = malloc(num_ctus_in_pic * sizeof(**g_alf_covariance[comp_idx][i]));
      for (int j = 0; j < num_ctus_in_pic; j++)
      {
        g_alf_covariance[comp_idx][i][j] = malloc(num_classes * sizeof(alf_covariance));
        for (int k = 0; k < num_classes; k++)
        {
          //g_alf_covariance[compIdx][i][j][k].create(m_filterShapes[chType][i].numCoeff);

          g_alf_covariance[comp_idx][i][j][k].num_coeff = num_coeffs;
          g_alf_covariance[comp_idx][i][j][k].y = malloc(num_coeffs * sizeof(*g_alf_covariance[comp_idx][i][j][k].y));
          g_alf_covariance[comp_idx][i][j][k].ee = malloc(num_coeffs * sizeof(**g_alf_covariance[comp_idx][i][j][k].ee));

          for (int ii = 0; ii < num_coeffs; ii++)
          {
            g_alf_covariance[comp_idx][i][j][k].ee[ii] = malloc(num_coeffs * sizeof(*g_alf_covariance[comp_idx][i][j][k].ee[ii]));
          }
        }
      }
    }
  }

  for (int i = 0; i != 1/*m_filterShapes[COMPONENT_Y].size()*/; i++)
  {
    for (int j = 0; j <= MAX_NUM_ALF_CLASSES + 1; j++)
    {
      //g_alf_covariance_merged[i][j].create(m_filterShapes[COMPONENT_Y][i].numCoeff);

      g_alf_covariance_merged[i][j].num_coeff = 13;

      g_alf_covariance_merged[i][j].y = malloc(13 * sizeof(*g_alf_covariance_merged[i][j].y));
      g_alf_covariance_merged[i][j].ee = malloc(13 * sizeof(**g_alf_covariance_merged[i][j].ee));

      for (int ii = 0; ii < 13; ii++)
      {
        g_alf_covariance_merged[i][j].ee[ii] = malloc(13 * sizeof(*g_alf_covariance_merged[i][j].ee[ii]));
      }
    }
  }
  g_filter_coeff_quant = malloc(MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
  g_filter_coeff_set = malloc(MAX_NUM_ALF_CLASSES * sizeof(int*));
  g_diff_filter_coeff = malloc(MAX_NUM_ALF_CLASSES * sizeof(int*));

  for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
  {
    g_filter_coeff_set[i] = malloc(MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
    g_diff_filter_coeff[i] = malloc(MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
  }
  g_aps_id_start = (int)MAX_NUM_APS;

  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    g_ctb_distortion_unfilter[comp] = malloc(num_ctus_in_pic * sizeof(double));
  }

  //g_alf_ctb_filter_set_index = malloc(num_ctus_in_pic * sizeof(short));
  g_alf_ctb_filter_index = malloc(num_ctus_in_pic * sizeof(*g_alf_ctb_filter_index));
  g_alf_ctb_filter_set_index_tmp = malloc(num_ctus_in_pic * sizeof(*g_alf_ctb_filter_set_index_tmp));
}


void kvz_alf_enc_destroy(encoder_state_t const *state)
{
  alf_info_t *alf = state->tile->frame->alf_info;
  int32_t pic_height = state->tile->frame->rec->height;

  if (!alf->created)
  {
    return;
  }

  for (int channel_idx = 0; channel_idx < MAX_NUM_CHANNEL_TYPE; channel_idx++)
  {
    if (g_alf_covariance_frame[channel_idx])
    {
      channel_type chType = channel_idx ? CHANNEL_TYPE_CHROMA : CHANNEL_TYPE_LUMA;
      int numClasses = channel_idx ? 1 : MAX_NUM_ALF_CLASSES;
      int num_coeff = channel_idx ? 7 : 13;
      for (int i = 0; i != 1/*m_filterShapes[chType].size()*/; i++)
      {
        for (int k = 0; k < numClasses; k++)
        {
          //g_alf_covariance_frame[channelIdx][i][k].destroy();
          for (int ii = 0; ii < num_coeff; ii++)
          {
            FREE_POINTER(g_alf_covariance_frame[channel_idx][i][k].ee[ii]);
          }

          FREE_POINTER(g_alf_covariance_frame[channel_idx][i][k].ee);
          FREE_POINTER(g_alf_covariance_frame[channel_idx][i][k].y);
        }
        FREE_POINTER(g_alf_covariance_frame[channel_idx][i]);
      }
      FREE_POINTER(g_alf_covariance_frame[channel_idx]);
    }
  }

  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    
    if (g_ctu_enable_flag && g_ctu_enable_flag[comp_idx])
    {
      FREE_POINTER(g_ctu_enable_flag[comp_idx]);
    }
    
    if (g_ctu_enable_flag_tmp && g_ctu_enable_flag_tmp[comp_idx] != NULL)
    {
      FREE_POINTER(g_ctu_enable_flag_tmp[comp_idx]);
    }

    if (g_alf_covariance[comp_idx])
    {
      channel_type chType = comp_idx ? CHANNEL_TYPE_CHROMA : CHANNEL_TYPE_LUMA;
      int numClasses = comp_idx ? 1 : MAX_NUM_ALF_CLASSES;

      for (int i = 0; i != 1/*m_filterShapes[chType].size()*/; i++)
      {
        for (int j = 0; j < num_ctus_in_pic; j++)
        {
          for (int k = 0; k < numClasses; k++)
          {
            //m_alfCovariance[compIdx][i][j][k].destroy(); 
            for (int ii = 0; ii < g_alf_covariance[comp_idx][i][j][k].num_coeff; ii++)
            {
              FREE_POINTER(g_alf_covariance[comp_idx][i][j][k].ee[ii]);
            }
            FREE_POINTER(g_alf_covariance[comp_idx][i][j][k].ee);
            FREE_POINTER(g_alf_covariance[comp_idx][i][j][k].y);
          }
          FREE_POINTER(g_alf_covariance[comp_idx][i][j]);
        }
        FREE_POINTER(g_alf_covariance[comp_idx][i]);
      }
      FREE_POINTER(g_alf_covariance[comp_idx]);
    }
  }

  for (int i = 0; i != 1/*m_filterShapes[COMPONENT_Y].size()*/; i++)
  {
    for (int j = 0; j <= MAX_NUM_ALF_CLASSES + 1; j++)
    {
      //m_alfCovarianceMerged[i][j].destroy();
      for (int ii = 0; ii < g_alf_covariance_merged[i][j].num_coeff; ii++)
      {
        FREE_POINTER(g_alf_covariance_merged[i][j].ee[ii]);
      }
      FREE_POINTER(g_alf_covariance_merged[i][j].ee);
      FREE_POINTER(g_alf_covariance_merged[i][j].y);
    }
  }

  if (g_filter_coeff_set)
  {
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
    {
      FREE_POINTER(g_filter_coeff_set[i]);
    }
    FREE_POINTER(g_filter_coeff_set);
  }

  if (g_diff_filter_coeff)
  {
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
    {
      FREE_POINTER(g_diff_filter_coeff[i]);
    }
    FREE_POINTER(g_diff_filter_coeff);
  }

  FREE_POINTER(g_filter_coeff_quant);

  /*
  if (g_ctb_distortion_fixed_filter != NULL) {
    FREE_POINTER(g_ctb_distortion_fixed_filter);
  }*/
  
  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    if (g_ctb_distortion_unfilter[comp] != NULL) {
      FREE_POINTER(g_ctb_distortion_unfilter[comp]);
    }
  }

  kvz_alf_destroy(state);

  /*
  for (int aps_idx = 0; aps_idx < MAX_NUM_APS; aps_idx++) {
    //FREE_POINTER(state->slice->param_set_map[aps_idx].p_nalu_data);
  }
  FREE_POINTER(state->slice->param_set_map);
  FREE_POINTER(state->slice->aps);
  FREE_POINTER(state->slice->tile_group_luma_aps_id);
  */
}

void kvz_alf_encoder(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  alf_aps *aps,
  channel_type channel)//CodingStructure& cs, AlfSliceParam& alfSliceParam, const PelUnitBuf& orgUnitBuf, const PelUnitBuf& recExtBuf, const PelUnitBuf& recBuf, const ChannelType channel)
{
  //const TempCtx  ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));
  //TempCtx        ctxBest(m_CtxCache);
  //cabac_data_t cabac = state->cabac;
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  kvz_config cfg = state->encoder_control->cfg;

  double cost_min = MAX_DOUBLE;
  
  g_bits_new_filter[channel] = 0;

  const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
  int ui_coeff_bits = 0;

  for (int i_shape_idx = 0; i_shape_idx < 1/*alfFilterShape.size()*/; i_shape_idx++)
  { 
    //m_alfSliceParamTemp = alfSliceParam;
    memcpy(g_alf_slice_aps_temp.enabled_flag, aps->enabled_flag, sizeof(g_alf_slice_aps_temp.enabled_flag));
    memcpy(g_alf_slice_aps_temp.luma_coeff, aps->luma_coeff, sizeof(g_alf_slice_aps_temp.luma_coeff));
    memcpy(g_alf_slice_aps_temp.chroma_coeff, aps->chroma_coeff, sizeof(g_alf_slice_aps_temp.chroma_coeff));
    memcpy(g_alf_slice_aps_temp.filter_coeff_delta_idx, aps->filter_coeff_delta_idx, sizeof(g_alf_slice_aps_temp.filter_coeff_delta_idx));
    memcpy(g_alf_slice_aps_temp.alf_luma_coeff_flag, aps->alf_luma_coeff_flag, sizeof(g_alf_slice_aps_temp.alf_luma_coeff_flag));
    g_alf_slice_aps_temp.num_luma_filters = aps->num_luma_filters;
    g_alf_slice_aps_temp.alf_luma_coeff_delta_flag = aps->alf_luma_coeff_delta_flag;
    g_alf_slice_aps_temp.alf_luma_coeff_delta_prediction_flag = aps->alf_luma_coeff_delta_prediction_flag;
    g_alf_slice_aps_temp.t_layer = aps->t_layer;
    memcpy(g_alf_slice_aps_temp.new_filter_flag, aps->new_filter_flag, sizeof(g_alf_slice_aps_temp.new_filter_flag));
    g_alf_slice_aps_temp.fixed_filter_pattern = aps->fixed_filter_pattern;
    memcpy(g_alf_slice_aps_temp.fixed_filter_idx, aps->fixed_filter_idx, sizeof(g_alf_slice_aps_temp.fixed_filter_idx));
    g_alf_slice_aps_temp.fixed_filter_set_index = aps->fixed_filter_set_index;

    //1. get unfiltered distortion
    double cost = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[channel][i_shape_idx], channel);
    cost /= 1.001; // slight preference for unfiltered choice

    if (cost < cost_min)
    {
      cost_min = cost;
      //setEnableFlag(alfSliceParam, channel, false);
      if (is_luma) {
        aps->enabled_flag[COMPONENT_Y] = 0;
      }
      else {
        aps->enabled_flag[COMPONENT_Cb] = 0;
        aps->enabled_flag[COMPONENT_Cr] = 0;
      }
      // no CABAC signalling
      //ctxBest = AlfCtx(ctxStart);
      //setCtuEnableFlag(m_ctuEnableFlagTmp, channel, 0);
      if (is_luma) {
        //g_ctu_enable_flag_tmp[COMPONENT_Y] = 0;
        memset(g_ctu_enable_flag_tmp[COMPONENT_Y], 0, sizeof(uint8_t) * num_ctus_in_pic);
      }
      else {
        //g_ctu_enable_flag_tmp[COMPONENT_Cb] = 0;
        //g_ctu_enable_flag_tmp[COMPONENT_Cr] = 0;
        memset(g_ctu_enable_flag_tmp[COMPONENT_Cb], 0, sizeof(uint8_t) * num_ctus_in_pic);
        memset(g_ctu_enable_flag_tmp[COMPONENT_Cr], 0, sizeof(uint8_t) * num_ctus_in_pic);
      }
    }

    //2. all CTUs are on
    //setEnableFlag(m_alfSliceParamTemp, channel, true);
    if (is_luma) {
      g_alf_slice_aps_temp.enabled_flag[COMPONENT_Y] = 1;
    }
    else {
      g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cb] = 1;
      g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cr] = 1;
    }

    //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
    //setCtuEnableFlag(m_ctuEnableFlag, channel, 1);
    if (is_luma) {
      memset(g_ctu_enable_flag[COMPONENT_Y], 1, sizeof(uint8_t) * num_ctus_in_pic);
    }
    else {
      memset(g_ctu_enable_flag[COMPONENT_Cb], 1, sizeof(uint8_t) * num_ctus_in_pic);
      memset(g_ctu_enable_flag[COMPONENT_Cr], 1, sizeof(uint8_t) * num_ctus_in_pic);
    }
    cost = kvz_alf_get_filter_coeff_and_cost(state, lcu, channel, 0, &ui_coeff_bits, i_shape_idx, false, false);

    if (cost < cost_min)
    {
      g_bits_new_filter[channel] = ui_coeff_bits;

      cost_min = cost;
      //copyAlfSliceParam(alfSliceParam, m_alfSliceParamTemp, channel);
      if (is_luma)
      {
        //memcpy(aps, &g_alf_slice_aps_temp, sizeof(alf_aps));
        memcpy(aps->enabled_flag, g_alf_slice_aps_temp.enabled_flag, sizeof(aps->enabled_flag));
        memcpy(aps->luma_coeff, g_alf_slice_aps_temp.luma_coeff, sizeof(aps->luma_coeff));
        memcpy(aps->chroma_coeff, g_alf_slice_aps_temp.chroma_coeff, sizeof(aps->chroma_coeff));
        memcpy(aps->filter_coeff_delta_idx, g_alf_slice_aps_temp.filter_coeff_delta_idx, sizeof(aps->filter_coeff_delta_idx));
        memcpy(aps->alf_luma_coeff_flag, g_alf_slice_aps_temp.alf_luma_coeff_flag, sizeof(aps->alf_luma_coeff_flag));
        aps->num_luma_filters = g_alf_slice_aps_temp.num_luma_filters;
        aps->alf_luma_coeff_delta_flag = g_alf_slice_aps_temp.alf_luma_coeff_delta_flag;
        aps->alf_luma_coeff_delta_prediction_flag = g_alf_slice_aps_temp.alf_luma_coeff_delta_prediction_flag;
        aps->t_layer = g_alf_slice_aps_temp.t_layer;
        memcpy(aps->new_filter_flag, g_alf_slice_aps_temp.new_filter_flag, sizeof(aps->new_filter_flag));
        aps->fixed_filter_pattern = g_alf_slice_aps_temp.fixed_filter_pattern;
        memcpy(aps->fixed_filter_idx, g_alf_slice_aps_temp.fixed_filter_idx, sizeof(aps->fixed_filter_idx));
        aps->fixed_filter_set_index = g_alf_slice_aps_temp.fixed_filter_set_index;
      }
      else
      {
        aps->enabled_flag[COMPONENT_Cb] = g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cb];
        aps->enabled_flag[COMPONENT_Cr] = g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cr];
        memcpy(aps->chroma_coeff, g_alf_slice_aps_temp.chroma_coeff, sizeof(aps->chroma_coeff));
      }

      //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
      //setCtuEnableFlag(m_ctuEnableFlagTmp, channel, 1);
      if (is_luma) {
        memset(g_ctu_enable_flag_tmp[COMPONENT_Y], 1, sizeof(uint8_t) * num_ctus_in_pic);
      }
      else {
        memset(g_ctu_enable_flag_tmp[COMPONENT_Cb], 1, sizeof(uint8_t) * num_ctus_in_pic);
        memset(g_ctu_enable_flag_tmp[COMPONENT_Cr], 1, sizeof(uint8_t) * num_ctus_in_pic);
      }
    }

    //3. CTU decision
    double dist_unfilter = 0;
    const int iter_num = is_luma ? (2 * 4 + 1) : (2 * 2 + 1);

    for (int iter = 0; iter < iter_num; iter++)
    {
      if ((iter & 0x01) == 0)
      {
        //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
        cost = g_lambda[channel] * ui_coeff_bits;
        cost += kvz_alf_derive_ctb_alf_enable_flags(state, lcu, channel, i_shape_idx, &dist_unfilter, num_classes);
        if (cost < cost_min)
        {
          g_bits_new_filter[channel] = ui_coeff_bits;

          cost_min = cost;
          //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
          //copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, channel);
          if (is_luma) {
            memcpy(g_ctu_enable_flag_tmp[COMPONENT_Y], g_ctu_enable_flag[COMPONENT_Y], sizeof(uint8_t) * num_ctus_in_pic);
          }
          else {
            memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cr], g_ctu_enable_flag[COMPONENT_Cr], sizeof(uint8_t) * num_ctus_in_pic);
            memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cb], g_ctu_enable_flag[COMPONENT_Cb], sizeof(uint8_t) * num_ctus_in_pic);
          }

          //copyAlfSliceParam(alfSliceParam, m_alfSliceParamTemp, channel);
          if (is_luma)
          {
            //memcpy(aps, &g_alf_slice_aps_temp, sizeof(alf_aps));
            memcpy(aps->enabled_flag, g_alf_slice_aps_temp.enabled_flag, sizeof(aps->enabled_flag));
            memcpy(aps->luma_coeff, g_alf_slice_aps_temp.luma_coeff, sizeof(aps->luma_coeff));
            memcpy(aps->chroma_coeff, g_alf_slice_aps_temp.chroma_coeff, sizeof(aps->chroma_coeff));
            memcpy(aps->filter_coeff_delta_idx, g_alf_slice_aps_temp.filter_coeff_delta_idx, sizeof(aps->filter_coeff_delta_idx));
            memcpy(aps->alf_luma_coeff_flag, g_alf_slice_aps_temp.alf_luma_coeff_flag, sizeof(aps->alf_luma_coeff_flag));
            aps->num_luma_filters = g_alf_slice_aps_temp.num_luma_filters;
            aps->alf_luma_coeff_delta_flag = g_alf_slice_aps_temp.alf_luma_coeff_delta_flag;
            aps->alf_luma_coeff_delta_prediction_flag = g_alf_slice_aps_temp.alf_luma_coeff_delta_prediction_flag;
            aps->t_layer = g_alf_slice_aps_temp.t_layer;
            memcpy(aps->new_filter_flag, g_alf_slice_aps_temp.new_filter_flag, sizeof(aps->new_filter_flag));
            aps->fixed_filter_pattern = g_alf_slice_aps_temp.fixed_filter_pattern;
            memcpy(aps->fixed_filter_idx, g_alf_slice_aps_temp.fixed_filter_idx, sizeof(aps->fixed_filter_idx));
            aps->fixed_filter_set_index = g_alf_slice_aps_temp.fixed_filter_set_index;
          }
          else
          {
            aps->enabled_flag[COMPONENT_Cb] = g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cb];
            aps->enabled_flag[COMPONENT_Cr] = g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cr];
            memcpy(aps->chroma_coeff, g_alf_slice_aps_temp.chroma_coeff, sizeof(aps->chroma_coeff));
          }

        }
      }
      else
      {
        // unfiltered distortion is added due to some CTBs may not use filter
        //cost = kvz_alf_get_filter_coeff_and_cost(cs, dist_unfilter, channel, true, i_shape_idx, ui_coeff_bits);
        cost = kvz_alf_get_filter_coeff_and_cost(state, lcu, channel, dist_unfilter, &ui_coeff_bits, i_shape_idx, true, false);
      }
    }//for iter
  }//for shapeIdx
  //m_CABACEstimator->getCtx() = AlfCtx(ctxBest);
}

void kvz_alf_encoder_ctb(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  alf_aps *aps)
{
  /*
  TempCtx        ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));
  TempCtx        ctxBest(m_CtxCache);
  TempCtx        ctxTempStart(m_CtxCache);
  TempCtx        ctxTempBest(m_CtxCache);
  AlfSliceParam  alfSliceParamNewFiltersBest = alfSliceParamNewFilters;
  APS**          apss = cs.slice->getAPSs();
  */

  //AlfSliceParam  alfSliceParamNewFiltersBest = alfSliceParamNewFilters;
  alf_aps alf_slice_param_new_filters_best;
  memcpy(alf_slice_param_new_filters_best.enabled_flag, aps->enabled_flag, sizeof(alf_slice_param_new_filters_best.enabled_flag));
  memcpy(alf_slice_param_new_filters_best.luma_coeff, aps->luma_coeff, sizeof(alf_slice_param_new_filters_best.luma_coeff));
  memcpy(alf_slice_param_new_filters_best.chroma_coeff, aps->chroma_coeff, sizeof(alf_slice_param_new_filters_best.chroma_coeff));
  memcpy(alf_slice_param_new_filters_best.filter_coeff_delta_idx, aps->filter_coeff_delta_idx, sizeof(alf_slice_param_new_filters_best.filter_coeff_delta_idx));
  memcpy(alf_slice_param_new_filters_best.alf_luma_coeff_flag, aps->alf_luma_coeff_flag, sizeof(alf_slice_param_new_filters_best.alf_luma_coeff_flag));
  alf_slice_param_new_filters_best.num_luma_filters = aps->num_luma_filters;
  alf_slice_param_new_filters_best.alf_luma_coeff_delta_flag = aps->alf_luma_coeff_delta_flag;
  alf_slice_param_new_filters_best.alf_luma_coeff_delta_prediction_flag = aps->alf_luma_coeff_delta_prediction_flag;
  alf_slice_param_new_filters_best.t_layer = aps->t_layer;
  memcpy(alf_slice_param_new_filters_best.new_filter_flag, aps->new_filter_flag, sizeof(alf_slice_param_new_filters_best.new_filter_flag));
  alf_slice_param_new_filters_best.fixed_filter_pattern = aps->fixed_filter_pattern;
  memcpy(alf_slice_param_new_filters_best.fixed_filter_idx, aps->fixed_filter_idx, sizeof(alf_slice_param_new_filters_best.fixed_filter_idx));
  alf_slice_param_new_filters_best.fixed_filter_set_index = aps->fixed_filter_set_index;
  
  alf_component_id component_id = state->tile->frame->alf_info->component_id;
  bool is_luma = component_id == COMPONENT_Y ? true : false;
  kvz_config cfg = state->encoder_control->cfg;
  //int32_t alf_ctu_flags[MAX_NUM_COMPONENT] = { cfg.alf_ctu_enable_flag };

  bool has_new_filters[2] = { aps->enabled_flag[COMPONENT_Y] , aps->enabled_flag[COMPONENT_Cb] || aps->enabled_flag[COMPONENT_Cr] };
  //initDistortion();
  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    for (int ctb_idx = 0; ctb_idx < num_ctus_in_pic; ctb_idx++)
    {
      g_ctb_distortion_unfilter[comp][ctb_idx] = get_unfiltered_distortion_cov_classes(g_alf_covariance[comp][0][ctb_idx], comp == 0 ? MAX_NUM_ALF_CLASSES : 1);
    }
  }

  //luma
  get_frame_stats(CHANNEL_TYPE_LUMA, 0, 1);
  double cost_off = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[CHANNEL_TYPE_LUMA][0], CHANNEL_TYPE_LUMA);
  //std::vector<int> apsIds = getAvaiApsIdsLuma(cs, newApsId);
  //int new_aps_id = 31;
  int new_aps_id;
  int aps_ids[6];
  int size_of_aps_ids = 0;

  kvz_alf_get_avai_aps_ids_luma(state, &new_aps_id, &aps_ids[6], &size_of_aps_ids);

  int best_aps_ids[7] = {-1, -1, -1, -1, -1, -1, -1};
  double cost_min = MAX_DOUBLE;
  kvz_alf_reconstruct_coeff_aps(state, true, false, true);

  int num_loops = has_new_filters[CHANNEL_TYPE_LUMA] ? 2 : 1;
  for (int use_new_filter = 0; use_new_filter < num_loops; use_new_filter++)
  {
    int bits_new_filter = 0;
    if (use_new_filter == 1)
    {
      if (!has_new_filters[CHANNEL_TYPE_LUMA])
      {
        continue;
      }
      else
      {
        bits_new_filter = g_bits_new_filter[CHANNEL_TYPE_LUMA];
        kvz_alf_reconstruct_coeff(state, aps, CHANNEL_TYPE_LUMA, true, true);
      }
    }
    int num_iter = use_new_filter ? 2 : 1;
    for (int num_temporal_aps = 0; num_temporal_aps <= 0/*apsIds.size()*/; num_temporal_aps++)
    {
      //cs.slice->setTileGroupNumAps(numTemporalAps + useNewFilter);
      //memset(state->slice->tile_group_num_aps, num_temporal_aps + use_new_filter, sizeof(int));
      state->slice->tile_group_num_aps = num_temporal_aps + use_new_filter;

      int num_filter_set = ALF_NUM_FIXED_FILTER_SETS + num_temporal_aps + use_new_filter;
      if (num_temporal_aps == size_of_aps_ids && num_temporal_aps > 0 && use_new_filter && new_aps_id == aps_ids[size_of_aps_ids-1] /*apsIds.back()*/) //last temporalAPS is occupied by new filter set and this temporal APS becomes unavailable
      {
        continue;
      }
      for (int iter = 0; iter < num_iter; iter++)
      {
        //g_alf_slice_aps_temp = aps;
        memcpy(g_alf_slice_aps_temp.enabled_flag, aps->enabled_flag, sizeof(g_alf_slice_aps_temp.enabled_flag));
        memcpy(g_alf_slice_aps_temp.luma_coeff, aps->luma_coeff, sizeof(g_alf_slice_aps_temp.luma_coeff));
        memcpy(g_alf_slice_aps_temp.chroma_coeff, aps->chroma_coeff, sizeof(g_alf_slice_aps_temp.chroma_coeff));
        memcpy(g_alf_slice_aps_temp.filter_coeff_delta_idx, aps->filter_coeff_delta_idx, sizeof(g_alf_slice_aps_temp.filter_coeff_delta_idx));
        memcpy(g_alf_slice_aps_temp.alf_luma_coeff_flag, aps->alf_luma_coeff_flag, sizeof(g_alf_slice_aps_temp.alf_luma_coeff_flag));
        g_alf_slice_aps_temp.num_luma_filters = aps->num_luma_filters;
        g_alf_slice_aps_temp.alf_luma_coeff_delta_flag = aps->alf_luma_coeff_delta_flag;
        g_alf_slice_aps_temp.alf_luma_coeff_delta_prediction_flag = aps->alf_luma_coeff_delta_prediction_flag;
        g_alf_slice_aps_temp.t_layer = aps->t_layer;
        memcpy(g_alf_slice_aps_temp.new_filter_flag, aps->new_filter_flag, sizeof(g_alf_slice_aps_temp.new_filter_flag));
        g_alf_slice_aps_temp.fixed_filter_pattern = aps->fixed_filter_pattern;
        memcpy(g_alf_slice_aps_temp.fixed_filter_idx, aps->fixed_filter_idx, sizeof(g_alf_slice_aps_temp.fixed_filter_idx));
        g_alf_slice_aps_temp.fixed_filter_set_index = aps->fixed_filter_set_index;

        g_alf_slice_aps_temp.enabled_flag[CHANNEL_TYPE_LUMA] = true;
        double cur_cost = get_tb_length(num_temporal_aps + use_new_filter, 6/*ALF_CTB_MAX_NUM_APS*/ + 1) * 981.62883931057581/*m_lambda[CHANNEL_TYPE_LUMA]*/;
        if (iter > 0)  //re-derive new filter-set
        {
          double d_dist_org_new_filter = 0;
          int blocks_using_new_filter = 0;
          for (int ctb_idx = 0; ctb_idx < num_ctus_in_pic; ctb_idx++)
          {
            if (g_ctu_enable_flag[COMPONENT_Y][ctb_idx] && g_alf_ctb_filter_index[ctb_idx] != ALF_NUM_FIXED_FILTER_SETS)
            {
              g_ctu_enable_flag[COMPONENT_Y][ctb_idx] = 0;
            }
            else if (g_ctu_enable_flag[COMPONENT_Y][ctb_idx] && g_alf_ctb_filter_index[ctb_idx] == ALF_NUM_FIXED_FILTER_SETS)
            {
              blocks_using_new_filter++;
              d_dist_org_new_filter += g_ctb_distortion_unfilter[COMPONENT_Y][ctb_idx];
              for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
              {
                short* p_coeff = g_coeff_final;
                for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
                {
                  g_filter_tmp[i] = p_coeff[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                }
                d_dist_org_new_filter += calc_error_for_coeffs(g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx].ee, g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx].y, g_filter_tmp, MAX_NUM_ALF_LUMA_COEFF, ALF_NUM_BITS);
              }
            }
          }
          if (blocks_using_new_filter > 0 && blocks_using_new_filter < num_ctus_in_pic)
          {
            int bit_nl[2] = { 0, 0 };
            double err_nl[2] = { 0.0, 0.0 };
            err_nl[1] = MAX_DOUBLE;

            //errNL[0] = getFilterCoeffAndCost(cs, 0, CHANNEL_TYPE_LUMA, true, 0, bitNL[0], true);
            err_nl[0] = kvz_alf_get_filter_coeff_and_cost(state, lcu, CHANNEL_TYPE_LUMA, 0, &bit_nl[0], 0, true, true);

            int bits_new_filter_temp_luma = bit_nl[0];
            double err = err_nl[0];
            if (err_nl[1]  < err_nl[0])
            {
              err = err_nl[1];
              bits_new_filter_temp_luma = bit_nl[1];
              //g_alf_slice_aps_temp = m_alfSliceParamTempNL;
            }
            if (d_dist_org_new_filter + 981.62883931057581/*m_lambda[CHANNEL_TYPE_LUMA]*/ * g_bits_new_filter[CHANNEL_TYPE_LUMA] < err) //re-derived filter is not good, skip
            {
              continue;
            }
            kvz_alf_reconstruct_coeff(state, &g_alf_slice_aps_temp, CHANNEL_TYPE_LUMA, true, true);
            bits_new_filter = bits_new_filter_temp_luma;
          }
          else //no blocks using new filter, skip
          {
            continue;
          }
        }

        //m_CABACEstimator->getCtx() = ctxStart;
        for (int ctb_idx = 0; ctb_idx < num_ctus_in_pic; ctb_idx++)
        {
          double dist_unfilter_ctb = g_ctb_distortion_unfilter[COMPONENT_Y][ctb_idx];
          //ctb on
          g_ctu_enable_flag[COMPONENT_Y][ctb_idx] = 1;
          double cost_on = MAX_DOUBLE;
          //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
          int i_best_filter_set_idx = 0;
          for (int filter_set_idx = 0; filter_set_idx < num_filter_set; filter_set_idx++)
          {
            //rate
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
            //m_CABACEstimator->resetBits();
            //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, COMPONENT_Y, &m_alfSliceParamTemp);
            //17
            //7 nyt
            code_alf_ctu_enable_flag(state, lcu, ctb_idx, COMPONENT_Y, &g_alf_slice_aps_temp);
            g_alf_ctb_filter_index[ctb_idx] = filter_set_idx;

            code_alf_ctu_filter_index(state, lcu, ctb_idx, g_alf_slice_aps_temp.enabled_flag[COMPONENT_Y]);

            double rate_on = 4.0255737304687500/*FracBitsScale *(double)m_CABACEstimator->getEstFracBits()*/;
            //distortion
            double dist = dist_unfilter_ctb;
            for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
            {
              if (filter_set_idx < ALF_NUM_FIXED_FILTER_SETS)
              {
                int filter_idx = g_class_to_filter_mapping[filter_set_idx][class_idx];
                dist += calc_error_for_coeffs(g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx].ee, g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx].y, g_fixed_filter_set_coeff[filter_idx], MAX_NUM_ALF_LUMA_COEFF, ALF_NUM_BITS);
              }
              else
              {
                short *p_coeff;
                if (use_new_filter && filter_set_idx == ALF_NUM_FIXED_FILTER_SETS)
                {
                  p_coeff = g_coeff_final;
                }
                else if (use_new_filter)
                {
                  p_coeff = g_coeff_aps_luma[filter_set_idx - 1 - ALF_NUM_FIXED_FILTER_SETS];
                }
                else
                {
                  p_coeff = g_coeff_aps_luma[filter_set_idx - ALF_NUM_FIXED_FILTER_SETS];
                }
                for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
                {
                  g_filter_tmp[i] = p_coeff[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                }
                dist += calc_error_for_coeffs(g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx].ee, g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx].y, g_filter_tmp, MAX_NUM_ALF_LUMA_COEFF, ALF_NUM_BITS);
              }
            }
            //cost
            double cost_on_tmp = dist + 981.62883931057581/*m_lambda[COMPONENT_Y]*/ * rate_on;
            if (cost_on_tmp < cost_on)
            {
              //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
              cost_on = cost_on_tmp;
              i_best_filter_set_idx = filter_set_idx;
            }
          }
          //ctb off
          g_ctu_enable_flag[COMPONENT_Y][ctb_idx] = 0;
          //rate
          //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
          //m_CABACEstimator->resetBits();

          //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, COMPONENT_Y, &m_alfSliceParamTemp);
          code_alf_ctu_enable_flag(state, lcu, ctb_idx, COMPONENT_Y, &g_alf_slice_aps_temp);

          //cost
          double cost_off = dist_unfilter_ctb + (1243132.9728229337-1237410.0000000000)/*m_lambda[COMPONENT_Y] * FracBitsScale*(double)m_CABACEstimator->getEstFracBits()*/;
          if (cost_on < cost_off)
          {
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
            g_ctu_enable_flag[COMPONENT_Y][ctb_idx] = 1;
            g_alf_ctb_filter_index[ctb_idx] = i_best_filter_set_idx;
            cur_cost += cost_on;
          }
          else
          {
            g_ctu_enable_flag[COMPONENT_Y][ctb_idx] = 0;
            cur_cost += cost_off;
          }
        } //for(ctbIdx)
        int tmp_bits = bits_new_filter + 5 * (num_filter_set - ALF_NUM_FIXED_FILTER_SETS) + (1/*cs.slice->isIntra()*/ ? 1 : get_tb_length(num_filter_set - ALF_NUM_FIXED_FILTER_SETS, 6/*ALF_CTB_MAX_NUM_APS*/ + 1));
        cur_cost += tmp_bits * 981.62883931057581/*m_lambda[COMPONENT_Y]*/;
        if (cur_cost < cost_min)
        {
          cost_min = cur_cost;
          //bestApsIds.resize(numFilterSet - alf_num_fixed_filter_sets);
          for (int i = 0; i < num_filter_set - ALF_NUM_FIXED_FILTER_SETS/*bestApsIds.size()*/; i++)
          {
            if (i == 0 && use_new_filter)
            {
              best_aps_ids[i] = new_aps_id;
            }
            else
            {
              best_aps_ids[i] = aps_ids[i - use_new_filter];
            }
          }
          //alfSliceParamNewFiltersBest = m_alfSliceParamTemp;
          memcpy(alf_slice_param_new_filters_best.enabled_flag, g_alf_slice_aps_temp.enabled_flag, sizeof(alf_slice_param_new_filters_best.enabled_flag));
          memcpy(alf_slice_param_new_filters_best.luma_coeff, g_alf_slice_aps_temp.luma_coeff, sizeof(alf_slice_param_new_filters_best.luma_coeff));
          memcpy(alf_slice_param_new_filters_best.chroma_coeff, g_alf_slice_aps_temp.chroma_coeff, sizeof(alf_slice_param_new_filters_best.chroma_coeff));
          memcpy(alf_slice_param_new_filters_best.filter_coeff_delta_idx, g_alf_slice_aps_temp.filter_coeff_delta_idx, sizeof(alf_slice_param_new_filters_best.filter_coeff_delta_idx));
          memcpy(alf_slice_param_new_filters_best.alf_luma_coeff_flag, g_alf_slice_aps_temp.alf_luma_coeff_flag, sizeof(alf_slice_param_new_filters_best.alf_luma_coeff_flag));
          alf_slice_param_new_filters_best.num_luma_filters = g_alf_slice_aps_temp.num_luma_filters;
          alf_slice_param_new_filters_best.alf_luma_coeff_delta_flag = g_alf_slice_aps_temp.alf_luma_coeff_delta_flag;
          alf_slice_param_new_filters_best.alf_luma_coeff_delta_prediction_flag = g_alf_slice_aps_temp.alf_luma_coeff_delta_prediction_flag;
          alf_slice_param_new_filters_best.t_layer = g_alf_slice_aps_temp.t_layer;
          memcpy(alf_slice_param_new_filters_best.new_filter_flag, g_alf_slice_aps_temp.new_filter_flag, sizeof(alf_slice_param_new_filters_best.new_filter_flag));
          alf_slice_param_new_filters_best.fixed_filter_pattern = g_alf_slice_aps_temp.fixed_filter_pattern;
          memcpy(alf_slice_param_new_filters_best.fixed_filter_idx, g_alf_slice_aps_temp.fixed_filter_idx, sizeof(alf_slice_param_new_filters_best.fixed_filter_idx));
          alf_slice_param_new_filters_best.fixed_filter_set_index = g_alf_slice_aps_temp.fixed_filter_set_index;

          //ctxBest = AlfCtx(m_CABACEstimator->getCtx());

          //copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, CHANNEL_TYPE_LUMA);
          memcpy(g_ctu_enable_flag_tmp[COMPONENT_Y], g_ctu_enable_flag[COMPONENT_Y], sizeof(uint8_t) * num_ctus_in_pic);

          for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
          {
            g_alf_ctb_filter_set_index_tmp[ctu_idx] = g_alf_ctb_filter_index[ctu_idx];
          }
          alf_slice_param_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA] = use_new_filter;
        }
      }//for (int iter = 0; iter < numIter; iter++)
    }// for (int numTemporalAps = 0; numTemporalAps < apsIds.size(); numTemporalAps++)
  }//for (int useNewFilter = 0; useNewFilter <= 1; useNewFilter++)

  if (cost_off <= cost_min)
  {
    //cs.slice->resetTileGroupAlfEnabledFlag();
    memset(state->slice->tile_group_alf_enabled_flag, 0, sizeof(state->slice->tile_group_alf_enabled_flag));

    //cs.slice->setTileGroupNumAps(0);
    state->slice->tile_group_num_aps = 0;

    //setCtuEnableFlag(m_ctuEnableFlag, CHANNEL_TYPE_LUMA, 0);
    memset(g_ctu_enable_flag[COMPONENT_Y], 0, sizeof(uint8_t) * num_ctus_in_pic);

    //setCtuEnableFlag(m_ctuEnableFlag, CHANNEL_TYPE_CHROMA, 0);
    memset(g_ctu_enable_flag[COMPONENT_Cb], 0, sizeof(uint8_t) * num_ctus_in_pic);
    memset(g_ctu_enable_flag[COMPONENT_Cr], 0, sizeof(uint8_t) * num_ctus_in_pic);

    return;
  }
  else
  {
    //alfSliceParamNewFiltersBest.tLayer = cs.slice->getTLayer();
    alf_slice_param_new_filters_best.t_layer = state->slice->id;
    //cs.slice->setTileGroupAlfEnabledFlag(COMPONENT_Y, true);
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Y] = 0;
    //cs.slice->setTileGroupNumAps((int)bestApsIds.size());
    int size_of_best_aps_ids = 0;
    for (int i = 0; i < 7; i++) {
      if (best_aps_ids[i] != -1) {
        size_of_best_aps_ids++;
      }
    }
    state->slice->tile_group_num_aps = size_of_best_aps_ids;
    //cs.slice->setAPSs(bestApsIds);
    //state->slice->tile_group_luma_aps_id = malloc(size_of_best_aps_ids * sizeof(int*));
    for (int i = 0; i < size_of_best_aps_ids/*state->slice->tile_group_num_aps*/; i++)
    {
      state->slice->tile_group_luma_aps_id[i] = best_aps_ids[i];
    }

    //copyCtuEnableFlag(m_ctuEnableFlag, m_ctuEnableFlagTmp, CHANNEL_TYPE_LUMA);
    memcpy(g_ctu_enable_flag[COMPONENT_Y], g_ctu_enable_flag_tmp[COMPONENT_Y], sizeof(uint8_t) * num_ctus_in_pic);
    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
    {
      g_alf_ctb_filter_index[ctu_idx] = g_alf_ctb_filter_set_index_tmp[ctu_idx];
    }

    if (alf_slice_param_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA])
    {
      //alf_aps* newAPS = &state->slice->aps[new_aps_id];//apss[new_aps_id]; //m_apsMap->getPS(newApsId);
      alf_aps* new_aps = &state->slice->param_set_map[new_aps_id].parameter_set;
      if (new_aps == NULL)
      {
        //newAPS = m_apsMap->allocatePS(newApsId);
        assert(new_aps_id < MAX_NUM_APS); //Invalid PS id
        bool found = false;
        for (int i = 0; i < (sizeof(state->slice->param_set_map) / sizeof(state->slice->param_set_map[0])); i++) {
          if (state->slice->param_set_map[i].parameter_set.aps_id == new_aps_id) {
            new_aps_id = true;
          }
        }
        if (!found) {
          state->slice->param_set_map[new_aps_id].b_changed = true;
          //state->slice->param_set_map[new_aps_id].p_nalu_data = 0;
          //state->slice->param_set_map[new_aps_id].parameter_set = malloc(sizeof(alf_aps));
          state->slice->param_set_map[new_aps_id].parameter_set.aps_id = new_aps_id;
        }
        new_aps->aps_id = state->slice->param_set_map[new_aps_id].parameter_set.aps_id;

      }
      //newAPS->setAlfAPSParam(alfSliceParamNewFiltersBest);
      memcpy(new_aps->enabled_flag, alf_slice_param_new_filters_best.enabled_flag, sizeof(new_aps->enabled_flag));
      memcpy(new_aps->luma_coeff, alf_slice_param_new_filters_best.luma_coeff, sizeof(new_aps->luma_coeff));
      memcpy(new_aps->chroma_coeff, alf_slice_param_new_filters_best.chroma_coeff, sizeof(new_aps->chroma_coeff));
      memcpy(new_aps->filter_coeff_delta_idx, alf_slice_param_new_filters_best.filter_coeff_delta_idx, sizeof(new_aps->filter_coeff_delta_idx));
      memcpy(new_aps->alf_luma_coeff_flag, alf_slice_param_new_filters_best.alf_luma_coeff_flag, sizeof(new_aps->alf_luma_coeff_flag));
      new_aps->num_luma_filters = alf_slice_param_new_filters_best.num_luma_filters;
      new_aps->alf_luma_coeff_delta_flag = alf_slice_param_new_filters_best.alf_luma_coeff_delta_flag;
      new_aps->alf_luma_coeff_delta_prediction_flag = alf_slice_param_new_filters_best.alf_luma_coeff_delta_prediction_flag;
      new_aps->t_layer = alf_slice_param_new_filters_best.t_layer;
      memcpy(new_aps->new_filter_flag, alf_slice_param_new_filters_best.new_filter_flag, sizeof(new_aps->new_filter_flag));
      new_aps->fixed_filter_pattern = alf_slice_param_new_filters_best.fixed_filter_pattern;
      memcpy(new_aps->fixed_filter_idx, alf_slice_param_new_filters_best.fixed_filter_idx, sizeof(new_aps->fixed_filter_idx));
      new_aps->fixed_filter_set_index = alf_slice_param_new_filters_best.fixed_filter_set_index;

      //m_apsMap->setChangedFlag(newApsId);
      state->slice->param_set_map[new_aps_id].b_changed = true;

      //m_apsIdStart = newApsId;
      g_aps_id_start = new_aps_id;
    }

    //std::vector<int> apsIds = cs.slice->getTileGroupApsIdLuma();
    int* aps_ids = state->slice->tile_group_luma_aps_id;
    for (int i = 0; i < state->slice->tile_group_num_aps; i++)
    {
      //apss[apsIds[i]] = m_apsMap->getPS(apsIds[i]);
      state->slice->aps[aps_ids[i]] = state->slice->param_set_map[aps_ids[i]].parameter_set;
    }
  }
  
  //chroma
  //setCtuEnableFlag(m_ctuEnableFlag, CHANNEL_TYPE_CHROMA, 1);
  memset(g_ctu_enable_flag[COMPONENT_Cb], 1, sizeof(uint8_t) * num_ctus_in_pic);
  memset(g_ctu_enable_flag[COMPONENT_Cr], 1, sizeof(uint8_t) * num_ctus_in_pic);

  get_frame_stats(CHANNEL_TYPE_CHROMA, 0, 1);

  cost_off = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA][0], CHANNEL_TYPE_CHROMA) + g_lambda[CHANNEL_TYPE_CHROMA]/*981.62883931057581*/ * 1.0;
  cost_min = MAX_DOUBLE;
  //m_CABACEstimator->getCtx() = AlfCtx(ctxBest);
  //ctxStart = AlfCtx(m_CABACEstimator->getCtx());
  int new_aps_id_chroma = -1;
  if (alf_slice_param_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA] && (alf_slice_param_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_slice_param_new_filters_best.enabled_flag[COMPONENT_Cr]))
  {
    new_aps_id_chroma = new_aps_id;
  }
  else if (alf_slice_param_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_slice_param_new_filters_best.enabled_flag[COMPONENT_Cr])
  {
    int cur_id = g_aps_id_start;
    while (new_aps_id_chroma < 0)
    {
      cur_id--;
      if (cur_id < 0)
      {
        cur_id = MAX_NUM_APS - 1;
      }
      
      bool found = false;
      for (int i = 0; i < 7/*sizeof(best_aps_ids)/sizeof(best_aps_ids[0])*/; i++) {
        if (cur_id == best_aps_ids[i]) {
          found = true;
        }
      }
      if (!found)//find(bestApsIds.begin(), bestApsIds.end(), curId) == bestApsIds.end())
      {
        new_aps_id_chroma = cur_id;
      }
    }
  }
  
  for (int cur_aps_id = 0; cur_aps_id < MAX_NUM_APS; cur_aps_id++)
  {
    if (1/*(cs.slice->getPendingRasInit() || cs.slice->isIDRorBLA() || cs.slice->isIntra())*/ && cur_aps_id != new_aps_id_chroma)
    {
      continue;
    }
    //APS* curAPS = m_apsMap->getPS(curApsId);
    alf_aps* cur_aps = &state->slice->param_set_map[cur_aps_id].parameter_set;
    double cur_cost = (1/*cs.slice->isIntra()*/ && state->slice->tile_group_num_aps == 1) ? 0 : (g_lambda[CHANNEL_TYPE_CHROMA]/*981.62883931057581*/ * 5);
    
    if (cur_aps_id == new_aps_id_chroma)
    {
      //g_alf_slice_aps_temp = aps;
      memcpy(g_alf_slice_aps_temp.enabled_flag, aps->enabled_flag, sizeof(g_alf_slice_aps_temp.enabled_flag));
      memcpy(g_alf_slice_aps_temp.luma_coeff, aps->luma_coeff, sizeof(g_alf_slice_aps_temp.luma_coeff));
      memcpy(g_alf_slice_aps_temp.chroma_coeff, aps->chroma_coeff, sizeof(g_alf_slice_aps_temp.chroma_coeff));
      memcpy(g_alf_slice_aps_temp.filter_coeff_delta_idx, aps->filter_coeff_delta_idx, sizeof(g_alf_slice_aps_temp.filter_coeff_delta_idx));
      memcpy(g_alf_slice_aps_temp.alf_luma_coeff_flag, aps->alf_luma_coeff_flag, sizeof(g_alf_slice_aps_temp.alf_luma_coeff_flag));
      g_alf_slice_aps_temp.num_luma_filters = aps->num_luma_filters;
      g_alf_slice_aps_temp.alf_luma_coeff_delta_flag = aps->alf_luma_coeff_delta_flag;
      g_alf_slice_aps_temp.alf_luma_coeff_delta_prediction_flag = aps->alf_luma_coeff_delta_prediction_flag;
      g_alf_slice_aps_temp.t_layer = aps->t_layer;
      memcpy(g_alf_slice_aps_temp.new_filter_flag, aps->new_filter_flag, sizeof(g_alf_slice_aps_temp.new_filter_flag));
      g_alf_slice_aps_temp.fixed_filter_pattern = aps->fixed_filter_pattern;
      memcpy(g_alf_slice_aps_temp.fixed_filter_idx, aps->fixed_filter_idx, sizeof(g_alf_slice_aps_temp.fixed_filter_idx));
      g_alf_slice_aps_temp.fixed_filter_set_index = aps->fixed_filter_set_index;

      cur_cost += g_lambda[CHANNEL_TYPE_CHROMA] * g_bits_new_filter[CHANNEL_TYPE_CHROMA];
    }
    else if (cur_aps && cur_aps->t_layer <= state->slice->id && cur_aps->new_filter_flag[CHANNEL_TYPE_CHROMA])
    {
      //g_alf_slice_aps_temp = curAPS;
      memcpy(g_alf_slice_aps_temp.enabled_flag, cur_aps->enabled_flag, sizeof(g_alf_slice_aps_temp.enabled_flag));
      memcpy(g_alf_slice_aps_temp.luma_coeff, cur_aps->luma_coeff, sizeof(g_alf_slice_aps_temp.luma_coeff));
      memcpy(g_alf_slice_aps_temp.chroma_coeff, cur_aps->chroma_coeff, sizeof(g_alf_slice_aps_temp.chroma_coeff));
      memcpy(g_alf_slice_aps_temp.filter_coeff_delta_idx, cur_aps->filter_coeff_delta_idx, sizeof(g_alf_slice_aps_temp.filter_coeff_delta_idx));
      memcpy(g_alf_slice_aps_temp.alf_luma_coeff_flag, cur_aps->alf_luma_coeff_flag, sizeof(g_alf_slice_aps_temp.alf_luma_coeff_flag));
      g_alf_slice_aps_temp.num_luma_filters = cur_aps->num_luma_filters;
      g_alf_slice_aps_temp.alf_luma_coeff_delta_flag = cur_aps->alf_luma_coeff_delta_flag;
      g_alf_slice_aps_temp.alf_luma_coeff_delta_prediction_flag = cur_aps->alf_luma_coeff_delta_prediction_flag;
      g_alf_slice_aps_temp.t_layer = cur_aps->t_layer;
      memcpy(g_alf_slice_aps_temp.new_filter_flag, cur_aps->new_filter_flag, sizeof(g_alf_slice_aps_temp.new_filter_flag));
      g_alf_slice_aps_temp.fixed_filter_pattern = cur_aps->fixed_filter_pattern;
      memcpy(g_alf_slice_aps_temp.fixed_filter_idx, cur_aps->fixed_filter_idx, sizeof(g_alf_slice_aps_temp.fixed_filter_idx));
      g_alf_slice_aps_temp.fixed_filter_set_index = cur_aps->fixed_filter_set_index;
    }
    else
    {
      continue;
    }
    kvz_alf_reconstruct_coeff(state, &g_alf_slice_aps_temp, CHANNEL_TYPE_CHROMA, true, true);
    //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
    for (int comp_id = 1; comp_id < MAX_NUM_COMPONENT; comp_id++)
    {
      g_alf_slice_aps_temp.enabled_flag[comp_id] = true;
      for (int ctb_idx = 0; ctb_idx < num_ctus_in_pic; ctb_idx++)
      {
        double dist_unfilter_ctu = g_ctb_distortion_unfilter[comp_id][ctb_idx];
        //cost on
        g_ctu_enable_flag[comp_id][ctb_idx] = 1;
        //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
        //rate
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
        //m_CABACEstimator->resetBits();
        //ctb flag
        //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, compId, &m_alfSliceParamTemp);
        code_alf_ctu_enable_flag(state, lcu, ctb_idx, comp_id, &g_alf_slice_aps_temp);
        double rate_on = frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;
        //distortion
        for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
        {
          g_filter_tmp[i] = g_chroma_coeff_final[i];
        }
        double dist = dist_unfilter_ctu + calc_error_for_coeffs(g_alf_covariance[comp_id][0][ctb_idx][0].ee, g_alf_covariance[comp_id][0][ctb_idx][0].y, g_filter_tmp, MAX_NUM_ALF_CHROMA_COEFF, ALF_NUM_BITS);
        double cost_on = dist + g_lambda[comp_id] * rate_on;
        //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
        //cost off
        g_ctu_enable_flag[comp_id][ctb_idx] = 0;
        //rate
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
        //m_CABACEstimator->resetBits();
        //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, compId, &m_alfSliceParamTemp);
        code_alf_ctu_enable_flag(state, lcu, ctb_idx, comp_id, &g_alf_slice_aps_temp);
        //cost
        double cost_off = dist_unfilter_ctu + g_lambda[comp_id] * frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;
        if (cost_on < cost_off)
        {
          //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
          g_ctu_enable_flag[comp_id][ctb_idx] = 1;
          cur_cost += cost_on;
        }
        else
        {
          g_ctu_enable_flag[comp_id][ctb_idx] = 0;
          cur_cost += cost_off;
        }
      }
    }
    //chroma idc
    //setEnableFlag(m_alfSliceParamTemp, CHANNEL_TYPE_CHROMA, m_ctuEnableFlag);
    const alf_component_id comp_id_first = COMPONENT_Cb;
    const alf_component_id comp_id_last = COMPONENT_Cr;
    for (int comp_id = comp_id_first; comp_id <= comp_id_last; comp_id++)
    {
      g_alf_slice_aps_temp.enabled_flag[comp_id] = false;
      for (int i = 0; i < num_ctus_in_pic; i++)
      {
        if (g_ctu_enable_flag[comp_id][i])
        {
          g_alf_slice_aps_temp.enabled_flag[comp_id] = true;
          break;
        }
      }
    }

    const int alf_chroma_idc = g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cb] * 2 + g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cr];
    cur_cost += length_truncated_unary(alf_chroma_idc, 3) * g_lambda[CHANNEL_TYPE_CHROMA];

    if (cur_cost < cost_min)
    {
      cost_min = cur_cost;
      //cs.slice->setTileGroupApsIdChroma(curApsId);
      state->slice->tile_group_chroma_aps_id = cur_aps_id;
      //cs.slice->setTileGroupAlfEnabledFlag(COMPONENT_Cb, m_alfSliceParamTemp.enabledFlag[COMPONENT_Cb]);
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cb];
      //cs.slice->setTileGroupAlfEnabledFlag(COMPONENT_Cr, m_alfSliceParamTemp.enabledFlag[COMPONENT_Cr]);
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cr];
      //copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, CHANNEL_TYPE_CHROMA);
      
      memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cb], g_ctu_enable_flag[COMPONENT_Cb], sizeof(uint8_t) * num_ctus_in_pic);
      memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cr], g_ctu_enable_flag[COMPONENT_Cr], sizeof(uint8_t) * num_ctus_in_pic);
    }
  }
  if (cost_off < cost_min)
  {
    //cs.slice->setTileGroupAlfEnabledFlag(COMPONENT_Cb, false);
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = false;
    //cs.slice->setTileGroupAlfEnabledFlag(COMPONENT_Cr, false);
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = false;
    //setCtuEnableFlag(m_ctuEnableFlag, CHANNEL_TYPE_CHROMA, 0);
    memset(g_ctu_enable_flag[COMPONENT_Cb], 0, sizeof(uint8_t) * num_ctus_in_pic);
    memset(g_ctu_enable_flag[COMPONENT_Cr], 0, sizeof(uint8_t) * num_ctus_in_pic);
  }
  else
  {
    //copyCtuEnableFlag(m_ctuEnableFlag, m_ctuEnableFlagTmp, CHANNEL_TYPE_CHROMA);
    memcpy(g_ctu_enable_flag[COMPONENT_Cb], g_ctu_enable_flag_tmp[COMPONENT_Cb], sizeof(uint8_t) * num_ctus_in_pic);
    memcpy(g_ctu_enable_flag[COMPONENT_Cr], g_ctu_enable_flag_tmp[COMPONENT_Cr], sizeof(uint8_t) * num_ctus_in_pic);
    if (state->slice->tile_group_chroma_aps_id == new_aps_id_chroma)  //new filter
    {
      //APS* newAPS = m_apsMap->getPS(new_aps_id_chroma);
      alf_aps* new_aps = &state->slice->param_set_map[new_aps_id_chroma].parameter_set;
      if (new_aps == NULL)
      {
        //newAPS = m_apsMap->allocatePS(new_aps_id_chroma);
        assert(new_aps_id_chroma < MAX_NUM_APS); //Invalid PS id
        bool found = false;
        for (int i = 0; i < (sizeof(state->slice->param_set_map)/sizeof(state->slice->param_set_map[0])); i++) {
          if (state->slice->param_set_map[i].parameter_set.aps_id == new_aps_id_chroma) {
            found = true;
          }
        }
        if (!found)
        {
          state->slice->param_set_map[new_aps_id_chroma].b_changed = true;
          //state->slice->param_set_map[new_aps_id_chroma].p_nalu_data = 0;
          //state->slice->param_set_map[new_aps_id_chroma].parameter_set = malloc(sizeof(alf_aps));
          //setID(m_paramsetMap[psId].parameterSet, psId);
          state->slice->param_set_map[new_aps_id_chroma].parameter_set.aps_id = new_aps_id_chroma;
        }

        new_aps->aps_id = state->slice->param_set_map[new_aps_id_chroma].parameter_set.aps_id;
        
        //newAPS->getAlfAPSParam().reset();
        memset(new_aps->enabled_flag, false, sizeof(new_aps->enabled_flag));
        memset(new_aps->luma_coeff, 0, sizeof(new_aps->luma_coeff));
        memset(new_aps->chroma_coeff, 0, sizeof(new_aps->chroma_coeff));
        memset(new_aps->filter_coeff_delta_idx, 0, sizeof(new_aps->filter_coeff_delta_idx));
        memset(new_aps->alf_luma_coeff_flag, true, sizeof(new_aps->alf_luma_coeff_flag));
        new_aps->num_luma_filters = 1;
        new_aps->alf_luma_coeff_delta_flag = false;
        new_aps->alf_luma_coeff_delta_prediction_flag = false;
        new_aps->t_layer = 0;
        memset(new_aps->new_filter_flag, 0, sizeof(new_aps->new_filter_flag));
        new_aps->fixed_filter_pattern = 0;
        memset(new_aps->fixed_filter_idx, 0, sizeof(new_aps->fixed_filter_idx));
        new_aps->fixed_filter_set_index = 0;;
      }
      //newAPS->getAlfAPSParam().newFilterFlag[CHANNEL_TYPE_CHROMA] = true;
      new_aps->new_filter_flag[CHANNEL_TYPE_CHROMA] = true;
      //newAPS->getAlfAPSParam().tLayer = cs.slice->getTLayer();
      new_aps->t_layer = state->slice->id;
      for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
      {
        //newAPS->getAlfAPSParam().chromaCoeff[i] = alfSliceParamNewFilters.chromaCoeff[i];
        new_aps->chroma_coeff[i] = aps->chroma_coeff[i];
      }
      //m_apsMap->setChangedFlag(new_aps_id_chroma);
      state->slice->param_set_map[new_aps_id_chroma].b_changed = true;
      g_aps_id_start = new_aps_id_chroma;
    }
    //apss[cs.slice->getTileGroupApsIdChroma()] = m_apsMap->getPS(cs.slice->getTileGroupApsIdChroma());
    //apss[state->slice->tile_group_chroma_aps_id] = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set;
    state->slice->aps[state->slice->tile_group_chroma_aps_id].aps_id = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.aps_id;
    memcpy(state->slice->aps[state->slice->tile_group_chroma_aps_id].enabled_flag, state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.enabled_flag, sizeof(state->slice->aps[state->slice->tile_group_chroma_aps_id].enabled_flag));
    memcpy(state->slice->aps[state->slice->tile_group_chroma_aps_id].luma_coeff, state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.luma_coeff, sizeof(state->slice->aps[state->slice->tile_group_chroma_aps_id].luma_coeff));
    memcpy(state->slice->aps[state->slice->tile_group_chroma_aps_id].chroma_coeff, state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.chroma_coeff, sizeof(state->slice->aps[state->slice->tile_group_chroma_aps_id].chroma_coeff));
    memcpy(state->slice->aps[state->slice->tile_group_chroma_aps_id].filter_coeff_delta_idx, state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.filter_coeff_delta_idx, sizeof(state->slice->aps[state->slice->tile_group_chroma_aps_id].filter_coeff_delta_idx));
    memcpy(state->slice->aps[state->slice->tile_group_chroma_aps_id].alf_luma_coeff_flag, state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.alf_luma_coeff_flag, sizeof(state->slice->aps[state->slice->tile_group_chroma_aps_id].alf_luma_coeff_flag));
    state->slice->aps[state->slice->tile_group_chroma_aps_id].num_luma_filters = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.num_luma_filters;
    state->slice->aps[state->slice->tile_group_chroma_aps_id].alf_luma_coeff_delta_flag = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.alf_luma_coeff_delta_flag;
    state->slice->aps[state->slice->tile_group_chroma_aps_id].alf_luma_coeff_delta_prediction_flag = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.alf_luma_coeff_delta_prediction_flag;
    state->slice->aps[state->slice->tile_group_chroma_aps_id].t_layer = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.t_layer;
    memcpy(state->slice->aps[state->slice->tile_group_chroma_aps_id].new_filter_flag, state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.new_filter_flag, sizeof(state->slice->aps[state->slice->tile_group_chroma_aps_id].new_filter_flag));
    state->slice->aps[state->slice->tile_group_chroma_aps_id].fixed_filter_pattern = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.fixed_filter_pattern;
    memcpy(state->slice->aps[state->slice->tile_group_chroma_aps_id].fixed_filter_idx, state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.fixed_filter_idx, sizeof(state->slice->aps[state->slice->tile_group_chroma_aps_id].fixed_filter_idx));
    state->slice->aps[state->slice->tile_group_chroma_aps_id].fixed_filter_set_index = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id].parameter_set.fixed_filter_set_index;
  }
}

void kvz_alf_get_avai_aps_ids_luma(encoder_state_t *const state, 
  int *newApsId, 
  int aps_ids[6],
  int *size_of_aps_ids)
{
  //APS** apss = cs.slice->getAPSs();
  //alf_aps** apss = &state->slice->aps;
  param_set_map *aps_set = state->slice->param_set_map;

  for (int i = 0; i < MAX_NUM_APS; i++)
  {
    //apss[i] = m_apsMap->getPS(i);
    state->slice->aps[i] = aps_set[i].parameter_set;
  }

  //std::vector<int> result;
  int apsIdChecked = 0, curApsId = g_aps_id_start;
  if (curApsId < MAX_NUM_APS)
  {
    while (apsIdChecked < MAX_NUM_APS && /*!cs.slice->isIntra() &&*/ *size_of_aps_ids < (6/*ALF_CTB_MAX_NUM_APS*/ - 1) /*&& /*!cs.slice->getPendingRasInit() &&*/ /*!cs.slice->isIDRorBLA()*/)
    {
      //APS* curAPS = cs.slice->getAPSs()[curApsId];
      alf_aps *curAPS = &state->slice->aps[curApsId];

      if (curAPS && curAPS->t_layer/*curAPS->getTemporalId()*/ <= state->slice->id/*cs.slice->getTLayer()*/ && curAPS->new_filter_flag[CHANNEL_TYPE_LUMA])
      {
        //result.push_back(curApsId);
        aps_ids[*size_of_aps_ids] = curApsId;
        *size_of_aps_ids++;
      }
      apsIdChecked++;
      curApsId = (curApsId + 1) % MAX_NUM_APS;
    }
  }
  //cs.slice->setTileGroupNumAps((int)result.size());
  //memset(state->slice[0].tile_group_num_aps, *size_of_aps_ids, sizeof(int));
  state->slice->tile_group_num_aps = *size_of_aps_ids;
  //cs.slice->setAPSs(result);
  for (int i = 0; i < state->slice->tile_group_num_aps; i++)
  {
    state->slice->tile_group_luma_aps_id[i] = aps_ids[i];
  }

  *newApsId = g_aps_id_start - 1;
  if (*newApsId < 0)
  {
    *newApsId = (int)MAX_NUM_APS - 1;
  }

  assert(*newApsId <= (int)MAX_NUM_APS); //Wrong APS index assignment in getAvaiApsIdsLuma
  //return result;
}

void kvz_alf_reconstructor(encoder_state_t const *state,
  const lcu_order_element_t *const lcu)
{
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;

  if (!state->slice->tile_group_alf_enabled_flag[COMPONENT_Y])
  {
    //return;
  }

  kvz_alf_reconstruct_coeff_aps(state, true, state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] || state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr], false);
  short* alf_ctu_filter_index = g_alf_ctb_filter_index;//cs.slice->getPic()->getAlfCtbFilterIndex();
  //PelUnitBuf& recBuf = cs.getRecoBufRef();
  int luma_height = state->tile->frame->height;
  int luma_width = state->tile->frame->width;
  //const PreCalcValues& pcv = *cs.pcv;

  int ctu_idx = 0;
  for (int y_pos = 0; y_pos < luma_height /*pcv.lumaHeight*/; y_pos += LCU_WIDTH)
  {
    for (int x_pos = 0; x_pos < luma_width; x_pos += LCU_WIDTH)
    {
      const int width = (x_pos + LCU_WIDTH > luma_width) ? (luma_width - x_pos) : LCU_WIDTH;
      const int height = (y_pos + LCU_WIDTH > luma_height) ? (luma_height - y_pos) : LCU_WIDTH;

      //const UnitArea area(cs.area.chromaFormat, Area(xPos, yPos, width, height));
      if (g_ctu_enable_flag[COMPONENT_Y][ctu_idx])
      {
        //Area blk(xPos, yPos, width, height);
        short filter_set_index = alf_ctu_filter_index[ctu_idx];
        short *coeff;

        if (filter_set_index >= ALF_NUM_FIXED_FILTER_SETS)
        {
          coeff = g_coeff_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
        }
        else
        {
          coeff = g_fixed_filter_set_coeff_dec[filter_set_index];
        }
        //m_filter7x7Blk(m_classifier, recBuf, recExtBuf, blk, COMPONENT_Y, coeff, m_clpRngs.comp[COMPONENT_Y], cs);
        kvz_alf_filter_block(state, lcu, coeff, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y, width, height, x_pos, y_pos);
      }

      for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
      {
        //ComponentID compID = ComponentID(compIdx);
        alf_component_id comp_id = comp_idx;
        const int chroma_scale_x = ((comp_id == COMPONENT_Y) || (chroma_fmt == KVZ_CSP_444)) ? 0 : 1;
        const int chroma_scale_y = ((comp_id == COMPONENT_Y) || (chroma_fmt != KVZ_CSP_420)) ? 0 : 1;
        if (g_ctu_enable_flag[comp_idx][ctu_idx])
        {
          //Area blk(xPos >> chromaScaleX, yPos >> chromaScaleY, width >> chromaScaleX, height >> chromaScaleY);
          //m_filter5x5Blk(m_classifier, recBuf, recExtBuf, blk, compID, m_chromaCoeffFinal, m_clpRngs.comp[compIdx], cs);
          kvz_alf_filter_block(state, lcu, g_chroma_coeff_final, g_clp_rngs.comp[comp_idx], comp_idx, 
                               width >> chroma_scale_x, height >> chroma_scale_y, 
                               x_pos >> chroma_scale_x, y_pos >> chroma_scale_y);
        }
      }
      ctu_idx++;
    }
  }
}

void kvz_alf_derive_stats_for_filtering(encoder_state_t *const state,
  const lcu_order_element_t *const lcu)
{
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;

  int32_t pic_width = state->tile->frame->rec->width;
  int32_t pic_height = state->tile->frame->rec->height;

  int ctu_rs_addr = 0;
  const int number_of_components = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_COMPONENT;

  // init CTU stats buffers
  for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
  {
    //const ComponentID compID = ComponentID(compIdx);
    bool is_luma = comp_idx == 0 ? 1 : 0;
    const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;

    for (int shape = 0; shape != 1 /*m_filterShapes[toChannelType(compID)].size()*/; shape++)
    {
      for (int class_idx = 0; class_idx < num_classes; class_idx++)
      {
        for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
        {
          alf_covariance *alf_cov = &g_alf_covariance[comp_idx][shape][ctu_idx][class_idx];
          alf_cov->pix_acc = 0;
          memset(alf_cov->y, 0, sizeof(*alf_cov->y) * alf_cov->num_coeff);
          for (int i = 0; i < alf_cov->num_coeff; i++)
          {
            memset(alf_cov->ee[i], 0, sizeof(*alf_cov->ee[i]) * alf_cov->num_coeff);
          }
        }
      }
    }
  }

  // init Frame stats buffers
  const int number_of_channels = 2 /*getNumberValidChannels(m_chromaFormat)*/;
  for (int channel_idx = 0; channel_idx < number_of_channels; channel_idx++)
  {
    //const ChannelType channelID = ChannelType(channelIdx);
    bool is_luma = channel_idx == 0 ? true : false;
    const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;

    for (int shape = 0; shape != 1/*m_filterShapes[channelIdx].size()*/; shape++)
    {
      for (int class_idx = 0; class_idx < num_classes; class_idx++)
      {
        alf_covariance *alf_cov_f = &g_alf_covariance_frame[channel_idx][shape][class_idx];
        alf_cov_f->pix_acc = 0;
        memset(alf_cov_f->y, 0, sizeof(alf_cov_f->y) * alf_cov_f->num_coeff);
        for (int i = 0; i < alf_cov_f->num_coeff; i++)
        {
          memset(alf_cov_f->ee[i], 0, sizeof(alf_cov_f->ee[i]) * alf_cov_f->num_coeff);
        }
      }
    }
  }

  int max_cu_width = LCU_WIDTH;
  int max_cu_height = LCU_WIDTH;
  for (int y_pos = 0; y_pos < pic_height; y_pos += max_cu_height)
  {
    for (int x_pos = 0; x_pos < pic_width; x_pos += max_cu_width)
    {
      int width = (x_pos + max_cu_width > pic_width) ? (pic_width - x_pos) : max_cu_width;
      int height = (y_pos + max_cu_height > pic_height) ? (pic_height - y_pos) : max_cu_height;
      //const UnitArea area(m_chromaFormat, Area(xPos, yPos, width, height));

      for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
      {
        const alf_component_id comp_id = comp_idx;

        width = comp_id            == COMPONENT_Y ? width  : width / 2;
        height = comp_id           == COMPONENT_Y ? height : height / 2;
        int pos_x = comp_id        == COMPONENT_Y ? x_pos  : x_pos / 2;
        int pos_y = comp_id        == COMPONENT_Y ? y_pos  : y_pos / 2;
 
        int32_t org_stride = comp_id == COMPONENT_Y ? state->tile->frame->source->stride : state->tile->frame->source->stride / 2;
        int32_t rec_stride = comp_id == COMPONENT_Y ? state->tile->frame->rec->stride    : state->tile->frame->rec->stride / 2;

        kvz_pixel *org = comp_id ? (comp_id - 1 ? &state->tile->frame->source->v[(y_pos + 2) * org_stride + x_pos + 2] : &state->tile->frame->source->u[(y_pos + 2) * org_stride + x_pos + 2]) 
                                 : &state->tile->frame->source->y[(y_pos + 3) * org_stride + x_pos + 3];
        kvz_pixel *rec = comp_id ? (comp_id - 1 ? &state->tile->frame->rec->v[(y_pos + 2) * rec_stride + x_pos + 2]    : &state->tile->frame->rec->u[(y_pos + 2) * rec_stride + x_pos + 2])
                                 : &state->tile->frame->rec->y[(y_pos + 3) * rec_stride + x_pos + 3];

        channel_type channel = comp_id == COMPONENT_Y ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;

        for (int shape = 0; shape != 1/*m_filterShapes[chType].size()*/; shape++)
        {
          kvz_alf_get_blk_stats(state, lcu, channel, g_alf_covariance[comp_idx][shape][ctu_rs_addr], comp_idx ? NULL : classifier, org, org_stride, rec, rec_stride, pos_x, pos_y, width, height);

          const int num_classes = channel ? 1 : MAX_NUM_ALF_CLASSES;

          for (int class_idx = 0; class_idx < num_classes; class_idx++)
          {
            //g_alf_covariance_frame[channel][shape][classIdx] += g_alf_covariance[compIdx][shape][ctuRsAddr][classIdx];
            int num_coeff = g_alf_covariance_frame[channel][shape][class_idx].num_coeff;
            for (int j = 0; j < num_coeff; j++)
            {
              for (int i = 0; i < num_coeff; i++)
              {
                g_alf_covariance_frame[channel][shape][class_idx].ee[j][i] += g_alf_covariance[comp_idx][shape][ctu_rs_addr][class_idx].ee[j][i];
              }
              g_alf_covariance_frame[channel][shape][class_idx].y[j] += g_alf_covariance[comp_idx][shape][ctu_rs_addr][class_idx].y[j];
            }
            g_alf_covariance_frame[channel][shape][class_idx].pix_acc += g_alf_covariance[comp_idx][shape][ctu_rs_addr][class_idx].pix_acc;
          }
        }
      }
      ctu_rs_addr++;
    }
  }
}

void kvz_alf_get_blk_stats(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  alf_covariance *alf_covariace,
  alf_classifier **classifier,
  kvz_pixel *org,
  int32_t org_stride,
  kvz_pixel *rec,
  int32_t rec_stride,
  const int x_pos,
  const int y_pos,
  const int width,
  const int height)
{
  //alf_classifier **classifier = state->tile->frame->alf_info->classifier;
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? true : false;
  const int num_coeff = is_luma ? 13 : 7;
  const int *pattern = is_luma ? alf_pattern_7 : alf_pattern_5;
  const int filter_length = is_luma ? 7 : 5;
  const int loop_adjust = is_luma ? 3 : 2;

  static int e_local[MAX_NUM_ALF_LUMA_COEFF];

  int transpose_idx = 0;
  int class_idx = 0;

  for (int i = 0; i < height; i++)
  {
    for (int j = 0; j < width; j++)
    {
      if (classifier && classifier[y_pos + i][x_pos + j].class_idx == ALF_UNUSED_CLASS_IDX && classifier[y_pos + i][x_pos + j].transpose_idx == ALF_UNUSED_TRANSPOSE_IDX)
      {
        continue;
      }
      //memset(e_local, 0, 13 * sizeof(int));
      e_local[0] = e_local[1] = e_local[2] = e_local[3]  = e_local[4]  = e_local[5]  = e_local[6] = 
      e_local[7] = e_local[8] = e_local[9] = e_local[10] = e_local[11] = e_local[12] = 0;
      if (classifier)
      {
        alf_classifier *cl = &classifier[y_pos + i][x_pos + j];
        transpose_idx = cl->transpose_idx;
        class_idx = cl->class_idx;
      }

      double weight = 1.0;
      if (0 /*m_alfWSSD*/)
      {
        //weight = m_lumaLevelToWeightPLUT[org[j]];
      }

      int y_local = org[j] - rec[j];
      kvz_alf_calc_covariance(e_local, rec + j, rec_stride, pattern, filter_length >> 1, transpose_idx);
      for (int k = 0; k < num_coeff; k++)
      {
        for (int l = k; l < num_coeff; l++)
        {
          if (0/*m_alfWSSD*/)
          {
            alf_covariace[class_idx].ee[k][l] += weight * (double)(e_local[k] * e_local[l]);
          }
          else
            alf_covariace[class_idx].ee[k][l] += (double)(e_local[k] * e_local[l]);
        }
        if (0/*m_alfWSSD*/)
        {
          alf_covariace[class_idx].y[k] += weight * (double)(e_local[k] * y_local);
        }
        else
          alf_covariace[class_idx].y[k] += (double)(e_local[k] * y_local);
      }
      if (0/*m_alfWSSD*/)
      {
        alf_covariace[class_idx].pix_acc += weight * (double)(y_local * y_local);
      }
      else
        alf_covariace[class_idx].pix_acc += (double)(y_local * y_local);
    }
    org += org_stride;
    rec += rec_stride;
  }

  int num_classes = classifier ? MAX_NUM_ALF_CLASSES : 1;
  for (class_idx = 0; class_idx < num_classes; class_idx++)
  {
    for (int k = 1; k < num_coeff; k++)
    {
      for (int l = 0; l < k; l++)
      {
        alf_covariace[class_idx].ee[k][l] = alf_covariace[class_idx].ee[l][k];
      }
    }
  }
}

void kvz_alf_calc_covariance(int *e_local,
  const kvz_pixel *rec,
  const int stride,
  const int *filter_pattern,
  const int half_filter_length,
  const int transpose_idx)
{
  int k = 0;

  if (transpose_idx == 0)
  {
    for (int i = -half_filter_length; i < 0; i++)
    {
      const kvz_pixel* rec0 = rec + i * stride;
      const kvz_pixel* rec1 = rec - i * stride;

      for (int j = -half_filter_length - i; j <= half_filter_length + i; j++)
      {
        e_local[filter_pattern[k++]] += rec0[j] + rec1[-j];
      }
    }
    for (int j = -half_filter_length; j < 0; j++)
    {
      e_local[filter_pattern[k++]] += rec[j] + rec[-j];
    }
  }
  else if (transpose_idx == 1)
  {
    for (int j = -half_filter_length; j < 0; j++)
    {
      const kvz_pixel* rec0 = rec + j;
      const kvz_pixel* rec1 = rec - j;

      for (int i = -half_filter_length - j; i <= half_filter_length + j; i++)
      {
        e_local[filter_pattern[k++]] += rec0[i * stride] + rec1[-i * stride];
      }
    }

    for (int i = -half_filter_length; i < 0; i++)
    {
      e_local[filter_pattern[k++]] += rec[i*stride] + rec[-i * stride];
    }
  }
  else if (transpose_idx == 2)
  {
    for (int i = -half_filter_length; i < 0; i++)
    {
      const kvz_pixel* rec0 = rec + i * stride;
      const kvz_pixel* rec1 = rec - i * stride;

      for (int j = half_filter_length + i; j >= -half_filter_length - i; j--)
      {
        e_local[filter_pattern[k++]] += rec0[j] + rec1[-j];
      }
    }
    for (int j = -half_filter_length; j < 0; j++)
    {
      e_local[filter_pattern[k++]] += rec[j] + rec[-j];
    }
  }
  else
  {
    for (int j = -half_filter_length; j < 0; j++)
    {
      const kvz_pixel* rec0 = rec + j;
      const kvz_pixel* rec1 = rec - j;

      for (int i = half_filter_length + j; i >= -half_filter_length - j; i--)
      {
        e_local[filter_pattern[k++]] += rec0[i * stride] + rec1[-i * stride];
      }
    }
    for (int i = -half_filter_length; i < 0; i++)
    {
      e_local[filter_pattern[k++]] += rec[i*stride] + rec[-i * stride];
    }
  }
  e_local[filter_pattern[k++]] += rec[0];
}

double kvz_alf_get_filter_coeff_and_cost(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  double dist_unfilter,
  int *ui_coeff_bits,
  int i_shape_idx,
  bool b_re_collect_stat,
  bool only_filter_cost)
{
  alf_aps *aps = state->slice->aps;
  short *chroma_coeff = aps->chroma_coeff;

  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int32_t slice_flags[3] = { 0,0,0 };
  slice_flags[0] = state->encoder_control->cfg.alf_slice_enable_flag[0];
  slice_flags[1] = state->encoder_control->cfg.alf_slice_enable_flag[1];
  slice_flags[2] = state->encoder_control->cfg.alf_slice_enable_flag[2];

  //collect stat based on CTU decision
  if (b_re_collect_stat)
  {
    //getFrameStats(channel, iShapeIdx);
  }

  double dist = dist_unfilter;
  *ui_coeff_bits = 0;
  int ui_slice_flag = 0;
  
  //AlfFilterShape& alfFilterShape = m_alfSliceParamTemp.filterShapes[channel][iShapeIdx];
  //get filter coeff
  if (is_luma)
  {
   dist += kvz_alf_merge_filters_and_cost(state, lcu, channel, ui_coeff_bits, g_alf_covariance_frame[channel][i_shape_idx], g_alf_covariance_merged[i_shape_idx]);
  }
  else
  {
    //distortion
    dist += g_alf_covariance_frame[channel][i_shape_idx][0].pix_acc + kvz_alf_derive_coeff_quant(channel, g_filter_coeff_quant, g_alf_covariance_frame[channel][i_shape_idx][0].ee, g_alf_covariance_frame[channel][i_shape_idx][0].y);
    memcpy(g_filter_coeff_set[0], g_filter_coeff_quant, sizeof(*g_filter_coeff_quant) * MAX_NUM_ALF_LUMA_COEFF);

    //setEnableFlag( m_alfSliceParamTemp, channel, m_ctuEnableFlag );
    // Slice Flag

    //const int alfChromaIdc = m_alfSliceParamTemp.enabledFlag[COMPONENT_Cb] * 2 + m_alfSliceParamTemp.enabledFlag[COMPONENT_Cr];
    const int alf_chroma_idc = slice_flags[COMPONENT_Cb] * 2 + slice_flags[COMPONENT_Cr];

    for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
    {
      chroma_coeff[i] = g_filter_coeff_quant[i];
    }
    ui_coeff_bits += get_coeff_rate(aps, true);

    ui_slice_flag = length_truncated_unary(alf_chroma_idc, 3);
  }

  if (only_filter_cost)
  {
    return dist + 981.62883931057581/*m_lambda[channel]*/ * *ui_coeff_bits;
  }

  double rate = *ui_coeff_bits + ui_slice_flag;

  //N‰m‰ kaksi funktiota: miten pit‰isi tehd‰ vastaavat kvazaarilla?
  //m_CABACEstimator->resetBits();
  //m_CABACEstimator->codeAlfCtuEnableFlags(cs, channel, &m_alfSliceParamTemp);
  //Vastaa alimmaista: (?)
  if (is_luma)
  {
    if (g_alf_slice_aps_temp.enabled_flag[COMPONENT_Y]) {
      for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ++ctu_idx) {
        code_alf_ctu_enable_flag(state, lcu, ctu_idx, COMPONENT_Y, &g_alf_slice_aps_temp);
      }
    }
  }
  else
  {
    if (state->encoder_control->cfg.alf_enable) {
      for (int ctu_idx = 0; ctu_idx < state->lcu_order_count; ++ctu_idx) {
        code_alf_ctu_enable_flag(state, lcu, ctu_idx, COMPONENT_Cr, &g_alf_slice_aps_temp);
        code_alf_ctu_enable_flag(state, lcu, ctu_idx, COMPONENT_Cb, &g_alf_slice_aps_temp);
      }
    }
  }

  rate += frac_bits_scale * 0;/*(double)m_CABACEstimator->getEstFracBits();*/ // <-- Not supported
  return dist + 981.62883931057581/*m_lambda[channel]*/ * rate;
}

int kvz_alf_derive_filter_coefficients_prediction_mode(channel_type channel,
  int **filter_set,
  int** filter_coeff_diff,
  const int num_filters,
  int *pred_mode)
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;

  int rate_pred_mode0 = get_cost_filter_coeff(channel, filter_set, num_filters);

  for (int ind = 0; ind < num_filters; ++ind)
  {
    if (ind == 0)
    {
      memcpy(filter_coeff_diff[ind], filter_set[ind], sizeof(int) * num_coeff);
    }
    else
    {
      for (int i = 0; i < num_coeff; i++)
      {
        filter_coeff_diff[ind][i] = filter_set[ind][i] - filter_set[ind - 1][i];
      }
    }
  }

  int rate_pred_mode1 = get_cost_filter_coeff(channel, filter_coeff_diff, num_filters);

  *pred_mode = (rate_pred_mode1 < rate_pred_mode0 && num_filters > 1) ? 1 : 0;

  return (num_filters > 1 ? 1 : 0)        // coeff_delta_pred_mode_flag
    + (pred_mode ? rate_pred_mode1 : rate_pred_mode0); // min_golomb_order, golomb_order_increase_flag, alf_coeff_luma_delta
}

void kvz_alf_merge_classes(channel_type channel,
  alf_covariance* cov,
  alf_covariance* cov_merged, 
  const int num_classes, 
  short filter_indices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES])
{
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;

  static bool available_class[MAX_NUM_ALF_CLASSES];
  static uint8_t index_list[MAX_NUM_ALF_CLASSES];
  static uint8_t index_list_temp[MAX_NUM_ALF_CLASSES];
  int num_remaining = num_classes;

  memset(filter_indices, 0, sizeof(short) * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_CLASSES);

  for (int i = 0; i < num_classes; i++)
  {
    filter_indices[num_remaining - 1][i] = i;
    index_list[i] = i;
    available_class[i] = true;
    //cov_merged[i] = cov[i];
    int num_coeff = cov_merged[i].num_coeff;
    for (int j = 0; j < num_coeff; j++)
    {
      memcpy(cov_merged[i].ee[j], cov[i].ee[j], sizeof(*cov_merged[i].ee[j]) * num_coeff);
    }
    memcpy(cov_merged[i].y, cov[i].y, sizeof(*cov_merged[i].y) * num_coeff);
    cov_merged[i].pix_acc = cov[i].pix_acc;

  }

  // Try merging different covariance matrices

  // temporal AlfCovariance structure is allocated as the last element in covMerged array, the size of covMerged is MAX_NUM_ALF_CLASSES + 1
  alf_covariance tmp_cov = cov_merged[MAX_NUM_ALF_CLASSES];

  while (num_remaining > 2)
  {
    double error_min = MAX_DOUBLE; //std::numeric_limits<double>::max();
    int best_to_merge_idx1 = 0, best_to_merge_idx2 = 1;

    for (int i = 0; i < num_classes - 1; i++)
    {
      if (available_class[i])
      {
        for (int j = i + 1; j < num_classes; j++)
        {
          if (available_class[j])
          {
            double error1 = calculate_error(cov_merged[i]);
            double error2 = calculate_error(cov_merged[j]);

            //tmp_cov.add(cov_merged[i], cov_merged[j]);
            for (int j = 0; j < num_coeff; j++)
            {
              for (int i = 0; i < num_coeff; i++)
              {
                tmp_cov.ee[j][i] = cov_merged[i].ee[j][i] + cov_merged[j].ee[j][i];
              }
              tmp_cov.y[j] = cov_merged[i].y[j] + cov_merged[j].y[j];
            }

            tmp_cov.pix_acc = cov_merged[i].pix_acc + cov_merged[j].pix_acc;

            double error = calculate_error(tmp_cov) - error1 - error2;

            if (error < error_min)
            {
              error_min = error;
              best_to_merge_idx1 = i;
              best_to_merge_idx2 = j;
            }
          }
        }
      }
    }

    //cov_merged[best_to_merge_idx1] += cov_merged[best_to_merge_idx2];
    for (int j = 0; j < num_coeff; j++)
    {
      for (int i = 0; i < num_coeff; i++)
      {
        cov_merged[best_to_merge_idx1].ee[j][i] += cov_merged[best_to_merge_idx2].ee[j][i];
      }
      cov_merged[best_to_merge_idx1].y[j] += cov_merged[best_to_merge_idx2].y[j];
    }
    cov_merged[best_to_merge_idx1].pix_acc += cov_merged[best_to_merge_idx2].pix_acc;


    available_class[best_to_merge_idx2] = false;

    for (int i = 0; i < num_classes; i++)
    {
      if (index_list[i] == best_to_merge_idx2)
      {
        index_list[i] = best_to_merge_idx1;
      }
    }

    num_remaining--;
    if (num_remaining <= num_classes)
    {
      memcpy(index_list_temp, index_list, sizeof(uint8_t) * num_classes);

      bool exist = false;
      int ind = 0;

      for (int j = 0; j < num_classes; j++)
      {
        exist = false;
        for (int i = 0; i < num_classes; i++)
        {
          if (index_list_temp[i] == j)
          {
            exist = true;
            break;
          }
        }

        if (exist)
        {
          for (int i = 0; i < num_classes; i++)
          {
            if (index_list_temp[i] == j)
            {
              filter_indices[num_remaining - 1][i] = ind;
              index_list_temp[i] = -1;
            }
          }
          ind++;
        }
      }
    }
  }
}

double kvz_alf_merge_filters_and_cost(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  int *ui_coeff_bits,
  alf_covariance *cov_frame,
  alf_covariance *cov_merged)
{
  alf_aps *aps = state->slice->aps;
  short *luma_coeff = aps->luma_coeff;
  short *filter_coeff_delta_idx = aps->filter_coeff_delta_idx;
  bool *alf_luma_coeff_flag = aps->alf_luma_coeff_flag;
  int *num_luma_filters = &aps->num_luma_filters;
  bool *alf_luma_coeff_delta_flag = &aps->alf_luma_coeff_delta_flag;
  bool *alf_luma_coeff_delta_prediction_flag = &aps->alf_luma_coeff_delta_prediction_flag;
  int *fixed_filter_pattern = &aps->fixed_filter_pattern;
  int *fixed_filter_set_index = &aps->fixed_filter_set_index;

  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int num_filters_best = 0;
  int num_filters = MAX_NUM_ALF_CLASSES;
  static bool coded_var_bins[MAX_NUM_ALF_CLASSES];
  static double error_force_0_coeff_tab[MAX_NUM_ALF_CLASSES][2];

  double cost, cost0, dist, dist_force0, cost_min = MAX_DOUBLE;
  int pred_mode = 0, best_pred_mode = 0, coeff_bits, coeff_bits_force0;

  kvz_alf_merge_classes(channel, cov_frame, cov_merged, MAX_NUM_ALF_CLASSES, g_filter_indices);

  while (num_filters >= 1)
  {
    dist = kvz_alf_derive_filter_coeffs(aps, channel, cov_frame, cov_merged, g_filter_indices[num_filters-1], num_filters,error_force_0_coeff_tab);

    // filter coeffs are stored in m_filterCoeffSet
    dist_force0 = get_dist_force_0(channel, num_filters, error_force_0_coeff_tab, coded_var_bins);
    coeff_bits = kvz_alf_derive_filter_coefficients_prediction_mode(channel, g_filter_coeff_set, g_diff_filter_coeff, num_filters, &pred_mode);
    coeff_bits_force0 = get_cost_filter_coeff_force_0(channel, g_filter_coeff_set, num_filters, coded_var_bins);

    cost = dist + 981.62883931057581/*m_lambda[COMPONENT_Y]*/ * coeff_bits;
    cost0 = dist_force0 + 981.62883931057581/*m_lambda[COMPONENT_Y]*/ * coeff_bits_force0;

    if (cost0 < cost)
    {
      cost = cost0;
    }

    if (*fixed_filter_set_index > 0)
    {
      int len = 0;
      len += get_tb_length(*fixed_filter_set_index - 1, ALF_NUM_FIXED_FILTER_SETS);
      len += 1; //fixed filter flag pattern
      if (*fixed_filter_pattern > 0)
      {
        len += MAX_NUM_ALF_CLASSES;  //"fixed_filter_flag" for each class
      }
      cost += 981.62883931057581/*m_lambda[COMPONENT_Y]*/ * len;
    }

    if (cost <= cost_min)
    {
      cost_min = cost;
      num_filters_best = num_filters;
      best_pred_mode = pred_mode;
    }
    num_filters--;
  }

  dist = kvz_alf_derive_filter_coeffs(aps, channel, cov_frame, cov_merged, g_filter_indices[num_filters_best - 1], num_filters_best, error_force_0_coeff_tab);

  coeff_bits = kvz_alf_derive_filter_coefficients_prediction_mode(channel, g_filter_coeff_set, g_diff_filter_coeff, num_filters_best, &pred_mode);
  dist_force0 = get_dist_force_0(channel, num_filters_best, error_force_0_coeff_tab, coded_var_bins);
  coeff_bits_force0 = get_cost_filter_coeff_force_0(channel, g_filter_coeff_set, num_filters_best, coded_var_bins);

  cost = dist + 981.62883931057581/*m_lambda[COMPONENT_Y]*/ * coeff_bits;
  cost0 = dist_force0 + 981.62883931057581/*m_lambda[COMPONENT_Y]*/ * coeff_bits_force0;

  *num_luma_filters = num_filters_best;
  double dist_return;
  if (cost <= cost0)
  {
    dist_return = dist;
    *alf_luma_coeff_delta_flag = 0;
    *ui_coeff_bits = coeff_bits;
    *alf_luma_coeff_delta_prediction_flag = best_pred_mode;
  }
  else
  {
    
    dist_return = dist_force0;
    *alf_luma_coeff_delta_flag = 1;
    *ui_coeff_bits = coeff_bits_force0;
    memcpy(alf_luma_coeff_flag, coded_var_bins, sizeof(coded_var_bins));
    *alf_luma_coeff_delta_prediction_flag = 0;

    for (int var_ind = 0; var_ind < num_filters_best; var_ind++)
    {
      if (coded_var_bins[var_ind] == 0)
      {
        memset(g_filter_coeff_set[var_ind], 0, sizeof(int) * MAX_NUM_ALF_LUMA_COEFF);
      }
    }
  }

  for (int ind = 0; ind < *num_luma_filters; ++ind)
  {
    for (int i = 0; i < num_coeff; i++)
    {
      if (*alf_luma_coeff_delta_prediction_flag)
      {
        luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_diff_filter_coeff[ind][i];
      }
      else
      {
        luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_filter_coeff_set[ind][i];
      }
    }
  }

  memcpy(filter_coeff_delta_idx, g_filter_indices[num_filters_best - 1], sizeof(short) * MAX_NUM_ALF_CLASSES);
  ui_coeff_bits += get_non_filter_coeff_rate(aps);
  return dist_return;
}

double kvz_alf_derive_filter_coeffs(alf_aps *aps,
  channel_type channel,
  alf_covariance *cov,
  alf_covariance *covMerged,
  short* filter_indices,
  int numFilters,
  double errorTabForce0Coeff[MAX_NUM_ALF_CLASSES][2])
{
  int *fixed_filter_pattern = &aps->fixed_filter_pattern;
  int *fixed_filter_idx = aps->fixed_filter_idx;
  int *fixed_filter_set_index = &aps->fixed_filter_set_index;

  int num_filters = MAX_NUM_ALF_CLASSES;
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int *weights = channel == CHANNEL_TYPE_LUMA ? alf_weights_7 : alf_weights_5;

  double error = 0.0;
  
  alf_covariance tmp_cov = covMerged[MAX_NUM_ALF_CLASSES];
  *fixed_filter_set_index = 0;
  alf_covariance tmp_cov_ff = covMerged[MAX_NUM_ALF_CLASSES + 1];
  double factor = 1 << (ALF_NUM_BITS - 1);
  double error_min = 0;
  double error_min_per_class[MAX_NUM_ALF_CLASSES] = { 0 };
  double error_cur_set_per_class[MAX_NUM_ALF_CLASSES] = { 0 };
  int    fixed_filter_flag_per_class[MAX_NUM_ALF_CLASSES] = { 0 };
 
  for (int filter_set_idx = 0; filter_set_idx < ALF_NUM_FIXED_FILTER_SETS; filter_set_idx++)
  {
    double error_cur = 0;
    for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
    {
      int fixed_filter_idx = g_class_to_filter_mapping[filter_set_idx][class_idx];
      error_cur_set_per_class[class_idx] = calc_error_for_coeffs(cov[class_idx].ee, cov[class_idx].y, g_fixed_filter_set_coeff[fixed_filter_idx], MAX_NUM_ALF_LUMA_COEFF, ALF_NUM_BITS);

      if (error_cur_set_per_class[class_idx] >= 0)
      {
        error_cur_set_per_class[class_idx] = 0;
        fixed_filter_flag_per_class[class_idx] = 0;
      }
      else
      {
        error_cur += error_cur_set_per_class[class_idx];
        fixed_filter_flag_per_class[class_idx] = 1;
      }
    }

    if (error_cur < error_min)
    {
      memcpy(fixed_filter_idx, fixed_filter_flag_per_class, sizeof(fixed_filter_flag_per_class));
      *fixed_filter_set_index = filter_set_idx + 1;
      error_min = error_cur;
      memcpy(error_min_per_class, error_cur_set_per_class, sizeof(error_min_per_class));
    }
  }

  *fixed_filter_pattern = 0;
  if (*fixed_filter_set_index > 0)
  {
    for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
    {
      if (fixed_filter_idx[class_idx] == 0)
      {
        *fixed_filter_pattern = 1;
        break;
      }
    }
  }

  for( int filt_idx = 0; filt_idx < num_filters; filt_idx++ )
  {
    //tmp_cov.reset();
    tmp_cov.pix_acc = 0;
    memset(tmp_cov.y, 0, sizeof(*tmp_cov.y) * tmp_cov.num_coeff);
    for (int i = 0; i < tmp_cov.num_coeff; i++)
    {
      memset(tmp_cov.ee[i], 0, sizeof(*tmp_cov.ee[i]) * tmp_cov.num_coeff);
    }

    for( int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++ )
    {
      if( filter_indices[class_idx] == filt_idx )
      {
        //adjust stat
        tmp_cov_ff = cov[class_idx];
        if (*fixed_filter_set_index > 0 && fixed_filter_idx[class_idx] > 0)
        {
          int fixed_filter_idx = g_class_to_filter_mapping[*fixed_filter_set_index - 1][class_idx];
          tmp_cov_ff.pix_acc += error_min_per_class[class_idx];
          for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
          {
            double sum = 0;
            for (int j = 0; j < MAX_NUM_ALF_LUMA_COEFF; j++)
            {
              sum += tmp_cov_ff.ee[i][j] * g_fixed_filter_set_coeff[fixed_filter_idx][j];
            }
            sum /= factor;
            tmp_cov_ff.y[i] -= sum;
          }
        }

        //tmp_cov += tmp_cov_ff;
        for (int j = 0; j < tmp_cov.num_coeff; j++)
        {
          for (int i = 0; i < tmp_cov.num_coeff; i++)
          {
            tmp_cov.ee[j][i] += tmp_cov_ff.ee[j][i];
          }
          tmp_cov.y[j] += tmp_cov_ff.y[j];
        }
        tmp_cov.pix_acc += tmp_cov_ff.pix_acc;
      }
    }

    // Find coeffcients
    errorTabForce0Coeff[filt_idx][1] = tmp_cov.pix_acc + kvz_alf_derive_coeff_quant(channel, g_filter_coeff_quant, tmp_cov.ee, tmp_cov.y);

    errorTabForce0Coeff[filt_idx][0] = tmp_cov.pix_acc;
    error += errorTabForce0Coeff[filt_idx][1];

    // store coeff
    memcpy( g_filter_coeff_set[filt_idx], g_filter_coeff_quant, sizeof(*g_filter_coeff_quant) * MAX_NUM_ALF_LUMA_COEFF );
  }
  return error;
}

double kvz_alf_derive_coeff_quant(channel_type channel,
  int *filter_coeff_quant,
  double **ee,
  double *y)
{
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int *weights = channel == CHANNEL_TYPE_LUMA ? alf_weights_7 : alf_weights_5;
  const int factor = 1 << (ALF_NUM_BITS - 1);
  const int max_value = factor - 1;
  const int min_value = -factor + 1;

  static int filter_coeff_quant_mod[MAX_NUM_ALF_LUMA_COEFF];

  static double filter_coeff[MAX_NUM_ALF_LUMA_COEFF];

  gns_solve_by_chol(ee, y, filter_coeff, num_coeff);

  //roundFiltCoeff(filterCoeffQuant, filter_coeff, numCoeff, factor);
  for (int i = 0; i < num_coeff; i++)
  {
    int sign = filter_coeff[i] > 0 ? 1 : -1;
    filter_coeff_quant[i] = (int)(filter_coeff[i] * sign * factor + 0.5) * sign;
  }

  const int target_coeff_sum_int = 0;
  int quant_coeff_sum = 0;
  for (int i = 0; i < num_coeff; i++)
  {
    quant_coeff_sum += weights[i] * filter_coeff_quant[i];
  }
  int count = 0;
  while (quant_coeff_sum != target_coeff_sum_int && count < 10)
  {
    int sign = quant_coeff_sum > target_coeff_sum_int ? 1 : -1;
    int diff = (quant_coeff_sum - target_coeff_sum_int) * sign;

    double err_min = MAX_DOUBLE;
    int min_ind = -1;

    for (int k = 0; k < num_coeff; k++)
    {
      if (weights[k] <= diff)
      {
        memcpy(filter_coeff_quant_mod, filter_coeff_quant, sizeof(*filter_coeff_quant) * num_coeff);

        filter_coeff_quant_mod[k] -= sign;
        double error = calc_error_for_coeffs(ee, y, filter_coeff_quant_mod, num_coeff, ALF_NUM_BITS);

        if (error < err_min)
        {
          err_min = error;
          min_ind = k;
        }
      }
    }

    if (min_ind != -1)
    {
      filter_coeff_quant[min_ind] -= sign;
    }

    quant_coeff_sum = 0;
    for (int i = 0; i < num_coeff; i++)
    {
      quant_coeff_sum += weights[i] * filter_coeff_quant[i];
    }
    ++count;
  }
  if (count == 10)
  {
    memset(filter_coeff_quant, 0, sizeof(int) * num_coeff);
  }

  for (int i = 0; i < num_coeff - 1; i++)
  {
    filter_coeff_quant[i] = MIN(max_value, MAX(min_value, filter_coeff_quant[i]));
    filter_coeff[i] = filter_coeff_quant[i] / (double)(factor);
  }

  quant_coeff_sum = 0;
  for (int i = 0; i < num_coeff - 1; i++)
  {
    quant_coeff_sum += weights[i] * filter_coeff_quant[i];
    filter_coeff[i] = filter_coeff_quant[i] / (double)(factor);
  }
  filter_coeff_quant[num_coeff - 1] = -quant_coeff_sum;
  filter_coeff[num_coeff - 1] = filter_coeff_quant[num_coeff - 1] / (double)(factor);


  //Restrict the range of the center coefficient
  int max_value_center = (2 * factor - 1) - factor;
  int min_value_center = 0 - factor;

  filter_coeff_quant[num_coeff - 1] = MIN(max_value_center, MAX(min_value_center, filter_coeff_quant[num_coeff - 1]));
  filter_coeff[num_coeff - 1] = filter_coeff_quant[num_coeff - 1] / (double)(factor);

  int coeff_quant_adjust[MAX_NUM_ALF_LUMA_COEFF];
  int adjusted_total_coeff = (num_coeff - 1) << 1;

  count = 0;
  quant_coeff_sum += filter_coeff_quant[num_coeff - 1];
  while (quant_coeff_sum != target_coeff_sum_int && count < 15)
  {
    int sign = quant_coeff_sum > target_coeff_sum_int ? 1 : -1;
    int diff = (quant_coeff_sum - target_coeff_sum_int) * sign;

    if (diff > 4 * adjusted_total_coeff)     sign = sign * 8;
    else if (diff > 2 * adjusted_total_coeff)     sign = sign * 4;
    else if (diff >     adjusted_total_coeff)     sign = sign * 2;

    double err_min = MAX_DOUBLE;
    int    min_ind = -1;

    for (int k = 0; k < num_coeff - 1; k++)
    {
      memcpy(coeff_quant_adjust, filter_coeff_quant, sizeof(*filter_coeff_quant) * num_coeff);

      coeff_quant_adjust[k] -= sign;

      if (coeff_quant_adjust[k] <= max_value && coeff_quant_adjust[k] >= min_value)
      {
        double error = calc_error_for_coeffs(ee, y, coeff_quant_adjust, num_coeff, ALF_NUM_BITS);

        if (error < err_min)
        {
          err_min = error;
          min_ind = k;
        }
      }
    }

    if (min_ind != -1)
    {
      filter_coeff_quant[min_ind] -= sign;
      quant_coeff_sum -= (weights[min_ind] * sign);
    }

    ++count;
  }

  if (quant_coeff_sum != target_coeff_sum_int)
  {
    memset(filter_coeff_quant, 0, sizeof(int) * num_coeff);
  }

  for (int i = 0; i < num_coeff - 1; i++)
  {
    //CHECK(filter_coeff_quant[i] > max_value || filter_coeff_quant[i] < min_value, "filter_coeff_quant[i]>max_value || filter_coeff_quant[i]<min_value");
    filter_coeff[i] = filter_coeff_quant[i] / (double)(factor);
  }
  //CHECK(filter_coeff_quant[num_coeff - 1] > max_value_center || filter_coeff_quant[num_coeff - 1] < min_value_center, "filter_coeff_quant[num_coeff-1]>max_value_center || filter_coeff_quant[num_coeff-1]<min_value_center");
  filter_coeff[num_coeff - 1] = filter_coeff_quant[num_coeff - 1] / (double)(factor);


  double error = calc_error_for_coeffs(ee, y, filter_coeff_quant, num_coeff, ALF_NUM_BITS);
  return error;
}

double kvz_alf_derive_ctb_alf_enable_flags(encoder_state_t * const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  const int i_shape_idx,
  double *dist_unfilter,
  const int num_classes)
{
  kvz_config cfg = state->encoder_control->cfg;
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  int ctu_index = lcu->index;

  const kvz_pixel comp_id_first = is_luma ? COMPONENT_Y : COMPONENT_Cb;
  const kvz_pixel comp_id_last = is_luma ? COMPONENT_Y : COMPONENT_Cr;

  int num_coeff = is_luma ? 13 : 7;

  double cost = 0;
  *dist_unfilter = 0;

  //setEnableFlag(m_alfSliceParamTemp, channel, true);
  if (is_luma) {
    g_alf_slice_aps_temp.enabled_flag[COMPONENT_Y] = 1;
  }
  else {
    g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cb] = 1;
    g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cr] = 1;
  }

  // Tehd‰‰n erikseen kaikille CTU:lle kerran
  kvz_alf_reconstruct_coeff(state, &g_alf_slice_aps_temp, channel, true, is_luma);

  //Samoin t‰m‰
  for (int class_idx = 0; class_idx < (is_luma ? MAX_NUM_ALF_CLASSES : 1); class_idx++)
  {
    for (int i = 0; i < (is_luma ? MAX_NUM_ALF_LUMA_COEFF : MAX_NUM_ALF_CHROMA_COEFF); i++)
    {
      g_filter_coeff_set[class_idx][i] = is_luma ? g_coeff_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + i] : g_chroma_coeff_final[i];
    }
  }

  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    for (int comp_id = comp_id_first; comp_id <= comp_id_last; comp_id++)
    {

      //double distUnfilterCtu = getUnfilteredDistortion(g_alf_covariance[comp_id][i_shape_idx][ctu_idx], numClasses);

      double dist_unfilter_ctu = 0;
      alf_covariance *cov = g_alf_covariance[comp_id][i_shape_idx][ctu_index];
      for (int classIdx = 0; classIdx < num_classes; classIdx++)
      {
        dist_unfilter_ctu += cov[classIdx].pix_acc;
      }

      //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
      //m_CABACEstimator->resetBits();
      //m_ctuEnableFlag[compID][ctu_idx] = 1;
      g_ctu_enable_flag[comp_id][ctu_idx] = 1;

      //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctu_idx, comp_id, &m_alfSliceParamTemp);
      code_alf_ctu_enable_flag(state, lcu, ctu_idx, comp_id, &g_alf_slice_aps_temp);

      //double costOn = distUnfilterCtu + getFilteredDistortion(g_alf_covariance[comp_id][iShapeIdx][ctu_idx], numClasses, m_alfSliceParamTemp.numLumaFilters - 1, numCoeff);
      double dist = 0;
      for (int class_idx = 0; class_idx < num_classes; class_idx++)
      {
        double **ee = cov[class_idx].ee;
        double *y = cov[class_idx].y;
        int *coeff = g_filter_coeff_set[class_idx];

        double factor = 1 << (ALF_NUM_BITS/*m_NUM_BITS*/ - 1);
        double error = 0;

        for (int i = 0; i < num_coeff; i++)   //diagonal
        {
          double sum = 0;
          for (int j = i + 1; j < num_coeff; j++)
          {
            // E[j][i] = E[i][j], sum will be multiplied by 2 later
            sum += ee[i][j] * coeff[j];
          }
          error += ((ee[i][i] * coeff[i] + sum * 2) / factor - 2 * y[i]) * coeff[i];
        }

        dist += error / factor; //calcErrorForCoeffs(cov[classIdx].E, cov[classIdx].y, m_filterCoeffSet[classIdx], numCoeff, m_NUM_BITS);
      }

      double cost_on = dist_unfilter_ctu + dist;

      //const double ctuLambda = m_lambda[comp_id];
      double ctu_lambda = g_lambda[comp_id];

      cost_on += ctu_lambda * frac_bits_scale * 0/*(double)m_CABACEstimator->getEstFracBits()*/;
      //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());

      //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
      //m_CABACEstimator->resetBits();

      //m_ctuEnableFlag[comp_id][ctu_idx] = 0;
      //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctu_idx, comp_id, &m_alfSliceParamTemp);

      double cost_off = dist_unfilter_ctu + ctu_lambda * frac_bits_scale * 0/*(double)m_CABACEstimator->getEstFracBits()*/;

      if (cost_on < cost_off)
      {
        cost += cost_on;
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
        g_ctu_enable_flag[comp_id][ctu_idx] = 1;
      }
      else
      {
        cost += cost_off;
        g_ctu_enable_flag[comp_id][ctu_idx] = 0;
        *dist_unfilter += dist_unfilter_ctu;
      }
    }
  }

  if (!is_luma)
  {
    //Slice flag
    //setEnableFlag(m_alfSliceParamTemp, channel, m_ctuEnableFlag);
    const alf_component_id compIDFirst = COMPONENT_Cb;
    const alf_component_id compIDLast = COMPONENT_Cr;
    for (int compId = compIDFirst; compId <= compIDLast; compId++)
    {
      g_alf_slice_aps_temp.enabled_flag[compId] = false;
      for (int i = 0; i < num_ctus_in_pic; i++)
      {
        if (g_ctu_enable_flag[compId][i])
        {
          g_alf_slice_aps_temp.enabled_flag[compId] = true;
          break;
        }
      }
    }

    //const int alfChromaIdc = m_alfSliceParamTemp.enabledFlag[COMPONENT_Cb] * 2 + m_alfSliceParamTemp.enabledFlag[COMPONENT_Cr];
    int alf_chroma_idc = g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cb] * 2 + g_alf_slice_aps_temp.enabled_flag[COMPONENT_Cr];
    cost += length_truncated_unary(alf_chroma_idc, 3) * g_lambda[channel];
  }

  return cost;
}

//----------------------------------------------------------------------

//-------------------------cabac writer functions------------------------

void code_alf_ctu_enable_flag(encoder_state_t * const state,
  const lcu_order_element_t *const lcu,
  uint32_t ctu_rs_addr,
  int8_t component_id,
  alf_aps *aps)
{
  const encoder_control_t * const encoder = state->encoder_control;
  cabac_data_t * const cabac = &state->cabac;
  
  const bool alf_component_enabled = (aps != NULL) ? aps->enabled_flag[component_id] : state->slice->tile_group_alf_enabled_flag[component_id];

  if (encoder->cfg.alf_enable && alf_component_enabled)
  {
    int frame_width_in_ctus = state->tile->frame->width_in_lcu;
    //int ry = ctu_rs_addr / frame_width_in_ctus;
    //int rx = ctu_rs_addr - ry * frame_width_in_ctus;
    const uint32_t curSliceIdx = state->slice->id;
    //const Position      pos(rx * cs.pcv->maxCUWidth, ry * cs.pcv->maxCUHeight);
    //const uint32_t          curSliceIdx = cs.slice->getIndependentSliceIdx();
    //const uint32_t      curTileIdx = cs.picture->brickMap->getBrickIdxRsMap(pos);

    //bool leftAvail = cs.getCURestricted(pos.offset(-(int)pcv.maxCUWidth, 0), curSliceIdx, curTileIdx, CH_L) ? true : false;
    //bool aboveAvail = cs.getCURestricted(pos.offset(0, -(int)pcv.maxCUHeight), curSliceIdx, curTileIdx, CH_L) ? true : false;
    bool left_avail = lcu->left ? 1 : 0;
    bool above_avail = lcu->above ? 1 : 0;

    int left_ctu_addr = left_avail ? ctu_rs_addr - 1 : -1;
    int above_ctu_addr = above_avail ? ctu_rs_addr - frame_width_in_ctus : -1;

    uint8_t* ctb_alf_flag = g_ctu_enable_flag[component_id];
    uint8_t flag = ctb_alf_flag[ctu_rs_addr];

    int ctx = 0;
    ctx += left_ctu_addr > -1 ? (ctb_alf_flag[left_ctu_addr] ? 1 : 0) : 0;
    ctx += above_ctu_addr > -1 ? (ctb_alf_flag[above_ctu_addr] ? 1 : 0) : 0;

    //m_BinEncoder.encodeBin(ctbAlfFlag[ctuRsAddr], Ctx::ctbAlfFlag(compIdx * 3 + ctx));
    cabac->cur_ctx = &(cabac->ctx.alf_ctb_flag_model[component_id * 3 + ctx]);
    CABAC_BIN(cabac, ctb_alf_flag[ctu_rs_addr], "alf_ctb_flag");
  }
}

void code_alf_ctu_filter_index(encoder_state_t * const state,
  const lcu_order_element_t *const lcu,
  uint32_t ctu_rs_addr,
  bool alf_enable_luma)
{
  bitstream_t *stream = &state->stream;

  const encoder_control_t * const encoder = state->encoder_control;
  cabac_data_t * const cabac = &state->cabac;

  if (!encoder->cfg.alf_enable || !alf_enable_luma)//(!cs.sps->getALFEnabledFlag()) || (!alfEnableLuma))
  {
    return;
  }

  //uint8_t* ctbAlfFlag = cs.slice->getPic()->getAlfCtuEnableFlag(COMPONENT_Y);
  if (!g_ctu_enable_flag[COMPONENT_Y][ctu_rs_addr])
  {
    return;
  }

  const unsigned filter_set_idx = g_alf_ctb_filter_index[ctu_rs_addr];
  unsigned num_aps = state->slice->tile_group_num_aps; /*cs.slice->getTileGroupNumAps()*/;
  unsigned num_available_filt_sets = num_aps + ALF_NUM_FIXED_FILTER_SETS;
  if (num_available_filt_sets > ALF_NUM_FIXED_FILTER_SETS)
  {
    int use_latest_filt = (filter_set_idx == ALF_NUM_FIXED_FILTER_SETS) ? 1 : 0;
    //m_BinEncoder.encodeBin(use_latest_filt, Ctx::AlfUseLatestFilt());
    cabac->cur_ctx = &(cabac->ctx.alf_latest_filt);
    CABAC_BIN(cabac, use_latest_filt, "use_latest_filt");  

    if (!use_latest_filt)
    {

      if (num_aps == 1)
      {
        //CHECK(filter_set_idx >= alf_num_fixed_filter_sets, "fixed set numavail < num_fixed");
        assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //Fixed set numavail

        //xWriteTruncBinCode(filter_set_idx, NUM_FIXED_FILTER_SETS);
        kvz_cabac_encode_trunc_bin(cabac, filter_set_idx, ALF_NUM_FIXED_FILTER_SETS);
      }
      else
      {
        int use_temporal_filt = (filter_set_idx > ALF_NUM_FIXED_FILTER_SETS) ? 1 : 0;
        //m_BinEncoder.encodeBin(use_temporal_filt, Ctx::AlfUseTemporalFilt());
        cabac->cur_ctx = &(cabac->ctx.alf_temporal_filt);
        CABAC_BIN(cabac, use_temporal_filt, "use_temporal_filt");

        if (use_temporal_filt)
        {
          assert((filter_set_idx - (ALF_NUM_FIXED_FILTER_SETS + 1)) < (num_aps - 1)); //Temporal non-latest set
          //xWriteTruncBinCode(filter_set_idx - (NUM_FIXED_FILTER_SETS + 1), num_available_filt_sets - (NUM_FIXED_FILTER_SETS + 1));
          kvz_cabac_encode_trunc_bin(cabac, filter_set_idx - (ALF_NUM_FIXED_FILTER_SETS + 1), num_available_filt_sets - (ALF_NUM_FIXED_FILTER_SETS + 1));
        }
        else
        {
          assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //Fixed set larger than temporal
          //xWriteTruncBinCode(filter_set_idx, NUM_FIXED_FILTER_SETS);
          kvz_cabac_encode_trunc_bin(cabac, filter_set_idx, ALF_NUM_FIXED_FILTER_SETS);
        }
      }
    }
  }
  else
  {
    assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //Fixed set numavail < num_fixed
    //xWriteTruncBinCode(filterSetIdx, NUM_FIXED_FILTER_SETS);
    kvz_cabac_encode_trunc_bin(cabac, filter_set_idx, ALF_NUM_FIXED_FILTER_SETS);
  }
}

//---------------------------------------------------------------------

//-------------------------CTU functions--------------------------------

void kvz_alf_process(encoder_state_t const *state,
  const lcu_order_element_t const *lcu)
{
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;
  cabac_data_t * const cabac = &state->cabac;
  //cabac_data_t ctx = &state->cabac.ctx.alf_ctb_enable_flag;
  bool alf_ctb_flag = &state->cabac.ctx.alf_ctb_flag_model[COMPONENT_Y];
  //if (!cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Y) && !cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Cb) && !cs.slice->getTileGroupAlfEnabledFlag(COMPONENT_Cr))
  if (!state->slice->tile_group_alf_enabled_flag[COMPONENT_Y] && !state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] && !state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr])
  {
    return;
  }

  // set clipping range
  // done in INIT
  //m_clp_rngs = cs.slice->getClpRngs();

  // set CTU enable flags
  /* Ei tarpeellinen
  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++) {
      g_ctu_enable_flag[comp_idx][ctu_idx] = //cs.picture->getAlfCtuEnableFlag( compIdx );
    }
  }*/

  kvz_alf_reconstruct_coeff_aps(state,
    true,
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] || state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr],
    false);

  //PelUnitBuf recYuv = cs.getRecoBuf();
  kvz_picture *rec_yuv = state->tile->frame->rec;
  //m_tempBuf.copyFrom(recYuv);
  //PelUnitBuf tmpYuv = m_tempBuf.getBuf(cs.area);
  //tmpYuv.extendBorderPel(MAX_ALF_FILTER_LENGTH >> 1);

  //const PreCalcValues& pcv = *cs.pcv;

  int luma_height = state->tile->frame->height;
  int luma_width = state->tile->frame->width;

  int ctu_idx = 0;
  int max_cu_width = LCU_WIDTH;
  //for (int yPos = 0; yPos < pcv.lumaHeight; yPos += pcv.maxCUHeight)
  for (int y_pos = 0; y_pos < luma_height; y_pos += max_cu_width)
  {
    //for( int xPos = 0; xPos < pcv.lumaWidth; xPos += pcv.maxCUWidth )
    for (int x_pos = 0; x_pos < luma_width; x_pos += max_cu_width)
    {
      //const int width = (xPos + pcv.maxCUWidth > pcv.lumaWidth) ? (pcv.lumaWidth - xPos) : pcv.maxCUWidth;
      //const int height = (yPos + pcv.maxCUHeight > pcv.lumaHeight) ? (pcv.lumaHeight - yPos) : pcv.maxCUHeight;
      const int width = (x_pos + max_cu_width > luma_width) ? (luma_width - x_pos) : max_cu_width;
      const int height = (y_pos + max_cu_width > luma_height) ? (luma_height - y_pos) : max_cu_width;
      
      //const UnitArea area(cs.area.chromaFormat, Area(xPos, yPos, width, height));
      //x_pos, y_pos, width, height

      if(g_ctu_enable_flag[COMPONENT_Y][ctu_idx])
      {
        //Area blk(xPos, yPos, width, height);
        //deriveClassification(m_classifier, tmpYuv.get(COMPONENT_Y), blk);
        kvz_alf_derive_classification(state, lcu, width, height, x_pos, y_pos);

        //Area blkPCM(xPos, yPos, width, height);
        kvz_alf_reset_pcm_blk_class_info(state, lcu, width, height, x_pos, y_pos);

        short filter_set_index = g_alf_ctb_filter_index[ctu_idx];
        short *coeff;
        if (filter_set_index >= ALF_NUM_FIXED_FILTER_SETS)
        {
          coeff = g_coeff_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
        }
        else
        {
          coeff = g_fixed_filter_set_coeff_dec[filter_set_index];
        }
        
        //filter_blk(m_classifier, recYuv, tmpYuv, blk, COMPONENT_Y, coeff, clp_rngs.comp[COMPONENT_Y], cs);
        kvz_alf_filter_block(state, lcu, coeff, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y, width, height, x_pos, y_pos);

      }

      for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
      {
        alf_component_id comp_id = comp_idx;
        //const int chromaScaleX = getComponentScaleX(comp_id, tmpYuv.chromaFormat);
        //const int chromaScaleY = getComponentScaleY(comp_id, tmpYuv.chromaFormat);
        const int chroma_scale_x = ((comp_id == COMPONENT_Y) || (chroma_fmt == KVZ_CSP_444)) ? 0 : 1;
        const int chroma_scale_y = ((comp_id == COMPONENT_Y) || (chroma_fmt != KVZ_CSP_420)) ? 0 : 1;
        
        if (g_ctu_enable_flag[comp_idx][ctu_idx])
        {
          //Area blk(xPos >> chroma_scale_x, yPos >> chroma_scale_y, width >> chroma_scale_x, height >> chroma_scale_y);
          //m_filter5x5Blk(m_classifier, recYuv, tmpYuv, blk, comp_id, m_chromaCoeffFinal, clp_rngs.comp[comp_idx], cs);
          kvz_alf_filter_block(state, lcu, g_chroma_coeff_final, g_clp_rngs.comp[comp_idx], comp_idx, 
                               width >> chroma_scale_x, height >> chroma_scale_y, 
                               x_pos >> chroma_scale_x, y_pos >> chroma_scale_y);
        }
      }
      ctu_idx++;
    }
  }
}

void kvz_alf_reconstruct_coeff_aps(encoder_state_t *const state, bool luma, bool chroma, bool is_rdo)
{
  //luma
  //APS** aps = cs.slice->getAPSs();
  alf_aps** aps = &state->slice->aps;
  //AlfSliceParam alfSliceParamTmp;

  //APS* curAPS;
  alf_aps* cur_aps;

  if (luma)
  {
    for (int i = 0; i < state->slice->tile_group_num_aps /* 1,  cs.slice->getTileGroupNumAps()*/; i++) {
      //int apsIdx = cs.slice->getTileGroupApsIdLuma()[i];
      int aps_idx = state->slice->tile_group_luma_aps_id[i];

      //curAPS = aps[apsIdx];
      cur_aps = aps[aps_idx];

      //CHECK(curAPS == NULL, "invalid APS");
      //alfSliceParamTmp = curAPS->getAlfAPSParam();

      kvz_alf_reconstruct_coeff(state, cur_aps, CHANNEL_TYPE_LUMA, is_rdo, true);
      memcpy(g_coeff_aps_luma[i], g_coeff_final, sizeof(g_coeff_final));
    }
  }

  //chroma
  if (chroma)
  {
    //int apsIdxChroma = cs.slice->getTileGroupApsIdChroma();
    int aps_idx_chroma = state->slice->tile_group_chroma_aps_id;
    //curAPS = aps[apsIdxChroma];
    cur_aps = aps[aps_idx_chroma];
    //alfSliceParamTmp = curAPS->getAlfAPSParam();
    kvz_alf_reconstruct_coeff(state, cur_aps, CHANNEL_TYPE_CHROMA, is_rdo, true);
  }
}

//void reconstructCoeff(AlfSliceParam& alfSliceParam, ChannelType channel, const bool isRdo, const bool isRedo)
void kvz_alf_reconstruct_coeff(encoder_state_t *const state,
  alf_aps *aps,
  channel_type channel,
  const bool is_rdo,
  const bool is_redo)
{
  //alf_aps *aps = state->slice->aps;
  int *num_luma_filters = &aps->num_luma_filters;
  bool *alf_luma_coeff_delta_prediction_flag = &aps->alf_luma_coeff_delta_prediction_flag;
  short *luma_coeff = aps->luma_coeff;
  short *chroma_coeff = aps->chroma_coeff;
  short *filter_coeff_delta_idx = aps->filter_coeff_delta_idx;
  int *fixed_filter_set_index = &aps->fixed_filter_set_index;

  int factor = is_rdo ? 0 : (1 << (ALF_NUM_BITS - 1));

  alf_filter_type filter_type = channel == CHANNEL_TYPE_LUMA ? ALF_FILTER_7X7 : ALF_FILTER_5X5;
  int num_classes = channel == CHANNEL_TYPE_LUMA ? MAX_NUM_ALF_CLASSES : 1;
  int num_coeff = filter_type == ALF_FILTER_5X5 ? 7 : 13;
  int num_coeff_minus1 = num_coeff - 1;
  int num_filters = channel == CHANNEL_TYPE_LUMA ? *num_luma_filters : 1;
  short* coeff = channel == CHANNEL_TYPE_LUMA ? luma_coeff : chroma_coeff;

  if (*alf_luma_coeff_delta_prediction_flag && channel == CHANNEL_TYPE_LUMA)
  {
    for (int i = 1; i < num_filters; i++)
    {
      for (int j = 0; j < num_coeff_minus1; j++)
      {
        coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] += coeff[(i - 1) * MAX_NUM_ALF_LUMA_COEFF + j];
      }
    }
  }

  for (int filter_idx = 0; filter_idx < num_filters; filter_idx++)
  {
    int sum = 0;
    for (int i = 0; i < num_coeff_minus1; i++)
    {
      sum += (coeff[filter_idx* MAX_NUM_ALF_LUMA_COEFF + i] << 1);
    }
    coeff[filter_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor - sum;
  }

  if ( channel == CHANNEL_TYPE_CHROMA )
  {
    int sum = 0;
    for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
    {
      g_chroma_coeff_final[coeff_idx] = chroma_coeff[coeff_idx];
      sum += (g_chroma_coeff_final[coeff_idx] << 1);
    }
    g_chroma_coeff_final[num_coeff_minus1] = factor - sum;

    return;
  }

  /*
  g_fixed_filter_idx[22] = 1;
  g_fixed_filter_idx[23] = 1;
  */

  for (int class_idx = 0; class_idx < num_classes; class_idx++)
  {
    int filter_idx = filter_coeff_delta_idx[class_idx];
    int sum = 0;
    int fixed_filter_idx = *fixed_filter_set_index; //13
    if (fixed_filter_idx > 0 && aps->fixed_filter_idx[class_idx] > 0)
    {
      fixed_filter_idx = g_class_to_filter_mapping[fixed_filter_idx - 1][class_idx];
    }
    else
    {
      fixed_filter_idx = -1;
    }
    for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
    {
      g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] = coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx];
      //fixed filter
      if (fixed_filter_idx >= 0)
      {
        g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] += g_fixed_filter_set_coeff[fixed_filter_idx][coeff_idx];
      }
      sum += (g_coeff_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + coeff_idx] << 1);
    }
    g_coeff_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor - sum;
  }

  if (is_redo && state->cabac.ctx.alf_luma_coeff_delta_prediction_flag.state[0])
  {
    for (int i = num_filters - 1; i > 0; i--)
    {
      for (int j = 0; j < num_coeff_minus1; j++)
      {
        coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] = coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] - coeff[(i - 1) * MAX_NUM_ALF_LUMA_COEFF + j];
      }
    }
  }
}



//(const int picWidth, const int picHeight, const ChromaFormat format, const int maxCUWidth, const int maxCUHeight, const int maxCUDepth, const int inputBitDepth[MAX_NUM_CHANNEL_TYPE])
void kvz_alf_create(encoder_state_t const *state,
  const lcu_order_element_t const *lcu) //,
  //const int format)
{
  //memcpy(m_inputBitDepth, inputBitDepth, sizeof(m_inputBitDepth));

  int ***laplacian = state->tile->frame->alf_info->laplacian;
  alf_classifier ***classifier = &state->tile->frame->alf_info->classifier;
  bool *created = &state->tile->frame->alf_info->created;

  const int pic_width = state->tile->frame->width;
  const int pic_height = state->tile->frame->height;
  const int max_cu_width = LCU_WIDTH; //128
  const int max_cu_height = LCU_WIDTH; //128

  int32_t heigth = state->tile->frame->height_in_lcu;
  int32_t width = state->tile->frame->width_in_lcu;
  num_ctus_in_pic = heigth * width;

  /*
  int num_ctus_in_width = (pic_width / max_cu_width) + ((pic_width % max_cu_width) ? 1 : 0);
  int num_ctus_in_height = (pic_height / max_cu_height) + ((pic_height % max_cu_height) ? 1 : 0);
  int num_ctus_in_pic = num_ctus_in_height * num_ctus_in_width;
  */
  //int filter_shapes[CHANNEL_TYPE_LUMA].push_back(AlfFilterShape(7));
  //filter_shapes[CHANNEL_TYPE_CHROMA].push_back(AlfFilterShape(5));

  if (*created)
  {
    return;
  }

  //m_tempBuf.destroy();
  //m_tempBuf.create(format, Area(0, 0, pic_width, pic_height), max_cu_width, MAX_ALF_FILTER_LENGTH >> 1, 0, false);

  // Laplacian based activity
  for (int i = 0; i < NUM_DIRECTIONS; i++)
  {
    //if (laplacian[i] == NULL)
    //{
      laplacian[i] = malloc((CLASSIFICATION_BLK_SIZE + 5) * sizeof(int*));
      for (int y = 0; y < CLASSIFICATION_BLK_SIZE + 5; y++)
      {
        laplacian[i][y] = malloc((CLASSIFICATION_BLK_SIZE + 5) * sizeof(int));
      }
    //}
  }

  // Classification
  //if (*classifier == NULL)
  //{
    *classifier = malloc(pic_height * sizeof(**classifier));
    for (int i = 0; i < pic_height; i++)
    {
      (*classifier)[i] = malloc(pic_width * sizeof(alf_classifier));
    }
  //}

  for (int filter_set_index = 0; filter_set_index < ALF_NUM_FIXED_FILTER_SETS; filter_set_index++)
  {
    for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
    {
      int fixed_filter_idx = g_class_to_filter_mapping[filter_set_index][class_idx];
      int sum = 0;
      for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF - 1; i++)
      {
        sum += (g_fixed_filter_set_coeff[fixed_filter_idx][i] << 1);
        g_fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + i] = g_fixed_filter_set_coeff[fixed_filter_idx][i];
      }
      g_fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + MAX_NUM_ALF_LUMA_COEFF - 1] = (1 << (ALF_NUM_BITS - 1)) - sum;
    }
  }
  *created = true;
}


void kvz_alf_destroy(encoder_state_t const *state)
{
  alf_info_t *alf = state->tile->frame->alf_info;
  int ***laplacian = state->tile->frame->alf_info->laplacian;
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;
  int32_t pic_height = state->tile->frame->rec->height;

  if (!alf->created)
  {
    return;
  }

  for (int i = 0; i < NUM_DIRECTIONS; i++)
  {
    if (laplacian[i])
    {
      for (int y = 0; y < CLASSIFICATION_BLK_SIZE + 5; y++)
      {
        if (laplacian[i][y] != NULL) {
          FREE_POINTER(laplacian[i][y]);
        }
      }
      FREE_POINTER(laplacian[i]);
    }
  }

  if (classifier)
  {
    for (int i = 0; i < pic_height; i++)
    {
      FREE_POINTER(classifier[i]);
    }
    FREE_POINTER(classifier);
  }

  alf->created = false;
}


void kvz_alf_derive_classification(encoder_state_t *const state,
  const lcu_order_element_t * const lcu,
  const int width,
  const int height,
  int x_pos,
  int y_pos)//,
  //alf_classifier** classifier)
{
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;

  int max_height = y_pos + height;
  int max_width = x_pos + width;

  for (int i = y_pos; i < max_height; i += CLASSIFICATION_BLK_SIZE)
  {
    int n_height = MIN(i + CLASSIFICATION_BLK_SIZE, max_height) - i;

    for (int j = x_pos; j < max_width; j += CLASSIFICATION_BLK_SIZE)
    {
      int n_width = MIN(j + CLASSIFICATION_BLK_SIZE, max_width) - j;
      kvz_alf_derive_classification_blk(state, lcu, classifier, alf_input_bit_depth[CHANNEL_TYPE_LUMA] + 4, n_height, n_width, j, i);
    }
  }
}

//Turha jos PCM on pois p‰‰lt‰.
void kvz_alf_reset_pcm_blk_class_info(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  const int width,
  const int height,
  int x_pos,
  int y_pos)
{
  if (!state->encoder_control->cfg.implicit_rdpcm) //!cs.sps->getPCMFilterDisableFlag()
  {
    return;
  }

  alf_classifier **classifier = state->tile->frame->alf_info->classifier;

  int max_height = y_pos + height;
  int max_width = x_pos + width;
  const int cls_size_y = 4;
  const int cls_size_x = 4;
  int class_idx = ALF_UNUSED_CLASS_IDX;
  int transpose_idx = ALF_UNUSED_TRANSPOSE_IDX;

  for (int i = y_pos; i < max_height; i += CLASSIFICATION_BLK_SIZE)
  {
    int n_height = MIN(i + CLASSIFICATION_BLK_SIZE, max_height) - i;

    for (int j = x_pos; j < max_width; j += CLASSIFICATION_BLK_SIZE)
    {
      int n_width = MIN(j + CLASSIFICATION_BLK_SIZE, max_width) - j;
      int pos_x = j;
      int pos_y = i;

      for (int subi = 0; subi < n_height; subi += cls_size_y)
      {
        for (int subj = 0; subj < n_width; subj += cls_size_x)
        {
          int y_offset = subi + pos_y;
          int x_offset = subj + pos_x;
          //Position pos(xOffset, yOffset);

          //const CodingUnit* cu = cs.getCU(pos, CH_L);
          if (1) //cu->ipcm ///ipcm_flag = 0
          {
            alf_classifier *cl0 = classifier[y_offset] + x_offset;
            alf_classifier *cl1 = classifier[y_offset + 1] + x_offset;
            alf_classifier *cl2 = classifier[y_offset + 2] + x_offset;
            alf_classifier *cl3 = classifier[y_offset + 3] + x_offset;

            cl0[0].class_idx = cl0[1].class_idx = cl0[2].class_idx = cl0[3].class_idx = 
            cl1[0].class_idx = cl1[1].class_idx = cl1[2].class_idx = cl1[3].class_idx = 
            cl2[0].class_idx = cl2[1].class_idx = cl2[2].class_idx = cl2[3].class_idx = 
            cl3[0].class_idx = cl3[1].class_idx = cl3[2].class_idx = cl3[3].class_idx = class_idx;

            cl0[0].transpose_idx = cl0[1].transpose_idx = cl0[2].transpose_idx = cl0[3].transpose_idx = 
            cl1[0].transpose_idx = cl1[1].transpose_idx = cl1[2].transpose_idx = cl1[3].transpose_idx = 
            cl2[0].transpose_idx = cl2[1].transpose_idx = cl2[2].transpose_idx = cl2[3].transpose_idx = 
            cl3[0].transpose_idx = cl3[1].transpose_idx = cl3[2].transpose_idx = cl3[3].transpose_idx = transpose_idx;
          }
        }
      }
    }
  }
}

void kvz_alf_derive_classification_blk(encoder_state_t * const state,
  const lcu_order_element_t * const lcu,
  alf_classifier** classifier,
  const int shift,
  const int n_height,
  const int n_width,
  const int x,
  const int y)
{
  videoframe_t* const frame = state->tile->frame;
  int ***laplacian = state->tile->frame->alf_info->laplacian;

  static const int th[16] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
  const int stride = frame->rec->stride;
  kvz_pixel *src = state->tile->frame->rec->y;
  const int max_activity = 15;

  int fl = 2;
  int fl_p1 = fl + 1;
  int fl2 = 2 * fl;

  int main_direction, secondary_direction, dir_temp_hv, dir_temp_d;
  int pix_y;

  kvz_pixel height = n_height + fl2;
  kvz_pixel width = n_width + fl2;
  int pos_x = x;
  int pos_y = y;
  int start_height = pos_y - fl_p1;

  for (int i = 0; i < height; i += 2)
  {
    int yoffset = (i + 1 + start_height) * stride - fl_p1;
    const kvz_pixel *src0 = &src[yoffset - stride];
    const kvz_pixel *src1 = &src[yoffset];
    const kvz_pixel *src2 = &src[yoffset + stride];
    const kvz_pixel *src3 = &src[yoffset + stride * 2];

    int *p_y_ver = laplacian[ALF_VER][i];
    int *p_y_hor = laplacian[ALF_HOR][i];
    int *p_y_dig0 = laplacian[ALF_DIAG0][i];
    int *p_y_dig1 = laplacian[ALF_DIAG1][i];

    for (int j = 0; j < width; j += 2)
    {
      pix_y = j + 1 + pos_x;
      const kvz_pixel *p_y = src1 + pix_y;
      const kvz_pixel *p_y_down = src0 + pix_y;
      const kvz_pixel *p_y_up = src2 + pix_y;
      const kvz_pixel *p_y_up2 = src3 + pix_y;

      const int16_t y0 = p_y[0] << 1;
      const int16_t y_up1 = p_y_up[1] << 1;

      p_y_ver[j] = abs(y0 - p_y_down[0] - p_y_up[0]) + abs(y_up1 - p_y[1] - p_y_up2[1]);
      p_y_hor[j] = abs(y0 - p_y[1] - p_y[-1]) + abs(y_up1 - p_y_up[2] - p_y_up[0]);
      p_y_dig0[j] = abs(y0 - p_y_down[-1] - p_y_up[1]) + abs(y_up1 - p_y[0] - p_y_up2[2]);
      p_y_dig1[j] = abs(y0 - p_y_up[-1] - p_y_down[1]) + abs(y_up1 - p_y_up2[0] - p_y[2]);

      if (j > 4 && (j - 6) % 4 == 0)
      {
        int j_m6 = j - 6;
        int j_m4 = j - 4;
        int j_m2 = j - 2;

        p_y_ver[j_m6] += p_y_ver[j_m4] + p_y_ver[j_m2] + p_y_ver[j];
        p_y_hor[j_m6] += p_y_hor[j_m4] + p_y_hor[j_m2] + p_y_hor[j];
        p_y_dig0[j_m6] += p_y_dig0[j_m4] + p_y_dig0[j_m2] + p_y_dig0[j];
        p_y_dig1[j_m6] += p_y_dig1[j_m4] + p_y_dig1[j_m2] + p_y_dig1[j];
      }
    }
  }

  // classification block size
  const int cls_size_y = 4;
  const int cls_size_x = 4;

  //for (int i = 0; i < blk.height; i += cls_size_y)
  for (int i = 0; i < n_height; i += cls_size_y)
  {
    int* p_y_ver = laplacian[ALF_VER][i];
    int* p_y_ver2 = laplacian[ALF_VER][i + 2];
    int* p_y_ver4 = laplacian[ALF_VER][i + 4];
    int* p_y_ver6 = laplacian[ALF_VER][i + 6];

    int* p_y_hor = laplacian[ALF_HOR][i];
    int* p_y_hor2 = laplacian[ALF_HOR][i + 2];
    int* p_y_hor4 = laplacian[ALF_HOR][i + 4];
    int* p_y_hor6 = laplacian[ALF_HOR][i + 6];

    int* p_y_dig0 = laplacian[ALF_DIAG0][i];
    int* p_y_dig02 = laplacian[ALF_DIAG0][i + 2];
    int* p_y_dig04 = laplacian[ALF_DIAG0][i + 4];
    int* p_y_dig06 = laplacian[ALF_DIAG0][i + 6];

    int* p_y_dig1 = laplacian[ALF_DIAG1][i];
    int* p_y_dig12 = laplacian[ALF_DIAG1][i + 2];
    int* p_y_dig14 = laplacian[ALF_DIAG1][i + 4];
    int* p_y_dig16 = laplacian[ALF_DIAG1][i + 6];

    //for (int j = 0; j < blk.width; j += cls_size_x)
    for (int j = 0; j < n_width; j += cls_size_x)
    {
      int sum_v = p_y_ver[j] + p_y_ver2[j] + p_y_ver4[j] + p_y_ver6[j];
      int sum_h = p_y_hor[j] + p_y_hor2[j] + p_y_hor4[j] + p_y_hor6[j];
      int sum_d0 = p_y_dig0[j] + p_y_dig02[j] + p_y_dig04[j] + p_y_dig06[j];
      int sum_d1 = p_y_dig1[j] + p_y_dig12[j] + p_y_dig14[j] + p_y_dig16[j];

      int temp_act = sum_v + sum_h;
      int activity = (kvz_pixel)alf_clip3(0, max_activity, (temp_act * 64) >> shift);
      int class_idx = th[activity];

      int hv1, hv0, d1, d0, hvd1, hvd0;

      if (sum_v > sum_h)
      {
        hv1 = sum_v;
        hv0 = sum_h;
        dir_temp_hv = 1;
      }
      else
      {
        hv1 = sum_h;
        hv0 = sum_v;
        dir_temp_hv = 3;
      }
      if (sum_d0 > sum_d1)
      {
        d1 = sum_d0;
        d0 = sum_d1;
        dir_temp_d = 0;
      }
      else
      {
        d1 = sum_d1;
        d0 = sum_d0;
        dir_temp_d = 2;
      }
      if (d1*hv0 > hv1*d0)
      {
        hvd1 = d1;
        hvd0 = d0;
        main_direction = dir_temp_d;
        secondary_direction = dir_temp_hv;
      }
      else
      {
        hvd1 = hv1;
        hvd0 = hv0;
        main_direction = dir_temp_hv;
        secondary_direction = dir_temp_d;
      }

      int direction_strength = 0;
      if (hvd1 > 2 * hvd0)
      {
        direction_strength = 1;
      }
      if (hvd1 * 2 > 9 * hvd0)
      {
        direction_strength = 2;
      }

      if (direction_strength)
      {
        class_idx += (((main_direction & 0x1) << 1) + direction_strength) * 5;
      }

      static const int transpose_table[8] = { 0, 1, 0, 2, 2, 3, 1, 3 };
      int transpose_idx = transpose_table[main_direction * 2 + (secondary_direction >> 1)];

      int y_offset = i + pos_y;
      int x_offset = j + pos_x;

      alf_classifier *cl0 = classifier[y_offset] + x_offset;
      alf_classifier *cl1 = classifier[y_offset + 1] + x_offset;
      alf_classifier *cl2 = classifier[y_offset + 2] + x_offset;
      alf_classifier *cl3 = classifier[y_offset + 3] + x_offset;

      cl0[0].class_idx = cl0[1].class_idx = cl0[2].class_idx = cl0[3].class_idx = 
      cl1[0].class_idx = cl1[1].class_idx = cl1[2].class_idx = cl1[3].class_idx = 
      cl2[0].class_idx = cl2[1].class_idx = cl2[2].class_idx = cl2[3].class_idx = 
      cl3[0].class_idx = cl3[1].class_idx = cl3[2].class_idx = cl3[3].class_idx = class_idx;

      cl0[0].transpose_idx = cl0[1].transpose_idx = cl0[2].transpose_idx = cl0[3].transpose_idx =
      cl1[0].transpose_idx = cl1[1].transpose_idx = cl1[2].transpose_idx = cl1[3].transpose_idx =
      cl2[0].transpose_idx = cl2[1].transpose_idx = cl2[2].transpose_idx = cl2[3].transpose_idx =
      cl3[0].transpose_idx = cl3[1].transpose_idx = cl3[2].transpose_idx = cl3[3].transpose_idx = transpose_idx;

    }
  }
}

void kvz_alf_filter_block(encoder_state_t * const state,
  const lcu_order_element_t * const lcu,
  short* filter_set,
  clp_rng clp_rng,
  alf_component_id component_id,
  const int width,
  const int height,
  int x_pos,
  int y_pos)
{
  videoframe_t* const frame = state->tile->frame;
  alf_filter_type const filter_type = component_id == COMPONENT_Y ? ALF_FILTER_7X7 : ALF_FILTER_5X5;
  const bool chroma = component_id == COMPONENT_Y ? 0 : 1;
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;

  if (chroma)
  {
    assert((int)filter_type == 0); //Chroma needs to have filtType == 0
  }
  
  //bool isDualTree = CS::isDualITree(cs);
  bool is_dual_tree = false;

  //bool isPCMFilterDisabled = sps->getPCMFilterDisableFlag();
  bool is_pcm_filter_disabled = state->encoder_control->cfg.implicit_rdpcm;

  //ChromaFormat nChromaFormat = sps->getChromaFormatIdc();
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;

  //const int srcStride = srcLuma.stride;
  const int src_stride = frame->rec->stride;
  //const int dstStride = dstLuma.stride;
  const int dst_stride = frame->rec->stride;

  const int start_height = y_pos;
  const int end_height = start_height + height;
  const int start_width = x_pos;
  const int end_width = start_width + width;

  const kvz_pixel *src = state->tile->frame->rec->y;
  kvz_pixel *dst = state->tile->frame->rec->y + start_height * dst_stride;

  const kvz_pixel *p_img_y_pad_0, *p_img_y_pad_1, *p_img_y_pad_2, *p_img_y_pad_3, *p_img_y_pad_4, *p_img_y_pad_5, *p_img_y_pad_6;
  const kvz_pixel *p_img_0, *p_img_1, *p_img_2, *p_img_3, *p_img_4, *p_img_5, *p_img_6;

  short *coef = filter_set;

  const int shift = ALF_NUM_BITS - 1;

  const int offset = 1 << (shift - 1);

  int transpose_idx = 0;
  const int cls_size_y = 4;
  const int cls_size_x = 4;

  bool pcm_flags_2x2[4] = { 0,0,0,0 };

  assert((start_height % cls_size_y) == 0); //Wrong startHeight in filtering
  assert((start_width % cls_size_x) == 0); //Wrong startWidth in filtering
  assert(((end_height - start_height) % cls_size_y) == 0); //Wrong endHeight in filtering
  assert(((end_width - start_width) % cls_size_x) == 0); //Wrong endWidth in filtering

  alf_classifier *p_class = NULL;

  int dst_stride2 = dst_stride * cls_size_y;
  int src_stride2 = src_stride * cls_size_y;

  //std::vector<Pel> filterCoeff(MAX_NUM_ALF_LUMA_COEFF);
  kvz_pixel filter_coeff[MAX_NUM_ALF_LUMA_COEFF];

  p_img_y_pad_0 = src + start_height * src_stride + start_width;
  p_img_y_pad_1 = p_img_y_pad_0 + src_stride;
  p_img_y_pad_2 = p_img_y_pad_0 - src_stride;
  p_img_y_pad_3 = p_img_y_pad_1 + src_stride;
  p_img_y_pad_4 = p_img_y_pad_2 - src_stride;
  p_img_y_pad_5 = p_img_y_pad_3 + src_stride;
  p_img_y_pad_6 = p_img_y_pad_4 - src_stride;

  kvz_pixel* p_rec_0 = dst + start_width;
  kvz_pixel* p_rec_1 = p_rec_0 + dst_stride;

  for (int i = 0; i < end_height - start_height; i += cls_size_y)
  {
    if (!chroma)
    {
      p_class = classifier[start_height + i] + start_width;
    }

    for (int j = 0; j < end_width - start_width; j += cls_size_x)
    {
      if (!chroma)
      {
        alf_classifier cl = p_class[j];
        transpose_idx = cl.transpose_idx;
        if (is_pcm_filter_disabled && cl.class_idx == ALF_UNUSED_CLASS_IDX && transpose_idx == ALF_UNUSED_TRANSPOSE_IDX)
        {
          continue;
        }
        coef = filter_set + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
      }
      else if (is_pcm_filter_disabled)
      {
        int  blk_x, blk_y;
        bool *flags = pcm_flags_2x2;

        // check which chroma 2x2 blocks use PCM
        // chroma PCM may not be aligned with 4x4 ALF processing grid
        for (blk_y = 0; blk_y < 4; blk_y += 2)
        {
          for (blk_x = 0; blk_x < 4; blk_x += 2)
          {
            //Position pos(j + startWidth + blkX, i + startHeight + blkY);
            //CodingUnit* cu = is_dual_tree ? cs.getCU(pos, CH_C) : cs.getCU(recalcPosition(nChromaFormat, CH_C, CH_L, pos), CH_L);
            
            *flags++ = 1; //cu->ipcm ? 1 : 0;
          }
        }

        // skip entire 4x4 if all chroma 2x2 blocks use PCM
        if (pcm_flags_2x2[0] && pcm_flags_2x2[1] && pcm_flags_2x2[2] && pcm_flags_2x2[3])
        {
          continue;
        }
      }


      if (filter_type == ALF_FILTER_7X7)
      {
        if (transpose_idx == 1)
        {
          filter_coeff[0] = coef[9];
          filter_coeff[1] = coef[4];
          filter_coeff[2] = coef[10];
          filter_coeff[3] = coef[8];
          filter_coeff[4] = coef[1];
          filter_coeff[5] = coef[5];
          filter_coeff[6] = coef[11];
          filter_coeff[7] = coef[7];
          filter_coeff[8] = coef[3];
          filter_coeff[9] = coef[0];
          filter_coeff[10] = coef[2];
          filter_coeff[11] = coef[6];
          filter_coeff[12] = coef[12];
        }
        else if (transpose_idx == 2)
        {
          filter_coeff[0] = coef[0];
          filter_coeff[1] = coef[3];
          filter_coeff[2] = coef[2];
          filter_coeff[3] = coef[1];
          filter_coeff[4] = coef[8];
          filter_coeff[5] = coef[7];
          filter_coeff[6] = coef[6];
          filter_coeff[7] = coef[5];
          filter_coeff[8] = coef[4];
          filter_coeff[9] = coef[9];
          filter_coeff[10] = coef[10];
          filter_coeff[11] = coef[11];
          filter_coeff[12] = coef[12];
        }
        else if (transpose_idx == 3)
        {
          filter_coeff[0] = coef[9];
          filter_coeff[1] = coef[8];
          filter_coeff[2] = coef[10];
          filter_coeff[3] = coef[4];
          filter_coeff[4] = coef[3];
          filter_coeff[5] = coef[7];
          filter_coeff[6] = coef[11];
          filter_coeff[7] = coef[5];
          filter_coeff[8] = coef[1];
          filter_coeff[9] = coef[0];
          filter_coeff[10] = coef[2];
          filter_coeff[11] = coef[6];
          filter_coeff[12] = coef[12];
        }
        else
        {
          filter_coeff[0] = coef[0];
          filter_coeff[1] = coef[1];
          filter_coeff[2] = coef[2];
          filter_coeff[3] = coef[3];
          filter_coeff[4] = coef[4];
          filter_coeff[5] = coef[5];
          filter_coeff[6] = coef[6];
          filter_coeff[7] = coef[7];
          filter_coeff[8] = coef[8];
          filter_coeff[9] = coef[9];
          filter_coeff[10] = coef[10];
          filter_coeff[11] = coef[11];
          filter_coeff[12] = coef[12];
        }
      }
      else
      {
        if (transpose_idx == 1)
        {
          filter_coeff[0] = coef[4];
          filter_coeff[1] = coef[1];
          filter_coeff[2] = coef[5];
          filter_coeff[3] = coef[3];
          filter_coeff[4] = coef[0];
          filter_coeff[5] = coef[2];
          filter_coeff[6] = coef[6];
        }
        else if (transpose_idx == 2)
        {
          filter_coeff[0] = coef[0];
          filter_coeff[1] = coef[3];
          filter_coeff[2] = coef[2];
          filter_coeff[3] = coef[1];
          filter_coeff[4] = coef[4];
          filter_coeff[5] = coef[5];
          filter_coeff[6] = coef[6];
        }
        else if (transpose_idx == 3)
        {
          filter_coeff[0] = coef[4];
          filter_coeff[1] = coef[3];
          filter_coeff[2] = coef[5];
          filter_coeff[3] = coef[1];
          filter_coeff[4] = coef[0];
          filter_coeff[5] = coef[2];
          filter_coeff[6] = coef[6];
        }
        else
        {
          filter_coeff[0] = coef[0];
          filter_coeff[1] = coef[1];
          filter_coeff[2] = coef[2];
          filter_coeff[3] = coef[3];
          filter_coeff[4] = coef[4];
          filter_coeff[5] = coef[5];
          filter_coeff[6] = coef[6];
        }
      }

      for (int ii = 0; ii < cls_size_y; ii++)
      {
        p_img_0 = p_img_y_pad_0 + j + ii * src_stride;
        p_img_1 = p_img_y_pad_1 + j + ii * src_stride;
        p_img_2 = p_img_y_pad_2 + j + ii * src_stride;
        p_img_3 = p_img_y_pad_3 + j + ii * src_stride;
        p_img_4 = p_img_y_pad_4 + j + ii * src_stride;
        p_img_5 = p_img_y_pad_5 + j + ii * src_stride;
        p_img_6 = p_img_y_pad_6 + j + ii * src_stride;

        p_rec_1 = p_rec_0 + j + ii * dst_stride;

        for (int jj = 0; jj < cls_size_x; jj++)
        {

          // skip 2x2 PCM chroma blocks
          if (chroma && is_pcm_filter_disabled)
          {
            if (pcm_flags_2x2[2 * (ii >> 1) + (jj >> 1)])
            {
              p_img_0++;
              p_img_1++;
              p_img_2++;
              p_img_3++;
              p_img_4++;
              p_img_5++;
              p_img_6++;
              continue;
            }
          }

          int sum = 0;

          if (filter_type == ALF_FILTER_7X7)
          {
            sum += filter_coeff[0] * (p_img_5[0] + p_img_6[0]);
            sum += filter_coeff[1] * (p_img_3[+1] + p_img_4[-1]);
            sum += filter_coeff[2] * (p_img_3[+0] + p_img_4[+0]);
            sum += filter_coeff[3] * (p_img_3[-1] + p_img_4[+1]);
            sum += filter_coeff[4] * (p_img_1[+2] + p_img_2[-2]);
            sum += filter_coeff[5] * (p_img_1[+1] + p_img_2[-1]);
            sum += filter_coeff[6] * (p_img_1[+0] + p_img_2[+0]);
            sum += filter_coeff[7] * (p_img_1[-1] + p_img_2[+1]);
            sum += filter_coeff[8] * (p_img_1[-2] + p_img_2[+2]);
            sum += filter_coeff[9] * (p_img_0[+3] + p_img_0[-3]);
            sum += filter_coeff[10] * (p_img_0[+2] + p_img_0[-2]);
            sum += filter_coeff[11] * (p_img_0[+1] + p_img_0[-1]);
            sum += filter_coeff[12] * (p_img_0[+0]);
          }
          else
          {
            sum += filter_coeff[0] * (p_img_3[+0] + p_img_4[+0]);
            sum += filter_coeff[1] * (p_img_1[+1] + p_img_2[-1]);
            sum += filter_coeff[2] * (p_img_1[+0] + p_img_2[+0]);
            sum += filter_coeff[3] * (p_img_1[-1] + p_img_2[+1]);
            sum += filter_coeff[4] * (p_img_0[+2] + p_img_0[-2]);
            sum += filter_coeff[5] * (p_img_0[+1] + p_img_0[-1]);
            sum += filter_coeff[6] * (p_img_0[+0]);
          }

          sum = (sum + offset) >> shift;

          p_rec_1[jj] = alf_clip_pixel(sum, clp_rng);

          p_img_0++;
          p_img_1++;
          p_img_2++;
          p_img_3++;
          p_img_4++;
          p_img_5++;
          p_img_6++;
        }
      }
    }

    p_rec_0 += dst_stride2;
    p_rec_1 += dst_stride2;

    p_img_y_pad_0 += src_stride2;
    p_img_y_pad_1 += src_stride2;
    p_img_y_pad_2 += src_stride2;
    p_img_y_pad_3 += src_stride2;
    p_img_y_pad_4 += src_stride2;
    p_img_y_pad_5 += src_stride2;
    p_img_y_pad_6 += src_stride2;
  }
}