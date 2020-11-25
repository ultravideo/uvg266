

#include "alf.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cabac.h"
#include "rdo.h"
#include "strategies/strategies-sao.h"
#include "kvz_math.h"

#if MAX_NUM_CC_ALF_FILTERS>1
typedef struct filter_idx_count
{
  uint64_t count;
  uint8_t filter_idx;
} filter_idx_count;

int comparator(const void *v1, const void *v2)
{
  const filter_idx_count *p1 = (filter_idx_count *)v1;
  const filter_idx_count *p2 = (filter_idx_count *)v2;
  return (p1->count < p2->count);
}
#endif

//-------------------------help functions---------------------------

void set_aps_map(kvz_config *const cfg)
{
  if (g_frame_count == state->frame->num) {
    return;
  }
  g_frame_count = state->frame->num;

  reset_alf_param(&alf_param);


  //int layerIdx = cs.vps == nullptr ? 0 : cs.vps->getGeneralLayerIdx(cs.slice->getPic()->layerId);
  int layer_idx = state->frame->num;

  if (layer_idx && (false/*cs.slice->getPendingRasInit()*/ || (state->frame->pictype == KVZ_NAL_IDR_W_RADL || state->frame->pictype == KVZ_NAL_IDR_N_LP)))
  {
    for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++) {
      //state->slice->apss[i].aps_id = 0;
      //state->slice->apss[i].aps_type = 0;
      state->slice->apss[i].temporal_id = 0;
      state->slice->apss[i].layer_id = 0;
      reset_alf_param(&state->slice->apss[i]);
      state->slice->apss[i].num_luma_filters = 0;
      state->slice->apss[i].num_alternatives_chroma = 0;

    }
    g_aps_id_start = ALF_CTB_MAX_NUM_APS;
  }

  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  chroma_scale_x = (chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  chroma_scale_y = (chroma_fmt != KVZ_CSP_420) ? 0 : 1;

  //Default clp_rng for a slice
  g_clp_rngs.comp[COMPONENT_Y].min = g_clp_rngs.comp[COMPONENT_Cb].min = g_clp_rngs.comp[COMPONENT_Cr].min = 0;
  g_clp_rngs.comp[COMPONENT_Y].max = (1 << kvz_bit_depth) - 1;
  g_clp_rngs.comp[COMPONENT_Y].bd = kvz_bit_depth;
  g_clp_rngs.comp[COMPONENT_Y].n = 0;
  g_clp_rngs.comp[COMPONENT_Cb].max = g_clp_rngs.comp[COMPONENT_Cr].max = (1 << kvz_bit_depth) - 1;
  g_clp_rngs.comp[COMPONENT_Cb].bd = g_clp_rngs.comp[COMPONENT_Cr].bd = kvz_bit_depth;
  g_clp_rngs.comp[COMPONENT_Cb].n = g_clp_rngs.comp[COMPONENT_Cr].n = 0;
  g_clp_rngs.used = g_clp_rngs.chroma = false;

  //int shiftLuma = 2 * 0;// DISTORTION_PRECISION_ADJUSTMENT(g_input_bit_depth[CHANNEL_TYPE_LUMA]);
  //int shiftChroma = 2 * 0;// DISTORTION_PRECISION_ADJUSTMENT(m_inputBitDepth[CHANNEL_TYPE_CHROMA]);
  g_lambda[COMPONENT_Y] = state->frame->lambda;// *double(1 << shiftLuma);
  g_lambda[COMPONENT_Cb] = state->frame->lambda;// *double(1 << shiftChroma);
  g_lambda[COMPONENT_Cr] = state->frame->lambda;// *double(1 << shiftChroma);

                                                //g_alf_covariance_cc_alf[0] = 0;
                                                //g_alf_covariance_cc_alf[1] = 0;
                                                //g_alf_covariance_frame_cc_alf[0] = 0;
                                                //g_alf_covariance_frame_cc_alf[1] = 0;
}
#endif // !FULL_FRAME

void init_ctu_alternative_chroma(uint8_t* ctu_alts[MAX_NUM_COMPONENT], const int32_t num_ctus)
{
  uint8_t alt_idx = 0;
  for (int ctu_idx = 0; ctu_idx < num_ctus; ++ctu_idx)
  {
    ctu_alts[COMPONENT_Cb][ctu_idx] = alt_idx;
    if ((ctu_idx + 1) * g_alf_aps_temp.num_alternatives_chroma >= (alt_idx + 1) * num_ctus)
      ++alt_idx;
  }
  alt_idx = 0;
  for (int ctu_idx = 0; ctu_idx < num_ctus; ++ctu_idx)
  {
    ctu_alts[COMPONENT_Cr][ctu_idx] = alt_idx;
    if ((ctu_idx + 1) * g_alf_aps_temp.num_alternatives_chroma >= (alt_idx + 1) * num_ctus)
      ++alt_idx;
  }
}

int16_t clip_alf(const int16_t clip, const int16_t ref, const int16_t val0, const int16_t val1)
{
  return alf_clip3(-clip, +clip, val0 - ref) + alf_clip3(-clip, +clip, val1 - ref);
}

int alf_clip_pixel(const int a, const clp_rng clp_rng)
{
  return MIN(MAX(clp_rng.min, a), clp_rng.max);
}

int16_t alf_clip3(const int16_t minVal, const int16_t maxVal, const int16_t a)
{
  return MIN(MAX(minVal, a), maxVal);
}

void get_clip_max(const alf_covariance *cov, int *clip_max)
{
  const int num_coeff = cov->num_coeff;
  for (int k = 0; k < num_coeff - 1; ++k)
  {
    clip_max[k] = 0;

    bool inc = true;
    while (inc && clip_max[k] + 1 < cov->num_bins && cov->y[clip_max[k] + 1][k] == cov->y[clip_max[k]][k])
    {
      for (int l = 0; inc && l < num_coeff; ++l)
      {
        if (cov->ee[clip_max[k]][0][k][l] != cov->ee[clip_max[k] + 1][0][k][l])
        {
          inc = false;
        }
      }
      if (inc)
      {
        ++clip_max[k];
      }
    }
  }
  clip_max[num_coeff - 1] = 0;
}

void reduce_clip_cost(const alf_covariance *cov, int *clip)
{
  for (int k = 0; k < cov->num_coeff - 1; ++k)
  {
    bool dec = true;
    while (dec && clip[k] > 0 && cov->y[clip[k] - 1][k] == cov->y[clip[k]][k])
    {
      for (int l = 0; dec && l < cov->num_coeff; ++l)
      {
        if (cov->ee[clip[k]][clip[l]][k][l] != cov->ee[clip[k] - 1][clip[l]][k][l])
        {
          dec = false;
        }
      }
      if (dec)
      {
        --clip[k];
      }
    }
  }
}

void set_ey_from_clip(const alf_covariance *cov,const int* clip, double ee[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double y[MAX_NUM_ALF_LUMA_COEFF], int size)
{
  for (int k = 0; k<size; k++)
  {
    y[k] = cov->y[clip[k]][k];
    for (int l = 0; l<size; l++)
    {
      ee[k][l] = cov->ee[clip[k]][clip[l]][k][l];
    }
  }
}

double optimize_filter(const alf_covariance *cov, int* clip, double *f, bool optimize_clip)
{
  const int size = cov->num_coeff;
  int clip_max[MAX_NUM_ALF_LUMA_COEFF];

  double err_best, err_last;

  double ke[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];
  double ky[MAX_NUM_ALF_LUMA_COEFF];

  if (optimize_clip)
  {
    // Start by looking for min clipping that has no impact => max_clipping
    get_clip_max(cov, clip_max);
    for (int k = 0; k<size; ++k)
    {
      clip[k] = MAX(clip_max[k], clip[k]);
      clip[k] = MIN(clip[k], cov->num_bins - 1);
    }
  }

  set_ey_from_clip(cov, clip, ke, ky, size);

  gns_solve_by_chol(ke, ky, f, size);
  err_best = calculate_error(cov, clip, f);

  int step = optimize_clip ? (cov->num_bins + 1) / 2 : 0;

  while (step > 0)
  {
    double err_min = err_best;
    int idx_min = -1;
    int inc_min = 0;

    for (int k = 0; k < size - 1; ++k)
    {
      if (clip[k] - step >= clip_max[k])
      {
        clip[k] -= step;
        ky[k] = cov->y[clip[k]][k];
        for (int l = 0; l < size; l++)
        {
          ke[k][l] = cov->ee[clip[k]][clip[l]][k][l];
          ke[l][k] = cov->ee[clip[l]][clip[k]][l][k];
        }

        gns_solve_by_chol(ke, ky, f, size);
        err_last = calculate_error(cov, clip, f);

        if (err_last < err_min)
        {
          err_min = err_last;
          idx_min = k;
          inc_min = -step;
        }
        clip[k] += step;
      }
      if (clip[k] + step < cov->num_bins)
      {
        clip[k] += step;
        ky[k] = cov->y[clip[k]][k];
        for (int l = 0; l < size; l++)
        {
          ke[k][l] = cov->ee[clip[k]][clip[l]][k][l];
          ke[l][k] = cov->ee[clip[l]][clip[k]][l][k];
        }

        gns_solve_by_chol(ke, ky, f, size);
        err_last = calculate_error(cov, clip, f);

        if (err_last < err_min)
        {
          err_min = err_last;
          idx_min = k;
          inc_min = step;
        }
        clip[k] -= step;

      }
      ky[k] = cov->y[clip[k]][k];
      for (int l = 0; l < size; l++)
      {
        ke[k][l] = cov->ee[clip[k]][clip[l]][k][l];
        ke[l][k] = cov->ee[clip[l]][clip[k]][l][k];
      }
    }

    if (idx_min >= 0)
    {
      err_best = err_min;
      clip[idx_min] += inc_min;
      ky[idx_min] = cov->y[clip[idx_min]][idx_min];
      for (int l = 0; l < size; l++)
      {
        ke[idx_min][l] = cov->ee[clip[idx_min]][clip[l]][idx_min][l];
        ke[l][idx_min] = cov->ee[clip[l]][clip[idx_min]][l][idx_min];
      }
    }
    else
    {
      --step;
    }
  }

  if (optimize_clip) {
    // test all max
    for (int k = 0; k < size - 1; ++k)
    {
      clip_max[k] = 0;
    }
    double ke_max[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];
    double ky_max[MAX_NUM_ALF_LUMA_COEFF];
    set_ey_from_clip(cov, clip_max, ke_max, ky_max, size);

    gns_solve_by_chol(ke_max, ky_max, f, size);
    err_last = calculate_error(cov, clip_max, f);
    if (err_last < err_best)
    {
      err_best = err_last;
      for (int k = 0; k<size; ++k)
      {
        clip[k] = clip_max[k];
      }
    }
    else
    {
      // update clip to reduce coding cost
      reduce_clip_cost(cov, clip);

      // update f with best solution
      gns_solve_by_chol(ke, ky, f, size);
    }
  }

  return err_best;
}

double optimize_filter_clip(alf_covariance *cov, int* clip)
{
  double f[MAX_NUM_ALF_LUMA_COEFF];
  return optimize_filter(cov, clip, f, true);
}

double optimize_filter_gns_calc(alf_covariance *cov, const int* clip, double *f, int size)
{
  gns_solve_by_chol_clip_gns(cov, clip, f, size);
  return calculate_error(cov, clip, f);
}

void gns_backsubstitution(double r[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double* z, int size, double* a)
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

void gns_transpose_backsubstitution(double u[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double* rhs, double* x, int order)
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

int gns_cholesky_dec(double inp_matr[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double out_matr[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], int num_eq)
{
  double inv_diag[MAX_NUM_ALF_LUMA_COEFF];  /* Vector of the inverse of diagonal entries of outMatr */
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
        {
          inv_diag[i] = 1.0 / (out_matr[i][i] = sqrt(scale));
        }
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

int gns_solve_by_chol(double lhs[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double rhs[MAX_NUM_ALF_LUMA_COEFF], double *x, int num_eq)
{
  double aux[MAX_NUM_ALF_LUMA_COEFF];     /* Auxiliary vector */
  double u[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];    /* Upper triangular Cholesky factor of lhs */
  int res = 1;  // Signal that Cholesky factorization is successfully performed

                /* The equation to be solved is LHSx = rhs */

                /* Compute upper triangular U such that U'*U = lhs */
  if (gns_cholesky_dec(lhs, u, num_eq)) /* If Cholesky decomposition has been successful */
  {
    /* Now, the equation is  U'*U*x = rhs, where U is upper triangular
    * Solve U'*aux = rhs for aux
    */
    gns_transpose_backsubstitution(u, rhs, aux, num_eq);

    /* The equation is now U*x = aux, solve it for x (new motion coefficients) */
    gns_backsubstitution(u, aux, num_eq, x);

  }
  else /* lhs was singular */
  {
    res = 0;

    /* Regularize lhs */
    for (int i = 0; i < num_eq; i++)
    {
      lhs[i][i] += REG;
    }

    /* Compute upper triangular U such that U'*U = regularized lhs */
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

int gns_solve_by_chol_clip_gns(alf_covariance *cov, const int *clip, double *x, int num_eq)
{
  double lhs[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];
  double rhs[MAX_NUM_ALF_LUMA_COEFF];

  set_ey_from_clip(cov, clip, lhs, rhs, num_eq);
  return gns_solve_by_chol(lhs, rhs, x, num_eq);
}

double calc_error_for_coeffs(const alf_covariance *cov, const int *clip, const int *coeff, const int num_coeff, const int bit_depth)
{
  double factor = 1 << (bit_depth - 1);
  double error = 0;

  for (int i = 0; i < num_coeff; i++)   //diagonal
  {
    double sum = 0;
    for (int j = i + 1; j < num_coeff; j++)
    {
      sum += cov->ee[clip[i]][clip[j]][i][j] * coeff[j];
    }
    error += ((cov->ee[clip[i]][clip[i]][i][i] * coeff[i] + sum * 2) / factor - 2 * cov->y[clip[i]][i]) * coeff[i];
  }

  return error / factor;
}

double calc_error_for_cc_alf_coeffs(const alf_covariance *cov, const int16_t* coeff, const int num_coeff, const int bit_depth)
{
  double factor = 1 << (bit_depth - 1);
  double error = 0;

  for (int i = 0; i < num_coeff; i++)   // diagonal
  {
    double sum = 0;
    for (int j = i + 1; j < num_coeff; j++)
    {
      // E[j][i] = E[i][j], sum will be multiplied by 2 later
      sum += cov->ee[0][0][i][j] * coeff[j];
    }
    error += ((cov->ee[0][0][i][i] * coeff[i] + sum * 2) / factor - 2 *cov->y[0][i]) * coeff[i];
  }

  return error / factor;
}

int length_uvlc(int ui_code)
{
  int ui_length = 1;
  int ui_temp = ++ui_code;

  assert(ui_temp); // "Integer overflow"

  while (1 != ui_temp)
  {
    ui_temp >>= 1;
    ui_length += 2;
  }
  // Take care of cases where ui_length > 32
  return (ui_length >> 1) + ((ui_length + 1) >> 1);
}

double get_dist_coeff_force_0(bool* coded_var_bins, double error_force_0_coeff_tab[MAX_NUM_ALF_CLASSES][2], int* bits_var_bin, int zero_bits_var_bin, const int num_filters, double lambda)
{
  double dist_force_0 = 0;
  memset(coded_var_bins, 0, sizeof(*coded_var_bins) * MAX_NUM_ALF_CLASSES);

  for (int filt_idx = 0; filt_idx < num_filters; filt_idx++)
  {
    double cost_diff = error_force_0_coeff_tab[filt_idx][0] - (error_force_0_coeff_tab[filt_idx][1] + g_lambda[COMPONENT_Y] * bits_var_bin[filt_idx]);
    coded_var_bins[filt_idx] = cost_diff > 0 ? true : false;
    dist_force_0 += error_force_0_coeff_tab[filt_idx][coded_var_bins[filt_idx] ? 1 : 0];
  }

  return dist_force_0;
}

double get_dist_force_0(channel_type channel, const int num_filters, double error_tab_force_0_coeff[MAX_NUM_ALF_CLASSES][2], bool* coded_var_bins, double lambda)
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int bits_var_bin[MAX_NUM_ALF_CLASSES];

  for (int ind = 0; ind < num_filters; ++ind)
  {
    bits_var_bin[ind] = 0;
    for (int i = 0; i < num_coeff - 1; i++)
    {
      bits_var_bin[ind] += length_uvlc(abs(g_filter_coeff_set[ind][i]));
      if (abs(g_filter_coeff_set[ind][i]) != 0)
        bits_var_bin[ind] += 1;
    }
  }

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
/*#else
  if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
#endif*/
  {
    for (int ind = 0; ind < num_filters; ++ind)
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        if (!abs(g_filter_coeff_set[ind][i]))
        {
          g_filter_clipp_set[ind][i] = 0;
        }
      }
    }
  }

  double dist_force_0 = get_dist_coeff_force_0(coded_var_bins, error_tab_force_0_coeff, bits_var_bin, num_filters);

  return dist_force_0;
}

int get_cost_filter_coeff_force_0(channel_type channel, int **p_diff_q_filter_coeff_int_pp, const int num_filters, bool* coded_var_bins)
{
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int len = num_filters; //filter_coefficient_flag[i]

  // Filter coefficients
  for (int ind = 0; ind < num_filters; ++ind)
  {
    if (coded_var_bins[ind])
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        len += length_uvlc(abs(p_diff_q_filter_coeff_int_pp[ind][i])); // alf_coeff_luma_delta[i][j]
        if ((abs(p_diff_q_filter_coeff_int_pp[ind][i]) != 0))
          len += 1;
      }
    }
  }
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
/*#else
  if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
#endif*/
  {
    for (int ind = 0; ind < num_filters; ++ind)
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        if (!abs(p_diff_q_filter_coeff_int_pp[ind][i]))
        {
          g_filter_clipp_set[ind][i] = 0;
        }
        len += 2;
      }
    }
  }

  return len;
}

int get_cost_filter_coeff(channel_type channel, int **p_diff_q_filter_coeff_int_pp, const int num_filters)
{
  // #if JVET_O0216_ALF_COEFF_EG3 
  return length_filter_coeffs(channel, num_filters, p_diff_q_filter_coeff_int_pp);  // alf_coeff_luma_delta[i][j];
  /* #else
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

  //len += lengthFilterCoeffs( alfShape, num_filters, p_diff_q_filter_coeff_int_pp, m_kMinTab );  // alf_coeff_luma_delta[i][j]
  for (int ind = 0; ind < num_filters; ++ind)
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      len += length_golomb(abs(p_diff_q_filter_coeff_int_pp[ind][i]), g_k_min_tab[golomb_idx[i]]);
    }
  }

  return len;
  */
}

int get_cost_filter_clipp(channel_type channel, int **p_diff_q_filter_coeff_int_pp, const int num_filters)
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  for (int filter_idx = 0; filter_idx < num_filters; ++filter_idx)
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      if (!abs(p_diff_q_filter_coeff_int_pp[filter_idx][i]))
      {
        g_filter_clipp_set[filter_idx][i] = 0;
      }
    }
  }
  return (num_filters * (num_coeff - 1)) << 1;
}

/*#if !JVET_O0491_HLS_CLEANUP
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
}*/



int get_non_filter_coeff_rate(alf_aps *aps)
{
  //short* filter_coeff_delta_idx = aps->filter_coeff_delta_idx;
  //int fixed_filter_pattern = aps->fixed_filter_pattern;
  //int fixed_filter_set_index = aps->fixed_filter_set_index;

  int len = 1   // alf_coefficients_delta_flag
          /*#if !JVET_O0491_HLS_CLEANUP
          + length_truncated_unary(0, 3)    // chroma_idc = 0, it is signalled when ALF is enabled for luma
          + get_tb_length(num_luma_filters - 1, MAX_NUM_ALF_CLASSES);   //numLumaFilters*/
          + 2                                          // slice_alf_chroma_idc                     u(2)
          + length_uvlc(aps->num_luma_filters - 1);  // alf_luma_num_filters_signalled_minus1   ue(v)

  if (aps->num_luma_filters > 1)
  {
    const int coeff_length = kvz_math_ceil_log2(aps->num_luma_filters); //#if JVET_O0491_HLS_CLEANUP
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
    {
      //len += get_tb_length((int)filter_coeff_delta_idx[i], num_luma_filters);  //#if !JVET_O0491_HLS_CLEANUP
      len += coeff_length;
    }
  }
  /*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  len++; //fixed filter set flag
  if (*fixed_filter_set_index > 0)
  {
    len += get_tb_length(*fixed_filter_set_index - 1, ALF_NUM_FIXED_FILTER_SETS);
    len += 1; //fixed filter flag pattern
    if (*fixed_filter_pattern > 0)
      len += MAX_NUM_ALF_CLASSES;  //"fixed_filter_flag" for each class
  }*/

  return len;
}

int length_filter_coeffs(channel_type channel, const int num_filters, int **filter_coeff)
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int bit_cnt = 0;

  for (int ind = 0; ind < num_filters; ++ind)
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      bit_cnt += length_uvlc(abs(filter_coeff[ind][i]));
      if (abs(filter_coeff[ind][i]) != 0)
        bit_cnt += 1;
    }
  }
  return bit_cnt;
}

double calculate_error(const alf_covariance *cov, const int *clip, const double *coeff)
{
  double sum = 0;
  for (int i = 0; i < cov->num_coeff; i++)
  {
    sum += coeff[i] * cov->y[clip[i]][i];
  }

  return cov->pix_acc - sum;
}

double calculate_error_opt_filt(alf_covariance *cov, const int *clip)
{
  double c[MAX_NUM_ALF_LUMA_COEFF];
  return optimize_filter_gns_calc(cov, clip, c, cov->num_coeff);
}

int get_chroma_coeff_rate(alf_aps* aps, int alt_idx)
{
  int i_bits = 0;

  const int num_coeff = 7;
  // Filter coefficients
  for (int i = 0; i < num_coeff - 1; i++)
  {
    i_bits += length_uvlc(abs(aps->chroma_coeff[alt_idx][i]));  // alf_coeff_chroma[alt_idx][i]
    if ((aps->chroma_coeff[alt_idx][i]) != 0)
      i_bits += 1;
  }
  if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_CHROMA])
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      if (!abs(aps->chroma_coeff[alt_idx][i]))
      {
        aps->chroma_clipp[alt_idx][i] = 0;
      }
    }
    i_bits += ((num_coeff - 1) << 1);
  }
  return i_bits;
}

/*#if !JVET_O0491_HLS_CLEANUP
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
}*/

double get_filtered_distortion(alf_covariance* cov, const int num_classes, const int num_filters_minus1, const int num_coeff)
{
  double dist = 0;

  for (int class_idx = 0; class_idx < num_classes; class_idx++)
  {
    dist += calc_error_for_coeffs(&cov[class_idx], g_filter_clipp_set[class_idx], g_filter_coeff_set[class_idx], num_coeff, bit_depth);
  }

  return dist;
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
    /*#if !JVET_O0491_HLS_CLEANUP
    dist = get_unfiltered_distortion_cov_classes(cov, 1) + length_truncated_unary(0, 3) * g_lambda[COMPONENT_Cb];*/
    dist = get_unfiltered_distortion_cov_classes(cov, 1);
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

void get_frame_stats(channel_type channel, const int32_t num_ctus)
{
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? true : false;
  int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  int num_alternatives = is_luma ? 1 : g_alf_aps_temp.num_alternatives_chroma;
  // When calling this function m_ctuEnableFlag shall be set to 0 for CTUs using alternative APS
  // Here we compute frame stats for building new alternative filters
  for (int alt_idx = 0; alt_idx < num_alternatives; ++alt_idx)
  {
    for (int i = 0; i < num_classes; i++)
    {
      reset_alf_covariance(&g_alf_covariance_frame[channel][is_luma ? i : alt_idx], MAX_ALF_NUM_CLIPPING_VALUES);
    }
    if (is_luma)
    {
      get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_LUMA], g_alf_covariance[COMPONENT_Y], g_ctu_enable_flag[COMPONENT_Y], NULL, num_classes, alt_idx, num_ctus);
    }
    else
    {
      get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA], g_alf_covariance[COMPONENT_Cb], g_ctu_enable_flag[COMPONENT_Cb], g_ctu_alternative[COMPONENT_Cb], num_classes, alt_idx, num_ctus);
      get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA], g_alf_covariance[COMPONENT_Cr], g_ctu_enable_flag[COMPONENT_Cr], g_ctu_alternative[COMPONENT_Cr], num_classes, alt_idx, num_ctus);
    }
  }
/*#else
  for (int i = 0; i < num_classes; i++)
  {
    reset_alf_covariance(&g_alf_covariance_frame[channel][i_shape_idx][i], g_alf_num_clipping_values[channel]);
  }
  if (channel == CHANNEL_TYPE_LUMA)
  {
    get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_LUMA][i_shape_idx], g_alf_covariance[COMPONENT_Y][i_shape_idx], g_ctu_enable_flag[COMPONENT_Y], num_classes);
  }
  else
  {
    get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA][i_shape_idx], g_alf_covariance[COMPONENT_Cb][i_shape_idx], g_ctu_enable_flag[COMPONENT_Cb], num_classes);
    get_frame_stat(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA][i_shape_idx], g_alf_covariance[COMPONENT_Cr][i_shape_idx], g_ctu_enable_flag[COMPONENT_Cr], num_classes);
  }
#endif*/
}

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
void get_frame_stat(alf_covariance* frame_cov, alf_covariance** ctb_cov, uint8_t* ctb_enable_flags, uint8_t* ctb_alt_idx, const int num_classes, int alt_idx, int ctu_idx)
/*#else
void get_frame_stat(alf_covariance* frame_cov, alf_covariance** ctb_cov, uint8_t* ctb_enable_flags, const int num_classes)
#endif*/
{
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  const channel_type channel = (!ctb_alt_idx ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA);
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? true : false;
  for (int ctu_idx = 0; ctu_idx < num_ctus; ctu_idx++)
  {
    if (ctb_enable_flags[ctu_idx])
    {
      for (int class_idx = 0; class_idx < num_classes; class_idx++)
      {
        if (is_luma || alt_idx == ctb_alt_idx[ctu_idx])
        {
          add_alf_cov(&frame_cov[is_luma ? class_idx : alt_idx], &ctb_cov[ctu_idx][class_idx]);
        }
      }
    }
  }
/*#else
  for (int i = 0; i < g_num_ctus_in_pic; i++)
  {
    if (ctb_enable_flags[i])
    {
      for (int j = 0; j < num_classes; j++)
      {
        add_alf_cov(&frame_cov[j], &ctb_cov[i][j]);
      }
    }
  }
#endif*/
}

void copy_cov(alf_covariance *dst, alf_covariance *src)
{
  dst->num_coeff = src->num_coeff;
  dst->num_bins = src->num_bins;
  memcpy(&dst->ee, &src->ee, sizeof(dst->ee));
  memcpy(&dst->y, &src->y, sizeof(dst->y));
  dst->pix_acc = src->pix_acc;
}

void copy_alf_param(alf_aps *dst, alf_aps *src)
{
  memcpy(dst->enabled_flag, src->enabled_flag, sizeof(dst->enabled_flag));
  memcpy(dst->non_linear_flag, src->non_linear_flag, sizeof(dst->non_linear_flag));
  memcpy(dst->luma_coeff, src->luma_coeff, sizeof(dst->luma_coeff));
  memcpy(dst->luma_clipp, src->luma_clipp, sizeof(dst->luma_clipp));
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  dst->num_alternatives_chroma = src->num_alternatives_chroma;
//#endif
  memcpy(dst->chroma_coeff, src->chroma_coeff, sizeof(dst->chroma_coeff));
  memcpy(dst->chroma_clipp, src->chroma_clipp, sizeof(dst->chroma_clipp));
  memcpy(dst->filter_coeff_delta_idx, src->filter_coeff_delta_idx, sizeof(dst->filter_coeff_delta_idx));
  memcpy(dst->alf_luma_coeff_flag, src->alf_luma_coeff_flag, sizeof(dst->alf_luma_coeff_flag));
  dst->num_luma_filters = src->num_luma_filters;
  dst->alf_luma_coeff_delta_flag = src->alf_luma_coeff_delta_flag;
//#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  //dst->alf_luma_coeff_delta_prediction_flag = src->alf_luma_coeff_delta_prediction_flag;
//#endif
  dst->t_layer = src->t_layer;
  memcpy(dst->new_filter_flag, src->new_filter_flag, sizeof(dst->new_filter_flag));
//#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  //dst->fixed_filter_pattern = src->fixed_filter_pattern;
  //memcpy(dst->fixed_filter_idx, src->fixed_filter_idx, sizeof(dst->fixed_filter_idx));
  //dst->fixed_filter_set_index = src->fixed_filter_set_index;
//#endif
}

void copy_alf_param_w_channel(alf_aps* dst, alf_aps* src, channel_type channel)
{
  if (channel == CHANNEL_TYPE_LUMA)
  {
    copy_alf_param(dst, src);
  }
  else
  {
/*#if !JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    alfParamDst.nonLinearFlag[channel] = alfParamSrc.nonLinearFlag[channel];
#endif*/
    dst->enabled_flag[COMPONENT_Cb] = src->enabled_flag[COMPONENT_Cb];
    dst->enabled_flag[COMPONENT_Cr] = src->enabled_flag[COMPONENT_Cr];
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    dst->num_alternatives_chroma = src->num_alternatives_chroma;
    dst->non_linear_flag[CHANNEL_TYPE_CHROMA] = src->non_linear_flag[CHANNEL_TYPE_CHROMA];
//#endif
    memcpy(dst->chroma_coeff, src->chroma_coeff, sizeof(dst->chroma_coeff));
    memcpy(dst->chroma_clipp, src->chroma_clipp, sizeof(dst->chroma_clipp));
  }
}

void copy_aps(alf_aps *dst, alf_aps *src)
{
  dst->aps_id = src->aps_id;
  dst->temporal_id = src->temporal_id;
  dst->layer_id = src->layer_id;
  dst->aps_type = src->aps_type;
  copy_alf_param(dst, src);
}

void copy_aps_to_map(param_set_map *dst, alf_aps *src, int8_t aps_id)
{
  assert(0 <= aps_id && aps_id < ALF_CTB_MAX_NUM_APS);
  bool found = false;
  for (int id = 0; id < ALF_CTB_MAX_NUM_APS; id++) 
  {
    if (dst[aps_id + T_ALF_APS].parameter_set.aps_id == id) {
      found = true;
    }
  }
  if (!found)
  {
    dst[aps_id + T_ALF_APS].b_changed = true;
    //apsMap[apsId].p_nalu_data = 0;
    dst[aps_id + T_ALF_APS].parameter_set.aps_id = aps_id;
    dst[aps_id + T_ALF_APS].parameter_set.temporal_id = src->temporal_id;
    dst[aps_id + T_ALF_APS].parameter_set.layer_id = src->layer_id;
    dst[aps_id + T_ALF_APS].parameter_set.aps_type = src->aps_type;
    copy_alf_param(&dst[aps_id + T_ALF_APS].parameter_set, src);
  }
}

void reset_alf_param(alf_aps *src)
{
  memset(src->enabled_flag, false, sizeof(src->enabled_flag));
  memset(src->non_linear_flag, false, sizeof(src->non_linear_flag));
  memset(src->luma_coeff, 0, sizeof(src->luma_coeff));
  memset(src->luma_clipp, 0, sizeof(src->luma_clipp));
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  src->num_alternatives_chroma = 1;
//#endif
  memset(src->chroma_coeff, 0, sizeof(src->chroma_coeff));
  memset(src->chroma_clipp, 0, sizeof(src->chroma_clipp));
  memset(src->filter_coeff_delta_idx, 0, sizeof(src->filter_coeff_delta_idx));
  memset(src->alf_luma_coeff_flag, true, sizeof(src->alf_luma_coeff_flag));
  src->num_luma_filters = 1;
  src->alf_luma_coeff_delta_flag = false;
/*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  alfLumaCoeffDeltaPredictionFlag = false;
#endif*/
  src->t_layer = 0;
  memset(src->new_filter_flag, 0, sizeof(src->new_filter_flag));
/*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  fixedFilterPattern = 0;
  std::memset(fixedFilterIdx, 0, sizeof(fixedFilterIdx));
  fixedFilterSetIndex = 0;
#endif*/
}

void add_alf_cov(alf_covariance *dst, alf_covariance *src)
{
  int num_bins = dst->num_bins;
  int num_coeff = dst->num_coeff;
  for (int b0 = 0; b0 < num_bins; b0++)
  {
    for (int b1 = 0; b1 < num_bins; b1++)
    {
      for (int j = 0; j < num_coeff; j++)
      {
        for (int i = 0; i < num_coeff; i++)
        {
          dst->ee[b0][b1][j][i] += src->ee[b0][b1][j][i];
        }
      }
    }
  }
  for (int b = 0; b < num_bins; b++)
  {
    for (int j = 0; j < num_coeff; j++)
    {
      dst->y[b][j] += src->y[b][j];
    }
  }
  dst->pix_acc += src->pix_acc;
}

void add_alf_cov_lhs_rhs(alf_covariance *dst, alf_covariance *lhs, alf_covariance *rhs)
{
  int num_coeff = lhs->num_coeff;
  int num_bins = lhs->num_bins;
  for (int b0 = 0; b0 < num_bins; b0++)
  {
    for (int b1 = 0; b1 < num_bins; b1++)
    {
      for (int j = 0; j < num_coeff; j++)
      {
        for (int i = 0; i < num_coeff; i++)
        {
          dst->ee[b0][b1][j][i] = lhs->ee[b0][b1][j][i] + rhs->ee[b0][b1][j][i];
        }
      }
    }
  }
  for (int b = 0; b < num_bins; b++)
  {
    for (int j = 0; j < num_coeff; j++)
    {
      dst->y[b][j] = lhs->y[b][j] + rhs->y[b][j];
    }
  }
  dst->pix_acc = lhs->pix_acc + rhs->pix_acc;
}

void reset_alf_covariance(alf_covariance *alf, int num_bins) {
  if (num_bins > 0) { alf->num_bins = num_bins; }
  alf->pix_acc = 0;
  memset(alf->y, 0, sizeof(alf->y));
  memset(alf->ee, 0, sizeof(alf->ee));
}

void reset_cc_alf_aps_param(cc_alf_filter_param *cc_alf) {
  memset(cc_alf->cc_alf_filter_enabled, false, sizeof(cc_alf->cc_alf_filter_enabled));
  memset(cc_alf->cc_alf_filter_idx_enabled, false, sizeof(cc_alf->cc_alf_filter_idx_enabled));
  memset(cc_alf->cc_alf_coeff, 0, sizeof(cc_alf->cc_alf_coeff));
  cc_alf->cc_alf_filter_count[0] = cc_alf->cc_alf_filter_count[1] = MAX_NUM_CC_ALF_FILTERS;
  cc_alf->number_valid_components = 3;
  cc_alf->new_cc_alf_filter[0] = cc_alf->new_cc_alf_filter[1] = 0;
}

void copy_pixels(kvz_pixel *src, int x_src_start, int y_src_start, int src_stride, 
  kvz_pixel *dst, int x_dst_start, int y_dst_start, int dst_stride,
  int width, int height)
{
  for (int y = 0; y < height; y++)
  {
    int src_y = y_src_start + y;
    int dst_y = y_dst_start + y;
    for (int x = 0; x < width; x++)
    {
      int src_x = x_src_start + x;
      int dst_x = x_dst_start + x;
      dst[dst_y*dst_stride + dst_x] = src[src_y*src_stride + src_x];
    }
  }
}

void adjust_pixels(kvz_pixel *src, int x_start, int x_end, int y_start, int y_end, int stride, int pic_width, int pic_height)
{
  assert(x_start <= x_end);
  assert(y_start <= y_end);
  assert(x_end <= pic_width);
  assert(y_end <= pic_height);

  //not on any edge
  if (x_start != 0 && y_start != 0 && x_end != pic_width && y_end != pic_height) {
    return;
  }

  bool top_left = (x_start == 0 && y_start == 0);
  bool top_right = (x_end == pic_width && y_start == 0);
  bool bottom_left = (x_start == 0 && y_end == pic_height);
  bool bottom_right = (x_end == pic_width && y_end == pic_height);

  //left side
  if (x_start == 0) {
    for (int y = y_start; y < y_end; y++) {
      src[y * stride - 4] =
        src[y * stride - 3] =
        src[y * stride - 2] =
        src[y * stride - 1] = src[y * stride];
    }
  }
  //right side
  if (x_end == pic_width) {
    const int x_px = x_end - 1;
    for (int y = y_start; y < y_end; y++) {
      src[y * stride + x_px + 4] =
        src[y * stride + x_px + 3] =
        src[y * stride + x_px + 2] =
        src[y * stride + x_px + 1] = src[y * stride + x_px];
    }
  }
  //top
  if (y_start == 0) {
    for (int x = x_start; x < x_end; x++) {
      src[-4 * stride + x] =
        src[-3 * stride + x] =
        src[-2 * stride + x] =
        src[-1 * stride + x] = src[x];
    }
  }
  //bottom
  if (y_end == pic_height) {
    const int y_px = y_end - 1;
    for (int x = x_start; x < x_end; x++) {
      src[x + stride * (4 + y_px)] =
        src[x + stride * (3 + y_px)] =
        src[x + stride * (2 + y_px)] =
        src[x + stride * (1 + y_px)] = src[x + stride * y_px];
    }
  }
  //left top corner
  if (top_left) {
    for (int x = -4; x < 0; x++) {
      src[-4 * stride + x] =
        src[-3 * stride + x] =
        src[-2 * stride + x] =
        src[-1 * stride + x] = src[0];
    }
  }
  //right top corner
  if (top_right) {
    const int x_px = x_end - 1;
    for (int x = pic_width; x < pic_width + 4; x++) {
      src[-4 * stride + x] =
        src[-3 * stride + x] =
        src[-2 * stride + x] =
        src[-1 * stride + x] = src[x_px];
    }
  }

  //left or right bottom corner
  if (bottom_left) {
    const int y_px = y_end - 1;
    for (int x = -4; x < 0; x++) {
      src[(4 + y_px) * stride + x] =
        src[(3 + y_px) * stride + x] =
        src[(2 + y_px) * stride + x] =
        src[(1 + y_px) * stride + x] = src[stride * y_px];
    }
  }
  if (bottom_right) {
    const int x_px = x_end - 1;
    const int y_px = y_end - 1;
    for (int x = x_end; x < x_end + 4; x++) {
      src[(4 + y_px) * stride + x] =
        src[(3 + y_px) * stride + x] =
        src[(2 + y_px) * stride + x] =
        src[(1 + y_px) * stride + x] = src[stride * y_px + x_px];
    }
  }
}

void adjust_pixels_CTU_plus_4_pix(kvz_pixel *src, int x_start, int x_end, int y_start, int y_end, int stride, int pic_width, int pic_height)
{
  assert(x_start <= x_end);
  assert(y_start <= y_end);
  assert(x_end <= pic_width);
  assert(y_end <= pic_height);

  //not on any edge
  if (x_start != 0 && y_start != 0 && x_end != pic_width && y_end != pic_height) {
    return;
  }

  bool top_left = (x_start == 0 && y_start == 0);
  bool top_right = (x_end == pic_width && y_start == 0);
  bool bottom_left = (x_start == 0 && y_end == pic_height);
  bool bottom_right = (x_end == pic_width && y_end == pic_height);

  //left side
  if (top_left && !bottom_left) {
    for (int y = y_start; y < y_end + MAX_ALF_PADDING_SIZE; y++) {
      src[y * stride - 4] = 
      src[y * stride - 3] = 
      src[y * stride - 2] = 
      src[y * stride - 1] = src[y * stride];
    }
  }
  else if (!top_left && bottom_left) {
    for (int y = y_start + MAX_ALF_PADDING_SIZE; y < y_end; y++) {
      src[y * stride - 4] =
        src[y * stride - 3] =
        src[y * stride - 2] =
        src[y * stride - 1] = src[y * stride];
    }
  }
  else if (top_left && bottom_left) {
    for (int y = y_start; y < y_end; y++) {
      src[y * stride - 4] =
        src[y * stride - 3] =
        src[y * stride - 2] =
        src[y * stride - 1] = src[y * stride];
    }
  }
  else if (x_start == 0) {
    for (int y = y_start + MAX_ALF_PADDING_SIZE; y < y_end + MAX_ALF_PADDING_SIZE; y++) {
      src[y * stride - 4] =
        src[y * stride - 3] =
        src[y * stride - 2] =
        src[y * stride - 1] = src[y * stride];
    }
  }//left side 

  //right side
  if (top_right && !bottom_right) {
    const int x_px = x_end - 1;
    for (int y = y_start; y < y_end + MAX_ALF_PADDING_SIZE; y++) {
      src[y * stride + x_px + 4] =
      src[y * stride + x_px + 3] =
      src[y * stride + x_px + 2] =
      src[y * stride + x_px + 1] = src[y * stride + x_px];
    }
  }
  else if (!top_right && bottom_right) {
    const int x_px = x_end - 1;
    for (int y = y_start + MAX_ALF_PADDING_SIZE; y < y_end; y++) {
      src[y * stride + x_px + 4] =
        src[y * stride + x_px + 3] =
        src[y * stride + x_px + 2] =
        src[y * stride + x_px + 1] = src[y * stride + x_px];
    }
  }
  else if (top_right && bottom_right) {
    const int x_px = x_end - 1;
    for (int y = y_start; y < y_end; y++) {
      src[y * stride + x_px + 4] =
        src[y * stride + x_px + 3] =
        src[y * stride + x_px + 2] =
        src[y * stride + x_px + 1] = src[y * stride + x_px];
    }
  }
  else if (x_end == pic_width) {
    const int x_px = x_end - 1;
    for (int y = y_start + MAX_ALF_PADDING_SIZE; y < y_end + MAX_ALF_PADDING_SIZE; y++) {
      src[y * stride + x_px + 4] =
        src[y * stride + x_px + 3] =
        src[y * stride + x_px + 2] =
        src[y * stride + x_px + 1] = src[y * stride + x_px];
    }
  }//right side

  //top
  if (top_left && !top_right) {
    for (int x = x_start; x < x_end + MAX_ALF_PADDING_SIZE; x++) {
      src[-4 * stride + x] = 
      src[-3 * stride + x] = 
      src[-2 * stride + x] = 
      src[-1 * stride + x] = src[x];
    }
  }
  else if (!top_left && top_right) {
    for (int x = x_start + MAX_ALF_PADDING_SIZE; x < x_end; x++) {
      src[-4 * stride + x] =
        src[-3 * stride + x] =
        src[-2 * stride + x] =
        src[-1 * stride + x] = src[x];
    }
  }
  else if (top_left && top_right) {
    for (int x = x_start; x < x_end; x++) {
      src[-4 * stride + x] =
        src[-3 * stride + x] =
        src[-2 * stride + x] =
        src[-1 * stride + x] = src[x];
    }
  }
  else if (y_start == 0) {
    for (int x = x_start + MAX_ALF_PADDING_SIZE; x < x_end + MAX_ALF_PADDING_SIZE; x++) {
      src[-4 * stride + x] =
      src[-3 * stride + x] =
      src[-2 * stride + x] =
      src[-1 * stride + x] = src[x];
    }
  }//top

  //bottom
  if (bottom_left && !bottom_right) {
    const int y_px = y_end - 1;
    for (int x = x_start; x < x_end + MAX_ALF_PADDING_SIZE; x++) {
      src[x + stride * (4 + y_px)] =
      src[x + stride * (3 + y_px)] =
      src[x + stride * (2 + y_px)] =
      src[x + stride * (1 + y_px)] = src[x + stride * y_px];
    }
  }
  else if (!bottom_left && bottom_right) {
    const int y_px = y_end - 1;
    for (int x = x_start + MAX_ALF_PADDING_SIZE; x < x_end; x++) {
      src[x + stride * (4 + y_px)] =
        src[x + stride * (3 + y_px)] =
        src[x + stride * (2 + y_px)] =
        src[x + stride * (1 + y_px)] = src[x + stride * y_px];
    }
  }
  else if (bottom_left && bottom_right) {
    const int y_px = y_end - 1;
    for (int x = x_start; x < x_end; x++) {
      src[x + stride * (4 + y_px)] =
        src[x + stride * (3 + y_px)] =
        src[x + stride * (2 + y_px)] =
        src[x + stride * (1 + y_px)] = src[x + stride * y_px];
    }
  }
  else if (y_end == pic_height) {
    const int y_px = y_end - 1;
    for (int x = x_start + MAX_ALF_PADDING_SIZE; x < x_end + MAX_ALF_PADDING_SIZE; x++) {
      src[x + stride * (4 + y_px)] =
        src[x + stride * (3 + y_px)] =
        src[x + stride * (2 + y_px)] =
        src[x + stride * (1 + y_px)] = src[x + stride * y_px];
    }
  }//bottom


  //left top corner
  if (top_left) {
    for (int x = -4; x < 0; x++) {
      src[-4 * stride + x] = 
      src[-3 * stride + x] = 
      src[-2 * stride + x] = 
      src[-1 * stride + x] = src[0];
    }
  }

  //right top corner
  if (top_right) {
    const int x_px = x_end - 1;
    for (int x = pic_width; x < pic_width + 4; x++) {
      src[-4 * stride + x] = 
      src[-3 * stride + x] = 
      src[-2 * stride + x] = 
      src[-1 * stride + x] = src[x_px];
    }
  }

  //left bottom corner
  if (bottom_left) {
    const int y_px = y_end - 1;
    for (int x = -4; x < 0; x++) {
      src[(4 + y_px) * stride + x] =
      src[(3 + y_px) * stride + x] =
      src[(2 + y_px) * stride + x] =
      src[(1 + y_px) * stride + x] = src[stride * y_px];
    }
  }

  //right bottom corner
  if (bottom_right) {
    const int x_px = x_end - 1;
    const int y_px = y_end - 1;
    for (int x = x_end; x < x_end + 4; x++) {
      src[(4 + y_px) * stride + x] =
      src[(3 + y_px) * stride + x] =
      src[(2 + y_px) * stride + x] =
      src[(1 + y_px) * stride + x] = src[stride * y_px + x_px];
    }
  }
}

//Need to adjust
void adjust_pixels_chroma(kvz_pixel *src, int x_start, int x_end, int y_start, int y_end, int stride, int pic_width, int pic_height)
{
  assert(x_start <= x_end);
  assert(y_start <= y_end);
  assert(x_end <= pic_width);
  assert(y_end <= pic_height);

  //not on any edge
  if (x_start != 0 && y_start != 0 && x_end != pic_width && y_end != pic_height) {
    return;
  }

  bool top_left = (x_start == 0 && y_start == 0);
  bool top_right = (x_end == pic_width && y_start == 0);
  bool bottom_left = (x_start == 0 && y_end == pic_height);
  bool bottom_right = (x_end == pic_width && y_end == pic_height);

  //left side
  if (x_start == 0) {
    for (int y = y_start; y < y_end; y++) {
      src[y * stride - 2] =
      src[y * stride - 1] = src[y * stride];
    }
  }
  //right side
  if (x_end == pic_width) {
    const int x_px = x_end - 1;
    for (int y = y_start; y < y_end; y++) {
      src[y * stride + x_px + 2] =
      src[y * stride + x_px + 1] = src[y * stride + x_px];
    }
  }
  //top
  if (y_start == 0) {
    for (int x = x_start; x < x_end; x++) {
      src[-2 * stride + x] =
      src[-1 * stride + x] = src[x];
    }
  }
  //bottom
  if (y_end == pic_height) {
    const int y_px = y_end - 1;
    for (int x = x_start; x < x_end; x++) {
      src[x + stride * (2 + y_px)] =
      src[x + stride * (1 + y_px)] = src[x + stride * y_px];
    }
  }
  //left top corner
  if (top_left) {
    for (int x = -2; x < 0; x++) {
      src[-2 * stride + x] =
      src[-1 * stride + x] = src[0];
    }
  }
  //right top corner
  if (top_right) {
    const int x_px = x_end - 1;
    for (int x = pic_width; x < pic_width + 2; x++) {
      src[-2 * stride + x] =
      src[-1 * stride + x] = src[x_px];
    }
  }

  //left or right bottom corner
  if (bottom_left) {
    const int y_px = y_end - 1;
    for (int x = -2; x < 0; x++) {
      src[(2 + y_px) * stride + x] =
      src[(1 + y_px) * stride + x] = src[stride * y_px];
    }
  }
  if (bottom_right) {
    const int x_px = x_end - 1;
    const int y_px = y_end - 1;
    for (int x = x_end; x < x_end + 2; x++) {
      src[(2 + y_px) * stride + x] =
      src[(1 + y_px) * stride + x] = src[stride * y_px + x_px];
    }
  }
}

void set_ctu_enable_flag(uint8_t **flags, channel_type channel, uint8_t value, const int32_t num_ctus)
{
  if (channel == CHANNEL_TYPE_LUMA) {
    memset(flags[COMPONENT_Y], value, sizeof(uint8_t) * num_ctus);
  }
  else {
    memset(flags[COMPONENT_Cr], value, sizeof(uint8_t) * num_ctus);
    memset(flags[COMPONENT_Cb], value, sizeof(uint8_t) * num_ctus);
  }
}

void copy_ctu_enable_flag(uint8_t **flags_dst, uint8_t **flags_src, channel_type channel, const int32_t num_ctus)
{
  if (channel == CHANNEL_TYPE_LUMA) {
    memcpy(flags_dst[COMPONENT_Y], flags_src[COMPONENT_Y], sizeof(uint8_t) * num_ctus);
  }
  else {
    memcpy(flags_dst[COMPONENT_Cr], flags_src[COMPONENT_Cr], sizeof(uint8_t) * num_ctus);
    memcpy(flags_dst[COMPONENT_Cb], flags_src[COMPONENT_Cb], sizeof(uint8_t) * num_ctus);
  }
}

//-------------------------------------------------------------------

//-------------------------encoding functions------------------------

void kvz_alf_enc_process(encoder_state_t *const state)
{
  kvz_alf_enc_create(state);

  if (1 /*!layerIdx*/ && (false/*cs.slice->getPendingRasInit()*/ || (state->frame->pictype == KVZ_NAL_IDR_W_RADL || state->frame->pictype == KVZ_NAL_IDR_N_LP)))
  {
    for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++) {
      state->slice->apss[i].aps_id = -1;
      state->slice->apss[i].aps_type = 0;
      state->slice->apss[i].temporal_id = 0;
      state->slice->apss[i].layer_id = 0;
      reset_alf_param(&state->slice->apss[i]);
      state->slice->apss[i].num_luma_filters = 0;
      state->slice->apss[i].num_alternatives_chroma = 0;

      if (state->encoder_control->cfg.param_set_map[i + T_ALF_APS].b_changed)
      {
        alf_aps* alf_aps = &state->encoder_control->cfg.param_set_map[i + T_ALF_APS].parameter_set;
        cc_alf_filter_param* cc_alf_aps = &state->encoder_control->cfg.param_set_map[i + T_ALF_APS].parameter_set.cc_alf_aps_param;
        state->encoder_control->cfg.param_set_map[i + T_ALF_APS].b_changed = false;
        if (alf_aps)
        {
          alf_aps->aps_id = 0;
          alf_aps->aps_type = 0;
          alf_aps->temporal_id = 0;
          alf_aps->layer_id = 0;
          reset_alf_param(alf_aps);
          reset_cc_alf_aps_param(cc_alf_aps);
        }
      }
    }
    g_aps_id_start = ALF_CTB_MAX_NUM_APS;
  }

  alf_aps alf_param;
  reset_alf_param(&alf_param);

  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  bool chroma_scale_x = (chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool chroma_scale_y = (chroma_fmt != KVZ_CSP_420) ? 0 : 1;
  int8_t kvz_bit_depth = state->encoder_control->bitdepth;
  const int32_t num_ctus_in_pic = state->lcu_order_count;

  //Default clp_rng
  g_clp_rngs.comp[COMPONENT_Y].min = g_clp_rngs.comp[COMPONENT_Cb].min = g_clp_rngs.comp[COMPONENT_Cr].min = 0;
  g_clp_rngs.comp[COMPONENT_Y].max = (1 << kvz_bit_depth) - 1;
  g_clp_rngs.comp[COMPONENT_Y].bd = kvz_bit_depth;
  g_clp_rngs.comp[COMPONENT_Y].n = 0;
  g_clp_rngs.comp[COMPONENT_Cb].max = g_clp_rngs.comp[COMPONENT_Cr].max = (1 << kvz_bit_depth) - 1;
  g_clp_rngs.comp[COMPONENT_Cb].bd = g_clp_rngs.comp[COMPONENT_Cr].bd = kvz_bit_depth;
  g_clp_rngs.comp[COMPONENT_Cb].n = g_clp_rngs.comp[COMPONENT_Cr].n = 0;
  g_clp_rngs.used = g_clp_rngs.chroma = false;
  double lambda_chroma_weight = 0.0;

  cabac_data_t ctx_start;
  cabac_data_t ctx_start_cc_alf;
  memcpy(&cabac_estimator, &state->cabac, sizeof(cabac_estimator));
  memcpy(&ctx_start, &state->cabac, sizeof(ctx_start));
  memcpy(&ctx_start_cc_alf, &cabac_estimator, sizeof(ctx_start_cc_alf));
  cabac_estimator.only_count = 1;
  ctx_start.only_count = 1;
  ctx_start_cc_alf.only_count = 1;

  // derive classification
  const int luma_height = state->tile->frame->height;
  const int luma_width = state->tile->frame->width;

  for (int y_pos = 0; y_pos < luma_height; y_pos += LCU_WIDTH)
  {
    for (int x_pos = 0; x_pos < lumaWidth; x_pos += maxCUWidth)
    {*/
  const int y_pos = lcu->position_px.y;
  const int x_pos = lcu->position_px.x;
  const int width = lcu->size.x; //(x_pos + maxCUWidth > lumaWidth) ? (lumaWidth - x_pos) : maxCUWidth;
  const int height = lcu->size.y; //(y_pos + maxCUHeight > lumaHeight) ? (lumaHeight - y_pos) : maxCUHeight;
  int raster_slice_alf_pad = 0;
  //T‰t‰ if-lauseen sis‰ll‰ olevaa algoritmia pit‰‰ viel‰ viilata
  if (is_crossed_by_virtual_boundaries(x_pos, y_pos, width, height, &clip_top, &clip_bottom, &clip_left, &clip_right, &num_hor_vir_bndry, &num_ver_vir_bndry, hor_vir_bndry_pos, ver_vir_bndry_pos, state))
  {
    int y_start = y_pos;
    for (int i = 0; i <= num_hor_vir_bndry; i++)
    {
      const int width = (x_pos + LCU_WIDTH > luma_width) ? (luma_width - x_pos) : LCU_WIDTH;
      const int height = (y_pos + LCU_WIDTH > luma_height) ? (luma_height - y_pos) : LCU_WIDTH;
      {
        kvz_alf_derive_classification(state, width, height, x_pos, y_pos, x_pos, y_pos);
      }
    }
  } 

  // get CTB stats for filtering 
  kvz_alf_derive_stats_for_filtering(state);

  //#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  /*for (int ctbIdx = 0; ctbIdx < m_numCTUsInPic; ctbIdx++)
  {
    g_alf_ctb_filter_index[ctb_iIdx] = ALF_NUM_FIXED_FILTER_SETS;
  }

  // consider using new filter (only)
  alf_param.new_filter_flag[CHANNEL_TYPE_LUMA] = true;
  alf_param.new_filter_flag[CHANNEL_TYPE_CHROMA] = true;
  state->slice->tile_group_num_aps = 1; // Only new filter for RD cost optimization

  // derive filter (luma)
  kvz_alf_encoder(state,
    &alf_param, CHANNEL_TYPE_LUMA,
    lambda_chroma_weight
  );

  // derive filter (chroma)
  if (state->encoder_control->chroma_format != KVZ_CSP_400) {
    kvz_alf_encoder(state,
      &alf_param, CHANNEL_TYPE_CHROMA,
      lambda_chroma_weight
    );
  }
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  // let alfEncoderCtb decide now
  alf_param.new_filter_flag[CHANNEL_TYPE_LUMA] = false;
  alf_param.new_filter_flag[CHANNEL_TYPE_CHROMA] = false;
  state->slice->tile_group_num_aps = 0;
//#endif

  //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
  memcpy(&cabac_estimator, &ctx_start, sizeof(cabac_estimator));
  kvz_alf_encoder_ctb(state, &alf_param, lambda_chroma_weight);

  //for (int s = 0; s < state.; s++) //numSliceSegments
  {
    if (state->encoder_control->cfg.lossless)
    {
      for (uint32_t ctb_idx = 0; ctb_idx < num_ctus_in_pic; ctb_idx++) //pcPic->slices[s]->getNumCtuInSlice()
      {
        //uint32_t ctuRsAddr = pcPic->slices[s]->getCtuAddrInSlice(ctuIdx);
        g_ctu_enable_flag[COMPONENT_Y][ctb_idx] = 0;
        g_ctu_enable_flag[COMPONENT_Cb][ctb_idx] = 0;
        g_ctu_enable_flag[COMPONENT_Cr][ctb_idx] = 0;
      }
    }
  }

  kvz_alf_reconstruct(state);

  // Do not transmit CC ALF if it is unchanged
  if(state->slice->tile_group_alf_enabled_flag[COMPONENT_Y])
  {
    for (int32_t luma_alf_aps_id = 0; luma_alf_aps_id < state->slice->tile_group_num_aps; luma_alf_aps_id++ )
    {
      //APS* aps = (luma_alf_aps_id >= 0) ? m_apsMap->getPS((luma_alf_aps_id << NUM_APS_TYPE_LEN) + ALF_APS) : nullptr;
      int aps_id = state->slice->tile_group_luma_aps_id[luma_alf_aps_id];
      alf_aps* aps = (aps_id >= 0) ? &state->encoder_control->cfg.param_set_map[aps_id + T_ALF_APS].parameter_set : NULL;
      bool changed = state->encoder_control->cfg.param_set_map[aps_id + T_ALF_APS].b_changed;
      if (aps && changed)
      {
        aps->cc_alf_aps_param.new_cc_alf_filter[0] = false;
        aps->cc_alf_aps_param.new_cc_alf_filter[1] = false;
      }
    }
  }
  int chroma_alf_aps_id = (state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] || state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr]) ? state->slice->tile_group_chroma_aps_id : -1;
  alf_aps* aps = (chroma_alf_aps_id >= 0) ? &state->encoder_control->cfg.param_set_map[chroma_alf_aps_id + T_ALF_APS].parameter_set : NULL;
  bool changed = (chroma_alf_aps_id >= 0) ? state->encoder_control->cfg.param_set_map[chroma_alf_aps_id + T_ALF_APS].b_changed : 0;
  if (aps && changed)
  {
    aps->cc_alf_aps_param.new_cc_alf_filter[0] = false;
    aps->cc_alf_aps_param.new_cc_alf_filter[1] = false;
  }
  
  if (state->encoder_control->cfg.alf_type != 2)
  {
    return;
  }

  const kvz_picture *org_yuv = state->tile->frame->source;
  const kvz_picture *rec_yuv = state->tile->frame->rec;

  const int luma_stride = state->tile->frame->rec->stride;
  const int chroma_stride = luma_stride >> chroma_scale_x;
  const int chroma_height = luma_height >> chroma_scale_y;
  const int chroma_padding = MAX_ALF_PADDING_SIZE >> chroma_scale_x;

  const int index_chroma = -(chroma_stride * chroma_padding + chroma_padding);

  //Copy reconstructed samples to a buffer.
  memcpy(&alf_tmp_u[index_chroma], &state->tile->frame->rec->u[index_chroma],
    sizeof(kvz_pixel) * chroma_stride * (chroma_height + chroma_padding * 2));
  memcpy(&alf_tmp_v[index_chroma], &state->tile->frame->rec->v[index_chroma],
    sizeof(kvz_pixel) * chroma_stride * (chroma_height + chroma_padding * 2));

  adjust_pixels_chroma(alf_tmp_u,
    0,
    rec_yuv->width >> chroma_scale_x,
    0,
    rec_yuv->height >> chroma_scale_y,
    rec_yuv->stride >> chroma_scale_x,
    rec_yuv->width >> chroma_scale_x,
    rec_yuv->height >> chroma_scale_y);
  adjust_pixels_chroma(alf_tmp_v,
    0,
    rec_yuv->width >> chroma_scale_x,
    0,
    rec_yuv->height >> chroma_scale_y,
    rec_yuv->stride >> chroma_scale_x,
    rec_yuv->width >> chroma_scale_x,
    rec_yuv->height >> chroma_scale_y);

  const int num_ctus_in_width = state->tile->frame->width_in_lcu;
  derive_stats_for_cc_alf_filtering(state, org_yuv, COMPONENT_Cb, num_ctus_in_width, (0 + 1));
  derive_stats_for_cc_alf_filtering(state, org_yuv, COMPONENT_Cr, num_ctus_in_width, (0 + 1));
  init_distortion_cc_alf(num_ctus_in_pic);

  memcpy(&cabac_estimator, &ctx_start_cc_alf, sizeof(cabac_estimator));
  derive_cc_alf_filter(state, COMPONENT_Cb, org_yuv, rec_yuv);
  memcpy(&cabac_estimator, &ctx_start_cc_alf, sizeof(cabac_estimator));
  derive_cc_alf_filter(state, COMPONENT_Cr, org_yuv, rec_yuv);

  setup_cc_alf_aps(state);

  for (int comp_idx = 1; comp_idx < (state->encoder_control->chroma_format == KVZ_CSP_400 ? 1 : MAX_NUM_COMPONENT); comp_idx++)
  {
    if (g_cc_alf_filter_param.cc_alf_filter_enabled[comp_idx - 1])
    {
      const kvz_pixel* rec_uv = comp_idx == COMPONENT_Cb ? rec_yuv->u : rec_yuv->v;
      const int luma_stride = rec_yuv->stride;
      apply_cc_alf_filter(state, comp_idx, rec_uv, alf_tmp_y, luma_stride, g_cc_alf_filter_control[comp_idx - 1],
        g_cc_alf_filter_param.cc_alf_coeff[comp_idx - 1], -1);
    }
  }
}

double kvz_alf_derive_ctb_alf_enable_flags(encoder_state_t * const state,
  channel_type channel,
  double *dist_unfilter,
  const int num_classes,
  const double chroma_weight
  )
{
/*  TempCtx        ctxTempStart(m_CtxCache);
  TempCtx        ctxTempBest(m_CtxCache);*/
  cabac_data_t ctx_temp_start;
  cabac_data_t ctx_temp_best;

/*#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  TempCtx        ctxTempAltStart(m_CtxCache);
  TempCtx        ctxTempAltBest(m_CtxCache);*/
  cabac_data_t ctx_temp_alt_start;
  //cabac_data_t ctx_temp_alt_best;
//#endif

  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;

  const kvz_pixel comp_id_first = is_luma ? COMPONENT_Y : COMPONENT_Cb;
  const kvz_pixel comp_id_last = is_luma ? COMPONENT_Y : COMPONENT_Cr;
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  const int num_alts = is_luma ? 1 : g_alf_aps_temp.num_alternatives_chroma;
//#endif

  int num_coeff = is_luma ? 13 : 7;

  double cost = 0;
  double lambda = state->frame->lambda;
  *dist_unfilter = 0;

  if (is_luma) {
    g_alf_aps_temp.enabled_flag[COMPONENT_Y] = 1;
  }
  else {
    g_alf_aps_temp.enabled_flag[COMPONENT_Cb] = 1;
    g_alf_aps_temp.enabled_flag[COMPONENT_Cr] = 1;
  }

  assert((chroma_weight <= 0.0) && (state->slice->start_in_rs == 0)); //"incompatible start CTU address, must be 0"

  kvz_alf_reconstruct_coeff(state, &g_alf_aps_temp, channel, true, is_luma);

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  for (int alt_idx = 0; alt_idx < (is_luma ? 1 : MAX_NUM_ALF_ALTERNATIVES_CHROMA); alt_idx++)
  {
    for (int class_idx = 0; class_idx < (is_luma ? MAX_NUM_ALF_CLASSES : 1); class_idx++)
    {
      for (int i = 0; i < (is_luma ? MAX_NUM_ALF_LUMA_COEFF : MAX_NUM_ALF_CHROMA_COEFF); i++)
      {
        g_filter_coeff_set[is_luma ? class_idx : alt_idx][i] = is_luma ? g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + i] : g_chroma_coeff_final[alt_idx][i];
        g_filter_clipp_set[is_luma ? class_idx : alt_idx][i] = is_luma ? g_clipp_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + i] : g_chroma_clipp_final[alt_idx][i];
      }
    }
  }
/*#else
  for (int class_idx = 0; class_idx < (is_luma ? MAX_NUM_ALF_CLASSES : 1); class_idx++)
  {
    for (int i = 0; i < (is_luma ? MAX_NUM_ALF_LUMA_COEFF : MAX_NUM_ALF_CHROMA_COEFF); i++)
    {
      g_filter_coeff_set[class_idx][i] = is_luma ? g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + i] : g_chroma_coeff_final[i];
      g_filter_clipp_set[class_idx][i] = is_luma ? g_clipp_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + i] : g_chroma_clipp_final[i];
    }
  }
#endif*/

  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    for (int comp_id = comp_id_first; comp_id <= comp_id_last; comp_id++)
    {
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
//#if ENABLE_QPA
      const double ctu_lambda = chroma_weight > 0.0 ? (is_luma ? 0/*cs.picture->m_uEnerHpCtu[ctuIdx]*/ : 0/*cs.picture->m_uEnerHpCtu[ctuIdx]*/ / chroma_weight) : g_lambda[comp_id];
/*#else
      const double ctu_lambda = m_lambda[comp_id];
#endif
#endif*/

      double dist_unfilter_ctu = get_unfiltered_distortion_cov_classes(g_alf_covariance[comp_id][ctu_idx], num_classes);

      //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
      memcpy(&ctx_temp_start, &cabac_estimator, sizeof(ctx_temp_start));
      //m_CABACEstimator->resetBits();
      kvz_cabac_reset_bits(&cabac_estimator);
      cabac_estimator.only_count = 1;
      g_ctu_enable_flag[comp_id][ctu_idx] = 1;
      code_alf_ctu_enable_flag(state, &cabac_estimator, ctu_idx, comp_id, &g_alf_aps_temp);

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
        if (is_luma)
        {
          // Evaluate cost of signaling filter set index for convergence of filters enabled flag / filter derivation
          assert(g_alf_ctb_filter_index[ctu_idx] == ALF_NUM_FIXED_FILTER_SETS);
          assert(state->slice->tile_group_num_aps == 1);
          code_alf_ctu_filter_index(state, &cabac_estimator, ctu_idx, g_alf_aps_temp.enabled_flag[COMPONENT_Y]);
        }
      double cost_on = dist_unfilter_ctu + ctu_lambda * (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3);
/*#else
      double costOn = distUnfilterCtu + getFilteredDistortion(m_alfCovariance[compID][iShapeIdx][ctuIdx], numClasses, m_alfParamTemp.numLumaFilters - 1, numCoeff);
#if ENABLE_QPA
      const double ctu_lambda = chromaWeight > 0.0 ? (isLuma(channel) ? cs.picture->m_uEnerHpCtu[ctuIdx] : cs.picture->m_uEnerHpCtu[ctuIdx] / chromaWeight) : m_lambda[compID];
#else
      const double ctu_lambda = m_lambda[compID];
#endif
      costOn += ctu_lambda * FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
#endif*/

      //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
      memcpy(&ctx_temp_best, &cabac_estimator, sizeof(ctx_temp_best));

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
      if (is_luma)
      {
        cost_on += get_filtered_distortion(g_alf_covariance[comp_id][ctu_idx], num_classes, g_alf_aps_temp.num_luma_filters - 1, num_coeff, bit_depth);
      }
      else
      {
        double best_alt_cost = MAX_DOUBLE;
        int best_alt_idx = -1;
        //ctxTempAltStart = AlfCtx(ctxTempBest);
        memcpy(&ctx_temp_alt_start, &ctx_temp_best, sizeof(ctx_temp_alt_start));
        for (int alt_idx = 0; alt_idx < num_alts; ++alt_idx)
        {
          if (alt_idx) 
          {
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempAltStart);
            memcpy(&cabac_estimator, &ctx_temp_alt_start, sizeof(cabac_estimator));
          }
          //m_CABACEstimator->resetBits();
          kvz_cabac_reset_bits(&cabac_estimator);
          cabac_estimator.only_count = 1;
          g_ctu_alternative[comp_id][ctu_idx] = alt_idx;
          code_alf_ctu_alternative_ctu(state, &cabac_estimator, ctu_idx, comp_id, &g_alf_aps_temp);
          double r_altCost = ctu_lambda * (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale * 0/*m_CABACEstimator->getEstFracBits()*/;

          double alt_dist = 0.;
          alt_dist += calc_error_for_coeffs(&g_alf_covariance[comp_id][ctu_idx][0], g_filter_clipp_set[alt_idx], g_filter_coeff_set[alt_idx], num_coeff, bit_depth);

          double alt_cost = alt_dist + r_altCost;
          if (alt_cost < best_alt_cost)
          {
            best_alt_cost = alt_cost;
            best_alt_idx = alt_idx;
            //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
            memcpy(&ctx_temp_best, &cabac_estimator, sizeof(ctx_temp_best));
          }
        }
        g_ctu_alternative[comp_id][ctu_idx] = best_alt_idx;
        cost_on += best_alt_cost;
      }
//#endif

      //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
      memcpy(&cabac_estimator, &ctx_temp_start, sizeof(cabac_estimator));
      //m_CABACEstimator->resetBits();
      kvz_cabac_reset_bits(&cabac_estimator);
      cabac_estimator.only_count = 1;
      g_ctu_enable_flag[comp_id][ctu_idx] = 0;
      code_alf_ctu_enable_flag(state, &cabac_estimator, ctu_idx, comp_id, &g_alf_aps_temp);
      double cost_off = dist_unfilter_ctu + ctu_lambda * (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale * 0;// m_CABACEstimator->getEstFracBits();

      if (cost_on < cost_off)
      {
        cost += cost_on;
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
        memcpy(&cabac_estimator, &ctx_temp_best, sizeof(cabac_estimator));
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
    const alf_component_id compIDFirst = COMPONENT_Cb;
    const alf_component_id compIDLast = COMPONENT_Cr;
    for (int compId = compIDFirst; compId <= compIDLast; compId++)
    {
      g_alf_aps_temp.enabled_flag[compId] = false;
      for (int i = 0; i < num_ctus_in_pic; i++)
      {
        if (g_ctu_enable_flag[compId][i])
        {
          g_alf_aps_temp.enabled_flag[compId] = true;
          break;
        }
      }
    }
    /*#if !JVET_O0491_HLS_CLEANUP
    const int alf_chroma_idc = g_alf_aps_temp.enabled_flag[COMPONENT_Cb] * 2 + g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
    cost += length_truncated_unary(alf_chroma_idc, 3) * g_lambda[channel];*/
  }

  return cost;
}


void kvz_alf_enc_create(encoder_state_t * const state)
{
  if (g_created) {
    return;
  }

  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  const int pic_width = state->tile->frame->width;
  const int pic_height = state->tile->frame->height;
  const int8_t input_bitdepth = state->encoder_control->bitdepth;
  const int32_t num_ctus_in_pic = state->lcu_order_count;

  assert(MAX_ALF_NUM_CLIPPING_VALUES > 0); //"g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] must be at least one"
  g_alf_clipping_values[CHANNEL_TYPE_LUMA][0] = 1 << input_bitdepth;
  int shift_luma = input_bitdepth - 8;
  for (int i = 1; i < MAX_ALF_NUM_CLIPPING_VALUES; ++i)
  {
    g_alf_clipping_values[CHANNEL_TYPE_LUMA][i] = 1 << (7 - 2 * i + shift_luma);
  }

  assert(MAX_ALF_NUM_CLIPPING_VALUES > 0); //"g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] must be at least one"
  g_alf_clipping_values[CHANNEL_TYPE_CHROMA][0] = 1 << input_bitdepth;
  int shift_chroma = input_bitdepth - 8;
  for (int i = 1; i < MAX_ALF_NUM_CLIPPING_VALUES; ++i)
  {
    g_alf_clipping_values[CHANNEL_TYPE_CHROMA][i] = 1 << (7 - 2 * i + shift_chroma);
  }

  // Classification
  g_classifier = malloc(pic_height * sizeof(**g_classifier));
  g_classifier[0] = malloc(pic_height * pic_width * sizeof(*g_classifier));

  for (int i = 1; i < pic_height; i++)
  {
    g_classifier[i] = g_classifier[0] + i * pic_width;
  }

  for (int filter_set_index = 0; filter_set_index < ALF_NUM_FIXED_FILTER_SETS; filter_set_index++)
  {
    for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
    {
      int fixed_filter_idx = g_class_to_filter_mapping[filter_set_index][class_idx];
      for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF - 1; i++)
      {
        g_fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + i] = g_fixed_filter_set_coeff[fixed_filter_idx][i];
      }
      g_fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + MAX_NUM_ALF_LUMA_COEFF - 1] = (1 << (input_bitdepth - 1));
    }
  }

  for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF * MAX_NUM_ALF_CLASSES; i++)
  {
    g_clip_default[i] = g_alf_clipping_values[CHANNEL_TYPE_LUMA][0];
  }

  for (int j = 0; j <= MAX_NUM_ALF_CLASSES + 1; j++)
  {
    g_alf_covariance_merged[j].num_coeff = 13;
    g_alf_covariance_merged[j].num_bins = MAX_ALF_NUM_CLIPPING_VALUES;
    memset(g_alf_covariance_merged[j].y, 0, sizeof(g_alf_covariance_merged[j].y));
    memset(g_alf_covariance_merged[j].ee, 0, sizeof(g_alf_covariance_merged[j].ee));
  }

  g_cc_alf_filter_control[0] = malloc(num_ctus_in_pic * sizeof(*g_cc_alf_filter_control));
  g_cc_alf_filter_control[1] = malloc(num_ctus_in_pic * sizeof(*g_cc_alf_filter_control));

  for (int channel_idx = 0; channel_idx < MAX_NUM_CHANNEL_TYPE; channel_idx++)
  {
    channel_type ch_type = (channel_type)channel_idx;

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    int num_classes = channel_idx ? MAX_NUM_ALF_ALTERNATIVES_CHROMA : MAX_NUM_ALF_CLASSES;
/*#else
    int num_classes = channel_idx ? 1 : MAX_NUM_ALF_CLASSES;
#endif*/    
    int num_coeffs = channel_idx ? 7 : 13;
    //m_alfCovarianceFrame[ch_type] = new AlfCovariance*[m_filterShapes[ch_type].size()];
    g_alf_covariance_frame[ch_type] = malloc(num_classes * sizeof(alf_covariance));
    for (int k = 0; k < num_classes; k++)
    {
      g_alf_covariance_frame[ch_type][k].num_coeff = num_coeffs;
      g_alf_covariance_frame[ch_type][k].num_bins = MAX_ALF_NUM_CLIPPING_VALUES;
      g_alf_covariance_frame[ch_type][k].pix_acc = 0;
      memset(g_alf_covariance_frame[ch_type][k].y, 0, sizeof(g_alf_covariance_frame[ch_type][k].y));
      memset(g_alf_covariance_frame[ch_type][k].ee, 0, sizeof(g_alf_covariance_frame[ch_type][k].ee));
    }
  }

  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    g_ctu_enable_flag[comp_idx] = malloc(num_ctus_in_pic * sizeof(*g_ctu_enable_flag[comp_idx]));
    g_ctu_enable_flag_tmp[comp_idx] = malloc(num_ctus_in_pic * sizeof(*g_ctu_enable_flag_tmp[comp_idx]));
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
   // g_ctu_enable_flag_tmp2[comp_idx] = malloc(g_num_ctus_in_pic * sizeof(*g_ctu_enable_flag_tmp2[comp_idx]));
    if (comp_idx == COMPONENT_Y)
    {
      g_ctu_alternative_tmp[comp_idx] = NULL;
      g_ctu_alternative[comp_idx] = NULL;
    }
    else
    {
      g_ctu_alternative_tmp[comp_idx] = malloc(num_ctus_in_pic * sizeof(*g_ctu_alternative_tmp[comp_idx]));
      g_ctu_alternative[comp_idx] = malloc(num_ctus_in_pic * sizeof(*g_ctu_alternative[comp_idx]));

      for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++) {
        g_ctu_alternative_tmp[comp_idx][ctu_idx] = 0;
        g_ctu_alternative[comp_idx][ctu_idx] = 0;
      }
    }
//#endif
    int num_classes = comp_idx ? 1 : MAX_NUM_ALF_CLASSES;
    int num_coeffs = comp_idx ? 7 : 13;
    int alf_num_clipping_values = MAX_ALF_NUM_CLIPPING_VALUES;

    g_alf_covariance[comp_idx] = malloc(num_ctus_in_pic * sizeof(**g_alf_covariance[comp_idx]));
    for (int j = 0; j < num_ctus_in_pic; j++)
    {
      g_alf_covariance[comp_idx][j] = malloc(num_classes * sizeof(alf_covariance));
      for (int k = 0; k < num_classes; k++)
      {
        g_alf_covariance[comp_idx][j][k].num_coeff = num_coeffs;
        g_alf_covariance[comp_idx][j][k].num_bins = alf_num_clipping_values;
        g_alf_covariance[comp_idx][j][k].pix_acc = 0;
        memset(g_alf_covariance[comp_idx][j][k].y, 0, sizeof(g_alf_covariance[comp_idx][j][k].y));
        memset(g_alf_covariance[comp_idx][j][k].ee, 0, sizeof(g_alf_covariance[comp_idx][j][k].ee));
      }
    }
  }

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  g_filter_coeff_set = malloc(/*MAX(*/MAX_NUM_ALF_CLASSES/*, MAX_NUM_ALF_ALTERNATIVES_CHROMA)*/ * sizeof(int*));
  g_filter_clipp_set = malloc(/*MAX(*/MAX_NUM_ALF_CLASSES/*, MAX_NUM_ALF_ALTERNATIVES_CHROMA)*/ * sizeof(int*));
/*#else
  g_filter_coeff_set = malloc(MAX_NUM_ALF_CLASSES * sizeof(int*));
  g_filter_clipp_set = malloc(MAX_NUM_ALF_CLASSES * sizeof(int*));
#endif*/

  for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
  {
    g_filter_coeff_set[i] = malloc(MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
    g_filter_clipp_set[i] = malloc(MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
  }

  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    g_ctb_distortion_unfilter[comp] = malloc(num_ctus_in_pic * sizeof(double));
  }
  
  g_alf_ctb_filter_index = malloc(num_ctus_in_pic * sizeof(*g_alf_ctb_filter_index));
  g_alf_ctb_filter_set_index_tmp = malloc(num_ctus_in_pic * sizeof(*g_alf_ctb_filter_set_index_tmp));

  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  const int number_of_components = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_COMPONENT;

  /*Turha, alustetaan jo ylemp‰n‰
  // init CTU stats buffers
  for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
  {
    bool is_luma = comp_idx == 0 ? 1 : 0;
    const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;

    for (int shape = 0; shape != 1 /*m_filterShapes[toChannelType(comp_id)].size()*//*; shape++)
    {
      for (int class_idx = 0; class_idx < num_classes; class_idx++)
      {
        for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++) {
          reset_alf_covariance(&g_alf_covariance[comp_idx][shape][ctu_idx][class_idx],
            g_alf_num_clipping_values[comp_idx == COMPONENT_Y ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA]);
        }
      }
    }
  }*/

  /*// Turha, alustetaan jo ylmep‰n‰ 
  // init Frame stats buffers
  const int number_of_channels = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_CHANNEL_TYPE;
  for (int channel_idx = 0; channel_idx < number_of_channels; channel_idx++)
  {
    const channel_type channel_id = channel_idx;
    bool is_luma = channel_id == CHANNEL_TYPE_LUMA ? true : false;
    //#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    const int num_alts = is_luma ? 1 : MAX_NUM_ALF_ALTERNATIVES_CHROMA;
    //#endif
    const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
    //#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    for (int alt_idx = 0; alt_idx < num_alts; ++alt_idx)
    {
      //#endif
      for (int shape = 0; shape != 1/*m_filterShapes[channel_idx].size()*/; shape++)
      {
        for (int class_idx = 0; class_idx < num_classes; class_idx++)
        {
          //#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
          reset_alf_covariance(&g_alf_covariance_frame[channel_idx][shape][is_luma ? class_idx : alt_idx], g_alf_num_clipping_values[channel_id]);
          /*#else
          reset_alf_covariance(&g_alf_covariance_frame[channel_idx][shape][class_idx], g_alf_num_clipping_values[channel_id]);
          #endif*/
        }
      }
    }
  }*/
  // init alf enable flags
  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++) {
      g_ctu_enable_flag[comp_idx][ctu_idx] = 0; //cs.picture->getAlfCtuEnableFlag( comp_idx );
      if (comp_idx != 0) {
        g_ctu_alternative[comp_idx][ctu_idx] = 0; //cs.picture->getAlfCtuAlternativeData(comp_idx);
      }
    }
  }

  const size_t simd_padding_width = 64;
  int width = state->tile->frame->width;
  int height = state->tile->frame->height;
  int stride = state->tile->frame->source->stride;
  unsigned int luma_size = (width + 8) * (height + 8);
  unsigned chroma_sizes[] = { 0, luma_size / 4, luma_size / 2, luma_size };
  unsigned chroma_size = chroma_sizes[chroma_fmt];

  alf_fulldata_buf = MALLOC_SIMD_PADDED(kvz_pixel, (luma_size + 2 * chroma_size), simd_padding_width * 2);
  alf_fulldata = &alf_fulldata_buf[4 * (width + 8) + 4] + simd_padding_width / sizeof(kvz_pixel);
  alf_tmp_y = &alf_fulldata[0];

  if (chroma_fmt == KVZ_CSP_400) {
    alf_tmp_u = NULL;
    alf_tmp_v = NULL;
  }
  else {
    alf_tmp_u = &alf_fulldata[luma_size - (4 * (width + 8) + 4) + (2 * (stride / 2) + 2)];
    alf_tmp_v = &alf_fulldata[luma_size - (4 * (width + 8) + 4) + chroma_size + (2 * (stride / 2) + 2)];
  }

  g_aps_id_cc_alf_start[0] = (int)MAX_NUM_APS;
  g_aps_id_cc_alf_start[1] = (int)MAX_NUM_APS;
  for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    int num_filters = MAX_NUM_CC_ALF_FILTERS;

    g_alf_covariance_frame_cc_alf[comp_idx - 1] = malloc(num_filters * sizeof(*g_alf_covariance_frame_cc_alf[comp_idx - 1]));
    for (int k = 0; k < num_filters; k++)
    {
      g_alf_covariance_frame_cc_alf[comp_idx - 1][k].num_coeff = 8;
      g_alf_covariance_frame_cc_alf[comp_idx - 1][k].num_bins = MAX_ALF_NUM_CLIPPING_VALUES;
      g_alf_covariance_frame_cc_alf[comp_idx - 1][k].pix_acc = 0;
      memset(g_alf_covariance_frame_cc_alf[comp_idx - 1][k].y, 0, sizeof(g_alf_covariance_frame_cc_alf[comp_idx][k].y));
      memset(g_alf_covariance_frame_cc_alf[comp_idx - 1][k].ee, 0, sizeof(g_alf_covariance_frame_cc_alf[comp_idx][k].ee));
    }

    g_alf_covariance_cc_alf[comp_idx - 1] = malloc(num_filters * sizeof(**g_alf_covariance_cc_alf[comp_idx - 1]));
    for (int j = 0; j < num_filters; j++)
    {
      g_alf_covariance_cc_alf[comp_idx - 1][j] = malloc(num_ctus_in_pic * sizeof(*g_alf_covariance_cc_alf[comp_idx - 1][j]));
      for (int k = 0; k < num_ctus_in_pic; k++)
      {
        //g_alf_covariance_cc_alf[comp_idx - 1][i][j] = malloc(g_num_ctus_in_pic * sizeof(*g_alf_covariance_cc_alf[comp_idx - 1][i][j]));
        g_alf_covariance_cc_alf[comp_idx - 1][j][k].num_coeff = 8;
        g_alf_covariance_cc_alf[comp_idx - 1][j][k].num_bins = MAX_ALF_NUM_CLIPPING_VALUES;
        g_alf_covariance_cc_alf[comp_idx - 1][j][k].pix_acc = 0;
        memset(g_alf_covariance_cc_alf[comp_idx - 1][j][k].y, 0, sizeof(g_alf_covariance_cc_alf[comp_idx - 1][j][k].y));
        memset(g_alf_covariance_cc_alf[comp_idx - 1][j][k].ee, 0, sizeof(g_alf_covariance_cc_alf[comp_idx - 1][j][k].ee));
      }
    }
  }

  g_training_cov_control = malloc(num_ctus_in_pic * sizeof(*g_training_cov_control));

  for (int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++)
  {
    g_training_distortion[i] = malloc(num_ctus_in_pic * sizeof(*g_training_distortion[i]));
  }
  g_filter_control = malloc(num_ctus_in_pic * sizeof(*g_filter_control));
  g_luma_swing_greater_than_threshold_count = malloc(num_ctus_in_pic * sizeof(*g_luma_swing_greater_than_threshold_count));
  g_chroma_sample_count_near_mid_point = malloc(num_ctus_in_pic * sizeof(*g_chroma_sample_count_near_mid_point));

  g_best_filter_control = malloc(sizeof(*g_best_filter_control) * num_ctus_in_pic);;

  g_created = true;
}

void kvz_alf_reconstruct(encoder_state_t * const state)
{
  if (g_created)
  {
    kvz_alf_reconstructor(state);
  }
}

void kvz_alf_enc_destroy(videoframe_t * const frame)
{
  if (!g_created)
  {
    return;
  }

  const int32_t num_ctus_in_pic = frame->height_in_lcu * frame->width_in_lcu;
  
  for (int channel_idx = 0; channel_idx < MAX_NUM_CHANNEL_TYPE; channel_idx++)
  {
    if (g_alf_covariance_frame[channel_idx])
    {
      FREE_POINTER(g_alf_covariance_frame[channel_idx]);
    }
  }

  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    
    if (g_ctu_enable_flag[comp_idx])
    {
      FREE_POINTER(g_ctu_enable_flag[comp_idx]);
    }
    
    if (g_ctu_enable_flag_tmp[comp_idx])
    {
      FREE_POINTER(g_ctu_enable_flag_tmp[comp_idx]);
    }

    if (g_ctu_alternative_tmp[comp_idx])
    {
      FREE_POINTER(g_ctu_alternative_tmp[comp_idx]);
    }
    if (g_ctu_alternative[comp_idx])
    {
      FREE_POINTER(g_ctu_alternative[comp_idx]);
    }
//#endif

    if (g_alf_covariance[comp_idx])
    {
      for (int k = 0; k < num_ctus_in_pic; k++)
      {
        FREE_POINTER(g_alf_covariance[comp_idx][k]);
      }
      FREE_POINTER(g_alf_covariance[comp_idx]);
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

  if (g_filter_clipp_set)
  {
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
    {
      FREE_POINTER(g_filter_clipp_set[i]);
    }
    FREE_POINTER(g_filter_clipp_set);
  }
  
  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    if (g_ctb_distortion_unfilter[comp] != NULL) {
      FREE_POINTER(g_ctb_distortion_unfilter[comp]);
    }
  }

  if (g_alf_ctb_filter_index)
  {
    FREE_POINTER(g_alf_ctb_filter_index);
  }

  if (g_alf_ctb_filter_set_index_tmp)
  {
    FREE_POINTER(g_alf_ctb_filter_set_index_tmp);
  }

  if (alf_tmp_y)
  {
    alf_tmp_y = NULL;
  }
  if (alf_tmp_u)
  {
    alf_tmp_u = NULL;
  }
  if (alf_tmp_v)
  {
    alf_tmp_v = NULL;
  }
  if (alf_fulldata)
  {
    alf_fulldata = NULL;
  }
  if (alf_fulldata_buf)
  {
    FREE_POINTER(alf_fulldata_buf);
  }

  for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    int num_filters = MAX_NUM_CC_ALF_FILTERS;
    if (g_alf_covariance_frame_cc_alf[comp_idx - 1])
    {
      FREE_POINTER(g_alf_covariance_frame_cc_alf[comp_idx - 1]);
    }

    if (g_alf_covariance_cc_alf[comp_idx - 1])
    {
      for (int j = 0; j < num_filters; j++)
      {
        FREE_POINTER(g_alf_covariance_cc_alf[comp_idx - 1][j]);
      }
      FREE_POINTER(g_alf_covariance_cc_alf[comp_idx - 1]);
    }
  }

  if (g_training_cov_control)
  {
    FREE_POINTER(g_training_cov_control);
  }

  for (int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++)
  {
    if (g_training_distortion[i])
    {
      FREE_POINTER(g_training_distortion[i]);
    }
  }

  if (g_filter_control)
  {
    FREE_POINTER(g_filter_control);
  }

  if (g_best_filter_control)
  {
    FREE_POINTER(g_best_filter_control);
  }

  if (g_luma_swing_greater_than_threshold_count)
  {
    FREE_POINTER(g_luma_swing_greater_than_threshold_count);
  }
  if (g_chroma_sample_count_near_mid_point)
  {
    FREE_POINTER(g_chroma_sample_count_near_mid_point);
  }

  if (g_classifier)
  {
    FREE_POINTER(g_classifier[0]);
    FREE_POINTER(g_classifier);
  }

  g_created = false;

  if (g_cc_alf_filter_control[0])
  {
    FREE_POINTER(g_cc_alf_filter_control[0])
  }

  if (g_cc_alf_filter_control[1])
  {
    FREE_POINTER(g_cc_alf_filter_control[1])
  }
}


void kvz_alf_encoder(encoder_state_t * const state,
  alf_aps *aps,
  channel_type channel,
  const double lambda_chroma_weight // = 0.0
  )
{
  //const TempCtx  ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));
  cabac_data_t ctx_start;
  memcpy(&ctx_start, &cabac_estimator, sizeof(ctx_start));
  //TempCtx        ctxBest(m_CtxCache);
  cabac_data_t ctx_best;

  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  kvz_config cfg = state->encoder_control->cfg;

  double cost_min = MAX_DOUBLE;
  double lambda = state->frame->lambda;
  g_bits_new_filter[channel] = 0;
  const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
  const int32_t num_ctus_in_pic = state->lcu_order_count;
  int ui_coeff_bits = 0;

  //m_alfSliceParamTemp = alfSliceParam;
  copy_alf_param(&g_alf_aps_temp, aps);

  //1. get unfiltered distortion
  if (!is_luma)
  {
    g_alf_aps_temp.num_alternatives_chroma = 1;
  }
  double cost = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[channel], channel);
  cost /= 1.001; // slight preference for unfiltered choice

  if (cost < cost_min)
  {
    cost_min = cost;
    if (is_luma) 
    {
      aps->enabled_flag[COMPONENT_Y] = 0;
    }
    else 
    {
      aps->enabled_flag[COMPONENT_Cb] = 0;
      aps->enabled_flag[COMPONENT_Cr] = 0;
    }
    // no CABAC signalling
    //ctxBest = AlfCtx(ctxStart);
    memcpy(&ctx_best, &ctx_start, sizeof(ctx_best));
    //setCtuEnableFlag(m_ctuEnableFlagTmp, channel, 0);
    set_ctu_enable_flag(g_ctu_enable_flag_tmp, channel, 0, num_ctus_in_pic);
    if (!is_luma) 
    {
      for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++) {
        g_ctu_alternative_tmp[COMPONENT_Cb][ctu_idx] = 0;
        g_ctu_alternative_tmp[COMPONENT_Cr][ctu_idx] = 0;
      }
    }
  }

  const int non_linear_flag_max =
    (is_luma ? cfg.alf_non_linear_luma : cfg.alf_non_linear_chroma) // For Chroma non linear flag is check for each alternative filter
    ? 2 : 1;

  for (int non_linear_flag = 0; non_linear_flag < non_linear_flag_max; non_linear_flag++)
  {
    for (int num_alternatives = is_luma ? 1 : MIN(num_ctus_in_pic * 2, MAX_NUM_ALF_ALTERNATIVES_CHROMA); num_alternatives > 0; num_alternatives--)
    {
      if (!is_luma)
      {
        g_alf_aps_temp.num_alternatives_chroma = num_alternatives;
      }
      //2. all CTUs are on
      if (is_luma) 
      {
        g_alf_aps_temp.enabled_flag[COMPONENT_Y] = 1;
      }
      else 
      {
        g_alf_aps_temp.enabled_flag[COMPONENT_Cb] = 1;
        g_alf_aps_temp.enabled_flag[COMPONENT_Cr] = 1;
      }
      g_alf_aps_temp.non_linear_flag[channel] = non_linear_flag;

      //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
      memcpy(&cabac_estimator, &ctx_start, sizeof(cabac_estimator));
      //setCtuEnableFlag(m_ctuEnableFlag, channel, 1);
      set_ctu_enable_flag(g_ctu_enable_flag, channel, 1, num_ctus_in_pic);
      // all alternatives are on
      if (!is_luma)
      {
        init_ctu_alternative_chroma(g_ctu_alternative, num_ctus_in_pic);
      }
      cost = kvz_alf_get_filter_coeff_and_cost(state, channel, 0, &ui_coeff_bits, true, false);

      if (cost < cost_min)
      {
        g_bits_new_filter[channel] = ui_coeff_bits;
        cost_min = cost;
        copy_alf_param_w_channel(aps, &g_alf_aps_temp, channel);
        //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
        memcpy(&ctx_best, &cabac_estimator, sizeof(ctx_best));
        //setCtuEnableFlag(m_ctuEnableFlagTmp, channel, 1);
        set_ctu_enable_flag(g_ctu_enable_flag_tmp, channel, 1, num_ctus_in_pic);
        if (!is_luma) 
        {
          memcpy(g_ctu_alternative_tmp[COMPONENT_Cb], g_ctu_alternative[COMPONENT_Cb], sizeof(uint8_t) * num_ctus_in_pic);
          memcpy(g_ctu_alternative_tmp[COMPONENT_Cr], g_ctu_alternative[COMPONENT_Cr], sizeof(uint8_t) * num_ctus_in_pic);
//#endif

      }

      //3. CTU decision
      double dist_unfilter = 0;
      double prev_it_cost = MAX_DOUBLE;
      const int iter_num = is_luma ? (2 * 4 + 1) : (2 * (2 + g_alf_aps_temp.num_alternatives_chroma - 1) + 1);
  /*#else
        cost = kvz_alf_get_filter_coeff_and_cost(state, channel, 0, &ui_coeff_bits, i_shape_idx, non_linear_flag != 0, false);
  #endif*/

      for (int iter = 0; iter < iter_num; iter++)
      {
        if ((iter & 0x01) == 0)
        {
          //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
          memcpy(&cabac_estimator, &ctx_start, sizeof(cabac_estimator));
          cost = lambda * ui_coeff_bits;
          cost += kvz_alf_derive_ctb_alf_enable_flags(state, channel, &dist_unfilter, num_classes, lambda_chroma_weight);
          if (cost < cost_min)
          {
            g_bits_new_filter[channel] = ui_coeff_bits;
            cost_min = cost;
            //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
            memcpy(&ctx_best, &cabac_estimator, sizeof(ctx_best));
            //copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, channel);
            copy_ctu_enable_flag(g_ctu_enable_flag_tmp, g_ctu_enable_flag, channel, num_ctus_in_pic);
            if (!is_luma) 
            {
              for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++) {
                g_ctu_alternative_tmp[COMPONENT_Cb][ctu_idx] = g_ctu_alternative[COMPONENT_Cb][ctu_idx];
                g_ctu_alternative_tmp[COMPONENT_Cr][ctu_idx] = g_ctu_alternative[COMPONENT_Cr][ctu_idx];
              }
            }
            copy_alf_param_w_channel(aps, &g_alf_aps_temp, channel);
  //#endif
          }
          else if (cost >= prev_it_cost)
          {
            // High probability that we have converged or we are diverging
            break;
  /*#else
            cost = kvz_alf_get_filter_coeff_and_cost(state, channel, dist_unfilter, &ui_coeff_bits, i_shape_idx, true, false);
  #endif*/
          }
          prev_it_cost = cost;
        }
        else
        {
          // unfiltered distortion is added due to some CTBs may not use filter
          // no need to reset CABAC here, since uiCoeffBits is not affected
          /*cost = */kvz_alf_get_filter_coeff_and_cost(state, channel, dist_unfilter, &ui_coeff_bits, true, false);
        }
      }//for iter
    // Decrease number of alternatives and reset ctu params and filters
    }
  }//for non_linea_flag
  //m_CABACEstimator->getCtx() = AlfCtx(ctxBest);
  memcpy(&cabac_estimator, &ctx_best, sizeof(cabac_estimator));
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  if (!is_luma) {
    memcpy(g_ctu_alternative[COMPONENT_Cb], g_ctu_alternative_tmp[COMPONENT_Cb], sizeof(uint8_t) * num_ctus_in_pic);
    memcpy(g_ctu_alternative[COMPONENT_Cr], g_ctu_alternative_tmp[COMPONENT_Cr], sizeof(uint8_t) * num_ctus_in_pic);
  }
  copy_ctu_enable_flag(g_ctu_enable_flag, g_ctu_enable_flag_tmp, channel, num_ctus_in_pic);
//#endif
}

#if !FULL_FRAME
void kvz_alf_get_avai_aps_ids_luma(encoder_state_t *const state, 
  int *new_aps_id, 
  int *aps_ids,
  int *size_of_aps_ids)
{
  param_set_map *aps_set = state->encoder_control->cfg.param_set_map;
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    state->slice->apss[i] = aps_set[i + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
  }

  //std::vector<int> result;
  int aps_id_checked = 0, cur_aps_id = g_aps_id_start;
  if (cur_aps_id < ALF_CTB_MAX_NUM_APS)
  {
    while (aps_id_checked < ALF_CTB_MAX_NUM_APS && (state->frame->slicetype == KVZ_SLICE_I) && *size_of_aps_ids < ALF_CTB_MAX_NUM_APS /*&& /*!cs.slice->getPendingRasInit()*/ && (state->frame->pictype == KVZ_NAL_IDR_W_RADL || state->frame->pictype == KVZ_NAL_IDR_N_LP))
    {
      alf_aps *cur_aps = &state->slice->apss[cur_aps_id];
      bool aps_found = aps_set[cur_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed;

      if (aps_found/*cur_aps*/ && cur_aps->t_layer/*cur_aps->getTemporalId()*/ <= state->slice->id/*cs.slice->getTLayer()*/ && cur_aps->new_filter_flag[CHANNEL_TYPE_LUMA])
      {
        //result.push_back(cur_aps_id);
        bool add_aps = true;
        for (int aps_idx = 0; aps_idx < (*size_of_aps_ids); aps_idx++) 
        {
          if (aps_ids[aps_idx] == cur_aps_id)
          {
            add_aps = false;
            continue;
          }
        }
        if (add_aps)
        {
          aps_ids[*size_of_aps_ids] = cur_aps_id;
          (*size_of_aps_ids)++;
        }
      }
      aps_id_checked++;
      cur_aps_id = (cur_aps_id + 1) % ALF_CTB_MAX_NUM_APS;
    }
  }
  state->slice->tile_group_num_aps = *size_of_aps_ids;
  for (int i = 0; i < state->slice->tile_group_num_aps; i++)
  {
    state->slice->tile_group_luma_aps_id[i] = aps_ids[i];
  }

  //*new_aps_id = g_aps_id_start - 1;
  *new_aps_id = ALF_CTB_MAX_NUM_APS - *size_of_aps_ids - 1;
  if (*new_aps_id < 0)
  {
    *new_aps_id = (int)ALF_CTB_MAX_NUM_APS - 1;
  }
  assert(*new_aps_id < (int)MAX_NUM_APS); //Wrong APS index assignment in getAvaiApsIdsLuma
}
#else
void kvz_alf_get_avai_aps_ids_luma(encoder_state_t * const state,
  int *new_aps_id,
  int *aps_ids,
  int *size_of_aps_ids)
{
  //alf_aps *apss = state->slice->apss;
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    copy_aps(&state->slice->apss[i], &state->encoder_control->cfg.param_set_map[i + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set);
  }

  //std::vector<int> result;
  int aps_id_checked = 0, cur_aps_id = g_aps_id_start;
  if (cur_aps_id < ALF_CTB_MAX_NUM_APS)
  {
    while ((aps_id_checked < ALF_CTB_MAX_NUM_APS) && state->frame->slicetype != KVZ_SLICE_I && *size_of_aps_ids < ALF_CTB_MAX_NUM_APS /*&& /*!cs.slice->getPendingRasInit()*/ && !(state->frame->pictype == KVZ_NAL_IDR_W_RADL || state->frame->pictype == KVZ_NAL_IDR_N_LP))
    {
      alf_aps *cur_aps = &state->slice->apss[cur_aps_id];
      bool aps_found = (0 <= cur_aps->aps_id && cur_aps->aps_id < ALF_CTB_MAX_NUM_APS);

      if (aps_found/*cur_aps*/ && cur_aps->layer_id == 0/*cs.slice->getPic()->layerId*/ && cur_aps->temporal_id <= state->slice->id /*cs.slice->getTLayer(*/ && cur_aps->new_filter_flag[CHANNEL_TYPE_LUMA])
      {
        for (int id = 0; id < ALF_CTB_MAX_NUM_APS; id++)
        {
          if (aps_ids[id] == -1)
          {
            aps_ids[id] = cur_aps_id;
            (*size_of_aps_ids)++;
            break;
          }
        }
      }
      aps_id_checked++;
      cur_aps_id = (cur_aps_id + 1) % ALF_CTB_MAX_NUM_APS;
    }
  }
  state->slice->tile_group_num_aps = *size_of_aps_ids;
  for (int i = 0; i < state->slice->tile_group_num_aps; i++)
  {
    state->slice->tile_group_luma_aps_id[i] = aps_ids[i];
  }

  //*new_aps_id = g_aps_id_start - 1;
  *new_aps_id =  g_aps_id_start - 1;
  if (*new_aps_id < 0)
  {
    *new_aps_id = (int)ALF_CTB_MAX_NUM_APS - 1;
  }
  assert(*new_aps_id < (int)ALF_CTB_MAX_NUM_APS); //Wrong APS index assignment in getAvaiApsIdsLuma
}

void kvz_alf_derive_stats_for_filtering(encoder_state_t * const state)
{
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  bool chroma_scale_x = (chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool chroma_scale_y = (chroma_fmt != KVZ_CSP_420) ? 0 : 1;

  const int32_t num_ctus_in_pic = state->lcu_order_count;
  const int alf_vb_luma_ctu_height = LCU_WIDTH;
  const int alf_vb_chma_ctu_height = (LCU_WIDTH >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0));
  const int alf_vb_luma_pos = LCU_WIDTH - ALF_VB_POS_ABOVE_CTUROW_LUMA;
  const int alf_vb_chma_pos = (LCU_WIDTH >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0)) - ALF_VB_POS_ABOVE_CTUROW_CHMA;
  int32_t pic_width = state->tile->frame->width;
  int32_t pic_height = state->tile->frame->height;
  int ctu_rs_addr = 0;

  const int number_of_components = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_COMPONENT;

  // init CTU stats buffers
  for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
  {
    bool is_luma = comp_idx == 0 ? 1 : 0;
    const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;

    for (int class_idx = 0; class_idx < num_classes; class_idx++)
    {
      for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
      {
        reset_alf_covariance(&g_alf_covariance[comp_idx][ctu_idx][class_idx], MAX_ALF_NUM_CLIPPING_VALUES);
      }
    }
  }

  //kerran jossain muualla (kai?)
  // init Frame stats buffers
  const int number_of_channels = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_CHANNEL_TYPE;
  for (int channel_idx = 0; channel_idx < number_of_channels; channel_idx++)
  {
    const channel_type channel_id = channel_idx;
    const int num_classes = channel_id == CHANNEL_TYPE_LUMA ? MAX_NUM_ALF_CLASSES : 1;
    for (int class_idx = 0; class_idx < num_classes; class_idx++)
    {
      reset_alf_covariance(&g_alf_covariance_frame[channel_idx][class_idx], MAX_ALF_NUM_CLIPPING_VALUES);
    }
}

  int max_cu_height = LCU_WIDTH;

#if FULL_FRAME
  for (int y_pos = 0; y_pos < pic_height; y_pos += LCU_WIDTH)
  {
    for (int x_pos = 0; x_pos < pic_width; x_pos += LCU_WIDTH)
    {
      const int width = (x_pos + LCU_WIDTH > pic_width) ? (pic_width - x_pos) : LCU_WIDTH;
      const int height = (y_pos + LCU_WIDTH > pic_height) ? (pic_height - y_pos) : LCU_WIDTH;
      {
        for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
        {
          const bool is_luma = comp_idx == COMPONENT_Y ? 1 : 0;
          channel_type ch_type = is_luma ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;

          int blk_w = is_luma ? width : width >> chroma_scale_x;
          int blk_h = is_luma ? height : height >> chroma_scale_y;
          int pos_x = is_luma ? x_pos : x_pos >> chroma_scale_x;
          int pos_y = is_luma ? y_pos : y_pos >> chroma_scale_y;

          int32_t org_stride = is_luma ? state->tile->frame->source->stride : state->tile->frame->source->stride >> chroma_scale_x;
          int32_t rec_stride = is_luma ? state->tile->frame->rec->stride    : state->tile->frame->rec->stride >> chroma_scale_x;

          kvz_pixel *org = comp_idx ? (comp_idx - 1 ? &state->tile->frame->source->v[pos_x + pos_y * org_stride] : &state->tile->frame->source->u[pos_x + pos_y * org_stride]) : &state->tile->frame->source->y[pos_x + pos_y * org_stride];
          kvz_pixel *rec = comp_idx ? (comp_idx - 1 ? &state->tile->frame->rec->v[pos_x + pos_y * rec_stride]    : &state->tile->frame->rec->u[pos_x + pos_y * rec_stride])    : &state->tile->frame->rec->y[pos_x + pos_y * rec_stride];

          kvz_alf_get_blk_stats(state, ch_type,
            &g_alf_covariance[comp_idx][ctu_rs_addr],
            comp_idx ? NULL : g_classifier,
            org, org_stride, rec, rec_stride, pos_x, pos_y, pos_x, pos_y, blk_w, blk_h,
            (is_luma ? alf_vb_luma_ctu_height : alf_vb_chma_ctu_height),
            (is_luma) ? alf_vb_luma_pos : alf_vb_chma_pos
          );

          const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;

          for (int class_idx = 0; class_idx < num_classes; class_idx++)
          {
            add_alf_cov(&g_alf_covariance_frame[ch_type][is_luma ? class_idx : 0],
              &g_alf_covariance[comp_idx][ctu_rs_addr][class_idx]
            );
          }
      }
      }
      ctu_rs_addr++;
    }
  }
}

void kvz_alf_get_blk_stats(encoder_state_t * const state,
  channel_type channel,
  alf_covariance **alf_covariance,
  alf_classifier **g_classifier,
  kvz_pixel *org,
  int32_t org_stride,
  kvz_pixel *rec,
  int32_t rec_stride,
  const int x_pos,
  const int y_pos,
  const int x_dst,
  const int y_dst,
  const int width,
  const int height,
  int vb_ctu_height,
  int vb_pos)
{
  int16_t e_local[MAX_NUM_ALF_LUMA_COEFF][MAX_ALF_NUM_CLIPPING_VALUES];

  const int num_bins = MAX_ALF_NUM_CLIPPING_VALUES;

  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int transpose_idx = 0;
  int class_idx = 0;

  for (int i = 0; i < height; i++)
  {
    int vb_distance = ((y_dst + i) % vb_ctu_height) - vb_pos;
    for (int j = 0; j < width; j++)
    {
      if (g_classifier && g_classifier[y_dst + i][x_dst + j].class_idx == ALF_UNUSED_CLASS_IDX && g_classifier[y_dst + i][x_dst + j].transpose_idx == ALF_UNUSED_TRANSPOSE_IDX)
      {
        continue;
      }
      memset(e_local, 0, sizeof(e_local));
      if (g_classifier)
      {
        alf_classifier* cl = &g_classifier[y_dst + i][x_dst + j];
        transpose_idx = cl->transpose_idx;
        class_idx = cl->class_idx;
      }

      double weight = 1.0;
      if (0/*m_alfWSSD*/)
      {
        //weight = g_luma_level_to_weight_plut[org[j]];
      }
      int16_t y_local = org[j] - rec[j];
      kvz_alf_calc_covariance(e_local, rec + j, rec_stride, channel, transpose_idx, vb_distance);
      for (int k = 0; k < num_coeff; k++)
      {
        for (int l = k; l < num_coeff; l++)
        {
          for (int b0 = 0; b0 < num_bins; b0++)
          {
            for (int b1 = 0; b1 < num_bins; b1++)
            {
              if (0/*m_alfWSSD*/)
              {
                (*alf_covariance)[class_idx].ee[b0][b1][k][l] += weight * (e_local[k][b0] * (double)e_local[l][b1]);
              }
              else
              {
                (*alf_covariance)[class_idx].ee[b0][b1][k][l] += e_local[k][b0] * (double)e_local[l][b1];
              }
            }
          }
        }
        for (int b = 0; b < num_bins; b++)
        {
          if (0/*m_alfWSSD*/)
          {
            (*alf_covariance)[class_idx].y[b][k] += weight * (e_local[k][b] * (double)y_local);
          }
          else
          {
            (*alf_covariance)[class_idx].y[b][k] += e_local[k][b] * (double)y_local;
          }
        }
      }
      if (0/*m_alfWSSD*/)
      {
        (*alf_covariance)[class_idx].pix_acc += weight * (y_local * (double)y_local);
      }
      else
      {
        (*alf_covariance)[class_idx].pix_acc += y_local * (double)y_local;
      }
    }
    org += org_stride;
    rec += rec_stride;
  }

  int num_classes = g_classifier ? MAX_NUM_ALF_CLASSES : 1;
  for (class_idx = 0; class_idx < num_classes; class_idx++)
  {
    for (int k = 1; k < num_coeff; k++)
    {
      for (int l = 0; l < k; l++)
      {
        for (int b0 = 0; b0 < num_bins; b0++)
        {
          for (int b1 = 0; b1 < num_bins; b1++)
          {
            (*alf_covariance)[class_idx].ee[b0][b1][k][l] = (*alf_covariance)[class_idx].ee[b1][b0][l][k];
          }
        }
      }
    }
  }
}

void kvz_alf_calc_covariance(int16_t e_local[MAX_NUM_ALF_LUMA_COEFF][MAX_ALF_NUM_CLIPPING_VALUES],
  const kvz_pixel *rec,
  const int stride,
  const channel_type channel,
  const int transpose_idx,
  int vb_distance)
{
  int clip_top_row = -4;
  int clip_bot_row = 4;
  if (vb_distance >= -3 && vb_distance < 0)
  {
    clip_bot_row = -vb_distance - 1;
    clip_top_row = -clip_bot_row; // symmetric
  }
  else if (vb_distance >= 0 && vb_distance < 3)
  {
    clip_top_row = -vb_distance;
    clip_bot_row = -clip_top_row; // symmetric
  }

  const bool is_luma = channel == CHANNEL_TYPE_LUMA;
  const int *filter_pattern = is_luma ? alf_pattern_7 : alf_pattern_5;
  const int half_filter_length = (is_luma ? 7 : 5) >> 1;
  const short* clip = g_alf_clipping_values[channel];
  const int num_bins = MAX_ALF_NUM_CLIPPING_VALUES;

  int k = 0;

  const int16_t curr = rec[0];

  if (transpose_idx == 0)
  {
    for (int i = -half_filter_length; i < 0; i++)
    {
      const kvz_pixel* rec0 = rec + MAX(i, clip_top_row) * stride;
      const kvz_pixel* rec1 = rec - MAX(i, -clip_bot_row) * stride;
      for (int j = -half_filter_length - i; j <= half_filter_length + i; j++, k++)
      {
        for (int b = 0; b < num_bins; b++)
        {
          e_local[filter_pattern[k]][b] += clip_alf(clip[b], curr, rec0[j], rec1[-j]);
        }
      }
    }
    for (int j = -half_filter_length; j < 0; j++, k++)
    {
      for (int b = 0; b < num_bins; b++)
      {
        e_local[filter_pattern[k]][b] += clip_alf(clip[b], curr, rec[j], rec[-j]);
      }
    }
  }
  else if (transpose_idx == 1)
  {
    for (int j = -half_filter_length; j < 0; j++)
    {
      const kvz_pixel* rec0 = rec + j;
      const kvz_pixel* rec1 = rec - j;

      for (int i = -half_filter_length - j; i <= half_filter_length + j; i++, k++)
      {
        for (int b = 0; b < num_bins; b++)
        {
          e_local[filter_pattern[k]][b] += clip_alf(clip[b], curr, rec0[MAX(i, clip_top_row) * stride], rec1[-MAX(i, -clip_bot_row) * stride]);
        }
      }
    }
    for (int i = -half_filter_length; i < 0; i++, k++)
    {
      for (int b = 0; b < num_bins; b++)
      {
        e_local[filter_pattern[k]][b] += clip_alf(clip[b], curr, rec[MAX(i, clip_top_row) * stride], rec[-MAX(i, -clip_bot_row) * stride]);
      }
    }
  }
  else if (transpose_idx == 2)
  {
    for (int i = -half_filter_length; i < 0; i++)
    {
      const kvz_pixel* rec0 = rec + MAX(i, clip_top_row) * stride;
      const kvz_pixel* rec1 = rec - MAX(i, -clip_bot_row) * stride;

      for (int j = half_filter_length + i; j >= -half_filter_length - i; j--, k++)
      {
        for (int b = 0; b < num_bins; b++)
        {
          e_local[filter_pattern[k]][b] += clip_alf(clip[b], curr, rec0[j], rec1[-j]);
        }
      }
    }
    for (int j = -half_filter_length; j < 0; j++, k++)
    {
      for (int b = 0; b < num_bins; b++)
      {
        e_local[filter_pattern[k]][b] += clip_alf(clip[b], curr, rec[j], rec[-j]);
      }
    }
  }
  else
  {
    for (int j = -half_filter_length; j < 0; j++)
    {
      const kvz_pixel* rec0 = rec + j;
      const kvz_pixel* rec1 = rec - j;

      for (int i = half_filter_length + j; i >= -half_filter_length - j; i--, k++)
      {
        for (int b = 0; b < num_bins; b++)
        {
          e_local[filter_pattern[k]][b] += clip_alf(clip[b], curr, rec0[MAX(i, clip_top_row) * stride], rec1[-MAX(i, -clip_bot_row) * stride]);
        }
      }
    }
    for (int i = -half_filter_length; i < 0; i++, k++)
    {
      for (int b = 0; b < num_bins; b++)
      {
        e_local[filter_pattern[k]][b] += clip_alf(clip[b], curr, rec[MAX(i, clip_top_row) * stride], rec[-MAX(i, -clip_bot_row) * stride]);
      }
    }

  }
  for (int b = 0; b < num_bins; b++)
  {
    e_local[filter_pattern[k]][b] += curr;
  }
}

double kvz_alf_get_filter_coeff_and_cost(encoder_state_t * const state,
  channel_type channel,
  double dist_unfilter,
  int *ui_coeff_bits,
  bool b_re_collect_stat,
  bool only_filter_cost)
{
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  double lambda = state->frame->lambda;
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
  const int8_t bit_depth = state->encoder_control->bitdepth;
  const int32_t num_ctus_in_pic = state->lcu_order_count;

  //collect stat based on CTU decision
  if (b_re_collect_stat)
  {
    get_frame_stats(channel, num_ctus_in_pic);
  }

  double dist = dist_unfilter;
  (*ui_coeff_bits) = 0;
  /*#if !JVET_O0491_HLS_CLEANUP
  int ui_slice_flag = 0;*/
  
  //get filter coeff
  if (is_luma)
  {
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    //Tarvitaanko t‰t‰ alustusta ollenkaan?
    const int fill_val = g_alf_aps_temp.non_linear_flag[channel][0] ? g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] / 2 : 0;
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++) {
      for (int j = 0; j < MAX_NUM_ALF_CLASSES; j++) {
        for (int k = 0; k < MAX_NUM_ALF_LUMA_COEFF; k++) {
          clip_merged[i][j][k] = fill_val;
        }
      }
    }

/*#else
    std::fill_n(m_alfClipMerged[iShapeIdx][0][0], MAX_NUM_ALF_LUMA_COEFF*MAX_NUM_ALF_CLASSES*MAX_NUM_ALF_CLASSES, m_alfParamTemp.nonLinearFlag[channel] ? AlfNumClippingValues[CHANNEL_TYPE_LUMA] / 2 : 0);
#endif*/

    // Reset Merge Tmp Cov
    reset_alf_covariance(&g_alf_covariance_merged[MAX_NUM_ALF_CLASSES], MAX_ALF_NUM_CLIPPING_VALUES);
    reset_alf_covariance(&g_alf_covariance_merged[MAX_NUM_ALF_CLASSES + 1], MAX_ALF_NUM_CLIPPING_VALUES);
    //distortion
    dist += kvz_alf_merge_filters_and_cost(state, &g_alf_aps_temp, channel, ui_coeff_bits, g_alf_covariance_frame[channel], g_alf_covariance_merged, clip_merged);
  }
  else
  {
    //distortion
/*#if !JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    assert(num_coeff == g_alf_covariance_frame[channel][i_shape_idx][0].num_coeff);
    //std::fill_n(m_filterClippSet[0], MAX_NUM_ALF_CHROMA_COEFF, m_alfParamTemp.non_linear_flag[channel] ? AlfNumClippingValues[CHANNEL_TYPE_CHROMA] / 2 : 0);
    const int fill_val = g_alf_aps_temp.non_linear_flag[channel] ? g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] / 2 : 0;
    for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++) {
      g_filter_clipp_set[0][i] = fill_val;
    }
    dist += g_alf_covariance_frame[channel][i_shape_idx][0].pix_acc + kvz_alf_derive_coeff_quant(channel, g_filter_clipp_set[0], g_filter_coeff_set[0], &g_alf_covariance_frame[channel][i_shape_idx][0], ALF_NUM_BITS, g_alf_aps_temp.non_linear_flag[channel]);
#endif*/

/*#if !JVET_O0491_HLS_CLEANUP
    //setEnableFlag( m_alfSliceParamTemp, channel, m_ctuEnableFlag );
    const int alf_chroma_idc = g_alf_aps_temp.enabled_flag[COMPONENT_Cb] * 2 + g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
#endif*/

    for (int alt_idx = 0; alt_idx < g_alf_aps_temp.num_alternatives_chroma; ++alt_idx)
    {
      assert(num_coeff == g_alf_covariance_frame[channel][alt_idx].num_coeff);
      alf_aps best_slice_param;
      double best_cost = MAX_DOUBLE;
      double best_dist = MAX_DOUBLE;
      int best_coeff_bits = 0;
      const int non_linear_flag_max = state->encoder_control->cfg.alf_non_linear_chroma ? 2 : 1;

      for (int non_linear_flag = 0; non_linear_flag < non_linear_flag_max; non_linear_flag++)
      {
        int current_non_linear_flag = g_alf_aps_temp.non_linear_flag[channel] ? 1 : 0;
        if (non_linear_flag != current_non_linear_flag)
        {
          continue;
        }

        int fill_val = non_linear_flag ? MAX_ALF_NUM_CLIPPING_VALUES / 2 : 0;
        for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++) {
          g_filter_clipp_set[alt_idx][i] = fill_val;
        }

        double dist = g_alf_covariance_frame[channel][alt_idx].pix_acc + kvz_alf_derive_coeff_quant(channel, g_filter_clipp_set[alt_idx], g_filter_coeff_set[alt_idx], &g_alf_covariance_frame[channel][alt_idx], bit_depth, non_linear_flag);
        for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
        {
          g_alf_aps_temp.chroma_coeff[alt_idx][i] = g_filter_coeff_set[alt_idx][i];
          g_alf_aps_temp.chroma_clipp[alt_idx][i] = g_filter_clipp_set[alt_idx][i];
        }
        int coeff_bits = get_chroma_coeff_rate(&g_alf_aps_temp, alt_idx);
        double cost = dist + lambda * coeff_bits;
        if (cost < best_cost)
        {
          best_cost = cost;
          best_dist = dist;
          best_coeff_bits = coeff_bits;
          copy_alf_param(&best_slice_param, &g_alf_aps_temp);
        }
      }
      *ui_coeff_bits += best_coeff_bits;
      dist += best_dist;
      copy_alf_param(&g_alf_aps_temp, &best_slice_param);
    }
    (*ui_coeff_bits) += length_uvlc(g_alf_aps_temp.num_alternatives_chroma - 1);
    (*ui_coeff_bits)++;
/*#if !JVET_O0491_HLS_CLEANUP
    uiSliceFlag = lengthTruncatedUnary(alfChromaIdc, 3)
      - lengthTruncatedUnary(0, 3);  // rate already put on Luma
#endif*/
/*#else
    for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
    {
      g_alf_aps_temp.chroma_coeff[i] = g_filter_coeff_set[0][i];
      g_alf_aps_temp.chroma_clipp[i] = g_filter_clipp_set[0][i];
    }
#endif*/
  }

  if (only_filter_cost)
  {
    return dist + lambda * (*ui_coeff_bits);
  }
  /*#if !JVET_O0491_HLS_CLEANUP
  double rate = *ui_coeff_bits + ui_slice_flag;*/
  double rate = *ui_coeff_bits;
  //m_CABACEstimator->resetBits();
  kvz_cabac_reset_bits(&cabac_estimator);
  //m_CABACEstimator->codeAlfCtuEnableFlags(cs, channel, &m_alfParamTemp);
  code_alf_ctu_enable_flags_channel(state, &cabac_estimator, channel, &g_alf_aps_temp);

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  //for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++)
  {
    if (is_luma)
    {
      // Evaluate cost of signaling filter set index for convergence of filters enabled flag / filter derivation
      assert(g_alf_ctb_filter_index[ctu_idx] == ALF_NUM_FIXED_FILTER_SETS);
      assert(state->slice->tile_group_num_aps == 1);
      //m_CABACEstimator->codeAlfCtuFilterIndex(cs, ctu_idx, &m_alfParamTemp.enabledFlag[COMPONENT_Y]);
      code_alf_ctu_filter_index(state, &cabac_estimator, ctu_idx, g_alf_aps_temp.enabled_flag[COMPONENT_Y]);
    }
  }
  //m_CABACEstimator->codeAlfCtuAlternatives(cs, channel, &m_alfParamTemp);
  code_alf_ctu_alternatives_channel(state, &cabac_estimator, channel, &g_alf_aps_temp);
//#endif

  rate += (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale * 0;/*(double)m_CABACEstimator->getEstFracBits();*/
  return dist + lambda * rate;
}

int kvz_alf_derive_filter_coefficients_prediction_mode(channel_type channel,
  int **filter_set,
  const int num_filters)
{
  return (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? get_cost_filter_clipp(channel, filter_set, num_filters) : 0) + get_cost_filter_coeff(channel, filter_set, num_filters);
/* #else
int kvz_alf_derive_filter_coefficients_prediction_mode(channel_type channel,
  int **filter_set,
  int** filter_coeff_diff,
  const int num_filters)
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
    */
}

void kvz_alf_merge_classes(channel_type channel,
  alf_covariance* cov,
  alf_covariance* cov_merged,
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  const int num_classes, 
  short filter_indices[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES])
{
  int tmp_clip[MAX_NUM_ALF_LUMA_COEFF];
  int best_merge_clip[MAX_NUM_ALF_LUMA_COEFF];
  double err[MAX_NUM_ALF_CLASSES];
  double best_merge_err = MAX_DOUBLE;
  bool available_class[MAX_NUM_ALF_CLASSES];
  int8_t index_list[MAX_NUM_ALF_CLASSES];
  int8_t index_list_temp[MAX_NUM_ALF_CLASSES];
  int num_remaining = num_classes;

  memset(filter_indices, 0, sizeof(short) * MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_CLASSES);

  for (int i = 0; i < num_classes; i++)
  {
    filter_indices[num_remaining - 1][i] = i;
    index_list[i] = i;
    available_class[i] = true;
    //cov_merged[i] = cov[i];
    copy_cov(&cov_merged[i], &cov[i]);
    
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    cov_merged[i].num_bins = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA][0] ? g_alf_num_clipping_values[COMPONENT_Y] : 1;
/*#else
    cov_merged[i].num_bins = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? g_alf_num_clipping_values[COMPONENT_Y] : 1;
#endif*/
  }

  // Try merging different covariance matrices

  // temporal AlfCovariance structure is allocated as the last element in covMerged array, the size of covMerged is MAX_NUM_ALF_CLASSES + 1
  alf_covariance* tmp_cov = &cov_merged[MAX_NUM_ALF_CLASSES];

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  tmp_cov->num_bins = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA][0] ? g_alf_num_clipping_values[COMPONENT_Y] : 1;
/*#else
  tmp_cov->num_bins = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? g_alf_num_clipping_values[COMPONENT_Y] : 1;
#endif*/

  // init Clip
  for (int i = 0; i < num_classes; i++)
  {
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    for (int val = 0; val < MAX_NUM_ALF_LUMA_COEFF; val++) {
      clip_merged[num_remaining - 1][i][val] = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? MAX_ALF_NUM_CLIPPING_VALUES / 2 : 0;
    }
    if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
/*#else
    for (int val = 0; val < MAX_NUM_ALF_LUMA_COEFF; val++) {
      clip_merged[num_remaining - 1][i][val] = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] / 2 : 0;
    }
    if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
#endif*/
    {
      err[i] = optimize_filter_clip(&cov_merged[i], clip_merged[num_remaining - 1][i]);
    }
    else
    {
      err[i] = calculate_error_opt_filt(&cov_merged[i], clip_merged[num_remaining - 1][i]);
    }
  }

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
            double error1 = err[i];
            double error2 = err[j];

            add_alf_cov_lhs_rhs(tmp_cov, &cov_merged[i], &cov_merged[j]);
            for (int l = 0; l < MAX_NUM_ALF_LUMA_COEFF; ++l)
            {
              tmp_clip[l] = (clip_merged[num_remaining - 1][i][l] + clip_merged[num_remaining - 1][j][l] + 1) >> 1;
            }

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
            double error_merged = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA][0] ? optimize_filter_clip(tmp_cov, tmp_clip) : calculate_error_opt_filt(tmp_cov, tmp_clip);
/*#else
            double error_merged = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? optimize_filter_clip(tmp_cov, tmp_clip) : calculate_error_opt_filt(tmp_cov, tmp_clip);
#endif*/
            double error = error_merged - error1 - error2;

            if (error < error_min)
            {
              best_merge_err = error_merged;
              memcpy(best_merge_clip, tmp_clip, sizeof(best_merge_clip));
              error_min = error;
              best_to_merge_idx1 = i;
              best_to_merge_idx2 = j;
            }
          }
        }
      }
    }

    add_alf_cov(&cov_merged[best_to_merge_idx1], &cov_merged[best_to_merge_idx2]);
    memcpy(clip_merged[num_remaining - 2], clip_merged[num_remaining - 1], sizeof(int[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF]));
    memcpy(clip_merged[num_remaining - 2][best_to_merge_idx1], best_merge_clip, sizeof(best_merge_clip));
    err[best_to_merge_idx1] = best_merge_err;
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
      memcpy(index_list_temp, index_list, sizeof(int8_t) * num_classes);

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

double kvz_alf_merge_filters_and_cost(encoder_state_t * const state,
  alf_aps *alf_aps,
  channel_type channel,
  int *ui_coeff_bits,
  alf_covariance *cov_frame,
  alf_covariance *cov_merged, 
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF])
{
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int num_filters_best = 0;
  int num_filters = MAX_NUM_ALF_CLASSES;
  bool coded_var_bins[MAX_NUM_ALF_CLASSES];
  double error_force_0_coeff_tab[MAX_NUM_ALF_CLASSES][2];
  double lambda = state->frame->lambda;
  const int8_t bit_depth = state->encoder_control->bitdepth;

  double cost, cost0, dist, dist_force0, cost_min = MAX_DOUBLE;
  int coeff_bits, coeff_bits_force0;

  //clip_merged:i‰ ei tarvitse nollata ennen
  kvz_alf_merge_classes(channel, cov_frame, cov_merged, clip_merged, MAX_NUM_ALF_CLASSES, g_filter_indices);

  while (num_filters >= 1)
  {
    dist = kvz_alf_derive_filter_coeffs(alf_aps, channel, cov_frame, cov_merged, g_filter_indices[num_filters-1], num_filters, error_force_0_coeff_tab, clip_merged, bit_depth);

    // filter coeffs are stored in m_filterCoeffSet
    dist_force0 = get_dist_force_0(channel, num_filters, error_force_0_coeff_tab, coded_var_bins, lambda);
    coeff_bits = kvz_alf_derive_filter_coefficients_prediction_mode(channel, g_filter_coeff_set, num_filters);
    coeff_bits_force0 = get_cost_filter_coeff_force_0(channel, g_filter_coeff_set, num_filters, coded_var_bins);

    cost = dist + lambda * coeff_bits;
    cost0 = dist_force0 + lambda * coeff_bits_force0;

    if (cost0 < cost)
    {
      cost = cost0;
    }
    /* #if !JVET_O0669_REMOVE_ALF_COEFF_PRED
    if (*fixed_filter_set_index > 0)
    {
      int len = 0;
      len += get_tb_length(*fixed_filter_set_index - 1, ALF_NUM_FIXED_FILTER_SETS);
      len += 1; //fixed filter flag pattern
      if (*fixed_filter_pattern > 0)
      {
        len += MAX_NUM_ALF_CLASSES;  //"fixed_filter_flag" for each class
      }
      cost += g_lambda[COMPONENT_Y] * len;
    }*/

    if (cost <= cost_min)
    {
      cost_min = cost;
      num_filters_best = num_filters;
      //best_pred_mode = pred_mode;   #if !JVET_O0669_REMOVE_ALF_COEFF_PRED
    }
    num_filters--;
  }

  dist = kvz_alf_derive_filter_coeffs(alf_aps, channel, cov_frame, cov_merged, g_filter_indices[num_filters_best - 1], num_filters_best, error_force_0_coeff_tab, clip_merged, bit_depth);

  coeff_bits = kvz_alf_derive_filter_coefficients_prediction_mode(channel, g_filter_coeff_set, num_filters_best);
  dist_force0 = get_dist_force_0(channel, num_filters_best, error_force_0_coeff_tab, coded_var_bins, lambda);
  coeff_bits_force0 = get_cost_filter_coeff_force_0(channel, g_filter_coeff_set, num_filters_best, coded_var_bins);

  cost = dist + lambda * coeff_bits;
  cost0 = dist_force0 + lambda * coeff_bits_force0;

  alf_aps->num_luma_filters = num_filters_best;
  double dist_return;
  if (cost <= cost0)
  {
    dist_return = dist;
    alf_aps->alf_luma_coeff_delta_flag = 0;
    *ui_coeff_bits = coeff_bits;
    //alf_aps->alf_luma_coeff_delta_prediction_flag = best_pred_mode;   #if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  }
  else
  {
    dist_return = dist_force0;
    alf_aps->alf_luma_coeff_delta_flag = 1;
    *ui_coeff_bits = coeff_bits_force0;
    memcpy(alf_aps->alf_luma_coeff_flag, coded_var_bins, sizeof(coded_var_bins));
    //*alf_luma_coeff_delta_prediction_flag = 0;    #if !JVET_O0669_REMOVE_ALF_COEFF_PRED

    for (int var_ind = 0; var_ind < num_filters_best; var_ind++)
    {
      if (coded_var_bins[var_ind] == 0)
      {
        memset(g_filter_coeff_set[var_ind], 0, sizeof(int) * MAX_NUM_ALF_LUMA_COEFF);
        memset(g_filter_clipp_set[var_ind], 0, sizeof(int) * MAX_NUM_ALF_LUMA_COEFF);
      }
    }
  }

  for (int ind = 0; ind < alf_aps->num_luma_filters; ++ind)
  {
    for (int i = 0; i < num_coeff; i++)
    {
      // #if JVET_O0669_REMOVE_ALF_COEFF_PRED
      alf_aps->luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_filter_coeff_set[ind][i];
      /* #else
      if (alf_aps->alf_luma_coeff_delta_prediction_flag)
      {
        alf_aps->luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_diff_filter_coeff[ind][i];
      }
      else
      {
        alf_aps->luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_filter_coeff_set[ind][i];
      }*/
      alf_aps->luma_clipp[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_filter_clipp_set[ind][i];
    }
  }

  memcpy(alf_aps->filter_coeff_delta_idx, g_filter_indices[num_filters_best - 1], sizeof(short) * MAX_NUM_ALF_CLASSES);
  *ui_coeff_bits += get_non_filter_coeff_rate(alf_aps);
  return dist_return;
}

double kvz_alf_derive_filter_coeffs(alf_aps *aps,
  channel_type channel,
  alf_covariance *cov,
  alf_covariance *covMerged,
  short* filter_indices,
  int num_filters,
  double error_tab_force_0_coeff[MAX_NUM_ALF_CLASSES][2],
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  const int bit_depth)
{
  // #if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  //int *fixed_filter_pattern = &aps->fixed_filter_pattern;
  //int *fixed_filter_idx = aps->fixed_filter_idx;
  //int *fixed_filter_set_index = &aps->fixed_filter_set_index;

  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
   
  double error = 0.0;
  alf_covariance *tmp_cov = &covMerged[MAX_NUM_ALF_CLASSES];

  /* #if !JVET_O0669_REMOVE_ALF_COEFF_PRED

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
  }*/

  for( int filt_idx = 0; filt_idx < num_filters; filt_idx++ )
  {
    reset_alf_covariance(tmp_cov, -1);

    bool found_clip = false;
    for( int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++ )
    {
      if( filter_indices[class_idx] == filt_idx )
      {
        //tmp_cov += cov[class_idx];
        add_alf_cov(tmp_cov, &cov[class_idx]);

        /* #if !JVET_O0669_REMOVE_ALF_COEFF_PRED  
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
        */
        if (!found_clip)
        {
          found_clip = true; // clip should be at the adress of shortest one
          memcpy(g_filter_clipp_set[filt_idx], clip_merged[num_filters - 1][class_idx], sizeof(int[MAX_NUM_ALF_LUMA_COEFF]));
        }
      }
    }

    // Find coeffcients
    assert(num_coeff == tmp_cov->num_coeff);
    error_tab_force_0_coeff[filt_idx][1] = tmp_cov->pix_acc + kvz_alf_derive_coeff_quant(channel, g_filter_clipp_set[filt_idx], g_filter_coeff_set[filt_idx], tmp_cov, bit_depth, false);
    error_tab_force_0_coeff[filt_idx][0] = tmp_cov->pix_acc;
    error += error_tab_force_0_coeff[filt_idx][1];
  }
  return error;
}

double kvz_alf_derive_coeff_quant(channel_type channel,
  int *filter_clipp,
  int *filter_coeff_quant,
  const alf_covariance* cov,
  const int bit_depth,
  const bool optimize_clip)
{
  const bool is_luma = channel == CHANNEL_TYPE_LUMA ? true : false;
  const int num_coeff = is_luma ? 13 : 7;
  const int factor = 1 << (bit_depth - 1);
  const int max_value = factor - 1;
  const int min_value = -factor + 1;
  double filter_coeff[MAX_NUM_ALF_LUMA_COEFF];

  optimize_filter(cov, filter_clipp, filter_coeff, optimize_clip);

  //roundFiltCoeff(filter_coeff_quant, filter_coeff, num_coeff, factor);
  for (int i = 0; i < num_coeff; i++)
  {
    int sign = filter_coeff[i] > 0 ? 1 : -1;
    filter_coeff_quant[i] = (int)(filter_coeff[i] * sign * factor + 0.5) * sign;
  }

  for (int i = 0; i < num_coeff - 1; i++)
  {
    filter_coeff_quant[i] = MIN(max_value, MAX(min_value, filter_coeff_quant[i]));
  }
  filter_coeff_quant[num_coeff - 1] = 0;

  int modified = 1;

  double err_ref = calc_error_for_coeffs(cov, filter_clipp, filter_coeff_quant, num_coeff, bit_depth);
  int sign;
  while (modified)
  {
    modified = 0;
    for (int sign_count = 0; sign_count <= 1; sign_count++)
    {
      sign = sign_count == 0 ? 1 : -1;
      double err_min = MAX_DOUBLE;
      int min_ind = -1;

      for (int k = 0; k < num_coeff - 1; k++)
      {
        if (filter_coeff_quant[k] - sign > max_value || filter_coeff_quant[k] - sign < min_value)
        {
          continue;
        }

        filter_coeff_quant[k] -= sign;

        double error = calc_error_for_coeffs(cov, filter_clipp, filter_coeff_quant, num_coeff, bit_depth);
        if (error < err_min)
        {
          err_min = error;
          min_ind = k;
        }
        filter_coeff_quant[k] += sign;
      }
      if (err_min < err_ref)
      {
        filter_coeff_quant[min_ind] -= sign;
        modified++;
        err_ref = err_min;
      }
    }
  }

  return err_ref;
}

void kvz_alf_encoder_ctb(encoder_state_t * const state,
  alf_aps *aps,
  const double lambda_chroma_weight)
{
  //TempCtx        ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));
  cabac_data_t ctx_start;
  memcpy(&ctx_start, &cabac_estimator, sizeof(ctx_start));
  //TempCtx        ctxBest(m_CtxCache);
  cabac_data_t ctx_best;
  //TempCtx        ctxTempStart(m_CtxCache);
  cabac_data_t ctx_temp_start;
  //TempCtx        ctxTempBest(m_CtxCache);*/
  cabac_data_t ctx_temp_best;
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  /*TempCtx        ctxTempAltStart(m_CtxCache);
  TempCtx        ctxTempAltBest(m_CtxCache);*/
  cabac_data_t ctx_temp_alt_start;
  //cabac_data_t ctx_temp_alt_best;
//#endif

  int best_aps_ids[ALF_CTB_MAX_NUM_APS] = { -1, -1, -1, -1, -1, -1, -1, -1 };
  int size_of_best_aps_ids = 0;
  int clip_default[13] = { 0,0,0,0,0,0,0,0,0,0,0,0,0 };
  const int8_t bit_depth = state->encoder_control->bitdepth;
  double lambda = state->frame->lambda;
  int size_of_aps_ids = 0;
  const int32_t num_ctus_in_pic = state->lcu_order_count;
  alf_aps alf_aps_temp_nl;

  //AlfSliceParam  alfSliceParamNewFiltersBest = alfSliceParamNewFilters;
  alf_aps alf_aps_new_filters_best;
  copy_alf_param(&alf_aps_new_filters_best, aps);
  alf_aps* apss = state->slice->apss;

  bool has_new_filters[2] = { aps->enabled_flag[COMPONENT_Y] , aps->enabled_flag[COMPONENT_Cb] || aps->enabled_flag[COMPONENT_Cr] };

  //initDistortion();
  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
    {
      g_ctb_distortion_unfilter[comp][ctu_idx] = get_unfiltered_distortion_cov_classes(g_alf_covariance[comp][ctu_idx], comp == 0 ? MAX_NUM_ALF_CLASSES : 1);
    }
  }

  //luma
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  copy_alf_param(&g_alf_aps_temp, aps);
  memset(g_ctu_enable_flag[COMPONENT_Y], 1, sizeof(uint8_t) * num_ctus_in_pic);
  get_frame_stats(CHANNEL_TYPE_LUMA, num_ctus_in_pic);
  memset(g_ctu_enable_flag[COMPONENT_Y], 0, sizeof(uint8_t) * num_ctus_in_pic);
  
  double cost_off = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[CHANNEL_TYPE_LUMA], CHANNEL_TYPE_LUMA);

  int new_aps_id;
  int aps_ids[ALF_CTB_MAX_NUM_APS];
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    aps_ids[i] = -1;
  }

  kvz_alf_get_avai_aps_ids_luma(state, &new_aps_id, aps_ids, &size_of_aps_ids);

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
    for (int num_temporal_aps = 0; num_temporal_aps <= size_of_aps_ids/*apsIds.size()*/; num_temporal_aps++)
    {
      if (num_temporal_aps + use_new_filter >= ALF_CTB_MAX_NUM_APS)
      {
        continue;
      }
      //cs.slice->setTileGroupNumAps(numTemporalAps + useNewFilter);
      state->slice->tile_group_num_aps = num_temporal_aps + use_new_filter;

      int num_filter_set = ALF_NUM_FIXED_FILTER_SETS + num_temporal_aps + use_new_filter;
      if (num_temporal_aps == size_of_aps_ids && num_temporal_aps > 0 && use_new_filter && new_aps_id == aps_ids[size_of_aps_ids - 1] /*apsIds.back()*/) //last temporalAPS is occupied by new filter set and this temporal APS becomes unavailable
      {
        continue;
      }
      for (int iter = 0; iter < num_iter; iter++)
      {
        //g_alf_aps_temp = aps;
        copy_alf_param(&g_alf_aps_temp, aps);
        g_alf_aps_temp.enabled_flag[CHANNEL_TYPE_LUMA] = true;
        double cur_cost = 3 * lambda;

        if (iter > 0)  //re-derive new filter-set
        {
          double d_dist_org_new_filter = 0;
          int blocks_using_new_filter = 0;
          for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
          {
            if (g_ctu_enable_flag[COMPONENT_Y][ctu_idx] && g_alf_ctb_filter_index[ctu_idx] != ALF_NUM_FIXED_FILTER_SETS)
            {
              g_ctu_enable_flag[COMPONENT_Y][ctu_idx] = 0;
            }
            else if (g_ctu_enable_flag[COMPONENT_Y][ctu_idx] && g_alf_ctb_filter_index[ctu_idx] == ALF_NUM_FIXED_FILTER_SETS)
            {
              blocks_using_new_filter++;
              d_dist_org_new_filter += g_ctb_distortion_unfilter[COMPONENT_Y][ctu_idx];
              for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
              {
                short* p_coeff = g_coeff_final;
                int16_t* p_clipp = g_clipp_final;
                for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
                {
                  g_filter_tmp[i] = p_coeff[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                  g_clip_tmp[i] = p_clipp[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                }
                d_dist_org_new_filter += calc_error_for_coeffs(&g_alf_covariance[COMPONENT_Y][ctu_idx][class_idx], g_clip_tmp, g_filter_tmp, MAX_NUM_ALF_LUMA_COEFF, bit_depth);
              }
            }
          } //for ctb
          if (blocks_using_new_filter > 0 && blocks_using_new_filter < num_ctus_in_pic)
          {
            int bit_nl[2] = { 0, 0 };
            double err_nl[2] = { 0.0, 0.0 };
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
            g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA][0] = 1;
/*#else
            g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] = 1;
#endif*/            
            if (state->encoder_control->cfg.alf_non_linear_luma)
            {
              err_nl[1] = kvz_alf_get_filter_coeff_and_cost(state, CHANNEL_TYPE_LUMA, 0, &bit_nl[1], true, true);
              copy_alf_param(&alf_aps_temp_nl, &g_alf_aps_temp);
            }
            else
            {
              err_nl[1] = MAX_DOUBLE;
            }
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
            g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA][0] = 0;
/*#else
            g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] = 0;
            err_nl[0] = kvz_alf_get_filter_coeff_and_cost(state, CHANNEL_TYPE_LUMA, 0, &bit_nl[0], true, true);

            int bits_new_filter_temp_luma = bit_nl[0];
            double err = err_nl[0];
            if (err_nl[1]  < err_nl[0])
            {
              err = err_nl[1];
              bits_new_filter_temp_luma = bit_nl[1];
              copy_alf_param(&g_alf_aps_temp, &alf_aps_temp_nl);
            }
            if (d_dist_org_new_filter + lambda * g_bits_new_filter[CHANNEL_TYPE_LUMA] < err) //re-derived filter is not good, skip
            {
              continue;
            }
            kvz_alf_reconstruct_coeff(state, &g_alf_aps_temp, CHANNEL_TYPE_LUMA, true, true);
            bits_new_filter = bits_new_filter_temp_luma;
          }
          else //no blocks using new filter, skip
          {
            continue;
          }
        }

        //m_CABACEstimator->getCtx() = ctxStart;
        memcpy(&cabac_estimator, &ctx_start, sizeof(cabac_estimator));
        for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
        {
          double dist_unfilter_ctb = g_ctb_distortion_unfilter[COMPONENT_Y][ctu_idx];
          //ctb on
          g_ctu_enable_flag[COMPONENT_Y][ctu_idx] = 1;
          double cost_on = MAX_DOUBLE;
          //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
          memcpy(&ctx_temp_start, &cabac_estimator, sizeof(ctx_temp_start));
          ctx_temp_start.only_count = 1;
          int i_best_filter_set_idx = 0;
          for (int filter_set_idx = 0; filter_set_idx < num_filter_set; filter_set_idx++)
          {
            //rate
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
            memcpy(&cabac_estimator, &ctx_temp_start, sizeof(cabac_estimator));
            //m_CABACEstimator->resetBits();
            kvz_cabac_reset_bits(&cabac_estimator);
            //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, COMPONENT_Y, &m_alfSliceParamTemp);
            code_alf_ctu_enable_flag(state, &cabac_estimator, ctu_idx, COMPONENT_Y, &g_alf_aps_temp);
            g_alf_ctb_filter_index[ctu_idx] = filter_set_idx;
            code_alf_ctu_filter_index(state, &cabac_estimator, ctu_idx, g_alf_aps_temp.enabled_flag[COMPONENT_Y]);

            double rate_on = (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale * 0; /*(double)m_CABACEstimator->getEstFracBits()*/ ;
            //distortion
            double dist = dist_unfilter_ctb;
            for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
            {
              if (filter_set_idx < ALF_NUM_FIXED_FILTER_SETS)
              {
                int filter_idx = g_class_to_filter_mapping[filter_set_idx][class_idx];
                dist += calc_error_for_coeffs(&g_alf_covariance[COMPONENT_Y][ctu_idx][class_idx], clip_default, g_fixed_filter_set_coeff[filter_idx], MAX_NUM_ALF_LUMA_COEFF, bit_depth);
              }
              else
              {
                short *p_coeff;
                int16_t *p_clipp;
                if (use_new_filter && filter_set_idx == ALF_NUM_FIXED_FILTER_SETS)
                {
                  p_coeff = g_coeff_final;
                  p_clipp = g_clipp_final;
                }
                else if (use_new_filter)
                {
                  p_coeff = g_coeff_aps_luma[filter_set_idx - 1 - ALF_NUM_FIXED_FILTER_SETS];
                  p_clipp = g_clipp_aps_luma[filter_set_idx - 1 - ALF_NUM_FIXED_FILTER_SETS];
                }
                else
                {
                  p_coeff = g_coeff_aps_luma[filter_set_idx - ALF_NUM_FIXED_FILTER_SETS];
                  p_clipp = g_clipp_aps_luma[filter_set_idx - ALF_NUM_FIXED_FILTER_SETS];
                }
                for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
                {
                  g_filter_tmp[i] = p_coeff[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                  g_clip_tmp[i] = p_clipp[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                }
                dist += calc_error_for_coeffs(&g_alf_covariance[COMPONENT_Y][ctu_idx][class_idx], g_clip_tmp, g_filter_tmp, MAX_NUM_ALF_LUMA_COEFF, bit_depth);
              }
            }
            //cost
            double cost_on_tmp = dist + lambda * rate_on;
            if (cost_on_tmp < cost_on)
            {
              //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
              memcpy(&ctx_temp_best, &cabac_estimator, sizeof(ctx_temp_best));
              ctx_temp_best.only_count = 1;
              cost_on = cost_on_tmp;
              i_best_filter_set_idx = filter_set_idx;
            }
          }
          //ctb off
          g_ctu_enable_flag[COMPONENT_Y][ctu_idx] = 0;
          //rate
          //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
          memcpy(&cabac_estimator, &ctx_temp_start, sizeof(cabac_estimator));
          //m_CABACEstimator->resetBits();
          kvz_cabac_reset_bits(&cabac_estimator);
          //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, COMPONENT_Y, &m_alfSliceParamTemp);
          code_alf_ctu_enable_flag(state, &cabac_estimator, ctu_idx, COMPONENT_Y, &g_alf_aps_temp);

          //cost
          double cost_off = dist_unfilter_ctb + lambda * (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3);// frac_bits_scale * 0; /* (double)m_CABACEstimator->getEstFracBits()*/ ;
          if (cost_on < cost_off)
          {
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
            memcpy(&cabac_estimator, &ctx_temp_best, sizeof(cabac_estimator));
            g_ctu_enable_flag[COMPONENT_Y][ctu_idx] = 1;
            g_alf_ctb_filter_index[ctu_idx] = i_best_filter_set_idx;
            cur_cost += cost_on;
          }
          else
          {
            g_ctu_enable_flag[COMPONENT_Y][ctu_idx] = 0;
            cur_cost += cost_off;
          }
        } //for(ctbIdx)

        int tmp_bits = bits_new_filter + 3 * (num_filter_set - ALF_NUM_FIXED_FILTER_SETS);
        cur_cost += tmp_bits * lambda;
        if (cur_cost < cost_min)
        {
          cost_min = cur_cost;
          size_of_best_aps_ids = num_filter_set - ALF_NUM_FIXED_FILTER_SETS;
          for (int i = 0; i < size_of_best_aps_ids; i++)
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
          copy_alf_param(&alf_aps_new_filters_best, &g_alf_aps_temp);

          //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
          memcpy(&ctx_best, &cabac_estimator, sizeof(ctx_best));
          //copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, CHANNEL_TYPE_LUMA);
          memcpy(g_ctu_enable_flag_tmp[COMPONENT_Y], g_ctu_enable_flag[COMPONENT_Y], sizeof(uint8_t) * num_ctus_in_pic);
          for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
          {
            g_alf_ctb_filter_set_index_tmp[ctu_idx] = g_alf_ctb_filter_index[ctu_idx];
          }
          alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA] = use_new_filter;
        }
      }//for (int iter = 0; iter < numIter; iter++)
    }// for (int numTemporalAps = 0; numTemporalAps < apsIds.size(); numTemporalAps++)
  }//for (int useNewFilter = 0; useNewFilter <= 1; useNewFilter++)

  state->slice->tile_group_cc_alf_cb_aps_id = new_aps_id;
  state->slice->tile_group_cc_alf_cr_aps_id = new_aps_id;

  if (cost_off <= cost_min)
  {
    memset(state->slice->tile_group_alf_enabled_flag, 0, sizeof(state->slice->tile_group_alf_enabled_flag));
    state->slice->tile_group_num_aps = 0;
    for (int i = 0; i < MAX_NUM_COMPONENT; i++) {
      memset(g_ctu_enable_flag[i], 0, sizeof(uint8_t) * num_ctus_in_pic);
    }
    return;
  }
  else
  {
    //alfSliceParamNewFiltersBest.tLayer = cs.slice->getTLayer();
    alf_aps_new_filters_best.t_layer = state->slice->id;
    //cs.slice->setTileGroupAlfEnabledFlag(COMPONENT_Y, true);
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Y] = true;
    //cs.slice->setTileGroupNumAps((int)bestApsIds.size());
    
    state->slice->tile_group_num_aps = size_of_best_aps_ids;
    //cs.slice->setAPSs(bestApsIds);
    for (int i = 0; i < size_of_best_aps_ids; i++)
    {
      state->slice->tile_group_luma_aps_id[i] = best_aps_ids[i];
    }

    //copyCtuEnableFlag(m_ctuEnableFlag, m_ctuEnableFlagTmp, CHANNEL_TYPE_LUMA);
    copy_ctu_enable_flag(g_ctu_enable_flag, g_ctu_enable_flag_tmp, CHANNEL_TYPE_LUMA, num_ctus_in_pic);
    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
    {
      g_alf_ctb_filter_index[ctu_idx] = g_alf_ctb_filter_set_index_tmp[ctu_idx];
    }

    if (alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA])
    {
      //APS* newAPS = m_apsMap->getPS((new_aps_id << NUM_APS_TYPE_LEN) + T_ALF_APS);
      alf_aps* new_aps = &state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
      if (new_aps->aps_id < 0 || new_aps->aps_id >= ALF_CTB_MAX_NUM_APS) // new_aps == NULL
      {
        //newAPS = m_apsMap->allocatePS(new_aps_id);
        assert(new_aps_id  + NUM_APS_TYPE_LEN + T_ALF_APS < MAX_NUM_APS); //Invalid PS id
        bool found = false;
        for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++) {
          if (state->encoder_control->cfg.param_set_map[i].parameter_set.aps_id == new_aps_id  + NUM_APS_TYPE_LEN + T_ALF_APS) {
            found = true;
          }
        }
        if (!found) {
          state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
          //state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN+ T_ALF_APS].p_nalu_data = 0;
          //state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set = malloc(sizeof(alf_aps));
          state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.aps_id = new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS;
        }
        copy_alf_param(new_aps, &state->slice->param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set);
        
        new_aps->aps_id = new_aps_id;
        new_aps->aps_type = T_ALF_APS;
      }
      copy_alf_param(new_aps, &alf_aps_new_filters_best);
      new_aps->new_filter_flag[CHANNEL_TYPE_CHROMA] = false;
      state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
      g_aps_id_start = new_aps_id;
    }

    int8_t *aps_ids = state->slice->tile_group_luma_aps_id;
    for (int i = 0; i < state->slice->tile_group_num_aps; i++)
    {
      copy_aps(&apss[aps_ids[i]], &state->encoder_control->cfg.param_set_map[aps_ids[i] + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set);
    }
  }

  //chroma
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  copy_alf_param(&g_alf_aps_temp, &alf_aps_new_filters_best);
  if (g_alf_aps_temp.num_alternatives_chroma < 1)
  {
    g_alf_aps_temp.num_alternatives_chroma = 1;
  }
  //set_ctu_alternative_chroma(m_ctuAlternative, 0);
  //for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++) 
  {
    g_ctu_alternative[COMPONENT_Cb][ctu_idx] = 0;
    g_ctu_alternative[COMPONENT_Cr][ctu_idx] = 0;
  }
//#endif
  //memset(g_ctu_enable_flag[COMPONENT_Cb], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
  //memset(g_ctu_enable_flag[COMPONENT_Cr], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
  set_ctu_enable_flag(g_ctu_enable_flag, CHANNEL_TYPE_CHROMA, ctu_idx, 1);
  get_frame_stats(CHANNEL_TYPE_CHROMA, 0, ctu_idx);
  cost_off = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA][0], CHANNEL_TYPE_CHROMA);
  cost_min = MAX_DOUBLE;
  //m_CABACEstimator->getCtx() = AlfCtx(ctxBest);
  memcpy(&cabac_estimator, &ctx_best, sizeof(cabac_estimator));
  //ctxStart = AlfCtx(m_CABACEstimator->getCtx());
  memcpy(&ctx_start, &cabac_estimator, sizeof(ctx_start));
  ctx_start.only_count = 1;
  int new_aps_id_chroma = -1;
  if (alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA] && (alf_aps_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_aps_new_filters_best.enabled_flag[COMPONENT_Cr]))
  {
    new_aps_id_chroma = new_aps_id;
  }
  else if (alf_aps_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_aps_new_filters_best.enabled_flag[COMPONENT_Cr])
  {
    int cur_id = g_aps_id_start;
    if (size_of_aps_ids < 8 || state->slice->tile_group_num_aps < 8)
    {
      g_alf_aps_temp.num_alternatives_chroma = 1;
    }
    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++) 
    {
      g_ctu_alternative[COMPONENT_Cb][ctu_idx] = 0;
      g_ctu_alternative[COMPONENT_Cr][ctu_idx] = 0;
    }
    set_ctu_enable_flag(g_ctu_enable_flag, CHANNEL_TYPE_CHROMA, 1, num_ctus_in_pic);
    get_frame_stats(CHANNEL_TYPE_CHROMA, num_ctus_in_pic);
    cost_off = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA], CHANNEL_TYPE_CHROMA);
    cost_min = MAX_DOUBLE;
    //m_CABACEstimator->getCtx() = AlfCtx(ctxBest);
    memcpy(&cabac_estimator, &ctx_best, sizeof(cabac_estimator));
    //ctxStart = AlfCtx(m_CABACEstimator->getCtx());
    memcpy(&ctx_start, &cabac_estimator, sizeof(ctx_start));
    ctx_start.only_count = 1;
    int new_aps_id_chroma = -1;
    if (alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA] && (alf_aps_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_aps_new_filters_best.enabled_flag[COMPONENT_Cr]))
    {
      new_aps_id_chroma = new_aps_id;
    }
    else if (alf_aps_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_aps_new_filters_best.enabled_flag[COMPONENT_Cr])
    {
      int cur_id = g_aps_id_start;
      if (size_of_aps_ids < 8 || state->slice->tile_group_num_aps < 8)
      {
        while (new_aps_id_chroma < 0)
        {
          cur_id--;
          if (cur_id < 0)
          {
            cur_id = ALF_CTB_MAX_NUM_APS - 1;
          }

          bool found = false;
          for (int i = 0; i < 8; i++) {
            if (cur_id == best_aps_ids[i]) {
              found = true;
            }
          }
          if (!found)
          {
            new_aps_id_chroma = cur_id;
          }
        }
      }
    }

    for (int cur_aps_id = 0; cur_aps_id < ALF_CTB_MAX_NUM_APS; cur_aps_id++)
    {
      if ((/*(cs.slice->getPendingRasInit() ||*/ (state->frame->pictype == KVZ_NAL_IDR_W_RADL || state->frame->pictype == KVZ_NAL_IDR_N_LP) || (state->frame->slicetype == KVZ_SLICE_I)) && cur_aps_id != new_aps_id_chroma)
      {
        continue;
      }
      //APS* cur_aps = m_apsMap->getPS(cur_aps_id);
      alf_aps* cur_aps = &state->encoder_control->cfg.param_set_map[cur_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
      if (cur_aps && cur_aps->layer_id != 0 /*cs.slice->getPic()->layerId*/)
      {
        continue;
      }
      double cur_cost = lambda * 3;

    if (cur_aps_id == new_aps_id_chroma)
    {
      copy_alf_param(&g_alf_aps_temp, aps);
      cur_cost += g_lambda[CHANNEL_TYPE_CHROMA] * g_bits_new_filter[CHANNEL_TYPE_CHROMA];
    }
    else if (cur_aps && cur_aps->t_layer <= state->slice->id && cur_aps->new_filter_flag[CHANNEL_TYPE_CHROMA])
    {
      //g_alf_slice_aps_temp = cur_aps;
      copy_alf_param(&g_alf_aps_temp, cur_aps);
    }
    else
    {
      continue;
    }
    kvz_alf_reconstruct_coeff(state, &g_alf_aps_temp, CHANNEL_TYPE_CHROMA, true, true);
    //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
    memcpy(&cabac_estimator, &ctx_start, sizeof(cabac_estimator));
    for (int comp_id = 1; comp_id < MAX_NUM_COMPONENT; comp_id++)
    {
      g_alf_aps_temp.enabled_flag[comp_id] = true;
      //for (int ctb_idx = 0; ctb_idx < g_num_ctus_in_pic; ctb_idx++)
      {
        double dist_unfilter_ctu = g_ctb_distortion_unfilter[comp_id][ctu_idx];
        //cost on
        g_ctu_enable_flag[comp_id][ctu_idx] = 1;
        //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
        memcpy(&ctx_temp_start, &cabac_estimator, sizeof(ctx_temp_start));
        ctx_temp_start.only_count = 1;
        //rate
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
        memcpy(&cabac_estimator, &ctx_temp_start, sizeof(cabac_estimator));
        //m_CABACEstimator->resetBits();
        kvz_cabac_reset_bits(&cabac_estimator);
        //ctb flag
        code_alf_ctu_enable_flag(state, &cabac_estimator, ctu_idx, comp_id, &g_alf_aps_temp);
        double rate_on = (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
//#if ENABLE_QPA
        const double ctu_lambda = lambda_chroma_weight > 0.0 ? 0/*cs.picture->m_uEnerHpCtu[ctbIdx]*/ / lambda_chroma_weight : g_lambda[comp_id];
/*#else
        const double ctu_lambda = m_lambda[compId];
#endif*/
        double dist = MAX_DOUBLE;
        int num_alts = g_alf_aps_temp.num_alternatives_chroma;
        //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
        memcpy(&ctx_temp_best, &cabac_estimator, sizeof(ctx_temp_best));
        ctx_temp_best.only_count = 1;
        double best_alt_rate = 0;
        double best_alt_cost = MAX_DOUBLE;
        int best_alt_idx = -1;
        //ctxTempAltStart = AlfCtx(ctxTempBest);
        memcpy(&ctx_temp_alt_start, &ctx_temp_best, sizeof(ctx_temp_alt_start));
        for (int alt_idx = 0; alt_idx < num_alts; ++alt_idx)
        {
          double dist_unfilter_ctu = g_ctb_distortion_unfilter[comp_id][ctu_idx];
          //cost on
          g_ctu_enable_flag[comp_id][ctu_idx] = 1;
          memcpy(&ctx_temp_start, &cabac_estimator, sizeof(ctx_temp_start));
          ctx_temp_start.only_count = 1;
          //rate
          //memcpy(&cabac_estimator, &ctx_temp_start, sizeof(cabac_estimator));
          kvz_cabac_reset_bits(&cabac_estimator);
          //ctb flag
          code_alf_ctu_enable_flag(state, &cabac_estimator, ctu_idx, comp_id, &g_alf_aps_temp);
          double rate_on = (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;
          const double ctu_lambda = lambda_chroma_weight > 0.0 ? 0/*cs.picture->m_uEnerHpCtu[ctbIdx]*/ / lambda_chroma_weight : lambda;
          double dist = MAX_DOUBLE;
          int num_alts = g_alf_aps_temp.num_alternatives_chroma;
          //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
          memcpy(&ctx_temp_best, &cabac_estimator, sizeof(ctx_temp_best));
          ctx_temp_best.only_count = 1;
          double best_alt_rate = 0;
          double best_alt_cost = MAX_DOUBLE;
          int best_alt_idx = -1;
          //ctxTempAltStart = AlfCtx(ctxTempBest);
          memcpy(&ctx_temp_alt_start, &ctx_temp_best, sizeof(ctx_temp_alt_start));
          for (int alt_idx = 0; alt_idx < num_alts; ++alt_idx)
          {
            if (alt_idx) {
              //m_CABACEstimator->getCtx() = AlfCtx(ctxTempAltStart);
              memcpy(&cabac_estimator, &ctx_temp_alt_start, sizeof(cabac_estimator));
            }
            //m_CABACEstimator->resetBits();
            kvz_cabac_reset_bits(&cabac_estimator);
            g_ctu_alternative[comp_id][ctu_idx] = alt_idx;
            //m_CABACEstimator->codeAlfCtuAlternative(cs, ctbIdx, compId, &m_alfParamTemp);
            code_alf_ctu_alternative_ctu(state, &cabac_estimator, ctu_idx, comp_id, &g_alf_aps_temp);
            double alt_rate = (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale * 0/*m_CABACEstimator->getEstFracBits()*/;
            double r_alt_cost = ctu_lambda * alt_rate;

            //distortion
            for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
            {
              g_filter_tmp[i] = g_chroma_coeff_final[alt_idx][i];
              g_clip_tmp[i] = g_chroma_clipp_final[alt_idx][i];
            }
            double alt_dist = calc_error_for_coeffs(&g_alf_covariance[comp_id][ctu_idx][0], g_clip_tmp, g_filter_tmp, MAX_NUM_ALF_CHROMA_COEFF, bit_depth);
            double alt_cost = alt_dist + r_alt_cost;
            if (alt_cost < best_alt_cost)
            {
              best_alt_cost = alt_cost;
              best_alt_idx = alt_idx;
              best_alt_rate = alt_rate;
              //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
              memcpy(&ctx_temp_best, &cabac_estimator, sizeof(ctx_temp_best));
              ctx_temp_best.only_count = 1;
              dist = alt_dist;
            }
          }
          g_ctu_alternative[comp_id][ctu_idx] = best_alt_idx;
          rate_on += best_alt_rate;
          dist += dist_unfilter_ctu;
          //cost
          double cost_on = dist + ctu_lambda * rate_on;
          //cost off
          g_ctu_enable_flag[comp_id][ctu_idx] = 0;
          //rate
          memcpy(&cabac_estimator, &ctx_temp_start, sizeof(cabac_estimator));
          kvz_cabac_reset_bits(&cabac_estimator);
          code_alf_ctu_enable_flag(state, &cabac_estimator, ctu_idx, comp_id, &g_alf_aps_temp);
          //cost
          double cost_off = dist_unfilter_ctu + lambda * (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;
          if (cost_on < cost_off)
          {
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
            memcpy(&cabac_estimator, &ctx_temp_best, sizeof(cabac_estimator));
            g_ctu_enable_flag[comp_id][ctu_idx] = 1;
            cur_cost += cost_on;
          }
        }
        g_ctu_alternative[comp_id][ctu_idx] = best_alt_idx;
        rate_on += best_alt_rate;
        dist += dist_unfilter_ctu;
        //cost
        double cost_on = dist + ctu_lambda * rate_on;

/*#else
        //distortion
        for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
        {
          g_filter_tmp[i] = g_chroma_coeff_final[i];
          g_clip_tmp[i] = g_chroma_clipp_final[i];
        }
        double dist = dist_unfilter_ctu + calc_error_for_coeffs(&g_alf_covariance[comp_id][0][ctb_idx][0], g_clip_tmp, g_filter_tmp, MAX_NUM_ALF_CHROMA_COEFF, ALF_NUM_BITS);
        double cost_on = dist + g_lambda[comp_id] * rate_on;
        //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
        memcpy(&ctx_temp_best, &cabac_estimator, sizeof(ctx_temp_best));
#endif*/

        //cost off
        g_ctu_enable_flag[comp_id][ctu_idx] = 0;
        //rate
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
        memcpy(&cabac_estimator, &ctx_temp_start, sizeof(cabac_estimator));
        //m_CABACEstimator->resetBits();
        kvz_cabac_reset_bits(&cabac_estimator);
        code_alf_ctu_enable_flag(state, &cabac_estimator, ctu_idx, comp_id, &g_alf_aps_temp);
        //cost
        double cost_off = dist_unfilter_ctu + g_lambda[comp_id] * (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3); //frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;
        if (cost_on < cost_off)
        {
          if (g_ctu_enable_flag[comp_id][ctu_idx])
          {
            g_alf_aps_temp.enabled_flag[comp_id] = true;
            break;
          }
        }
      }

      if (cur_cost < cost_min)
      {
        cost_min = cur_cost;
        state->slice->tile_group_chroma_aps_id = cur_aps_id;
        state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = g_alf_aps_temp.enabled_flag[COMPONENT_Cb];
        state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
        copy_ctu_enable_flag(g_ctu_enable_flag_tmp, g_ctu_enable_flag, CHANNEL_TYPE_CHROMA, num_ctus_in_pic);

        for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
        {
          g_ctu_alternative_tmp[COMPONENT_Cb][ctu_idx] = g_ctu_alternative[COMPONENT_Cb][ctu_idx];
          g_ctu_alternative_tmp[COMPONENT_Cr][ctu_idx] = g_ctu_alternative[COMPONENT_Cr][ctu_idx];
        }
      }
    }
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
/*#if !JVET_O0491_HLS_CLEANUP
    curCost += (lengthTruncatedUnary(alfChromaIdc, 3) - lengthTruncatedUnary(0, 3)) * m_lambda[CHANNEL_TYPE_CHROMA];
#endif
#else
    cur_cost += length_truncated_unary(alf_chroma_idc, 3) * g_lambda[CHANNEL_TYPE_CHROMA];
#endif*/
    if (new_aps_id_chroma >= 0)
    {
      cost_min = cur_cost;
      state->slice->tile_group_chroma_aps_id = cur_aps_id;
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = g_alf_aps_temp.enabled_flag[COMPONENT_Cb];
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
      //memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cb], g_ctu_enable_flag[COMPONENT_Cb], sizeof(uint8_t) * g_num_ctus_in_pic);
      //memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cr], g_ctu_enable_flag[COMPONENT_Cr], sizeof(uint8_t) * g_num_ctus_in_pic);
      copy_ctu_enable_flag(g_ctu_enable_flag_tmp, g_ctu_enable_flag, CHANNEL_TYPE_CHROMA, ctu_idx);
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
      //for (int idx = 0; idx < g_num_ctus_in_pic; idx++) 
      {
        g_ctu_alternative_tmp[COMPONENT_Cb][ctu_idx] = g_ctu_alternative[COMPONENT_Cb][ctu_idx];
        g_ctu_alternative_tmp[COMPONENT_Cr][ctu_idx] = g_ctu_alternative[COMPONENT_Cr][ctu_idx];
      }
//#endif
    }
  }
  if (cost_off < cost_min)
  {
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = false;
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = false;
    //memset(g_ctu_enable_flag[COMPONENT_Cb], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
    //memset(g_ctu_enable_flag[COMPONENT_Cr], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
    set_ctu_enable_flag(g_ctu_enable_flag, CHANNEL_TYPE_CHROMA, ctu_idx, 0);
  }
  else
  {
    //memcpy(g_ctu_enable_flag[COMPONENT_Cb], g_ctu_enable_flag_tmp[COMPONENT_Cb], sizeof(uint8_t) * g_num_ctus_in_pic);
    //memcpy(g_ctu_enable_flag[COMPONENT_Cr], g_ctu_enable_flag_tmp[COMPONENT_Cr], sizeof(uint8_t) * g_num_ctus_in_pic);
    copy_ctu_enable_flag(g_ctu_enable_flag, g_ctu_enable_flag_tmp, CHANNEL_TYPE_CHROMA, ctu_idx);
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    //for (int idx = 0; idx < g_num_ctus_in_pic; idx++) 
    {
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = false;
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = false;
      set_ctu_enable_flag(g_ctu_enable_flag, CHANNEL_TYPE_CHROMA, 0, num_ctus_in_pic);
    }
//#endif
    if (state->slice->tile_group_chroma_aps_id == new_aps_id_chroma)  //new filter
    {
      copy_ctu_enable_flag(g_ctu_enable_flag, g_ctu_enable_flag_tmp, CHANNEL_TYPE_CHROMA, num_ctus_in_pic);
      for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
      {
        //newAPS = m_apsMap->allocatePS(new_aps_id);
        assert(new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS < MAX_NUM_APS); //Invalid PS id
        bool found = false;
        for (int i = 0; i < (sizeof(state->slice->param_set_map) / sizeof(state->slice->param_set_map[0])); i++) {
          if (state->slice->param_set_map[i].parameter_set.aps_id == new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS) {
            found = true;
          }
          if (!found) {
            state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
            //state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].p_nalu_data = 0;
            //state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set = malloc(sizeof(alf_aps));
            state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.aps_id = new_aps_id + T_ALF_APS;
          }
          //copy_alf_param(new_aps, &state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set);
          new_aps->aps_id = new_aps_id;
          new_aps->aps_type = T_ALF_APS;
          reset_alf_param(new_aps);
        }
        new_aps->new_filter_flag[CHANNEL_TYPE_CHROMA] = true;
        if (!alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA]) {
          new_aps->new_filter_flag[CHANNEL_TYPE_LUMA] = false;
        }

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
      new_aps->num_alternatives_chroma = aps->num_alternatives_chroma;
      for (int alt_idx = 0; alt_idx < MAX_NUM_ALF_ALTERNATIVES_CHROMA; ++alt_idx)
        new_aps->non_linear_flag[CHANNEL_TYPE_CHROMA][alt_idx] = aps->non_linear_flag[CHANNEL_TYPE_CHROMA][alt_idx];
/*#else
      new_aps->non_linear_flag[CHANNEL_TYPE_CHROMA] = aps->non_linear_flag[CHANNEL_TYPE_CHROMA];
#endif*/
      new_aps->t_layer = state->slice->id;
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
      for (int alt_idx = 0; alt_idx < MAX_NUM_ALF_ALTERNATIVES_CHROMA; ++alt_idx) 
      {
        for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
        {
          for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
          {
            new_aps->chroma_coeff[alt_idx][i] = aps->chroma_coeff[alt_idx][i];
            new_aps->chroma_clipp[alt_idx][i] = aps->chroma_clipp[alt_idx][i];
          }
        }
        state->encoder_control->cfg.param_set_map[new_aps_id_chroma + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
        g_aps_id_start = new_aps_id_chroma;
      }
/*#else
      for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
      {
        new_aps->chroma_coeff[i] = aps->chroma_coeff[i];
        new_aps->chroma_clipp[i] = aps->chroma_clipp[i];
      }
#endif*/
      state->slice->param_set_map[new_aps_id_chroma + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
      g_aps_id_start = new_aps_id_chroma;
    }
  }
}

void kvz_alf_reconstructor(encoder_state_t * const state)
{
  if (!state->slice->tile_group_alf_enabled_flag[COMPONENT_Y])
  {
    return;
  }

  kvz_alf_reconstruct_coeff_aps(state, true, state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] || state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr], false);

  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  bool chroma_scale_x = (chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool chroma_scale_y = (chroma_fmt != KVZ_CSP_420) ? 0 : 1;

  const int alf_vb_luma_ctu_height = LCU_WIDTH;
  const int alf_vb_chma_ctu_height = (LCU_WIDTH >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0));
  const int alf_vb_luma_pos = LCU_WIDTH - ALF_VB_POS_ABOVE_CTUROW_LUMA;
  const int alf_vb_chma_pos = (LCU_WIDTH >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0)) - ALF_VB_POS_ABOVE_CTUROW_CHMA;
  const int luma_height = state->tile->frame->height;
  const int luma_width = state->tile->frame->width;
  const int max_cu_width = LCU_WIDTH;
  const int max_cu_height = LCU_WIDTH;

  int ctu_idx = 0;

  const int luma_stride = state->tile->frame->rec->stride;
  const int chroma_stride = luma_stride >> chroma_scale_x;
  const int chroma_height = luma_height >> chroma_scale_y;
  const int chroma_padding = MAX_ALF_PADDING_SIZE >> chroma_scale_x;

  const int index_luma = -(luma_stride * MAX_ALF_PADDING_SIZE + MAX_ALF_PADDING_SIZE);
  const int index_chroma = -(chroma_stride * chroma_padding + chroma_padding);

  //Copy reconstructed samples to a buffer.
  memcpy(&alf_tmp_y[index_luma], &state->tile->frame->rec->y[index_luma],
    sizeof(kvz_pixel) * luma_stride * (luma_height + MAX_ALF_PADDING_SIZE * 2));
  memcpy(&alf_tmp_u[index_chroma], &state->tile->frame->rec->u[index_chroma],
    sizeof(kvz_pixel) * chroma_stride * (chroma_height + chroma_padding * 2));
  memcpy(&alf_tmp_v[index_chroma], &state->tile->frame->rec->v[index_chroma],
    sizeof(kvz_pixel) * chroma_stride * (chroma_height + chroma_padding * 2));

  for (int y_pos = 0; y_pos < luma_height; y_pos += max_cu_height)
  {
    for (int x_pos = 0; x_pos < luma_width; x_pos += max_cu_width)
    {

      const int width = (x_pos + max_cu_width > luma_width) ? (luma_width - x_pos) : max_cu_width;
      const int height = (y_pos + max_cu_height > luma_height) ? (luma_height - y_pos) : max_cu_height;

      bool ctu_enable_flag = g_ctu_enable_flag[COMPONENT_Y][ctu_idx];
      for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
      {
        ctu_enable_flag |= g_ctu_enable_flag[comp_idx][ctu_idx] > 0;
      }

      {
        if (g_ctu_enable_flag[COMPONENT_Y][ctu_idx])
        {
          short filter_set_index = g_alf_ctb_filter_index[ctu_idx];
          short *coeff;
          int16_t *clip;
          if (filter_set_index >= ALF_NUM_FIXED_FILTER_SETS)
          {
            coeff = g_coeff_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
            clip = g_clipp_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
          }
          else
          {
            coeff = g_fixed_filter_set_coeff_dec[filter_set_index];
            clip = g_clip_default;
          }
          kvz_alf_filter_block(state,
            alf_tmp_y, state->tile->frame->rec->y,
            luma_stride, luma_stride,
            coeff, clip, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y,
            width, height, x_pos, y_pos, x_pos, y_pos,
            alf_vb_luma_pos, alf_vb_luma_ctu_height);
        }
        for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
        {
          coeff = g_fixed_filter_set_coeff_dec[filter_set_index];
          clip = g_clip_default;
        }
        kvz_alf_filter_block(state,
          state->tile->frame->rec->y, alf_tmp_y,
          state->tile->frame->rec->stride, state->tile->frame->rec->stride,
          coeff, clip, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y,
          width, height, x_pos, y_pos, x_pos, y_pos,
          ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_luma_pos),
          g_alf_vb_luma_ctu_height);
      }
      /*else
      {
        int stride = state->tile->frame->rec->stride;
        for (int h = y_pos; h < y_pos + height; h++) {
          for (int w = x_pos; w < x_pos + width; w++) {
            alf_tmp_y[h * stride + w] = state->tile->frame->rec->y[h * stride + w];
          }
        }
      }*/
      for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
      {
        alf_component_id comp_id = comp_idx;

          if (g_ctu_enable_flag[comp_idx][ctu_idx])
          {
            kvz_pixel *dst_pixels = comp_id - 1 ? state->tile->frame->rec->v : state->tile->frame->rec->u;
            const kvz_pixel *src_pixels = comp_id - 1 ? alf_tmp_v : alf_tmp_u;

          const int alt_num = g_ctu_alternative[comp_id][ctu_idx];
          kvz_alf_filter_block(state,
            src_pixels, dst_pixels,
            src_stride, dst_stride, 
            g_chroma_coeff_final[alt_num], g_chroma_clipp_final[alt_num], g_clp_rngs.comp[comp_idx], comp_idx,
            width >> chroma_scale_x, height >> chroma_scale_y,
            x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
            x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
            ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
            g_alf_vb_chma_ctu_height);
/*#else
          kvz_alf_filter_block(state, g_chroma_coeff_final, g_chroma_clipp_final, g_clp_rngs.comp[comp_idx], comp_idx,
            width >> chroma_scale_x, height >> chroma_scale_y,
            x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
            x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
            ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
            g_alf_vb_chma_ctu_height);
#endif*/
        }
        /*else
        {
          int stride = state->tile->frame->rec->stride >> chroma_scale_y;
          int h_start = y_pos >> chroma_scale_y;
          int w_start = x_pos >> chroma_scale_x;
          int c_width = width >> chroma_scale_x;
          int c_height = height >> chroma_scale_y;

          if (comp_idx == COMPONENT_Cb) 
          {
            for (int h = h_start; h < h_start + c_height; h++) {
              for (int w = w_start; w < w_start + c_width; w++) {
                alf_tmp_u[h * stride + w] = state->tile->frame->rec->u[h * stride + w];
              }
            }
          }
          if (comp_idx == COMPONENT_Cr)
          {
            for (int h = h_start; h < h_start + c_height; h++) {
              for (int w = w_start; w < w_start + c_width; w++) {
                alf_tmp_v[h * stride + w] = state->tile->frame->rec->v[h * stride + w];
              }
            }
          }
        }*/
      }
      ctu_idx++;
    }
  }
      //ctu_idx++;
    //}
  //}
//}
}

//----------------------------------------------------------------------

//-------------------------cabac writer functions------------------------

void kvz_cabac_reset_bits(cabac_data_t * const data)
{
  data->low = 0;
  data->bits_left = 23;
  data->num_buffered_bytes = 0;
  data->buffered_byte = 0xff;
}

void code_alf_ctu_enable_flags_channel(encoder_state_t * const state,
  cabac_data_t * const cabac,
  channel_type channel,
  alf_aps *aps)
{
  if (channel == CHANNEL_TYPE_LUMA)
  {
    if (aps->enabled_flag[COMPONENT_Y])
      code_alf_ctu_enable_flags_component(state, cabac, COMPONENT_Y, aps);
  }
  else
  {
    if (aps->enabled_flag[COMPONENT_Cb])
      code_alf_ctu_enable_flags_component(state, cabac, COMPONENT_Cb, aps);
    if (aps->enabled_flag[COMPONENT_Cr])
      code_alf_ctu_enable_flags_component(state, cabac, COMPONENT_Cr, aps);
  }
}

void code_alf_ctu_enable_flags_component(encoder_state_t * const state,
  cabac_data_t * const cabac,
  alf_component_id component_id,
  alf_aps *aps)
{
  const int32_t num_ctus_in_pic = state->lcu_order_count;
  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    code_alf_ctu_enable_flag(state, cabac, ctu_idx, component_id, aps);
  }
}

void code_alf_ctu_enable_flag(encoder_state_t * const state,
  cabac_data_t * const cabac,
  uint32_t ctu_rs_addr,
  alf_component_id component_id,
  alf_aps *aps)
{
  const encoder_control_t * const encoder = state->encoder_control;
  
  const bool alf_component_enabled = (aps != NULL) ? aps->enabled_flag[component_id] : state->slice->tile_group_alf_enabled_flag[component_id];

  if (encoder->cfg.alf_type && alf_component_enabled)
  {
    int frame_width_in_ctus = state->tile->frame->width_in_lcu;

    bool left_avail = state->lcu_order[ctu_rs_addr].left ? 1 : 0;
    bool above_avail = state->lcu_order[ctu_rs_addr].above ? 1 : 0;

    int left_ctu_addr = left_avail ? ctu_rs_addr - 1 : -1;
    int above_ctu_addr = above_avail ? ctu_rs_addr - frame_width_in_ctus : -1;

    uint8_t* ctb_alf_flag = g_ctu_enable_flag[component_id];

    int ctx = 0;
    ctx += left_ctu_addr > -1 ? (ctb_alf_flag[left_ctu_addr] ? 1 : 0) : 0;
    ctx += above_ctu_addr > -1 ? (ctb_alf_flag[above_ctu_addr] ? 1 : 0) : 0;

    cabac->cur_ctx = &(cabac->ctx.alf_ctb_flag_model[component_id * 3 + ctx]);
    CABAC_BIN(cabac, ctb_alf_flag[ctu_rs_addr], "alf_ctb_flag");
  }
}

void code_alf_ctu_filter_index(encoder_state_t * const state,
  cabac_data_t * const cabac,
  uint32_t ctu_rs_addr,
  bool alf_enable_luma)
{
  const encoder_control_t * const encoder = state->encoder_control;

  if (!encoder->cfg.alf_type || !alf_enable_luma)//(!cs.sps->getALFEnabledFlag()) || (!alfEnableLuma))
  {
    return;
  }

  if (!g_ctu_enable_flag[COMPONENT_Y][ctu_rs_addr])
  {
    return;
  }

  const unsigned filter_set_idx = g_alf_ctb_filter_index[ctu_rs_addr];
  unsigned num_aps = state->slice->tile_group_num_aps;
  unsigned num_available_filt_sets = num_aps + ALF_NUM_FIXED_FILTER_SETS;
  if (num_available_filt_sets > ALF_NUM_FIXED_FILTER_SETS)
  {
#if JVET_P0162_REMOVE_ALF_CTB_FIRST_USE_APS_FLAG
    int use_temporal_filt = (filter_set_idx >= ALF_NUM_FIXED_FILTER_SETS) ? 1 : 0;
    cabac->cur_ctx = &(cabac->ctx.alf_temporal_filt);
    CABAC_BIN(cabac, use_temporal_filt, "use_latest_filt");
    if (use_temporal_filt)
    {
      assert(filter_set_idx < num_available_filt_sets); //"temporal non-latest set"
      if (num_aps > 1)
      {
        kvz_cabac_encode_trunc_bin(cabac, filter_set_idx - ALF_NUM_FIXED_FILTER_SETS, num_available_filt_sets - ALF_NUM_FIXED_FILTER_SETS);
      }
    }
    else
    {
      assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //"fixed set larger than temporal"
      kvz_cabac_encode_trunc_bin(cabac, filter_set_idx, ALF_NUM_FIXED_FILTER_SETS);
    }
#else
    int use_latest_filt = (filter_set_idx == ALF_NUM_FIXED_FILTER_SETS) ? 1 : 0;
    /*if (num_aps == 0) {
      use_latest_filt = 1;
    }*/
    cabac->cur_ctx = &(cabac->ctx.alf_latest_filt);
    CABAC_BIN(cabac, use_latest_filt, "use_latest_filt");

    if (!use_latest_filt)
    {

      if (num_aps == 1)
      {
        assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //Fixed set numavail
        kvz_cabac_encode_trunc_bin(cabac, filter_set_idx, ALF_NUM_FIXED_FILTER_SETS);
      }
      else
      {
        int use_temporal_filt = (filter_set_idx > ALF_NUM_FIXED_FILTER_SETS) ? 1 : 0;
        cabac->cur_ctx = &(cabac->ctx.alf_temporal_filt);
        CABAC_BIN(cabac, use_temporal_filt, "use_temporal_filt");

        if (use_temporal_filt)
        {
          assert((filter_set_idx - (ALF_NUM_FIXED_FILTER_SETS + 1)) < (num_aps - 1)); //Temporal non-latest set
          if (num_aps > 2) 
          {
            kvz_cabac_encode_trunc_bin(cabac, filter_set_idx - (ALF_NUM_FIXED_FILTER_SETS + 1), num_available_filt_sets - (ALF_NUM_FIXED_FILTER_SETS + 1));
          }
        }
        else
        {
          assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //Fixed set larger than temporal
          kvz_cabac_encode_trunc_bin(cabac, filter_set_idx, ALF_NUM_FIXED_FILTER_SETS);
        }
      }
    }
#endif
  }
  else
  {
    assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //Fixed set numavail < num_fixed
    kvz_cabac_encode_trunc_bin(cabac, filter_set_idx, ALF_NUM_FIXED_FILTER_SETS);
  }
}

void code_alf_ctu_alternatives_channel(encoder_state_t * const state,
  cabac_data_t * const cabac,
  channel_type channel,
  alf_aps* aps)
{
  if (channel == CHANNEL_TYPE_CHROMA)
  {
    if (aps->enabled_flag[COMPONENT_Cb])
      code_alf_ctu_alternatives_component(state, cabac, COMPONENT_Cb, aps);
    if (aps->enabled_flag[COMPONENT_Cr])
      code_alf_ctu_alternatives_component(state, cabac, COMPONENT_Cr, aps);
  }
}
void code_alf_ctu_alternatives_component(encoder_state_t * const state,
  cabac_data_t * const cabac,
  alf_component_id comp_id,
  alf_aps* aps)
{
  if (comp_id == COMPONENT_Y)
    return;
  uint32_t num_ctus = state->lcu_order_count;
  uint8_t* ctb_alf_flag = g_ctu_enable_flag[comp_id];
  for (int ctu_idx = 0; ctu_idx < num_ctus; ctu_idx++)
  {
    if (ctb_alf_flag[ctu_idx])
    {
      code_alf_ctu_alternative_ctu(state, cabac, ctu_idx, comp_id, aps);
    }
  }
}

void code_alf_ctu_alternative_ctu(encoder_state_t * const state,
  cabac_data_t * const cabac,
  uint32_t ctu_rs_addr,
  const alf_component_id comp_idx,
  const alf_aps* aps)
{

  if (comp_idx == COMPONENT_Y)
    return;
  int aps_idx = aps ? 0 : state->slice->tile_group_chroma_aps_id;
  const alf_aps* alf_param_ref = aps ? (aps) : &state->slice->apss[aps_idx];

  if (aps || (state->encoder_control->cfg.alf_type && state->slice->tile_group_alf_enabled_flag[comp_idx]))
  {
    uint8_t* ctb_alf_flag = g_ctu_enable_flag[comp_idx];

    if (ctb_alf_flag[ctu_rs_addr])
    {
      const int num_alts = alf_param_ref->num_alternatives_chroma;
      uint8_t* ctb_alf_alternative = g_ctu_alternative[comp_idx];
      unsigned num_ones = ctb_alf_alternative[ctu_rs_addr];
      assert(ctb_alf_alternative[ctu_rs_addr] < num_alts);
      for (int i = 0; i < num_ones; ++i) {
        cabac->cur_ctx = &cabac->ctx.alf_ctb_alternatives[comp_idx - 1];
        CABAC_BIN(cabac, 1, "alf_ctb_alternatives");
      }
      if (num_ones < num_alts - 1) {
        cabac->cur_ctx = &cabac->ctx.alf_ctb_alternatives[comp_idx - 1];
        CABAC_BIN(cabac, 0, "alf_ctb_alternatives");
      }
    }
  }
}

void kvz_encode_alf_bits(encoder_state_t * const state, const int ctu_idx)
{
  if (state->encoder_control->cfg.alf_type)
  {
    for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
    {
      bool is_luma = comp_idx == COMPONENT_Y ? true : false;
      //Pit‰isi poistaa//
      /*if (!is_luma)
      {
        state->slice->tile_group_alf_enabled_flag[comp_idx] = false;
      }*/
      //---------------//
      code_alf_ctu_enable_flag(state, &state->cabac, ctu_idx, comp_idx, NULL);
      if (is_luma)
      {
        if (g_ctu_enable_flag[comp_idx][ctu_idx])
        {
          //int num_aps = state->slice->tile_group_num_aps;
          //state->slice->tile_group_num_aps = 0;
          code_alf_ctu_filter_index(state, &state->cabac, ctu_idx, state->slice->tile_group_alf_enabled_flag[COMPONENT_Y]);
          //state->slice->tile_group_num_aps = num_aps;
        }
      }
      if (!is_luma)
      {
        uint8_t* ctb_alf_flag = state->slice->tile_group_alf_enabled_flag[comp_idx] ? g_ctu_enable_flag[comp_idx] : NULL;
        if (ctb_alf_flag && ctb_alf_flag[ctu_idx])
        {
          code_alf_ctu_alternative_ctu(state, &state->cabac, ctu_idx, comp_idx, NULL);
        }
      }
    }

    int num_components = state->encoder_control->chroma_format == KVZ_CSP_400 ? 1 : MAX_NUM_COMPONENT;
    for (int comp_idx = 1; comp_idx < num_components; comp_idx++)
    {
      if (g_cc_alf_filter_param.cc_alf_filter_enabled[comp_idx - 1])
      {
        const int filter_count = g_cc_alf_filter_param.cc_alf_filter_count[comp_idx - 1];

        code_cc_alf_filter_control_idc(state, &state->cabac, g_cc_alf_filter_control[comp_idx - 1][ctu_idx], comp_idx,
          ctu_idx, g_cc_alf_filter_control[comp_idx - 1], filter_count);
      }
    }
  }
}

void encoder_state_write_adaptation_parameter_set(encoder_state_t * const state, alf_aps *aps)
{
#ifdef KVZ_DEBUG
  printf("=========== Adaptation Parameter Set  ===========\n");
#endif

  bitstream_t * const stream = &state->stream;

  WRITE_U(stream, (int)aps->aps_type, 3, "aps_params_type");
  WRITE_U(stream, aps->aps_id, 5, "adaptation_parameter_set_id");
  WRITE_U(stream, state->encoder_control->chroma_format != KVZ_CSP_400, "aps_chroma_present_flag");

  //WRITE_CODE(pcAPS->getAPSType(), 3, "aps_params_type");
  WRITE_U(stream, aps->aps_type, 3, "aps_params_type");

  if (aps->aps_type == T_ALF_APS)
  {
    encode_alf_aps_flags(state, aps);
  }
  /*else if (aps->aps_type == T_LMCS_APS)
  {
    codeLmcsAps(pcAPS);
  }*/
  WRITE_U(stream, 0, 1, "aps_extension_flag"); //Implementation when this flag is equal to 1 should be added when it is needed. Currently in the spec we don't have case when this flag is equal to 1
  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

void encode_alf_aps_flags(encoder_state_t * const state,
  alf_aps* aps)
{
  bitstream_t * const stream = &state->stream;

  WRITE_U(stream, aps->new_filter_flag[CHANNEL_TYPE_LUMA], 1, "alf_luma_new_filter");
  if (state->encoder_control->chroma_format != KVZ_CSP_400)
  {
    WRITE_U(stream, aps->new_filter_flag[CHANNEL_TYPE_CHROMA], 1, "alf_chroma_new_filter")
  }

  if (state->encoder_control->chroma_format != KVZ_CSP_400)
  {
    WRITE_U(stream, aps->cc_alf_aps_param.new_cc_alf_filter[COMPONENT_Cb - 1], 1, "alf_cc_cb_filter_signal_flag");
    WRITE_U(stream, aps->cc_alf_aps_param.new_cc_alf_filter[COMPONENT_Cr - 1], 1, "alf_cc_cr_filter_signal_flag");
  }

  if (aps->new_filter_flag[CHANNEL_TYPE_LUMA])
  {
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    //WRITE_FLAG(param.nonLinearFlag[CHANNEL_TYPE_LUMA][0], "alf_luma_clip");
    WRITE_U(stream, aps->non_linear_flag[CHANNEL_TYPE_LUMA], 1, "alf_luma_clip");
/*#else
    WRITE_FLAG(param.nonLinearFlag[CHANNEL_TYPE_LUMA], "alf_luma_clip");
#endif*/

//#if JVET_O0491_HLS_CLEANUP
    //WRITE_UVLC(param.numLumaFilters - 1, "alf_luma_num_filters_signalled_minus1");
    WRITE_UE(stream, aps->num_luma_filters - 1, "alf_luma_num_filters_signalled_minus1");
/*#else
    xWriteTruncBinCode(param.numLumaFilters - 1, MAX_NUM_ALF_CLASSES);  //number_of_filters_minus1
#endif*/

    if (aps->num_luma_filters > 1)
    {
//#if JVET_O0491_HLS_CLEANUP
      //const int length = ceilLog2(param.numLumaFilters);
      const int length = kvz_math_ceil_log2(aps->num_luma_filters);
//#endif
      for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
      {
//#if JVET_O0491_HLS_CLEANUP
        //WRITE_CODE(param.filterCoeffDeltaIdx[i], length, "alf_luma_coeff_delta_idx");
        WRITE_U(stream, aps->filter_coeff_delta_idx[i], length, "alf_luma_coeff_delta_idx");
/*#else
        xWriteTruncBinCode((uint32_t)param.filterCoeffDeltaIdx[i], param.numLumaFilters);  //filter_coeff_delta[i]
#endif*/
      }
    }
/*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
    WRITE_FLAG(param.fixedFilterSetIndex > 0 ? 1 : 0, "fixed_filter_set_flag");
    if (param.fixedFilterSetIndex > 0)
    {
      xWriteTruncBinCode(param.fixedFilterSetIndex - 1, NUM_FIXED_FILTER_SETS);
      WRITE_FLAG(param.fixedFilterPattern, "fixed_filter_flag_pattern");
      for (int classIdx = 0; classIdx < MAX_NUM_ALF_CLASSES; classIdx++)
      {
        if (param.fixedFilterPattern > 0)
        {
          WRITE_FLAG(param.fixedFilterIdx[classIdx], "fixed_filter_flag");
        }
        else
        {
          CHECK(param.fixedFilterIdx[classIdx] != 1, "Disabled fixed filter");
        }
      }
    }
#endif*/

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    alf_filter(state, aps, false, 0);
/*#else
    alfFilter(param, false);
#endif*/

  }
  if (aps->new_filter_flag[CHANNEL_TYPE_CHROMA])
  {
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
    if (MAX_NUM_ALF_ALTERNATIVES_CHROMA > 1)
    {
      //WRITE_UVLC(param.numAlternativesChroma - 1, "alf_chroma_num_alts_minus1");
      WRITE_UE(stream, aps->num_alternatives_chroma - 1, "alf_chroma_num_alts_minus1");
    }

    for (int alt_idx = 0; alt_idx < aps->num_alternatives_chroma; ++alt_idx)
    {
      encode_alf_aps_filter(state, aps, true, alt_idx);
    }
/*#else
    WRITE_FLAG(param.nonLinearFlag[CHANNEL_TYPE_CHROMA], "alf_chroma_clip");
    alfFilter(param, true);
#endif*/
  }

  for (int cc_idx = 0; cc_idx < 2; cc_idx++)
  {
    if (aps->cc_alf_aps_param.new_cc_alf_filter[cc_idx])
    {
      const int filter_count = aps->cc_alf_aps_param.cc_alf_filter_count[cc_idx];
      assert(filter_count <= MAX_NUM_CC_ALF_FILTERS); // "CC ALF Filter count is too large"
      assert(filter_count > 0); // "CC ALF Filter count is too small"

      if (MAX_NUM_CC_ALF_FILTERS > 1)
      {
        WRITE_UE(stream, filter_count - 1,
          cc_idx == 0 ? "alf_cc_cb_filters_signalled_minus1" : "alf_cc_cr_filters_signalled_minus1");
      }

      for (int filter_idx = 0; filter_idx < filter_count; filter_idx++)
      {
        int num_coeff = 8; //CC_ALF_FILTER

        const short *coeff = aps->cc_alf_aps_param.cc_alf_coeff[cc_idx][filter_idx];
        // Filter coefficients
        for (int i = 0; i < num_coeff - 1; i++)
        {
          if (coeff[i] == 0)
          {
            WRITE_U(stream, 0, CC_ALF_BITS_PER_COEFF_LEVEL,
              cc_idx == 0 ? "alf_cc_cb_mapped_coeff_abs" : "alf_cc_cr_mapped_coeff_abs");
          }
          else
          {
            WRITE_U(stream, 1 + kvz_math_floor_log2(abs(coeff[i])), CC_ALF_BITS_PER_COEFF_LEVEL,
              cc_idx == 0 ? "alf_cc_cb_mapped_coeff_abs" : "alf_cc_cr_mapped_coeff_abs");
            WRITE_U(stream, coeff[i] < 0 ? 1 : 0, 1, cc_idx == 0 ? "alf_cc_cb_coeff_sign" : "alf_cc_cr_coeff_sign");
          }
        }

        /*DTRACE(g_trace_ctx, D_SYNTAX, "%s coeff filter_idx %d: ", cc_idx == 0 ? "Cb" : "Cr", filter_idx);
        for (int i = 0; i < alfShape.numCoeff; i++)
        {
          DTRACE(g_trace_ctx, D_SYNTAX, "%d ", coeff[i]);
        }
        DTRACE(g_trace_ctx, D_SYNTAX, "\n");*/
      }
    }
  }
}

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
void alf_filter(encoder_state_t * const state,
  alf_aps* aps,
  const bool is_chroma,
  const int alt_idx)
/*#else
void HLSWriter::alfFilter(const AlfParam& alfParam, const bool is_chroma)
#endif*/
{
  bitstream_t * const stream = &state->stream;
  const int num_coeff = is_chroma ? 7 : 13;
/*#if !JVET_O0216_ALF_COEFF_EG3 || !JVET_O0064_SIMP_ALF_CLIP_CODING
  static int bitsCoeffScan[EncAdaptiveLoopFilter::m_MAX_SCAN_VAL][EncAdaptiveLoopFilter::m_MAX_EXP_GOLOMB];
  memset(bitsCoeffScan, 0, sizeof(bitsCoeffScan));
  const int maxGolombIdx = AdaptiveLoopFilter::getMaxGolombIdx(alfShape.filterType);
#endif*/
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  const short* coeff = is_chroma ? aps->chroma_coeff[alt_idx] : aps->luma_coeff;
  const int16_t* clipp = is_chroma ? aps->chroma_clipp[alt_idx] : aps->luma_clipp;
/*#else
  const short* coeff = is_chroma ? alfParam.chromaCoeff : alfParam.lumaCoeff;
  const short* clipp = is_chroma ? alfParam.chromaClipp : alfParam.lumaClipp;
#endif*/
  const int num_filters = is_chroma ? 1 : aps->num_luma_filters;

  // vlc for all
/*#if !JVET_O0216_ALF_COEFF_EG3
  for (int ind = 0; ind < num_filters; ++ind)
  {
    if (is_chroma || !alfParam.alfLumaCoeffDeltaFlag || alfParam.alfLumaCoeffFlag[ind])
    {
      for (int i = 0; i < alfShape.numCoeff - 1; i++)
      {
        int coeffVal = abs(coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i]);

        for (int k = 1; k < 15; k++)
        {
          bitsCoeffScan[alfShape.golombIdx[i]][k] += EncAdaptiveLoopFilter::lengthGolomb(coeffVal, k);
        }
      }
    }
  }
#endif*/
/*#if !JVET_O0216_ALF_COEFF_EG3 || !JVET_O0064_SIMP_ALF_CLIP_CODING
  static int kMinTab[MAX_NUM_ALF_COEFF];
#endif*/
/*#if !JVET_O0216_ALF_COEFF_EG3
  int kMin = EncAdaptiveLoopFilter::getGolombKMin(alfShape, num_filters, kMinTab, bitsCoeffScan);
  // Golomb parameters
  WRITE_UVLC(kMin - 1, is_chroma ? "alf_chroma_min_eg_order_minus1" : "alf_luma_min_eg_order_minus1");

  for (int idx = 0; idx < maxGolombIdx; idx++)
  {
    bool golombOrderIncreaseFlag = (kMinTab[idx] != kMin) ? true : false;
    CHECK(!(kMinTab[idx] <= kMin + 1), "ALF Golomb parameter not consistent");
    WRITE_FLAG(golombOrderIncreaseFlag, is_chroma ? "alf_chroma_eg_order_increase_flag" : "alf_luma_eg_order_increase_flag");
    kMin = kMinTab[idx];
  }
#endif*/
  if (!is_chroma)
  {
    if (aps->alf_luma_coeff_delta_flag)
    {
      for (int ind = 0; ind < num_filters; ++ind)
      {
        //WRITE_FLAG(alfParam.alfLumaCoeffFlag[ind], "alf_luma_coeff_flag[i]");
        WRITE_U(stream, aps->alf_luma_coeff_flag[ind], 1, "alf_luma_coeff_flag[i]");
      }
    }
  }

  // Filter coefficients
  for (int ind = 0; ind < num_filters; ++ind)
  {
    if (!is_chroma && !aps->alf_luma_coeff_flag[ind] && aps->alf_luma_coeff_delta_flag)
    {
      continue;
    }

    for (int i = 0; i < num_coeff - 1; i++)
    {
//#if JVET_O0216_ALF_COEFF_EG3
      alf_golomb_encode(state, coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i], 3, true);  // alf_coeff_chroma[i], alf_coeff_luma_delta[i][j]
/*#else
      alfGolombEncode(coeff[ind* MAX_NUM_ALF_LUMA_COEFF + i], kMinTab[alfShape.golombIdx[i]]);  // alf_coeff_chroma[i], alf_coeff_luma_delta[i][j]
#endif*/
    }
  }

  // Clipping values coding
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  if (aps->non_linear_flag[is_chroma][alt_idx])
/*#else
  if (alfParam.nonLinearFlag[is_chroma])
#endif*/
  {
//#if JVET_O0064_SIMP_ALF_CLIP_CODING
    for (int ind = 0; ind < num_filters; ++ind)
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        WRITE_U(stream, clipp[ind * MAX_NUM_ALF_LUMA_COEFF + i], 2, is_chroma ? "alf_chroma_clip_idx" : "alf_luma_clip_idx");
      }
    }
/*#else
    memset(bitsCoeffScan, 0, sizeof(bitsCoeffScan));

    short recCoeff[MAX_NUM_ALF_CLASSES * MAX_NUM_ALF_LUMA_COEFF];
    if (is_chroma)
    {
      memcpy(recCoeff, coeff, sizeof(short) * MAX_NUM_ALF_CHROMA_COEFF);
    }
    else
    {
      memcpy(recCoeff, coeff, sizeof(short) * num_filters * MAX_NUM_ALF_LUMA_COEFF);
#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
      if (alfParam.alfLumaCoeffDeltaPredictionFlag)
      {
        for (int i = 1; i < num_filters; i++)
        {
          for (int j = 0; j < alfShape.numCoeff - 1; j++)
          {
            recCoeff[i * MAX_NUM_ALF_LUMA_COEFF + j] += recCoeff[(i - 1) * MAX_NUM_ALF_LUMA_COEFF + j];
          }
        }
      }
#endif
    }
    // vlc for all
    for (int ind = 0; ind < num_filters; ++ind)
    {
      if (is_chroma || !alfParam.alfLumaCoeffDeltaFlag || alfParam.alfLumaCoeffFlag[ind])
      {
        for (int i = 0; i < alfShape.numCoeff - 1; i++)
        {
          if (!abs(recCoeff[ind * MAX_NUM_ALF_LUMA_COEFF + i]))
            continue;
          int coeffVal = abs(clipp[ind * MAX_NUM_ALF_LUMA_COEFF + i]);

          for (int k = 1; k < 15; k++)
          {
            bitsCoeffScan[alfShape.golombIdx[i]][k] += EncAdaptiveLoopFilter::lengthGolomb(coeffVal, k, false);
          }
        }
      }
    }
#if JVET_O0216_ALF_COEFF_EG3
    int kMin = EncAdaptiveLoopFilter::getGolombKMin(alfShape, num_filters, kMinTab, bitsCoeffScan);
#else
    kMin = EncAdaptiveLoopFilter::getGolombKMin(alfShape, num_filters, kMinTab, bitsCoeffScan);
#endif

    // Golomb parameters
    WRITE_UVLC(kMin - 1, "clip_min_golomb_order");

    for (int idx = 0; idx < maxGolombIdx; idx++)
    {
      bool golombOrderIncreaseFlag = (kMinTab[idx] != kMin) ? true : false;
      CHECK(!(kMinTab[idx] <= kMin + 1), "ALF Golomb parameter not consistent");
      WRITE_FLAG(golombOrderIncreaseFlag, "clip_golomb_order_increase_flag");
      kMin = kMinTab[idx];
    }

    // Filter coefficients
    for (int ind = 0; ind < num_filters; ++ind)
    {
      if (!is_chroma && !alfParam.alfLumaCoeffFlag[ind] && alfParam.alfLumaCoeffDeltaFlag)
      {
        continue;
      }

      for (int i = 0; i < alfShape.numCoeff - 1; i++)
      {
        if (!abs(recCoeff[ind * MAX_NUM_ALF_LUMA_COEFF + i]))
          continue;
        alfGolombEncode(clipp[ind* MAX_NUM_ALF_LUMA_COEFF + i], kMinTab[alfShape.golombIdx[i]], false);  // alf_coeff_chroma[i], alf_coeff_luma_delta[i][j]
      }
    }
#endif*/
  }
}

void encode_alf_adaptive_parameter_set(encoder_state_t * const state)
{
  //send LMCS APS when LMCSModel is updated. It can be updated even current slice does not enable reshaper.
  //For example, in RA, update is on intra slice, but intra slice may not use reshaper
  encode_alf_aps_lmcs(state);

  // only 1 SCALING LIST data for 1 picture
  encode_alf_aps_scaling_list(state);

  encode_alf_aps(state);
}


  // set CTU enable flags
  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++) {
      g_ctu_enable_flag[comp_idx][ctu_idx] = g_alf_ctu_enable_flag[comp_idx][ctu_idx];
      //#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
      g_ctu_alternative[comp_idx][ctu_idx] = g_alf_ctu_alternative[comp_idx][ctu_idx];
      //#endif
    }
  }
}

// only 1 SCALING LIST data for 1 picture
void encode_alf_aps_scaling_list(encoder_state_t * const state)
{
  if (0/*pcSlice->getSPS()->getScalingListFlag() && (m_pcCfg->getUseScalingListId() == SCALING_LIST_FILE_READ)*/)
  {/*
    int apsId = picHeader->getScalingListAPSId();
    ParameterSetMap<APS> *apsMap = m_pcEncLib->getApsMap();
    APS* aps = apsMap->getPS((apsId << NUM_APS_TYPE_LEN) + SCALING_LIST_APS);
    bool writeAPS = aps && apsMap->getChangedFlag((apsId << NUM_APS_TYPE_LEN) + SCALING_LIST_APS);
    if (writeAPS)
    {
#if JVET_R0433
      aps->chromaPresentFlag = pcSlice->getSPS()->getChromaFormatIdc() != CHROMA_400;
#endif
      actualTotalBits += xWriteAPS(accessUnit, aps, m_pcEncLib->getLayerId(), true);
      apsMap->clearChangedFlag((apsId << NUM_APS_TYPE_LEN) + SCALING_LIST_APS);
      CHECK(aps != picHeader->getScalingListAPS(), "Wrong SCALING LIST APS pointer in compressGOP");
    }*/
  }
}

void encode_alf_aps(encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
  bitstream_t * const stream = &state->stream;
  if (encoder->cfg.alf_type) // && (state->slice->tile_group_alf_enabled_flag[COMPONENT_Y] || state->slice->tile_group_cc_alf_cb_enabled_flag || state->slice->tile_group_cc_alf_cr_enabled_flag))
  {
    param_set_map *aps_map = state->encoder_control->cfg.param_set_map;
    for (int aps_id = 0; aps_id < ALF_CTB_MAX_NUM_APS; aps_id++)
    {
      alf_aps aps = aps_map[aps_id + T_ALF_APS + NUM_APS_TYPE_LEN].parameter_set;
      bool write_aps = aps_map[aps_id + T_ALF_APS + NUM_APS_TYPE_LEN].b_changed;
      /*if (!write_aps && state->slice->apss && state->slice->apss[aps_id].aps_id >= 0 && state->slice->apss[aps_id].aps_id < 8)
      {
        write_aps = true;
        aps = state->slice->apss[aps_id]; // use aps from slice header

        //*apsMap->allocatePS(apsId) = *aps; //allocate and cpy
        copy_aps_to_map(aps_map, &aps, aps_id + T_ALF_APS + NUM_APS_TYPE_LEN);
        //m_pcALF->setApsIdStart(apsId);
        g_aps_id_start = aps_id;
      }*/

      //else if
      /*if (state->slice->tile_group_cc_alf_cb_enabled_flag && !write_aps && aps_id == state->slice->tile_group_cc_alf_cb_aps_id)
      {
        write_aps = true;
        aps = aps_map[(state->slice->tile_group_cc_alf_cb_aps_id << NUM_APS_TYPE_LEN) + T_ALF_APS].parameter_set;
      }

      if(ctuEnableFlag && is_crossed_by_virtual_boundaries(x_pos, y_pos, width, height, &clip_top, &clip_bottom, &clip_left, &clip_right, &num_hor_vir_bndry, &num_ver_vir_bndry, hor_vir_bndry_pos, ver_vir_bndry_pos, state))
      {
        write_aps = true;
        aps = aps_map[(state->slice->tile_group_cc_alf_cr_aps_id << NUM_APS_TYPE_LEN) + T_ALF_APS].parameter_set;
      }*/
          const bool clipT = (i == 0 && clip_top) || (i > 0) || (y_start == 0);
          const bool clipB = (i == num_hor_vir_bndry && clip_bottom) || (i < num_hor_vir_bndry) || (y_end == luma_height);

          int x_start = x_pos;
          for (int j = 0; j <= num_ver_vir_bndry; j++)
          {
            const int x_end = j == num_ver_vir_bndry ? x_pos + width : ver_vir_bndry_pos[j];
            const int w = x_end - x_start;
            const bool clipL = (j == 0 && clip_left) || (j > 0) || (x_start == 0);
            const bool clipR = (j == num_ver_vir_bndry && clip_right) || (j < num_ver_vir_bndry) || (x_end == luma_width);

            const int wBuf = w + (clipL ? 0 : MAX_ALF_PADDING_SIZE) + (clipR ? 0 : MAX_ALF_PADDING_SIZE);
            const int hBuf = h + (clipT ? 0 : MAX_ALF_PADDING_SIZE) + (clipB ? 0 : MAX_ALF_PADDING_SIZE);
            /*
            PelUnitBuf buf = m_tempBuf2.subBuf(UnitArea(cs.area.chromaFormat, Area(0, 0, w_buf, h_buf)));
            buf.copyFrom(tmpYuv.subBuf(UnitArea(cs.area.chromaFormat, Area(x_start - (clip_l ? 0 : MAX_ALF_PADDING_SIZE), y_start - (clip_t ? 0 : MAX_ALF_PADDING_SIZE), w_buf, h_buf))));
            buf.extendBorderPel(MAX_ALF_PADDING_SIZE);
            buf = buf.subBuf(UnitArea(cs.area.chromaFormat, Area(clip_l ? 0 : MAX_ALF_PADDING_SIZE, clip_t ? 0 : MAX_ALF_PADDING_SIZE, w, h)));
            */
            if (g_ctu_enable_flag[COMPONENT_Y][ctu_idx])
            {
              //const Area blkSrc(0, 0, w, h);
              //const Area blkDst(x_start, y_start, w, h);
              //deriveClassification(m_classifier, buf.get(COMPONENT_Y), blkDst, blkSrc);
              kvz_alf_derive_classification(state, w, h, x_start, y_start, x_start, y_start);
              //const Area blkPCM(x_start, y_start, w, h);
              //#if !JVET_O0525_REMOVE_PCM
              //resetPCMBlkClassInfo(cs, m_classifier, buf.get(COMPONENT_Y), blkPCM);
              //kvz_alf_reset_pcm_blk_class_info(state, lcu, w, h, x_start, y_start);
              short filter_set_index = g_alf_ctb_filter_index[ctu_idx];
              short *coeff;
              int16_t *clip;
              if (filter_set_index >= ALF_NUM_FIXED_FILTER_SETS)
              {
                coeff = g_coeff_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
                clip = g_clipp_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
              }
              else
              {
                coeff = g_fixed_filter_set_coeff_dec[filter_set_index];
                clip = g_clip_default;
              }
              kvz_alf_filter_block(state,
                state->tile->frame->rec->y, alf_tmp_y,
                state->tile->frame->rec->stride, state->tile->frame->rec->stride,
                coeff, clip, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y,
                w, h, x_start, y_start, x_start, y_start, 
                ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_luma_pos),
                g_alf_vb_luma_ctu_height);
            }

            for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
            {
              alf_component_id comp_id = comp_idx;

              if (g_ctu_enable_flag[comp_idx][ctu_idx])
              {
                //const Area blkSrc(0, 0, w >> chromaScaleX, h >> chromaScaleY);
                //const Area blkDst(x_start >> chromaScaleX, y_start >> chromaScaleY, w >> chromaScaleX, h >> chromaScaleY);

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB

                const kvz_pixel *src_pixels = comp_id - 1 ? state->tile->frame->rec->v : state->tile->frame->rec->u;
                kvz_pixel *dst_pixels = comp_id - 1 ? alf_tmp_v : alf_tmp_u;
                const int src_stride = state->tile->frame->rec->stride >> 1;
                const int dst_stride = state->tile->frame->rec->stride >> 1;

                uint8_t alt_num = g_ctu_alternative[comp_idx][ctu_idx];
                kvz_alf_filter_block(state,
                  src_pixels, dst_pixels,
                  src_stride, dst_stride,
                  g_chroma_coeff_final[alt_num], g_chroma_clipp_final[alt_num], g_clp_rngs.comp[comp_idx],
                  comp_id, w >> chroma_scale_x, h >> chroma_scale_y,
                  x_start >> chroma_scale_x, y_start >> chroma_scale_y,
                  x_start >> chroma_scale_x, y_start >> chroma_scale_y,
                  ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
                  g_alf_vb_chma_ctu_height); 
/*#else
                kvz_alf_filter_block(state, g_chroma_coeff_final, g_chroma_clipp_final, g_clp_rngs.comp[comp_idx],
                  comp_id, w >> chroma_scale_x, h >> chroma_scale_y,
                  x_start >> chroma_scale_x, y_start >> chroma_scale_y,
                  x_start >> chroma_scale_x, y_start >> chroma_scale_y,
                  ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
                  g_alf_vb_chma_ctu_height);
#endif*/
              }
            }

            x_start = x_end;
          }

          y_start = y_end;
        }
      }
      else
      {
        //actualTotalBits += xWriteAPS(accessUnit, aps);
        kvz_nal_write(stream, NAL_UNIT_PREFIX_APS, 0, state->frame->first_nal);
        state->frame->first_nal = false;
        encoder_state_write_adaptation_parameter_set(state, &aps);
          //deriveClassification(m_classifier, tmpYuv.get(COMPONENT_Y), blk, blk);
        //apsMap->clearChangedFlag((apsId << NUM_APS_TYPE_LEN) + ALF_APS);
        aps_map[aps_id + T_ALF_APS].b_changed = false;
          //Area blkPCM(x_pos, y_pos, width, height);
          //#if !JVET_O0525_REMOVE_PCM
          //resetPCMBlkClassInfo(cs, m_classifier, tmpYuv.get(COMPONENT_Y), blkPCM);
          //kvz_alf_reset_pcm_blk_class_info(state, lcu, width, height, x_pos, y_pos);
          short filter_set_index = g_alf_ctb_filter_index[ctu_idx];
          short *coeff;
          int16_t *clip;
          if (filter_set_index >= ALF_NUM_FIXED_FILTER_SETS)
          {
            coeff = g_coeff_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
            clip = g_clipp_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
          }
          else
          {
            coeff = g_fixed_filter_set_coeff_dec[filter_set_index];
            clip = g_clip_default;
          }
          kvz_alf_filter_block(state,
            state->tile->frame->rec->y, alf_tmp_y,
            state->tile->frame->rec->stride, state->tile->frame->rec->stride,
            coeff, clip, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y,
            width, height, x_pos, y_pos, x_pos, y_pos,
            ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_luma_pos), 
            g_alf_vb_luma_ctu_height);
        }
        for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
        {
          alf_component_id comp_id = comp_idx;

          if (g_ctu_enable_flag[comp_idx][ctu_idx])
          {
            //Area blk(x_pos >> chroma_scale_x, y_pos >> chroma_scale_y, width >> chroma_scale_x, height >> chroma_scale_y);
            //m_filter5x5Blk(m_classifier, recYuv, tmpYuv, blk, comp_id, m_chromaCoeffFinal, clp_rngs.comp[comp_idx], cs);
//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB

            const kvz_pixel *src_pixels = comp_id - 1 ? state->tile->frame->rec->v : state->tile->frame->rec->u;
            kvz_pixel *dst_pixels = comp_id - 1 ? alf_tmp_v : alf_tmp_u;
            const int src_stride = state->tile->frame->rec->stride >> 1;
            const int dst_stride = state->tile->frame->rec->stride >> 1;

            uint8_t alt_num = g_ctu_alternative[comp_idx][ctu_idx];
            kvz_alf_filter_block(state,
              src_pixels, dst_pixels,
              src_stride, dst_stride, 
              g_chroma_coeff_final[alt_num], g_chroma_clipp_final[alt_num], g_clp_rngs.comp[comp_idx], comp_idx,
              width >> chroma_scale_x, height >> chroma_scale_y,
              x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
              x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
              ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
              g_alf_vb_chma_ctu_height); 
/*#else
            kvz_alf_filter_block(state, g_chroma_coeff_final, g_chroma_clipp_final, g_clp_rngs.comp[comp_idx], comp_idx,
              width >> chroma_scale_x, height >> chroma_scale_y,
              x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
              x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
              ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
              g_alf_vb_chma_ctu_height); 
#endif*/


          }
        }
      }
    }
  }
}

//------------------------- CC ALF cabac writer functions------------------------

void code_cc_alf_filter_control_idc(encoder_state_t * const state,
  cabac_data_t * const cabac, uint8_t idc_val, 
  const alf_component_id comp_id, const int ctu_idx, 
  const uint8_t *filter_control_idc,
  const int filter_count)
{
  assert(!(idc_val > filter_count)); //Filter index is too large

  bool left_avail = state->lcu_order[ctu_idx].left ? 1 : 0;
  bool above_avail = state->lcu_order[ctu_idx].above ? 1 : 0;
  int ctxt = 0;

  if (left_avail)
  {
    ctxt += (filter_control_idc[ctu_idx - 1]) ? 1 : 0;
  }
  if (above_avail)
  {
    ctxt += (filter_control_idc[ctu_idx - state->tile->frame->width_in_lcu]) ? 1 : 0;
  }
  ctxt += (comp_id == COMPONENT_Cr) ? 3 : 0;

  //m_BinEncoder.encodeBin((idc_val == 0) ? 0 : 1, Ctx::CcAlfFilterControlFlag(ctxt)); // ON/OFF flag is context coded
  cabac->cur_ctx = &(cabac->ctx.alf_cc_filter_control_flag[ctxt]);
  CABAC_BIN(cabac, (idc_val == 0) ? 0 : 1, "cc_alf_filter_control_flag");

  if (idc_val > 0)
  {
    int val = (idc_val - 1);
    while (val)
    {
      //m_BinEncoder.encodeBinEP(1);
      CABAC_BIN_EP(cabac, 1, "cc_alf_filter_control_flag");
      val--;
    }
    if (idc_val < filter_count)
    {
      //m_BinEncoder.encodeBinEP(0);
      CABAC_BIN_EP(cabac, 0, "cc_alf_filter_control_flag");
    }
  }
  //DTRACE(g_trace_ctx, D_SYNTAX, "ccAlfFilterControlIdc() comp_id=%d pos=(%d,%d) ctxt=%d, filter_count=%d, idc_val=%d\n", comp_id, lumaPos.x, lumaPos.y, ctxt, filter_count, idc_val);
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

//-------------------------CTU functions--------------------------------

void kvz_alf_reconstruct_coeff_aps(encoder_state_t * const state, bool luma, bool chroma, bool is_rdo)
{
  //luma
  alf_aps* apss = state->slice->apss;
  //AlfSliceParam alfSliceParamTmp;
  alf_aps alf_param_tmp;
  //APS* cur_aps;
  alf_aps* cur_aps;

  if (luma)
  {
    for (int i = 0; i < state->slice->tile_group_num_aps /* 1,  cs.slice->getTileGroupNumAps()*/; i++) {
      int aps_idx = state->slice->tile_group_luma_aps_id[i];
      cur_aps = &apss[aps_idx];

      assert(cur_aps != NULL); // "invalid APS"
      alf_param_tmp = *cur_aps;

      kvz_alf_reconstruct_coeff(state, &alf_param_tmp, CHANNEL_TYPE_LUMA, is_rdo, true);
      memcpy(g_coeff_aps_luma[i], g_coeff_final, sizeof(g_coeff_final));
      memcpy(g_clipp_aps_luma[i], g_clipp_final, sizeof(g_clipp_final));
    }
  }

  //chroma
  if (chroma)
  {
    int aps_idx_chroma = state->slice->tile_group_chroma_aps_id;
    cur_aps = &apss[aps_idx_chroma];
    //copy_alf_param(g_alf_aps_chroma, cur_aps);
    //copy_alf_param(&alf_param_tmp, g_alf_aps_chroma);
    copy_alf_param(&alf_param_tmp, cur_aps);
    kvz_alf_reconstruct_coeff(state, &alf_param_tmp, CHANNEL_TYPE_CHROMA, is_rdo, true);
  }
}

//void reconstructCoeff(AlfSliceParam& alfSliceParam, ChannelType channel, const bool isRdo, const bool isRedo)
void kvz_alf_reconstruct_coeff(encoder_state_t * const state,
  alf_aps *aps,
  channel_type channel,
  const bool is_rdo,
  const bool is_redo)
{
  const int8_t bit_depth = state->encoder_control->bitdepth;
  int factor = is_rdo ? 0 : (1 << (bit_depth - 1));
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  alf_filter_type filter_type = is_luma ? ALF_FILTER_7X7 : ALF_FILTER_5X5;
  int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
  int num_coeff = filter_type == ALF_FILTER_5X5 ? 7 : 13;
  int num_coeff_minus1 = num_coeff - 1;
/*#if !JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  int num_filters = is_luma ? num_luma_filters : 1;
  short* coeff = is_luma ? luma_coeff : chroma_coeff;
  short* clipp = is_luma ? luma_clipp : chroma_clipp;
#endif*/

//#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB
  const int num_alts = is_luma ? 1 : aps->num_alternatives_chroma;

  for (int alt_idx = 0; alt_idx < num_alts; ++alt_idx)
  {
    int num_filters = is_luma ? aps->num_luma_filters : 1;
    short* coeff = is_luma ? aps->luma_coeff : aps->chroma_coeff[alt_idx];
    int16_t* clipp = is_luma ? aps->luma_clipp : aps->chroma_clipp[alt_idx];

/*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
    if (alfParam.alfLumaCoeffDeltaPredictionFlag && isLuma(channel))
    {
      for (int i = 1; i < num_filters; i++)
      {
        for (int j = 0; j < numCoeffMinus1; j++)
        {
          coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] += coeff[(i - 1) * MAX_NUM_ALF_LUMA_COEFF + j];
        }
      }
    }
#endif*/
    for (int filter_idx = 0; filter_idx < num_filters; filter_idx++)
    {
      coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
    }

    if (!is_luma)
    {
      for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
      {
        g_chroma_coeff_final[alt_idx][coeff_idx] = coeff[coeff_idx];
        int clip_idx = aps->non_linear_flag[channel] ? clipp[coeff_idx] : 0;
        g_chroma_clipp_final[alt_idx][coeff_idx] = is_rdo ? clip_idx : g_alf_clipping_values[channel][clip_idx];
      }
      g_chroma_coeff_final[alt_idx][num_coeff_minus1] = factor;
      g_chroma_clipp_final[alt_idx][num_coeff_minus1] = is_rdo ? 0 : g_alf_clipping_values[channel][0];
      continue;
    }
    for (int class_idx = 0; class_idx < num_classes; class_idx++)
    {
      int filterIdx = aps->filter_coeff_delta_idx[class_idx];
/*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
      int fixedFilterIdx = alfParam.fixedFilterSetIndex;
      if (fixedFilterIdx > 0 && alfParam.fixedFilterIdx[class_idx] > 0)
      {
        fixedFilterIdx = m_classToFilterMapping[fixedFilterIdx - 1][class_idx];
      }
      else
      {
        fixedFilterIdx = -1;
      }
#endif*/
      for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
      {
        g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] = coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx];
/*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
        //fixed filter
        if (fixedFilterIdx >= 0)
        {
          m_coeffFinal[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] += m_fixedFilterSetCoeff[fixedFilterIdx][coeff_idx];
        }
#endif*/
      }
      g_coeff_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
      g_clipp_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = is_rdo ? 0 : g_alf_clipping_values[channel][0];
      for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
      {
        int clip_idx = aps->non_linear_flag[channel] ? clipp[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] : 0;
        assert((clip_idx >= 0 && clip_idx < MAX_ALF_NUM_CLIPPING_VALUES)); // "Bad clip idx in ALF"
        g_clipp_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] = is_rdo ? clip_idx : g_alf_clipping_values[channel][clip_idx];
      }
      g_clipp_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] =
        is_rdo ? 0 :
        g_alf_clipping_values[channel][0];
    }
  }
/*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED

  if (is_chroma(channel))
    return;

  if (isRedo && alfParam.alfLumaCoeffDeltaPredictionFlag)
  {
    int num_filters = alfParam.numLumaFilters;
    short* coeff = alfParam.lumaCoeff;

    for (int i = num_filters - 1; i > 0; i--)
    {
      for (int j = 0; j < numCoeffMinus1; j++)
      {
        coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] = coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] - coeff[(i - 1) * MAX_NUM_ALF_LUMA_COEFF + j];
      }
    }
  }
#endif*/
/*#else

  /*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  if (*alf_luma_coeff_delta_prediction_flag && is_luma)
  {
    for (int i = 1; i < num_filters; i++)
    {
      for (int j = 0; j < num_coeff_minus1; j++)
      {
        coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] += coeff[(i - 1) * MAX_NUM_ALF_LUMA_COEFF + j];
      }
    }
  }*//*

  for (int filter_idx = 0; filter_idx < num_filters; filter_idx++)
  {
    coeff[filter_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
  }

  if ( !is_luma )
  {
    for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
    {
      g_chroma_coeff_final[coeff_idx] = chroma_coeff[coeff_idx];

      g_chroma_coeff_final[coeff_idx] = chroma_coeff[coeff_idx];
      int clip_idx = aps->non_linear_flag[channel] ? clipp[coeff_idx] : 0;
      g_chroma_clipp_final[coeff_idx] = is_rdo ? clip_idx : g_alf_clipping_values[channel][clip_idx];

    }
    g_chroma_coeff_final[num_coeff_minus1] = factor;
    g_chroma_clipp_final[num_coeff_minus1] = is_rdo ? 0 : g_alf_clipping_values[channel][0];

    return;
  }

  for (int class_idx = 0; class_idx < num_classes; class_idx++)
  {
    int filter_idx = filter_coeff_delta_idx[class_idx];
    /*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
    int fixed_filter_idx = *fixed_filter_set_index; //13
    if (fixed_filter_idx > 0 && aps->fixed_filter_idx[class_idx] > 0)
    {
      fixed_filter_idx = g_class_to_filter_mapping[fixed_filter_idx - 1][class_idx];
    }
    else
    {
      fixed_filter_idx = -1;
    }*//*
    for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
    {
      g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] = coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx];
      /*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
      //fixed filter
      if (fixed_filter_idx >= 0)
      {
        g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] += g_fixed_filter_set_coeff[fixed_filter_idx][coeff_idx];
      }*//*
    }
    g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
    g_clipp_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = is_rdo ? 0 : g_alf_clipping_values[channel][0];
    for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
    {
      int clip_idx = aps->non_linear_flag[channel] ? (clipp + filter_idx * MAX_NUM_ALF_LUMA_COEFF)[coeff_idx] : 0;
      g_clipp_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] = is_rdo ? clip_idx : g_alf_clipping_values[channel][clip_idx];
    }
  }
  /*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  if (is_redo && state->cabac.ctx.alf_luma_coeff_delta_prediction_flag.state[0])
  {
    for (int i = num_filters - 1; i > 0; i--)
    {
      for (int j = 0; j < num_coeff_minus1; j++)
      {
        coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] = coeff[i * MAX_NUM_ALF_LUMA_COEFF + j] - coeff[(i - 1) * MAX_NUM_ALF_LUMA_COEFF + j];
      }
    }
  }*/
}

void kvz_alf_create(encoder_state_t const *state)
{
  const int pic_width = state->tile->frame->width;
  const int pic_height = state->tile->frame->height;
  const int max_cu_width = LCU_WIDTH; //128
  const int max_cu_height = LCU_WIDTH; //128
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;

  const int num_ctus_in_width = (pic_width / max_cu_width) + ((pic_width % max_cu_width) ? 1 : 0);
  const int num_ctus_in_height = (pic_height / max_cu_height) + ((pic_height % max_cu_height) ? 1 : 0);
  g_num_ctus_in_pic = num_ctus_in_width * num_ctus_in_height;

  g_alf_vb_luma_pos = max_cu_height - ALF_VB_POS_ABOVE_CTUROW_LUMA;
  g_alf_vb_chma_pos = (max_cu_height >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0)) - ALF_VB_POS_ABOVE_CTUROW_CHMA;

  g_alf_vb_luma_ctu_height = max_cu_height;
  g_alf_vb_chma_ctu_height = (max_cu_height >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0));

  assert(g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] > 0); //"g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] must be at least one"
  g_alf_clipping_values[CHANNEL_TYPE_LUMA][0] = 1 << g_input_bit_depth[CHANNEL_TYPE_LUMA];
  int shift_luma = g_input_bit_depth[CHANNEL_TYPE_LUMA] - 8;
  for (int i = 1; i < g_alf_num_clipping_values[CHANNEL_TYPE_LUMA]; ++i)
  {
    g_alf_clipping_values[CHANNEL_TYPE_LUMA][i] =
      (short)round(pow(2., g_input_bit_depth[CHANNEL_TYPE_LUMA] *
      (g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] - i) / g_alf_num_clipping_values[CHANNEL_TYPE_LUMA]));
  }

  assert(g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] > 0); //"g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] must be at least one"
  g_alf_clipping_values[CHANNEL_TYPE_CHROMA][0] = 1 << g_input_bit_depth[CHANNEL_TYPE_CHROMA];
  int shift_chroma = g_input_bit_depth[CHANNEL_TYPE_CHROMA] - 8;
  for (int i = 1; i < g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA]; ++i)
  {
    g_alf_clipping_values[CHANNEL_TYPE_CHROMA][i] =
      (short)round(pow(2., g_input_bit_depth[CHANNEL_TYPE_CHROMA] - 8 + 8. *
      (g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] - i - 1) / (g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] - 1)));
  }

  if (g_created)
  {
    return;
  }

  // Classification
  g_classifier = malloc(pic_height * sizeof(**g_classifier));
  g_classifier[0] = malloc(pic_height * pic_width * sizeof(*g_classifier));
  for (int i = 1; i < pic_height; i++)
  {
    g_classifier[i] = g_classifier[0] + i * pic_width;
  }

  for (int filter_set_index = 0; filter_set_index < ALF_NUM_FIXED_FILTER_SETS; filter_set_index++)
  {
    for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
    {
      int fixed_filter_idx = g_class_to_filter_mapping[filter_set_index][class_idx];
      for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF - 1; i++)
      {
        g_fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + i] = g_fixed_filter_set_coeff[fixed_filter_idx][i];
      }
      g_fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + MAX_NUM_ALF_LUMA_COEFF - 1] = (1 << (kvz_bit_depth - 1));
    }
  }

  for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF * MAX_NUM_ALF_CLASSES; i++)
  {
    g_clip_default[i] = g_alf_clipping_values[CHANNEL_TYPE_LUMA][0];
  }

  g_created = true;
  g_cc_alf_filter_control[0] = malloc(g_num_ctus_in_pic * sizeof(*g_cc_alf_filter_control));
  g_cc_alf_filter_control[1] = malloc(g_num_ctus_in_pic * sizeof(*g_cc_alf_filter_control));
}

void kvz_alf_destroy(videoframe_t * const frame)
{
  if (!g_created)
  {
    return;
  }

  if (g_classifier)
  {
    FREE_POINTER(g_classifier[0]);
    FREE_POINTER(g_classifier);
  }

  g_created = false;

  if (g_cc_alf_filter_control[0])
  {
    FREE_POINTER(g_cc_alf_filter_control[0])
  }

  if (g_cc_alf_filter_control[1])
  {
    FREE_POINTER(g_cc_alf_filter_control[1])
  }
}

void kvz_alf_derive_classification(encoder_state_t * const state,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  const int blk_dst_x,
  const int blk_dst_y)//,
  //alf_classifier** g_classifier)
{
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  bool chroma_scale_x = (chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool chroma_scale_y = (chroma_fmt != KVZ_CSP_420) ? 0 : 1;

  const int alf_vb_luma_ctu_height = LCU_WIDTH;
  const int alf_vb_luma_pos = LCU_WIDTH - ALF_VB_POS_ABOVE_CTUROW_LUMA;
  int32_t pic_height = state->tile->frame->rec->height;
  int32_t pic_width = state->tile->frame->rec->width;

  int max_height = y_pos + height;
  int max_width = x_pos + width;

  //Use if adjacent CTUs are not reconstructed 
  adjust_pixels(state->tile->frame->rec->y, x_pos, pic_width, y_pos, pic_height, state->tile->frame->rec->stride,
    pic_width, pic_height);
  //Use if adjacent CTUs are reconstructed 
  /*adjust_pixels_CTU_plus_4_pix(state->tile->frame->rec->y, x_pos, state->tile->frame->width, y_pos, state->tile->frame->height, state->tile->frame->rec->stride,
    state->tile->frame->width, state->tile->frame->height);*/

  adjust_pixels_chroma(state->tile->frame->rec->u,
    x_pos >> chroma_scale_x,
    pic_width >> chroma_scale_x,
    y_pos >> chroma_scale_y,
    pic_height >> chroma_scale_y,
    state->tile->frame->rec->stride >> chroma_scale_x,
    pic_width >> chroma_scale_x,
    pic_height >> chroma_scale_y);
  adjust_pixels_chroma(state->tile->frame->rec->v,
    x_pos >> chroma_scale_x,
    pic_width >> chroma_scale_x,
    y_pos >> chroma_scale_y,
    pic_height >> chroma_scale_y,
    state->tile->frame->rec->stride >> chroma_scale_x,
    pic_width >> chroma_scale_x,
    pic_height >> chroma_scale_y);

  for (int i = y_pos; i < max_height; i += CLASSIFICATION_BLK_SIZE)
  {
    int n_height = MIN(i + CLASSIFICATION_BLK_SIZE, max_height) - i;

    for (int j = x_pos; j < max_width; j += CLASSIFICATION_BLK_SIZE)
    {
      int n_width = MIN(j + CLASSIFICATION_BLK_SIZE, max_width) - j;
      kvz_alf_derive_classification_blk(state, state->encoder_control->cfg.input_bitdepth + 4, n_height, n_width, j, i,
        j - x_pos + blk_dst_x, i - y_pos + blk_dst_y,
        alf_vb_luma_ctu_height, 
        alf_vb_luma_pos);
    }
  }
}

void kvz_alf_derive_classification_blk(encoder_state_t * const state,
  const int shift,
  const int n_height,
  const int n_width,
  const int blk_pos_x,
  const int blk_pos_y,
  const int blk_dst_x,
  const int blk_dst_y,
  const int vb_ctu_height,
  int vb_pos)
{
  videoframe_t* const frame = state->tile->frame;
  //int ***g_laplacian = state->tile->frame->alf_info->g_laplacian;
  //alf_classifier **g_classifier = state->tile->frame->alf_info->g_classifier;
  //CHECK((vb_ctu_height & (vb_ctu_height - 1)) != 0, "vb_ctu_height must be a power of 2");

  static const int th[16] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
  int laplacian[NUM_DIRECTIONS][CLASSIFICATION_BLK_SIZE + 5][CLASSIFICATION_BLK_SIZE + 5];
  memset(laplacian, 0, sizeof(laplacian));

  const int stride = frame->rec->stride;
  kvz_pixel *src = state->tile->frame->rec->y;
  const int max_activity = 15;

  int fl = 2;
  int fl_p1 = fl + 1;
  int fl2 = 2 * fl;

  int main_direction, secondary_direction, dir_temp_hv, dir_temp_d;
  int pix_y;

  int height = n_height + fl2;
  int width = n_width + fl2;
  int pos_x = blk_pos_x;
  int pos_y = blk_pos_y;
  int start_height = pos_y - fl_p1;

  for (int i = 0; i < height; i += 2)
  {
    int yoffset = (i + 1 + start_height) * stride - fl_p1;
    const kvz_pixel *src0 = &src[yoffset - stride];
    const kvz_pixel *src1 = &src[yoffset];
    const kvz_pixel *src2 = &src[yoffset + stride];
    const kvz_pixel *src3 = &src[yoffset + stride * 2];

    const int y = blk_dst_y - 2 + i;
    if (y > 0 && (y & (vb_ctu_height - 1)) == vb_pos - 2)
    {
      src3 = &src[yoffset + stride];
    }
    else if (y > 0 && (y & (vb_ctu_height - 1)) == vb_pos)
    {
      src0 = &src[yoffset];
    }

    int *p_y_ver = g_laplacian[ALF_VER][i];
    int *p_y_hor = g_laplacian[ALF_HOR][i];
    int *p_y_dig0 = g_laplacian[ALF_DIAG0][i];
    int *p_y_dig1 = g_laplacian[ALF_DIAG1][i];

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

        p_y_ver[j_m_6] += p_y_ver[j_m_4] + p_y_ver[j_m_2] + p_y_ver[j];
        p_y_hor[j_m_6] += p_y_hor[j_m_4] + p_y_hor[j_m_2] + p_y_hor[j];
        p_y_dig_0[j_m_6] += p_y_dig_0[j_m_4] + p_y_dig_0[j_m_2] + p_y_dig_0[j];
        p_y_dig_1[j_m_6] += p_y_dig_1[j_m_4] + p_y_dig_1[j_m_2] + p_y_dig_1[j];
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
      int sum_v = 0; int sum_h = 0; int sum_d0 = 0; int sum_d1 = 0;

      if (((i + blk_dst_y) % vb_ctu_height) == (vb_pos - 4))
      {
        sum_v = p_y_ver[j] + p_y_ver2[j] + p_y_ver4[j];
        sum_h = p_y_hor[j] + p_y_hor2[j] + p_y_hor4[j];
        sum_d0 = p_y_dig0[j] + p_y_dig02[j] + p_y_dig04[j];
        sum_d1 = p_y_dig1[j] + p_y_dig12[j] + p_y_dig14[j];
      }
      else if (((i + blk_dst_y) % vb_ctu_height) == vb_pos)
      {
        sum_v = p_y_ver2[j] + p_y_ver4[j] + p_y_ver6[j];
        sum_h = p_y_hor2[j] + p_y_hor4[j] + p_y_hor6[j];
        sum_d0 = p_y_dig02[j] + p_y_dig04[j] + p_y_dig06[j];
        sum_d1 = p_y_dig12[j] + p_y_dig14[j] + p_y_dig16[j];
      }
      else
      {
        sum_v = p_y_ver[j] + p_y_ver2[j] + p_y_ver4[j] + p_y_ver6[j];
        sum_h = p_y_hor[j] + p_y_hor2[j] + p_y_hor4[j] + p_y_hor6[j];
        sum_d0 = p_y_dig0[j] + p_y_dig02[j] + p_y_dig04[j] + p_y_dig06[j];
        sum_d1 = p_y_dig1[j] + p_y_dig12[j] + p_y_dig14[j] + p_y_dig16[j];
      }

      int temp_act = sum_v + sum_h;
      int activity = 0;

      const int y = (i + blk_dst_y) & (vb_ctu_height - 1);
      if (y == vb_pos - 4 || y == vb_pos)
      {
        activity = alf_clip3(0, max_activity, (temp_act * 96) >> shift);
      }
      else
      {
        activity = alf_clip3(0, max_activity, (temp_act * 64) >> shift);
      }

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
      if((uint32_t)d1 * (uint32_t)hv0 > (uint32_t)hv1 * (uint32_t)d0)
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

      int y_offset = i + blk_dst_y;
      int x_offset = j + blk_dst_x;

      alf_classifier *cl0 = g_classifier[y_offset] + x_offset;
      alf_classifier *cl1 = g_classifier[y_offset + 1] + x_offset;
      alf_classifier *cl2 = g_classifier[y_offset + 2] + x_offset;
      alf_classifier *cl3 = g_classifier[y_offset + 3] + x_offset;

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
  const kvz_pixel *src_pixels,
  kvz_pixel *dst_pixels,
  const int src_stride,
  const int dst_stride,
  const short* filter_set,
  const int16_t *fClipSet,
  clp_rng clp_rng,
  alf_component_id component_id,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  int blk_dst_x,
  int blk_dst_y,
  int vb_pos,
  const int vb_ctu_height)
{
  alf_filter_type const filter_type = component_id == COMPONENT_Y ? ALF_FILTER_7X7 : ALF_FILTER_5X5;
  const bool chroma = component_id == COMPONENT_Y ? 0 : 1;
  const int8_t bit_depth = state->encoder_control->bitdepth;

  if (chroma)
  {
    assert((int)filter_type == 0); //Chroma needs to have filtType == 0
  }
  /*#if !JVET_O0525_REMOVE_PCM
  //bool isDualTree = CS::isDualITree(cs);
  bool is_dual_tree = false;
  bool is_pcm_filter_enabled = ENABLE_PCM;
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  */

  //const int srcStride = srcLuma.stride;
  //const int src_stride = frame->rec->stride;
  //const int dstStride = dstLuma.stride;
  //const int dst_stride = frame->rec->stride;

  const int start_height = y_pos;
  const int end_height = start_height + height;
  const int start_width = x_pos;
  const int end_width = start_width + width;

  const kvz_pixel *src = src_pixels;
  kvz_pixel *dst = dst_pixels + blk_dst_y * dst_stride;

  const kvz_pixel *p_img_y_pad_0, *p_img_y_pad_1, *p_img_y_pad_2, *p_img_y_pad_3, *p_img_y_pad_4, *p_img_y_pad_5, *p_img_y_pad_6;
  const kvz_pixel *p_img_0, *p_img_1, *p_img_2, *p_img_3, *p_img_4, *p_img_5, *p_img_6;

  const short *coef = filter_set;
  const int16_t *clip = fClipSet;

  const int shift = bit_depth - 1;

  const int offset = 1 << (shift - 1);

  int transpose_idx = 0;
  const int cls_size_y = 4;
  const int cls_size_x = 4;

  /*#if !JVET_O0525_REMOVE_PCM
  bool pcm_flags_2x2[4] = { 0,0,0,0 };*/

  assert((start_height % cls_size_y) == 0); //Wrong startHeight in filtering
  assert((start_width % cls_size_x) == 0); //Wrong startWidth in filtering
  assert(((end_height - start_height) % cls_size_y) == 0); //Wrong endHeight in filtering
  assert(((end_width - start_width) % cls_size_x) == 0); //Wrong endWidth in filtering

  alf_classifier *p_class = NULL;

  int dst_stride2 = dst_stride * cls_size_y;
  int src_stride2 = src_stride * cls_size_y;

  //std::vector<Pel> filter_coeff(MAX_NUM_ALF_LUMA_COEFF);
  int filter_coeff[MAX_NUM_ALF_LUMA_COEFF];
  memset(filter_coeff, 0, MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
  //std::array<int, MAX_NUM_ALF_LUMA_COEFF> filterClipp;
  int filter_clipp[MAX_NUM_ALF_LUMA_COEFF];
  memset(filter_clipp, 0, MAX_NUM_ALF_LUMA_COEFF * sizeof(int));

  p_img_y_pad_0 = src + start_height * src_stride + start_width;
  p_img_y_pad_1 = p_img_y_pad_0 + src_stride;
  p_img_y_pad_2 = p_img_y_pad_0 - src_stride;
  p_img_y_pad_3 = p_img_y_pad_1 + src_stride;
  p_img_y_pad_4 = p_img_y_pad_2 - src_stride;
  p_img_y_pad_5 = p_img_y_pad_3 + src_stride;
  p_img_y_pad_6 = p_img_y_pad_4 - src_stride;

  kvz_pixel* p_rec_0 = dst + blk_dst_x;//start_width;
  kvz_pixel* p_rec_1 = p_rec_0 + dst_stride;

  for (int i = 0; i < end_height - start_height; i += cls_size_y)
  {
    if (!chroma)
    {
      p_class = g_classifier[blk_dst_y + i] + blk_dst_x;
    }

    for (int j = 0; j < end_width - start_width; j += cls_size_x)
    {
      if (!chroma)
      {
        alf_classifier cl = p_class[j];
        transpose_idx = cl.transpose_idx;
        /*#if !JVET_O0525_REMOVE_PCM
        if (is_pcm_filter_enabled && cl.class_idx == ALF_UNUSED_CLASS_IDX && transpose_idx == ALF_UNUSED_TRANSPOSE_IDX)
        {
          continue;
        }*/
        coef = filter_set + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
        clip = fClipSet + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
      }
      /*#if !JVET_O0525_REMOVE_PCM
      else if (is_pcm_filter_enabled)
      {
        int  blk_x, blk_y;
        bool *flags = pcm_flags_2x2;

        // check which chroma 2x2 blocks use PCM
        // chroma PCM may not be aligned with 4x4 ALF processing grid
        for (blk_y = 0; blk_y < 4; blk_y += 2)
        {
          for (blk_x = 0; blk_x < 4; blk_x += 2)
          {
            //Position pos(j + blkDst.x + blkX, i + blkDst.y + blkY);
            //CodingUnit* cu = is_dual_tree ? cs.getCU(pos, CH_C) : cs.getCU(recalcPosition(nChromaFormat, CH_C, CH_L, pos), CH_L);
            *flags++ = 1; //cu->ipcm ? 1 : 0;
          }
        }

        // skip entire 4x4 if all chroma 2x2 blocks use PCM
        if (pcm_flags_2x2[0] && pcm_flags_2x2[1] && pcm_flags_2x2[2] && pcm_flags_2x2[3])
        {
          continue;
        }
      }*/


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

          filter_clipp[0] = clip[9];
          filter_clipp[1] = clip[4];
          filter_clipp[2] = clip[10];
          filter_clipp[3] = clip[8];
          filter_clipp[4] = clip[1];
          filter_clipp[5] = clip[5];
          filter_clipp[6] = clip[11];
          filter_clipp[7] = clip[7];
          filter_clipp[8] = clip[3];
          filter_clipp[9] = clip[0];
          filter_clipp[10] = clip[2];
          filter_clipp[11] = clip[6];
          filter_clipp[12] = clip[12];
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

          filter_clipp[0] = clip[0];
          filter_clipp[1] = clip[3];
          filter_clipp[2] = clip[2];
          filter_clipp[3] = clip[1];
          filter_clipp[4] = clip[8];
          filter_clipp[5] = clip[7];
          filter_clipp[6] = clip[6];
          filter_clipp[7] = clip[5];
          filter_clipp[8] = clip[4];
          filter_clipp[9] = clip[9];
          filter_clipp[10] = clip[10];
          filter_clipp[11] = clip[11];
          filter_clipp[12] = clip[12];

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
          
          filter_clipp[0] = clip[9];
          filter_clipp[1] = clip[8];
          filter_clipp[2] = clip[10];
          filter_clipp[3] = clip[4];
          filter_clipp[4] = clip[3];
          filter_clipp[5] = clip[7];
          filter_clipp[6] = clip[11];
          filter_clipp[7] = clip[5];
          filter_clipp[8] = clip[1];
          filter_clipp[9] = clip[0];
          filter_clipp[10] = clip[2];
          filter_clipp[11] = clip[6];
          filter_clipp[12] = clip[12];
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
          
          filter_clipp[0] = clip[0];
          filter_clipp[1] = clip[1];
          filter_clipp[2] = clip[2];
          filter_clipp[3] = clip[3];
          filter_clipp[4] = clip[4];
          filter_clipp[5] = clip[5];
          filter_clipp[6] = clip[6];
          filter_clipp[7] = clip[7];
          filter_clipp[8] = clip[8];
          filter_clipp[9] = clip[9];
          filter_clipp[10] = clip[10];
          filter_clipp[11] = clip[11];
          filter_clipp[12] = clip[12];
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
          
          filter_clipp[0] = clip[4];
          filter_clipp[1] = clip[1];
          filter_clipp[2] = clip[5];
          filter_clipp[3] = clip[3];
          filter_clipp[4] = clip[0];
          filter_clipp[5] = clip[2];
          filter_clipp[6] = clip[6];

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
          
          filter_clipp[0] = clip[0];
          filter_clipp[1] = clip[3];
          filter_clipp[2] = clip[2];
          filter_clipp[3] = clip[1];
          filter_clipp[4] = clip[4];
          filter_clipp[5] = clip[5];
          filter_clipp[6] = clip[6];

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
          
          filter_clipp[0] = clip[4];
          filter_clipp[1] = clip[3];
          filter_clipp[2] = clip[5];
          filter_clipp[3] = clip[1];
          filter_clipp[4] = clip[0];
          filter_clipp[5] = clip[2];
          filter_clipp[6] = clip[6];

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
          
          filter_clipp[0] = clip[0];
          filter_clipp[1] = clip[1];
          filter_clipp[2] = clip[2];
          filter_clipp[3] = clip[3];
          filter_clipp[4] = clip[4];
          filter_clipp[5] = clip[5];
          filter_clipp[6] = clip[6];

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

        const int y_vb = (blk_dst_y + i + ii) & (vb_ctu_height - 1);
        if (y_vb < vb_pos && (y_vb >= vb_pos - (chroma ? 2 : 4)))   // above
        {
          p_img_1 = (y_vb == vb_pos - 1) ? p_img_0 : p_img_1;
          p_img_3 = (y_vb >= vb_pos - 2) ? p_img_1 : p_img_3;
          p_img_5 = (y_vb >= vb_pos - 3) ? p_img_3 : p_img_5;

          p_img_2 = (y_vb == vb_pos - 1) ? p_img_0 : p_img_2;
          p_img_4 = (y_vb >= vb_pos - 2) ? p_img_2 : p_img_4;
          p_img_6 = (y_vb >= vb_pos - 3) ? p_img_4 : p_img_6;
        }
        else if (y_vb >= vb_pos && (y_vb <= vb_pos + (chroma ? 1 : 3)))   // bottom
        {
          p_img_2 = (y_vb == vb_pos) ? p_img_0 : p_img_2;
          p_img_4 = (y_vb <= vb_pos + 1) ? p_img_2 : p_img_4;
          p_img_6 = (y_vb <= vb_pos + 2) ? p_img_4 : p_img_6;

          p_img_1 = (y_vb == vb_pos) ? p_img_0 : p_img_1;
          p_img_3 = (y_vb <= vb_pos + 1) ? p_img_1 : p_img_3;
          p_img_5 = (y_vb <= vb_pos + 2) ? p_img_3 : p_img_5;
        }

        bool is_near_vb_above = y_vb < vb_pos && (y_vb >= vb_pos - 1);
        bool is_near_vb_below = y_vb >= vb_pos && (y_vb <= vb_pos);
        for (int jj = 0; jj < cls_size_x; jj++)
        {
          /*#if !JVET_O0525_REMOVE_PCM
          // skip 2x2 PCM chroma blocks
          if (chroma && is_pcm_filter_enabled)
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
          }*/

          int sum = 0;
          const kvz_pixel curr = p_img_0[+0];
          if (filter_type == ALF_FILTER_7X7)
          {
            sum += filter_coeff[0] * (clip_alf(filter_clipp[0], curr, p_img_5[+0], p_img_6[+0]));

            sum += filter_coeff[1] * (clip_alf(filter_clipp[1], curr, p_img_3[+1], p_img_4[-1]));
            sum += filter_coeff[2] * (clip_alf(filter_clipp[2], curr, p_img_3[+0], p_img_4[+0]));
            sum += filter_coeff[3] * (clip_alf(filter_clipp[3], curr, p_img_3[-1], p_img_4[+1]));

            sum += filter_coeff[4] * (clip_alf(filter_clipp[4], curr, p_img_1[+2], p_img_2[-2]));
            sum += filter_coeff[5] * (clip_alf(filter_clipp[5], curr, p_img_1[+1], p_img_2[-1]));
            sum += filter_coeff[6] * (clip_alf(filter_clipp[6], curr, p_img_1[+0], p_img_2[+0]));
            sum += filter_coeff[7] * (clip_alf(filter_clipp[7], curr, p_img_1[-1], p_img_2[+1]));
            sum += filter_coeff[8] * (clip_alf(filter_clipp[8], curr, p_img_1[-2], p_img_2[+2]));
            
            sum += filter_coeff[9] * (clip_alf(filter_clipp[9], curr, p_img_0[+3], p_img_0[-3]));
            sum += filter_coeff[10] * (clip_alf(filter_clipp[10], curr, p_img_0[+2], p_img_0[-2]));
            sum += filter_coeff[11] * (clip_alf(filter_clipp[11], curr, p_img_0[+1], p_img_0[-1]));
          }
          else
          {
            sum += filter_coeff[0] * (clip_alf(filter_clipp[0], curr, p_img_3[+0], p_img_4[+0]));

            sum += filter_coeff[1] * (clip_alf(filter_clipp[1], curr, p_img_1[+1], p_img_2[-1]));
            sum += filter_coeff[2] * (clip_alf(filter_clipp[2], curr, p_img_1[+0], p_img_2[+0]));
            sum += filter_coeff[3] * (clip_alf(filter_clipp[3], curr, p_img_1[-1], p_img_2[+1]));

            sum += filter_coeff[4] * (clip_alf(filter_clipp[4], curr, p_img_0[+2], p_img_0[-2]));
            sum += filter_coeff[5] * (clip_alf(filter_clipp[5], curr, p_img_0[+1], p_img_0[-1]));
          }

          sum = (sum + offset) >> shift;
          sum += curr;

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

//-------------------------CC ALF encoding functions------------------------

void apply_cc_alf_filter(encoder_state_t * const state, alf_component_id comp_id, const kvz_pixel *dst_buf,
  const kvz_pixel *rec_yuv_ext, const int luma_stride, uint8_t *filter_control,
  const short filter_set[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF],
  const int   selected_filter_idx)
{
  enum kvz_chroma_format chroma_format = state->encoder_control->chroma_format;
  uint8_t component_scale_y = (comp_id == COMPONENT_Y || chroma_format != KVZ_CSP_420) ? 0 : 1;
  uint8_t component_scale_x = (comp_id == COMPONENT_Y || chroma_format == KVZ_CSP_444) ? 0 : 1;
  const int pic_height = state->tile->frame->height;
  const int pic_width = state->tile->frame->width;
  const int max_ctu_height_log2 = kvz_math_floor_log2(LCU_WIDTH);
  const int max_ctu_width_log2 = kvz_math_floor_log2(LCU_WIDTH);
  const int width_in_ctus = state->tile->frame->width_in_lcu;
  const int alf_vb_luma_ctu_height = LCU_WIDTH;
  const int alf_vb_luma_pos = LCU_WIDTH - ALF_VB_POS_ABOVE_CTUROW_LUMA;

  int ctu_idx = 0;
  for (int y_pos = 0; y_pos < pic_height; y_pos += LCU_WIDTH)
  {
    for (int x_pos = 0; x_pos < pic_width; x_pos += LCU_WIDTH)
    {
      int filter_idx =
        (filter_control == NULL)
        ? selected_filter_idx
        : filter_control[(y_pos >> max_ctu_height_log2) * width_in_ctus + (x_pos >> max_ctu_width_log2)];
      bool skip_filtering = (filter_control != NULL && filter_idx == 0) ? true : false;
      if (!skip_filtering)
      {
        if (filter_control != NULL)
        {
          filter_idx--;
        }

        const int16_t *filter_coeff = filter_set[filter_idx];

        const int width = (x_pos + LCU_WIDTH > pic_width) ? (pic_width - x_pos) : LCU_WIDTH;
        const int height = (y_pos + LCU_WIDTH > pic_height) ? (pic_height - y_pos) : LCU_WIDTH;

        int raster_slice_alf_pad = 0;
        {
          filter_blk_cc_alf(state, dst_buf, rec_yuv_ext, luma_stride, comp_id, filter_coeff, g_clp_rngs, alf_vb_luma_ctu_height,
            alf_vb_luma_pos, x_pos >> component_scale_x, y_pos >> component_scale_y,
            width >> component_scale_x, height >> component_scale_y);
        }
      }
      ctu_idx++;
    }
  }
}

void setup_cc_alf_aps(encoder_state_t * const state)
{
  if (g_cc_alf_filter_param.cc_alf_filter_enabled[COMPONENT_Cb - 1])
  {
    int  cc_alf_cb_aps_id = state->slice->tile_group_cc_alf_cb_aps_id;
    alf_aps *aps = &state->encoder_control->cfg.param_set_map[cc_alf_cb_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
    if (aps->aps_id >= 0 && aps->aps_id < ALF_CTB_MAX_NUM_APS)
    {
      //aps = m_apsMap->allocatePS((cc_alf_cb_aps_id << NUM_APS_TYPE_LEN) + ALF_APS);
      aps->temporal_id = 0; // cs.slice->getTLayer()
    }
    aps->cc_alf_aps_param.cc_alf_filter_enabled[COMPONENT_Cb - 1] = 1;
    aps->cc_alf_aps_param.cc_alf_filter_count[COMPONENT_Cb - 1] = g_cc_alf_filter_param.cc_alf_filter_count[COMPONENT_Cb - 1];
    for ( int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++ )
    {
      aps->cc_alf_aps_param.cc_alf_filter_idx_enabled[COMPONENT_Cb - 1][filter_idx] =
        g_cc_alf_filter_param.cc_alf_filter_idx_enabled[COMPONENT_Cb - 1][filter_idx];
      memcpy(aps->cc_alf_aps_param.cc_alf_coeff[COMPONENT_Cb - 1][filter_idx],
        g_cc_alf_filter_param.cc_alf_coeff[COMPONENT_Cb - 1][filter_idx], sizeof(short) * MAX_NUM_CC_ALF_CHROMA_COEFF);
    }
    aps->aps_id = cc_alf_cb_aps_id;
    aps->aps_type = T_ALF_APS;
    if (g_reuse_aps_id[COMPONENT_Cb - 1] < 0)
    {
      aps->cc_alf_aps_param.new_cc_alf_filter[COMPONENT_Cb - 1] = 1;
      state->encoder_control->cfg.param_set_map[cc_alf_cb_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
      aps->temporal_id = 0; // cs.slice->getTLayer()
    }
    state->slice->tile_group_cc_alf_cb_enabled_flag = true;
  }
  else
  {
    state->slice->tile_group_cc_alf_cb_enabled_flag = false;
  }
  if (g_cc_alf_filter_param.cc_alf_filter_enabled[COMPONENT_Cr - 1])
  {
    int  cc_alf_cr_aps_id = state->slice->tile_group_cc_alf_cr_aps_id;
    alf_aps *aps = &state->encoder_control->cfg.param_set_map[cc_alf_cr_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
    if (aps->aps_id >= 0 && aps->aps_id < ALF_CTB_MAX_NUM_APS)
    {
      //aps = m_apsMap->allocatePS((cc_alf_cb_aps_id << NUM_APS_TYPE_LEN) + ALF_APS);
      aps->temporal_id = 0; // cs.slice->getTLayer()
    }
    aps->cc_alf_aps_param.cc_alf_filter_enabled[COMPONENT_Cr - 1] = 1;
    aps->cc_alf_aps_param.cc_alf_filter_count[COMPONENT_Cr - 1] = g_cc_alf_filter_param.cc_alf_filter_count[COMPONENT_Cr - 1];
    for (int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++)
    {
      aps->cc_alf_aps_param.cc_alf_filter_idx_enabled[COMPONENT_Cr - 1][filter_idx] =
        g_cc_alf_filter_param.cc_alf_filter_idx_enabled[COMPONENT_Cr - 1][filter_idx];
      memcpy(aps->cc_alf_aps_param.cc_alf_coeff[COMPONENT_Cr - 1][filter_idx],
        g_cc_alf_filter_param.cc_alf_coeff[COMPONENT_Cr - 1][filter_idx], sizeof(short) * MAX_NUM_CC_ALF_CHROMA_COEFF);
    }
    aps->aps_id = cc_alf_cr_aps_id;
    aps->aps_type = T_ALF_APS;
    if (g_reuse_aps_id[COMPONENT_Cr - 1] < 0)
    {
      aps->cc_alf_aps_param.new_cc_alf_filter[COMPONENT_Cr - 1] = 1;
      state->encoder_control->cfg.param_set_map[cc_alf_cr_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
      aps->temporal_id = 0; // cs.slice->getTLayer()
    }
    state->slice->tile_group_cc_alf_cr_enabled_flag = true;
  }
  else
  {
    state->slice->tile_group_cc_alf_cr_enabled_flag = false;
  }
}

void round_filt_coeff_cc_alf(int16_t *filter_coeff_quant, double *filter_coeff, const int num_coeff, const int factor )
{
  for( int i = 0; i < num_coeff; i++ )
  {
    int sign = filter_coeff[i] > 0 ? 1 : -1;
    double best_err = 128.0*128.0;
    int best_index = 0;
    for(int k = 0; k < CCALF_CANDS_COEFF_NR; k++)
    {
      double err = (filter_coeff[i] * sign * factor - cc_alf_small_tab[k]);
      err = err*err;
      if(err < best_err)
      {
        best_err = err;
        best_index = k;
      }
    }
    filter_coeff_quant[i] = cc_alf_small_tab[best_index] * sign;
  }
}

int get_coeff_rate_cc_alf(short chroma_coeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF], bool filter_enabled[MAX_NUM_CC_ALF_FILTERS], uint8_t filter_count, alf_component_id comp_id)
{
  int bits = 0;

  if ( filter_count > 0 )
  {
    bits += length_uvlc(filter_count - 1);
    int signaled_filter_count = 0;
    for ( int filterIdx=0; filterIdx<MAX_NUM_CC_ALF_FILTERS; filterIdx++ )
    {
      if (filter_enabled[filterIdx])
      {
        // Filter coefficients
        for (int i = 0; i < CC_ALF_NUM_COEFF - 1; i++)
        {
          bits += CCALF_BITS_PER_COEFF_LEVEL + (chroma_coeff[filterIdx][i] == 0 ? 0 : 1);
        }

        signaled_filter_count++;
      }
    }
    assert(signaled_filter_count == filter_count); //Number of filter signaled not same as indicated
  }

  return bits;
}

void derive_cc_alf_filter_coeff( alf_component_id comp_id, short filter_coeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF], const uint8_t filter_idx )
{
  int forward_tab[CCALF_CANDS_COEFF_NR * 2 - 1] = {0};
  for (int i = 0; i < CCALF_CANDS_COEFF_NR; i++)
  {
    forward_tab[CCALF_CANDS_COEFF_NR - 1 + i] = cc_alf_small_tab[i];
    forward_tab[CCALF_CANDS_COEFF_NR - 1 - i] = (-1) * cc_alf_small_tab[i];
  }

  double filter_coeff_dbl[MAX_NUM_CC_ALF_CHROMA_COEFF];
  int16_t filter_coeff_int[MAX_NUM_CC_ALF_CHROMA_COEFF];

  memset(filter_coeff_int, 0, sizeof(filter_coeff_int));

  double k_e[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];
  double ky[MAX_NUM_ALF_LUMA_COEFF];
  const int size = CC_ALF_NUM_COEFF - 1;

  for (int k = 0; k < size; k++)
  {
    ky[k] = g_alf_covariance_frame_cc_alf[comp_id - 1][filter_idx].y[0][k];
    for (int l = 0; l < size; l++)
    {
      k_e[k][l] = g_alf_covariance_frame_cc_alf[comp_id - 1][filter_idx].ee[0][0][k][l];
    }
  }

  //m_alfCovarianceFrameCcAlf[compID - 1][0][filterIdx].gnsSolveByChol(kE, ky, filterCoeffDbl, size);
  gns_solve_by_chol(k_e, ky, filter_coeff_dbl, size);
  round_filt_coeff_cc_alf(filter_coeff_int, filter_coeff_dbl, size, (1 << 7 /*m_scaleBits*/));

  for (int k = 0; k < size; k++)
  {
    assert(!(filter_coeff_int[k] < -(1 << CCALF_DYNAMIC_RANGE))); // this is not possible: filter_coeff_int[k] <  -(1 << CCALF_DYNAMIC_RANGE)
    assert(!(filter_coeff_int[k] > (1 << CCALF_DYNAMIC_RANGE))); // this is not possible: filter_coeff_int[k] >  (1 << CCALF_DYNAMIC_RANGE)
  }

  // Refine quanitzation
  int modified = 1;
  double err_ref = calc_error_for_cc_alf_coeffs(&g_alf_covariance_frame_cc_alf[comp_id - 1][filter_idx], filter_coeff_int, size, (7/*m_scaleBits*/+1));
  while (modified)
  {
    modified = 0;
    for (int i = 1; i > -2; i -= 2 )
    {
      int delta = i;
      double err_min = MAX_DOUBLE;
      int    idx_min = -1;
      int min_index = -1;

      for (int k = 0; k < size; k++)
      {
        int org_idx = -1;
        for (int i = 0; i < CCALF_CANDS_COEFF_NR * 2 - 1; i++)
        {
          if (forward_tab[i] == filter_coeff_int[k])
          {
            org_idx = i;
            break;
          }
        }
        assert(!(org_idx < 0)); //this is wrong, does not find coeff from forward_tab
        if ( (org_idx - delta < 0) || (org_idx - delta >= CCALF_CANDS_COEFF_NR * 2 - 1) )
          continue;

        filter_coeff_int[k] = forward_tab[org_idx - delta];
        double error = calc_error_for_cc_alf_coeffs(&g_alf_covariance_frame_cc_alf[comp_id - 1][filter_idx], filter_coeff_int, size, (7/*m_scaleBits*/+1));
        if( error < err_min )
        {
          err_min = error;
          idx_min = k;
          min_index = org_idx;
        }
        filter_coeff_int[k] = forward_tab[org_idx];
      }
      if (err_min < err_ref)
      {
        min_index -= delta;
        assert(!(min_index < 0));// this is wrong, index - delta < 0
        assert(!(min_index >= CCALF_CANDS_COEFF_NR * 2 - 1)); // this is wrong, index - delta >= CCALF_CANDS_COEFF_NR * 2 - 1
        filter_coeff_int[idx_min] = forward_tab[min_index];
        modified++;
        err_ref = err_min;
      }
    }
  }

  for (int k = 0; k < (size + 1); k++)
  {
    assert(!((filter_coeff_int[k] < -(1 << CCALF_DYNAMIC_RANGE)) || (filter_coeff_int[k] > (1 << CCALF_DYNAMIC_RANGE)))); //Exceeded valid range for CC ALF coefficient
    filter_coeff[filter_idx][k] = filter_coeff_int[k];
  }
}

void determine_control_idc_values(encoder_state_t *const state, const alf_component_id comp_id,
  const int ctu_width_c, const int ctu_height_c, const int pic_width_c,
  const int pic_height_c, double **unfiltered_distortion,
  uint64_t *training_distortion[MAX_NUM_CC_ALF_FILTERS],
  uint64_t *luma_swing_greater_than_threshold_count,
  uint64_t *chroma_sample_count_near_mid_point,
  bool reuse_temporal_filter_coeff, uint8_t *training_cov_control,
  uint8_t *filter_control, uint64_t *cur_total_distortion,
  double *cur_total_rate, bool filter_enabled[MAX_NUM_CC_ALF_FILTERS],
  uint8_t  map_filter_idx_to_filter_idc[MAX_NUM_CC_ALF_FILTERS + 1],
  uint8_t *cc_alf_filter_count)
{
  bool cur_filter_enabled[MAX_NUM_CC_ALF_FILTERS];
  //std::fill_n(cur_filter_enabled, MAX_NUM_CC_ALF_FILTERS, false);
  memset(cur_filter_enabled, false, sizeof(cur_filter_enabled));

#if MAX_NUM_CC_ALF_FILTERS>1
  filter_idx_count filter_idx_count[MAX_NUM_CC_ALF_FILTERS];
  for (int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++)
  {
    filter_idx_count[i].count = 0;
    filter_idx_count[i].filter_idx = i;
  }

  double prev_rate = (*cur_total_rate);
#endif

  cabac_data_t ctx_initial;
  cabac_data_t ctx_best;
  cabac_data_t ctx_start;
  memcpy(&ctx_initial, &cabac_estimator, sizeof(ctx_initial));
  memcpy(&ctx_best, &cabac_estimator, sizeof(ctx_best));
  cabac_estimator.only_count = 1;
  ctx_initial.only_count = 1;
  ctx_best.only_count = 1;
  ctx_start.only_count = 1;

  enum kvz_chroma_format chroma_format = state->encoder_control->chroma_format;
  uint8_t component_scale_y = (comp_id == COMPONENT_Y || chroma_format != KVZ_CSP_420) ? 0 : 1;
  uint8_t component_scale_x = (comp_id == COMPONENT_Y || chroma_format == KVZ_CSP_444) ? 0 : 1;
  double lambda = state->frame->lambda;
  bool limit_cc_alf = state->encoder_control->cfg.qp >= 37;
  
  int ctu_idx = 0;
  for (int y_ctu = 0; y_ctu < pic_height_c; y_ctu += ctu_height_c)
  {
    for (int x_ctu = 0; x_ctu < pic_width_c; x_ctu += ctu_width_c)
    {
      uint64_t ssd;
      double   rate;
      double   cost;

      uint64_t best_ssd = (MAX_INT64 >> 1);
      double   best_rate = MAX_DOUBLE;
      double   best_cost = MAX_DOUBLE;
      uint8_t  best_filter_idc = 0;
      uint8_t  best_filter_idx = 0;
      const uint32_t threshold_s = MIN(pic_height_c - y_ctu, ctu_height_c) << (chroma_format != KVZ_CSP_420 ? 0 : 1);
      const uint32_t number_of_chroma_samples = MIN(pic_height_c - y_ctu, ctu_height_c) * MIN(pic_height_c - x_ctu, ctu_width_c);
      const uint32_t threshold_c = (number_of_chroma_samples >> 2);

      memcpy(&cabac_estimator, &ctx_best, sizeof(cabac_estimator));
      memcpy(&ctx_start, &cabac_estimator, sizeof(ctx_start));

      for (int filter_idx = 0; filter_idx <= MAX_NUM_CC_ALF_FILTERS; filter_idx++)
      {
        uint8_t filter_idc = map_filter_idx_to_filter_idc[filter_idx];
        if (filter_idx < MAX_NUM_CC_ALF_FILTERS && !filter_enabled[filter_idx])
        {
          continue;
        }

        if (filter_idx == MAX_NUM_CC_ALF_FILTERS)
        {
          ssd = (uint64_t)unfiltered_distortion[comp_id][ctu_idx];   // restore saved distortion computation
        }
        else
        {
          ssd = training_distortion[filter_idx][ctu_idx];
        }

        memcpy(&cabac_estimator, &ctx_start, sizeof(cabac_estimator));
        kvz_cabac_reset_bits(&cabac_estimator);

        //const Position lumaPos = Position({ xCtu << getComponentScaleX(comp_id, cs.pcv->chrFormat),
        //  yCtu << getComponentScaleY(comp_id, cs.pcv->chrFormat) });
        code_cc_alf_filter_control_idc(state, &cabac_estimator, filter_idc, comp_id, ctu_idx,
          filter_control, *cc_alf_filter_count);
        //rate = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
        rate = (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3);
        cost = rate * lambda + ssd;

        bool limitation_exceeded = false;
        if (limit_cc_alf && filter_idx < MAX_NUM_CC_ALF_FILTERS)
        {
          limitation_exceeded = limitation_exceeded || (luma_swing_greater_than_threshold_count[ctu_idx] >= threshold_s);
          limitation_exceeded = limitation_exceeded || (chroma_sample_count_near_mid_point[ctu_idx] >= threshold_c);
        }
        if (cost < best_cost && !limitation_exceeded)
        {
          best_cost = cost;
          best_rate = rate;
          best_ssd = ssd;
          best_filter_idc = filter_idc;
          best_filter_idx = filter_idx;

          //ctx_best = SubCtx(Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx());
          memcpy(&ctx_best, &cabac_estimator, sizeof(ctx_best));

          training_cov_control[ctu_idx] = (filter_idx == MAX_NUM_CC_ALF_FILTERS) ? 0 : (filter_idx + 1);
          filter_control[ctu_idx] = (filter_idx == MAX_NUM_CC_ALF_FILTERS) ? 0 : (filter_idx + 1);
        }
      }
      if (best_filter_idc != 0)
      {
        cur_filter_enabled[best_filter_idx] = true;
#if MAX_NUM_CC_ALF_FILTERS>1
        filter_idx_count[best_filter_idx].count++;
#endif
      }
      (*cur_total_rate) += best_rate;
      (*cur_total_distortion) += best_ssd;
      ctu_idx++;
    }
  }

#if MAX_NUM_CC_ALF_FILTERS>1
  if (!reuse_temporal_filter_coeff)
  {
    memcpy(cur_filter_enabled, filter_enabled, sizeof(cur_filter_enabled));
    qsort(filter_idx_count, MAX_NUM_CC_ALF_FILTERS, sizeof(*filter_idx_count), comparator);

    int filter_idc = 1;
    (*cc_alf_filter_count) = 0;
    for (int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++)
    {
      const int filter_idx = filter_idx_count[i].filter_idx;
      if (filter_enabled[filter_idx])
      {
        map_filter_idx_to_filter_idc[filter_idx] = filter_idc;
        filter_idc++;
        (*cc_alf_filter_count)++;
      }
    }

    (*cur_total_rate) = prev_rate;
    //m_CABACEstimator->getCtx() = ctx_initial;
    memcpy(&cabac_estimator, &ctx_initial, sizeof(cabac_estimator));
    //m_CABACEstimator->resetBits();
    kvz_cabac_reset_bits(&cabac_estimator);

    int ctu_idx = 0;
    for (int y = 0; y < pic_height_c; y += ctu_height_c)
    {
      for (int x = 0; x < pic_width_c; x += ctu_width_c)
      {
        const int filter_idx_plus1 = filter_control[ctu_idx];
        code_cc_alf_filter_control_idc(state, &cabac_estimator, (filter_idx_plus1 == 0 ? 0
          : map_filter_idx_to_filter_idc[filter_idx_plus1 - 1]),
          comp_id, ctu_idx, filter_control, *cc_alf_filter_count);

        ctu_idx++;
      }
    }
    (*cur_total_rate) += (23 - cabac_estimator.bits_left) + (cabac_estimator.num_buffered_bytes << 3);
  }
#endif

  // restore for next iteration
  memcpy(&cabac_estimator, &ctx_initial, sizeof(cabac_estimator));
}

void get_available_cc_alf_aps_ids(encoder_state_t *const state, alf_component_id compID, 
  int *aps_ids_size, int *aps_ids)
{
  alf_aps* apss = state->slice->apss;
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    copy_aps(&apss[i], &state->encoder_control->cfg.param_set_map[i + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set);
  }

  int aps_id_checked = 0, cur_aps_id = g_aps_id_start;
  if (cur_aps_id < ALF_CTB_MAX_NUM_APS)
  {
    while (aps_id_checked < ALF_CTB_MAX_NUM_APS &&
      state->frame->slicetype != KVZ_SLICE_I &&
      (*aps_ids_size) < ALF_CTB_MAX_NUM_APS  &&
      !(state->frame->pictype == KVZ_NAL_IDR_W_RADL || state->frame->pictype == KVZ_NAL_IDR_N_LP)) //&& !cs.slice->getPendingRasInit()
    {
      alf_aps cur_aps = state->slice->apss[cur_aps_id];
      bool aps_found = (0 <= cur_aps.aps_id && cur_aps.aps_id < ALF_CTB_MAX_NUM_APS);
      if (aps_found && cur_aps.temporal_id <= state->slice->id && cur_aps.cc_alf_aps_param.new_cc_alf_filter[compID - 1])
      {
        aps_ids[(*aps_ids_size)] = cur_aps_id;
        (*aps_ids_size) += 1;
      }
      aps_id_checked++;
      cur_aps_id = (cur_aps_id + 1) % ALF_CTB_MAX_NUM_APS;
    }
  }
}

void derive_cc_alf_filter(encoder_state_t * const state, alf_component_id comp_id, 
  const kvz_picture *org_yuv, const kvz_picture *rec_dst_yuv )
{
  if (!state->slice->tile_group_alf_enabled_flag[COMPONENT_Y])
  {
    g_cc_alf_filter_param.cc_alf_filter_enabled[comp_id - 1] = false;
    return;
  }

  bool limit_cc_alf = state->encoder_control->cfg.qp >= 37; // m_encCfg->getCCALFQpThreshold();
  if (limit_cc_alf) // && state->slice. cs.slice->getSliceQp() <= m_encCfg->getBaseQP() + 1)
  {
    g_cc_alf_filter_param.cc_alf_filter_enabled[comp_id - 1] = false;
    return;
  }

  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  uint8_t best_map_filter_idx_to_filter_idc[MAX_NUM_CC_ALF_FILTERS+1];
  bool scale_x = (comp_id == COMPONENT_Y || chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool scale_y = (comp_id == COMPONENT_Y || chroma_fmt != KVZ_CSP_420) ? 0 : 1;
  const int ctu_width_c            = LCU_WIDTH >> scale_x;
  const int ctu_height_c           = LCU_WIDTH >> scale_y;
  const int pic_width_c            = state->tile->frame->width >> scale_x;
  const int pic_height_c           = state->tile->frame->height >> scale_y;
  const int pic_stride_c           = rec_dst_yuv->stride >> scale_x;
  const int8_t bit_depth = state->encoder_control->bitdepth;
  const int max_training_iter_count = 15;
  int max_ctu_height_log2 = kvz_math_floor_log2(LCU_WIDTH);
  int max_ctu_height_log2_chrma = kvz_math_floor_log2(LCU_WIDTH) >> scale_y;
  int max_ctu_width_log2 = kvz_math_floor_log2(LCU_WIDTH);
  int max_ctu_width_log2_chrma = kvz_math_floor_log2(LCU_WIDTH) >> scale_x;
  int32_t ctus_in_width = state->tile->frame->width_in_lcu;
  const uint32_t num_ctus_in_pic = state->lcu_order_count;
  short best_filter_coeff_set[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
  bool best_filter_idx_enabled[MAX_NUM_CC_ALF_FILTERS];
  uint8_t best_filter_count = 0;
  double lambda = state->frame->lambda;

  if (limit_cc_alf)
  {
    count_luma_swing_greater_than_threshold(rec_dst_yuv->y, rec_dst_yuv->stride, rec_dst_yuv->height, rec_dst_yuv->width,
      max_ctu_width_log2, max_ctu_height_log2, g_luma_swing_greater_than_threshold_count,
      ctus_in_width, bit_depth);
  }
  if (limit_cc_alf)
  {
    if (comp_id == COMPONENT_Cb)
    {
      count_luma_swing_greater_than_threshold(rec_dst_yuv->u, pic_stride_c, pic_height_c, pic_width_c,
        max_ctu_width_log2_chrma, max_ctu_height_log2_chrma, g_luma_swing_greater_than_threshold_count,
        ctus_in_width, bit_depth);

    }
    else if (comp_id == COMPONENT_Cr)
    {
      count_luma_swing_greater_than_threshold(rec_dst_yuv->v, pic_stride_c, pic_height_c, pic_width_c,
        max_ctu_width_log2_chrma, max_ctu_height_log2_chrma, g_luma_swing_greater_than_threshold_count,
        ctus_in_width, bit_depth);
    }
    else
    {
      assert(false); // Component ID not allowed.
    }
  }

  for ( int filter_idx = 0; filter_idx <= MAX_NUM_CC_ALF_FILTERS; filter_idx++ )
  {
    if ( filter_idx < MAX_NUM_CC_ALF_FILTERS)
    {
      memset(best_filter_coeff_set[filter_idx], 0, sizeof(best_filter_coeff_set[filter_idx]) );
      best_map_filter_idx_to_filter_idc[filter_idx] = filter_idx + 1;
    }
    else
    {
      best_map_filter_idx_to_filter_idc[filter_idx] = 0;
    }
  }
  
  memset(g_best_filter_control, 0, sizeof(uint8_t) * num_ctus_in_pic);
  int cc_alf_reuse_aps_id = -1;
  g_reuse_aps_id[comp_id - 1] = -1;

  //const TempCtx ctxStartCcAlfFilterControlFlag(m_CtxCache, SubCtx(Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx()));
  cabac_data_t ctx_start_cc_alf_filter_control_flag;
  memcpy(&ctx_start_cc_alf_filter_control_flag, &cabac_estimator, sizeof(ctx_start_cc_alf_filter_control_flag));
  ctx_start_cc_alf_filter_control_flag.only_count = 1;

  // compute cost of not filtering
  uint64_t unfiltered_distortion = 0;
  for (int ctb_idx = 0; ctb_idx < num_ctus_in_pic; ctb_idx++)
  {
    unfiltered_distortion += (uint64_t)g_alf_covariance_cc_alf[comp_id - 1][0][ctb_idx].pix_acc;
  }

  double best_unfiltered_total_cost = 1 * lambda + unfiltered_distortion;   // 1 bit is for gating flag

  bool cc_alf_filter_idx_enabled[MAX_NUM_CC_ALF_FILTERS];
  short cc_alf_filter_coeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
  uint8_t cc_alf_filter_count = MAX_NUM_CC_ALF_FILTERS;
  double best_filtered_total_cost = MAX_DOUBLE;
  bool best_reuse_temporal_filter_coeff = false;
  int    training_iter_count = 0;
  bool   keep_training = true;
  bool   improvement = false;
  double prev_total_cost = MAX_DOUBLE;
  const int num_coeff = CC_ALF_NUM_COEFF - 1;
  int log2_block_width = max_ctu_width_log2 - scale_x;
  int log2_block_height = max_ctu_height_log2 - scale_y;
  uint64_t cur_total_distortion = 0;
  double cur_total_rate = 0;
  int aps_ids_size = 0;
  int aps_ids[ALF_CTB_MAX_NUM_APS] = { -1, -1, -1, -1, -1, -1, -1, -1 };
  get_available_cc_alf_aps_ids(state, comp_id, &aps_ids_size, aps_ids);

  for (int test_filter_idx = 0; test_filter_idx < (aps_ids_size + 1 ); test_filter_idx++ )
  {
    bool referencing_existing_aps   = (test_filter_idx < aps_ids_size) ? true : false;
    int max_number_of_filters_being_tested = MAX_NUM_CC_ALF_FILTERS - (test_filter_idx - aps_ids_size);

    if (max_number_of_filters_being_tested < 0)
    {
      max_number_of_filters_being_tested = 1;
    }

    {
      // Instead of rewriting the control buffer for every training iteration just keep a mapping from filterIdx to filterIdc
      uint8_t map_filter_idx_to_filter_idc[MAX_NUM_CC_ALF_FILTERS + 1];
      for (int filter_idx = 0; filter_idx <= MAX_NUM_CC_ALF_FILTERS; filter_idx++)
      {
        if (filter_idx == MAX_NUM_CC_ALF_FILTERS)
        {
          map_filter_idx_to_filter_idc[filter_idx] = 0;
        }
        else
        {
          map_filter_idx_to_filter_idc[filter_idx] = filter_idx + 1;
        }
      }

      // initialize filters
      for ( int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++ )
      {
        cc_alf_filter_idx_enabled[filter_idx] = false;
        memset(cc_alf_filter_coeff[filter_idx], 0, sizeof(cc_alf_filter_coeff[filter_idx]));
      }
      if ( referencing_existing_aps )
      {
        max_number_of_filters_being_tested =
          state->encoder_control->cfg.param_set_map[aps_ids[test_filter_idx] + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.cc_alf_aps_param.cc_alf_filter_count[comp_id - 1];
        cc_alf_filter_count = max_number_of_filters_being_tested;
        for (int filter_idx = 0; filter_idx < max_number_of_filters_being_tested; filter_idx++)
        {
          cc_alf_filter_idx_enabled[filter_idx] = true;
          memcpy(cc_alf_filter_coeff[filter_idx], g_cc_alf_filter_param.cc_alf_coeff[comp_id - 1][filter_idx],
                 sizeof(cc_alf_filter_coeff[filter_idx]));
        }
        memcpy( cc_alf_filter_coeff, state->encoder_control->cfg.param_set_map[aps_ids[test_filter_idx] + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.cc_alf_aps_param.cc_alf_coeff[comp_id - 1], sizeof(cc_alf_filter_coeff) );
      }
      else
      {
        for (int i = 0; i < max_number_of_filters_being_tested; i++)
        {
          cc_alf_filter_idx_enabled[i] = true;
        }
        cc_alf_filter_count = max_number_of_filters_being_tested;
      }

      // initialize
      int control_idx = 0;
      const int column_size = (pic_width_c / max_number_of_filters_being_tested);
      for (int y = 0; y < pic_height_c; y += ctu_height_c)
      {
        for (int x = 0; x < pic_width_c; x += ctu_width_c)
        {
          g_training_cov_control[control_idx] = ( x / column_size ) + 1;
          control_idx++;
        }
      }

      // compute cost of filtering
      training_iter_count = 0;
      keep_training      = true;
      improvement       = false;
      prev_total_cost     = MAX_DOUBLE;
      while (keep_training)
      {
        improvement = false;
        for (int filter_idx = 0; filter_idx < max_number_of_filters_being_tested; filter_idx++)
        {
          if (cc_alf_filter_idx_enabled[filter_idx])
          {
            if (!referencing_existing_aps)
            {
              get_frame_stats_cc_alf(comp_id, (filter_idx + 1), state->lcu_order_count);
              derive_cc_alf_filter_coeff(comp_id, cc_alf_filter_coeff, filter_idx);
            }
            
            for (int y = 0; y < pic_height_c; y += (1 << log2_block_height))
            {
              for (int x = 0; x < pic_width_c; x += (1 << log2_block_width))
              {
                int ctu_idx = (y >> log2_block_height) * ctus_in_width + (x >> log2_block_width);
                g_training_distortion[filter_idx][ctu_idx] =
                  (int)(g_ctb_distortion_unfilter[comp_id][ctu_idx]
                      + calc_error_for_cc_alf_coeffs(&g_alf_covariance_cc_alf[comp_id - 1][0][ctu_idx],
                        cc_alf_filter_coeff[filter_idx], num_coeff, 7 + 1));
              }
            }
          }
        }

        memcpy(&cabac_estimator, &ctx_start_cc_alf_filter_control_flag, sizeof(cabac_estimator));

        cur_total_distortion = 0;
        cur_total_rate = 0;
        determine_control_idc_values(state, comp_id, ctu_width_c, ctu_height_c, pic_width_c, pic_height_c,
                                  g_ctb_distortion_unfilter, g_training_distortion,
                                  g_luma_swing_greater_than_threshold_count,
                                  g_chroma_sample_count_near_mid_point,
                                  (referencing_existing_aps == true),
                                  g_training_cov_control, g_filter_control, &cur_total_distortion, &cur_total_rate,
                                  cc_alf_filter_idx_enabled, map_filter_idx_to_filter_idc, &cc_alf_filter_count);

        // compute coefficient coding bit cost
        if (cc_alf_filter_count > 0)
        {
          if (referencing_existing_aps)
          {
            cur_total_rate += 1 + 3; // +1 for enable flag, +3 APS ID in slice header
          }
          else
          {
            cur_total_rate += get_coeff_rate_cc_alf(cc_alf_filter_coeff, cc_alf_filter_idx_enabled, cc_alf_filter_count, comp_id) + 1
            + 9;   // +1 for the enable flag, +9 3-bit for APS ID in slice header, 5-bit for APS ID in APS, a 1-bit
            // new filter flags (ignore shared cost such as other new-filter flags/NALU header/RBSP
            // terminating bit/byte alignment bits)*/
          }

          double cur_total_cost = cur_total_rate * lambda + cur_total_distortion;

          if (cur_total_cost < prev_total_cost)
          {
            prev_total_cost = cur_total_cost;
            improvement = true;
          }

          if (cur_total_cost < best_filtered_total_cost)
          {
            best_filtered_total_cost = cur_total_cost;
            memcpy(best_filter_idx_enabled, cc_alf_filter_idx_enabled, sizeof(cc_alf_filter_idx_enabled));
            memcpy(best_filter_coeff_set, cc_alf_filter_coeff, sizeof(cc_alf_filter_coeff));
            memcpy(g_best_filter_control, g_filter_control, sizeof(uint8_t) * num_ctus_in_pic);
            best_filter_count = cc_alf_filter_count;
            cc_alf_reuse_aps_id = referencing_existing_aps ? aps_ids[test_filter_idx] : -1;
            memcpy(best_map_filter_idx_to_filter_idc, map_filter_idx_to_filter_idc, sizeof(map_filter_idx_to_filter_idc));
          }
        }

        training_iter_count++;
        if (!improvement || training_iter_count > max_training_iter_count || referencing_existing_aps)
        {
          keep_training = false;
        }
      }
    }
  }

  if (best_unfiltered_total_cost < best_filtered_total_cost)
  {
    memset(g_best_filter_control, 0, sizeof(uint8_t) * num_ctus_in_pic);
  }

  // save best coeff and control
  bool atleast_one_block_undergoes_fitlering = false;
  for (int controlIdx = 0; best_filter_count > 0 && controlIdx < num_ctus_in_pic; controlIdx++)
  {
    if (g_best_filter_control[controlIdx])
    {
      atleast_one_block_undergoes_fitlering = true;
      break;
    }
  }
  g_cc_alf_filter_param.number_valid_components = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_COMPONENT;
  g_cc_alf_filter_param.cc_alf_filter_enabled[comp_id - 1] = atleast_one_block_undergoes_fitlering;
  if (atleast_one_block_undergoes_fitlering)
  {
    // update the filter control indicators
    if (best_reuse_temporal_filter_coeff!=1)
    {
      short stored_best_filter_coeff_set[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
      for (int filter_idx=0; filter_idx<MAX_NUM_CC_ALF_FILTERS; filter_idx++)
      {
        memcpy(stored_best_filter_coeff_set[filter_idx], best_filter_coeff_set[filter_idx], sizeof(best_filter_coeff_set[filter_idx]));
      }
      memcpy(g_filter_control, g_best_filter_control, sizeof(uint8_t) * num_ctus_in_pic);

      int filter_count = 0;
      for ( int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++ )
      {
        uint8_t cur_filter_idc = best_map_filter_idx_to_filter_idc[filter_idx];
        if (best_filter_idx_enabled[filter_idx])
        {
          for (int control_idx = 0; control_idx < num_ctus_in_pic; control_idx++)
          {
            if (g_filter_control[control_idx] == (filter_idx+1) )
            {
              g_best_filter_control[control_idx] = cur_filter_idc;
            }
          }
          memcpy( best_filter_coeff_set[cur_filter_idc-1], stored_best_filter_coeff_set[filter_idx], sizeof(stored_best_filter_coeff_set[filter_idx]) );
          filter_count++;
        }
        best_filter_idx_enabled[filter_idx] = ( filter_idx < best_filter_count ) ? true : false;
      }
      assert(filter_count == best_filter_count); //Number of filters enabled did not match the filter count
    }

    g_cc_alf_filter_param.cc_alf_filter_count[comp_id - 1] = best_filter_count;
    // cleanup before copying
    memset(g_cc_alf_filter_control[comp_id - 1], 0, sizeof(uint8_t) * num_ctus_in_pic);
    for ( int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++ )
    {
      memset(g_cc_alf_filter_param.cc_alf_coeff[comp_id - 1][filter_idx], 0,
             sizeof(g_cc_alf_filter_param.cc_alf_coeff[comp_id - 1][filter_idx]));
    }
    memset(g_cc_alf_filter_param.cc_alf_filter_idx_enabled[comp_id - 1], false,
           sizeof(g_cc_alf_filter_param.cc_alf_filter_idx_enabled[comp_id - 1]));
    for ( int filter_idx = 0; filter_idx < best_filter_count; filter_idx++ )
    {
      g_cc_alf_filter_param.cc_alf_filter_idx_enabled[comp_id - 1][filter_idx] = best_filter_idx_enabled[filter_idx];
      memcpy(g_cc_alf_filter_param.cc_alf_coeff[comp_id - 1][filter_idx], best_filter_coeff_set[filter_idx],
             sizeof(best_filter_coeff_set[filter_idx]));
    }
    memcpy(g_cc_alf_filter_control[comp_id - 1], g_best_filter_control, sizeof(uint8_t) * num_ctus_in_pic);
    if ( cc_alf_reuse_aps_id >= 0 )
    {
      g_reuse_aps_id[comp_id - 1] = cc_alf_reuse_aps_id;
      if (comp_id == COMPONENT_Cb)
      {
        state->slice->tile_group_cc_alf_cb_aps_id == cc_alf_reuse_aps_id;
      }
      else
      {
        state->slice->tile_group_cc_alf_cr_aps_id == cc_alf_reuse_aps_id;
      }
    }
  }
}


void derive_stats_for_cc_alf_filtering(encoder_state_t * const state,
  const kvz_picture *org_yuv,
  const int comp_idx, const int mask_stride,
  const uint8_t filter_idc)
{
  const int32_t num_ctus_in_pic = state->lcu_order_count;
  const int filter_idx = filter_idc - 1;

  // init CTU stats buffers
  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    reset_alf_covariance(&g_alf_covariance_cc_alf[comp_idx - 1][filter_idx][ctu_idx], -1);
  }

  // init Frame stats buffers
  reset_alf_covariance(&g_alf_covariance_frame_cc_alf[comp_idx - 1][filter_idx], -1);

  int ctu_rs_addr = 0;
  const int frame_height = state->tile->frame->height;
  const int frame_width = state->tile->frame->width;
  const int max_cu_width = LCU_WIDTH;
  const int max_cu_height = LCU_WIDTH;

  for (int y_pos = 0; y_pos < state->tile->frame->height; y_pos += max_cu_height)
  {
    for (int x_pos = 0; x_pos < state->tile->frame->width; x_pos += max_cu_width)
    {
      const int width = (x_pos + max_cu_width > frame_width) ? (frame_width - x_pos) : max_cu_width;
      const int height = (y_pos + max_cu_height > frame_height) ? (frame_height - y_pos) : max_cu_height;
      int       raster_slice_alf_pad = 0;
      get_blk_stats_cc_alf(state, &g_alf_covariance_cc_alf[comp_idx - 1][filter_idx][ctu_rs_addr],
                        org_yuv, comp_idx, x_pos, y_pos, width, height);
      add_alf_cov(&g_alf_covariance_frame_cc_alf[comp_idx - 1][filter_idx], &g_alf_covariance_cc_alf[comp_idx - 1][filter_idx][ctu_rs_addr]);
      ctu_rs_addr++;
    }
  }
}


void get_blk_stats_cc_alf(encoder_state_t * const state, 
  alf_covariance *alf_covariance,
  const kvz_picture *org_yuv,
  const alf_component_id comp_id,
  const int x_pos, const int y_pos,
  const int width, const int height)
{
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  bool chroma_scale_x = (chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool chroma_scale_y = (chroma_fmt != KVZ_CSP_420) ? 0 : 1;

  const int frame_height = state->tile->frame->height;
  const int alf_vb_luma_pos = LCU_WIDTH - ALF_VB_POS_ABOVE_CTUROW_LUMA;
  const int alf_vb_luma_ctu_height = LCU_WIDTH;
  const int max_cu_height = LCU_WIDTH;
  const int x_pos_c = x_pos >> chroma_scale_x;
  const int y_pos_c = y_pos >> chroma_scale_y;
  const int c_width = width >> chroma_scale_x;
  const int c_height = height >> chroma_scale_y;

  const int num_coeff = 8;
  const channel_type channel = (comp_id == COMPONENT_Y) ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;

  enum kvz_chroma_format chroma_format = state->encoder_control->chroma_format;
  const int number_of_components = (chroma_format == KVZ_CSP_400) ? 1 : MAX_NUM_COMPONENT;;
  int rec_stride[MAX_NUM_COMPONENT];
  int rec_pixel_idx[MAX_NUM_COMPONENT];
  const int luma_rec_pos = y_pos * state->tile->frame->rec->stride + x_pos;
  const int chroma_rec_pos = y_pos_c * (state->tile->frame->rec->stride >> chroma_scale_x) + x_pos_c;
  kvz_pixel *rec_y = &alf_tmp_y[luma_rec_pos];
  kvz_pixel *rec_u = &alf_tmp_u[chroma_rec_pos];
  kvz_pixel *rec_v = &alf_tmp_v[chroma_rec_pos];

  for (int c_idx = 0; c_idx < number_of_components; c_idx++)
  {
    bool is_luma = c_idx == COMPONENT_Y;
    rec_stride[c_idx] = state->tile->frame->rec->stride >> (is_luma ? 0 : chroma_scale_x);
    rec_pixel_idx[c_idx] = 0;
  }

  int org_stride = 0;
  const kvz_pixel *org = 0;
  if (comp_id == COMPONENT_Y)
  {
    org_stride = org_yuv->stride;
    org = &org_yuv->y[y_pos*org_stride + x_pos];
  }
  else if (comp_id == COMPONENT_Cb)
  {
    org_stride = org_yuv->stride >> chroma_scale_x;
    org = &org_yuv->u[y_pos_c*org_stride + x_pos_c];
  }
  else if (comp_id == COMPONENT_Cr)
  {
    org_stride = org_yuv->stride >> chroma_scale_x;
    org = &org_yuv->v[y_pos_c*org_stride + x_pos_c];
  }

  const int  num_bins = 1;
  int vb_ctu_height = alf_vb_luma_ctu_height;
  int vb_pos = alf_vb_luma_pos;
  if ((y_pos + max_cu_height) >= frame_height)
  {
    vb_pos = frame_height;
  }

  int32_t e_local[MAX_NUM_CC_ALF_CHROMA_COEFF][1];
  kvz_pixel *rec_pixels = (comp_id == COMPONENT_Y ? rec_y : (comp_id == COMPONENT_Cb ? rec_u : rec_v));
  uint8_t component_scale_y = (comp_id == COMPONENT_Y || chroma_format != KVZ_CSP_420) ? 0 : 1;
  uint8_t component_scale_x = (comp_id == COMPONENT_Y || chroma_format == KVZ_CSP_444) ? 0 : 1;
  int16_t y_local = 0;

  for (int i = 0; i < (comp_id == COMPONENT_Y ? height : c_height); i++)
  {
    int vb_distance = ((i << component_scale_y) % vb_ctu_height) - vb_pos;
    const bool skip_this_row = (component_scale_y == 0 && (vb_distance == 0 || vb_distance == 1));
    for (int j = 0; j < (comp_id == COMPONENT_Y ? width : c_width) && (!skip_this_row); j++)
    {
      memset(e_local, 0, sizeof(e_local));

      double weight = 1.0;
      if (0 /*g_alf_wssd*/)
      {
        //weight = m_lumaLevelToWeightPLUT[org[j]];
      }

      y_local = org[j] - rec_pixels[j + rec_pixel_idx[comp_id]];
      calc_covariance_cc_alf(e_local, rec_y + rec_pixel_idx[COMPONENT_Y] + (j << component_scale_x), rec_stride[COMPONENT_Y], vb_distance);

      for (int k = 0; k < (num_coeff - 1); k++)
      {
        for (int l = k; l < (num_coeff - 1); l++)
        {
          for (int b0 = 0; b0 < num_bins; b0++)
          {
            for (int b1 = 0; b1 < num_bins; b1++)
            {
              if (0 /*g_alf_wssd*/)
              {
                alf_covariance->ee[b0][b1][k][l] += weight * (e_local[k][b0] * (double)e_local[l][b1]);
              }
              else
              {
                alf_covariance->ee[b0][b1][k][l] += e_local[k][b0] * (double)e_local[l][b1];
              }
            }
          }
        }
        for (int b = 0; b < num_bins; b++)
        {
          if (0 /*g_alf_wssd*/)
          {
            alf_covariance->y[b][k] += weight * (e_local[k][b] * (double)y_local);
          }
          else
          {
            alf_covariance->y[b][k] += e_local[k][b] * (double)y_local;
          }
        }
      }
      if (0 /*g_alf_wssd*/)
      {
        alf_covariance->pix_acc += weight * (y_local * (double)y_local);
      }
      else
      {
        alf_covariance->pix_acc += y_local * (double)y_local;
      }
    }
    org += org_stride;
    for (int src_c_idx = 0; src_c_idx < number_of_components; src_c_idx++)
    {
      const channel_type c_channel = (src_c_idx == COMPONENT_Y) ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;
      if (c_channel == channel)
      {
        rec_pixel_idx[src_c_idx] += rec_stride[src_c_idx];
      }
      else
      {
        if (comp_id == COMPONENT_Y)
        {
          rec_pixel_idx[src_c_idx] += rec_stride[src_c_idx] >> ((src_c_idx == COMPONENT_Y || chroma_format != KVZ_CSP_420) ? 0 : 1);
        }
        else
        {
          rec_pixel_idx[src_c_idx] += rec_stride[src_c_idx] << ((comp_id == COMPONENT_Y || chroma_format != KVZ_CSP_420) ? 0 : 1);
        }
      }
    }
  }

  for (int k = 1; k < (MAX_NUM_CC_ALF_CHROMA_COEFF - 1); k++)
  {
    for (int l = 0; l < k; l++)
    {
      for (int b0 = 0; b0 < num_bins; b0++)
      {
        for (int b1 = 0; b1 < num_bins; b1++)
        {
          alf_covariance->ee[b0][b1][k][l] = alf_covariance->ee[b1][b0][l][k];
        }
      }
    }
  }
}

void calc_covariance_cc_alf(int32_t e_local[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const kvz_pixel *rec, const int stride, int vb_distance)
{
  const kvz_pixel *rec_y_m1 = rec - 1 * stride;
  const kvz_pixel *rec_y_0  = rec;
  const kvz_pixel *rec_y_p1 = rec + 1 * stride;
  const kvz_pixel *rec_y_p2 = rec + 2 * stride;

  if (vb_distance == -2 || vb_distance == +1)
  {
    rec_y_p2 = rec_y_p1;
  }
  else if (vb_distance == -1 || vb_distance == 0)
  {
    rec_y_m1 = rec_y_0;
    rec_y_p2 = rec_y_p1 = rec_y_0;
  }

  const kvz_pixel center_value = rec_y_0[+0];
  for (int b = 0; b < 1; b++)
  {
    e_local[0][b] += rec_y_m1[+0] - center_value;
    e_local[1][b] += rec_y_0[-1] - center_value;
    e_local[2][b] += rec_y_0[+1] - center_value;
    e_local[3][b] += rec_y_p1[-1] - center_value;
    e_local[4][b] += rec_y_p1[+0] - center_value;
    e_local[5][b] += rec_y_p1[+1] - center_value;
    e_local[6][b] += rec_y_p2[+0] - center_value;
  }
}

void count_luma_swing_greater_than_threshold(const kvz_pixel* luma, 
  int luma_stride, int height, int width,
  int log2_block_width, int log2_block_height, 
  uint64_t* luma_swing_greater_than_threshold_count, 
  int luma_count_stride,
  int8_t input_bit_depth)
{
  const int threshold = (1 << (input_bit_depth - 2 )) - 1;

  // 3x4 Diamond
  int x_support[] = {  0, -1, 0, 1, -1, 0, 1, 0 };
  int y_support[] = { -1,  0, 0, 0,  1, 1, 1, 2 };

  for (int y = 0; y < height; y += (1 << log2_block_height))
  {
    for (int x = 0; x < width; x += (1 << log2_block_width))
    {
      luma_swing_greater_than_threshold_count[(y >> log2_block_height) * luma_count_stride + (x >> log2_block_width)] = 0;

      for (int y_off = 0; y_off < (1 << log2_block_height); y_off++)
      {
        for (int x_off = 0; x_off < (1 << log2_block_width); x_off++)
        {
          if ((y + y_off) >= (height - 2) || (x + x_off) >= (width - 1) || (y + y_off) < 1 || (x + x_off) < 1) // only consider samples that are fully supported by picture
          {
            continue;
          }

          int min_val = ((1 << input_bit_depth) - 1);
          int max_val = 0;
          for (int i = 0; i < 8; i++)
          {
            kvz_pixel p = luma[(y_off + y_support[i]) * luma_stride + x + x_off + x_support[i]];

            if ( p < min_val )
            {
              min_val = p;
            }
            if ( p > max_val )
            {
              max_val = p;
            }
          }

          if ((max_val - min_val) > threshold)
          {
            luma_swing_greater_than_threshold_count[(y >> log2_block_height) * luma_count_stride + (x >> log2_block_width)]++;
          }
        }
      }
    }
    luma += (luma_stride << log2_block_height);
  }
}

/*
void EncAdaptiveLoopFilter::countChromaSampleValueNearMidPoint(const Pel* chroma, int chromaStride, int height, int width, int log2BlockWidth, int log2BlockHeight, uint64_t* chromaSampleCountNearMidPoint, int chromaSampleCountNearMidPointStride)
{
  const int midPoint  = (1 << m_inputBitDepth[CH_C]) >> 1;
  const int threshold = 16;

  for (int y = 0; y < height; y += (1 << log2BlockHeight))
  {
    for (int x = 0; x < width; x += (1 << log2BlockWidth))
    {
      chromaSampleCountNearMidPoint[(y >> log2BlockHeight)* chromaSampleCountNearMidPointStride + (x >> log2BlockWidth)] = 0;

      for (int yOff = 0; yOff < (1 << log2BlockHeight); yOff++)
      {
        for (int xOff = 0; xOff < (1 << log2BlockWidth); xOff++)
        {
          if ((y + yOff) >= height || (x + xOff) >= width)
          {
            continue;
          }

          int distanceToMidPoint = abs(chroma[yOff * chromaStride + x + xOff] - midPoint);
          if (distanceToMidPoint < threshold)
          {
            chromaSampleCountNearMidPoint[(y >> log2BlockHeight)* chromaSampleCountNearMidPointStride + (x >> log2BlockWidth)]++;
          }
        }
      }
    }
    chroma += (chromaStride << log2BlockHeight);
  }
}




*/

void init_distortion_cc_alf(const int num_ctus)
{
  for (int comp = 1; comp < MAX_NUM_COMPONENT; comp++)
  {
    for (int ctb_idx = 0; ctb_idx < num_ctus; ctb_idx++)
    {
      g_ctb_distortion_unfilter[comp][ctb_idx] = g_alf_covariance_cc_alf[comp - 1][0][ctb_idx].pix_acc;
    }
  }
}

void get_frame_stats_cc_alf(alf_component_id comp_idx, int filter_idc, const int num_ctus_in_frame)
{

  const int filter_idx = filter_idc - 1;

  // init Frame stats buffers
  reset_alf_covariance(&g_alf_covariance_frame_cc_alf[comp_idx - 1][filter_idx], -1);

  for (int ctu_rs_addr = 0; ctu_rs_addr < num_ctus_in_frame; ctu_rs_addr++)
  {
    if (g_training_cov_control[ctu_rs_addr] == filter_idc)
    {
      add_alf_cov(&g_alf_covariance_frame_cc_alf[comp_idx - 1][filter_idx],
        &g_alf_covariance_cc_alf[comp_idx - 1][0][ctu_rs_addr]);
    }
  }
}

void filter_blk_cc_alf(encoder_state_t * const state,
  const kvz_pixel *dst_buf, const kvz_pixel *rec_src,
  const int rec_luma_stride,
  const alf_component_id comp_id, const int16_t *filter_coeff,
  const clp_rngs clp_rngs, int vb_ctu_height, int vb_pos,
  const int x_pos, const int y_pos,
  const int blk_width,
  const int blk_height)
{
  
  assert(!(1 << kvz_math_floor_log2(vb_ctu_height) != vb_ctu_height)); //Not a power of 2

  assert(comp_id != COMPONENT_Y); //Must be chroma

  enum kvz_chroma_format chroma_format = state->encoder_control->chroma_format;
  uint8_t scale_y = (comp_id == COMPONENT_Y || chroma_format != KVZ_CSP_420) ? 0 : 1;
  uint8_t scale_x = (comp_id == COMPONENT_Y || chroma_format == KVZ_CSP_444) ? 0 : 1;
  const int cls_size_y = 4;
  const int cls_size_x = 4;
  const int start_height = y_pos;
  const int end_height = y_pos + blk_height;
  const int start_width = x_pos;
  const int end_width = x_pos + blk_width;
  const int luma_start_height = start_height << scale_y;
  const int luma_start_width = start_width << scale_x;

  assert(!(start_height % cls_size_y)); //Wrong start_height in filtering
  assert(!(start_width % cls_size_x)); //Wrong start_width in filtering
  assert(!((end_height - start_height) % cls_size_y)); //Wrong end_height in filtering
  assert(!((end_width - start_width) % cls_size_x)); //Wrong end_width in filtering

  kvz_pixel* src_buf = rec_src;
  const kvz_pixel* luma_ptr = src_buf + luma_start_height * rec_luma_stride + luma_start_width;

  const int chroma_stride = rec_luma_stride >> scale_x;
  kvz_pixel* chroma_ptr = dst_buf + start_height * chroma_stride + start_width;

  for (int i = 0; i < end_height - start_height; i += cls_size_y)
  {
    for (int j = 0; j < end_width - start_width; j += cls_size_x)
    {
      for (int ii = 0; ii < cls_size_y; ii++)
      {
        int row = ii;
        int col = j;
        kvz_pixel *src_self = chroma_ptr + col + row * chroma_stride;

        int offset1 = rec_luma_stride;
        int offset2 = -rec_luma_stride;
        int offset3 = 2 * rec_luma_stride;
        row <<= scale_y;
        col <<= scale_x;
        const kvz_pixel *src_cross = luma_ptr + col + row * rec_luma_stride;

        int pos = ((start_height + i + ii) << scale_y) & (vb_ctu_height - 1);
        if (scale_y == 0 && (pos == vb_pos || pos == vb_pos + 1))
        {
          continue;
        }
        if (pos == (vb_pos - 2) || pos == (vb_pos + 1))
        {
          offset3 = offset1;
        }
        else if (pos == (vb_pos - 1) || pos == vb_pos)
        {
          offset1 = 0;
          offset2 = 0;
          offset3 = 0;
        }

        for (int jj = 0; jj < cls_size_x; jj++)
        {
          const int jj2 = (jj << scale_x);
          const int offset0 = 0;

          int sum = 0;
          const kvz_pixel curr_src_cross = src_cross[offset0 + jj2];
          sum += filter_coeff[0] * (src_cross[offset2 + jj2] - curr_src_cross);
          sum += filter_coeff[1] * (src_cross[offset0 + jj2 - 1] - curr_src_cross);
          sum += filter_coeff[2] * (src_cross[offset0 + jj2 + 1] - curr_src_cross);
          sum += filter_coeff[3] * (src_cross[offset1 + jj2 - 1] - curr_src_cross);
          sum += filter_coeff[4] * (src_cross[offset1 + jj2] - curr_src_cross);
          sum += filter_coeff[5] * (src_cross[offset1 + jj2 + 1] - curr_src_cross);
          sum += filter_coeff[6] * (src_cross[offset3 + jj2] - curr_src_cross);

          sum = (sum + ((1 << 7/*m_scaleBits*/) >> 1)) >> 7/*m_scaleBits*/;
          const int offset = 1 << clp_rngs.comp[comp_id].bd >> 1;
          sum = alf_clip_pixel(sum + offset, clp_rngs.comp[comp_id]) - offset;
          sum += src_self[jj];
          src_self[jj] = alf_clip_pixel(sum, clp_rngs.comp[comp_id]);
        }
      }
    }

    chroma_ptr += chroma_stride * cls_size_y;

    luma_ptr += rec_luma_stride * cls_size_y << scale_y;
  }
}