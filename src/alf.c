

#include "alf.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cabac.h"
#include "rdo.h"
#include "strategies/strategies-sao.h"
#include "kvz_math.h"
#include "reshape.h"

#if MAX_NUM_CC_ALF_FILTERS>1
typedef struct filter_idx_count
{
  uint64_t count;
  uint8_t filter_idx;
} filter_idx_count;

static int comparator(const void *v1, const void *v2)
{
  const filter_idx_count *p1 = (filter_idx_count *)v1;
  const filter_idx_count *p2 = (filter_idx_count *)v2;
  return (p1->count < p2->count);
}
#endif

static void reset_alf_covariance(alf_covariance *alf, int num_bins) {
  if (num_bins > 0) { alf->num_bins = num_bins; }
  alf->pix_acc = 0;
  memset(alf->y, 0, sizeof(alf->y));
  memset(alf->ee, 0, sizeof(alf->ee));
}

void kvz_reset_cc_alf_aps_param(cc_alf_filter_param *cc_alf) {
  memset(cc_alf->cc_alf_filter_enabled, false, sizeof(cc_alf->cc_alf_filter_enabled));
  memset(cc_alf->cc_alf_filter_idx_enabled, false, sizeof(cc_alf->cc_alf_filter_idx_enabled));
  memset(cc_alf->cc_alf_coeff, 0, sizeof(cc_alf->cc_alf_coeff));
  cc_alf->cc_alf_filter_count[0] = cc_alf->cc_alf_filter_count[1] = MAX_NUM_CC_ALF_FILTERS;
  cc_alf->number_valid_components = 3;
  cc_alf->new_cc_alf_filter[0] = cc_alf->new_cc_alf_filter[1] = 0;
}

static void reset_alf_param(alf_aps *src)
{
  memset(src->enabled_flag, false, sizeof(src->enabled_flag));
  memset(src->non_linear_flag, false, sizeof(src->non_linear_flag));
  memset(src->luma_coeff, 0, sizeof(src->luma_coeff));
  memset(src->luma_clipp, 0, sizeof(src->luma_clipp));
  src->num_alternatives_chroma = 1;
  memset(src->chroma_coeff, 0, sizeof(src->chroma_coeff));
  memset(src->chroma_clipp, 0, sizeof(src->chroma_clipp));
  memset(src->filter_coeff_delta_idx, 0, sizeof(src->filter_coeff_delta_idx));
  memset(src->alf_luma_coeff_flag, true, sizeof(src->alf_luma_coeff_flag));
  src->num_luma_filters = 1;
  src->alf_luma_coeff_delta_flag = false;
  memset(src->new_filter_flag, 0, sizeof(src->new_filter_flag));
}

static void reset_aps(alf_aps *src, bool cc_alf_enabled)
{
  src->aps_type = 0;
  src->aps_id = -1;
  src->temporal_id = 0;
  src->layer_id = 0;
  reset_alf_param(src);
  if (cc_alf_enabled) {
    kvz_reset_cc_alf_aps_param(&src->cc_alf_aps_param);
  }
}

void kvz_set_aps_map(kvz_config *const cfg)
{
  cfg->param_set_map = malloc(ALF_CTB_MAX_NUM_APS * sizeof(param_set_map));
  for (int aps_idx = 0; aps_idx < ALF_CTB_MAX_NUM_APS; aps_idx++) {
    cfg->param_set_map[aps_idx + T_ALF_APS].b_changed = false;
    reset_aps(&cfg->param_set_map[aps_idx + T_ALF_APS].parameter_set, cfg->alf_type == KVZ_ALF_FULL);
  }
}

static void init_ctu_alternative_chroma(const alf_aps *alf_param, uint8_t* ctu_alts[MAX_NUM_COMPONENT], const int32_t num_ctus)
{
  uint8_t alt_idx = 0;
  for (int ctu_idx = 0; ctu_idx < num_ctus; ++ctu_idx)
  {
    ctu_alts[COMPONENT_Cb][ctu_idx] = alt_idx;
    if ((ctu_idx + 1) * alf_param->num_alternatives_chroma >= (alt_idx + 1) * num_ctus)
      ++alt_idx;
  }
  alt_idx = 0;
  for (int ctu_idx = 0; ctu_idx < num_ctus; ++ctu_idx)
  {
    ctu_alts[COMPONENT_Cr][ctu_idx] = alt_idx;
    if ((ctu_idx + 1) * alf_param->num_alternatives_chroma >= (alt_idx + 1) * num_ctus)
      ++alt_idx;
  }
}

static int16_t alf_clip3(const int16_t minVal, const int16_t maxVal, const int16_t a)
{
  return MIN(MAX(minVal, a), maxVal);
}

static int16_t clip_alf(const int16_t clip, const int16_t ref, const int16_t val0, const int16_t val1)
{
  return alf_clip3(-clip, +clip, val0 - ref) + alf_clip3(-clip, +clip, val1 - ref);
}

static int alf_clip_pixel(const int a, const clp_rng clp_rng)
{
  return MIN(MAX(clp_rng.min, a), clp_rng.max);
}

static int alf_clip3_int(const int minVal, const int maxVal, const int a)
{
  return MIN(MAX(minVal, a), maxVal);
}

static void get_clip_max(const alf_covariance *cov, int *clip_max)
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

static void reduce_clip_cost(const alf_covariance *cov, int *clip)
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

static void set_ey_from_clip(const alf_covariance *cov, const int* clip, double ee[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double y[MAX_NUM_ALF_LUMA_COEFF], int size)
{
  for (int k = 0; k < size; k++)
  {
    y[k] = cov->y[clip[k]][k];
    for (int l = 0; l < size; l++)
    {
      ee[k][l] = cov->ee[clip[k]][clip[l]][k][l];
    }
  }
}

static int gns_cholesky_dec(double inp_matr[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double out_matr[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], int num_eq)
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
        if (scale <= 0.0000001) // if(scale <= 0 )  /* If inpMatr is singular */
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

static void gns_transpose_backsubstitution(double u[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double* rhs, double* x, int order)
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

static void gns_backsubstitution(double r[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double* z, int size, double* a)
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

static int gns_solve_by_chol(double lhs[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double rhs[MAX_NUM_ALF_LUMA_COEFF], double *x, int num_eq)
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
      lhs[i][i] += 0.0001;
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

static int gns_solve_by_chol_clip_gns(alf_covariance *cov, const int *clip, double *x, int num_eq)
{
  double lhs[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];
  double rhs[MAX_NUM_ALF_LUMA_COEFF];

  set_ey_from_clip(cov, clip, lhs, rhs, num_eq);
  return gns_solve_by_chol(lhs, rhs, x, num_eq);
}

static double calculate_error(const alf_covariance *cov, const int *clip, const double *coeff)
{
  double sum = 0;
  for (int i = 0; i < cov->num_coeff; i++)
  {
    sum += coeff[i] * cov->y[clip[i]][i];
  }

  return cov->pix_acc - sum;
}

static double optimize_filter(const alf_covariance *cov, int* clip, double *f, bool optimize_clip)
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
    for (int k = 0; k < size; ++k)
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
      for (int k = 0; k < size; ++k)
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

static double optimize_filter_clip(alf_covariance *cov, int* clip)
{
  double f[MAX_NUM_ALF_LUMA_COEFF];
  return optimize_filter(cov, clip, f, true);
}

static double optimize_filter_gns_calc(alf_covariance *cov, const int* clip, double *f, int size)
{
  gns_solve_by_chol_clip_gns(cov, clip, f, size);
  return calculate_error(cov, clip, f);
}

static double calc_error_for_coeffs(const alf_covariance *cov, const int *clip, const int *coeff, const int num_coeff, const int bit_depth)
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

static double calc_error_for_cc_alf_coeffs(const alf_covariance *cov, const int16_t* coeff, const int num_coeff, const int bit_depth)
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
    error += ((cov->ee[0][0][i][i] * coeff[i] + sum * 2) / factor - 2 * cov->y[0][i]) * coeff[i];
  }

  return error / factor;
}

static int length_uvlc(int ui_code)
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

static double get_dist_coeff_force_0(bool* coded_var_bins, double error_force_0_coeff_tab[MAX_NUM_ALF_CLASSES][2], int* bits_var_bin, int zero_bits_var_bin, const int num_filters, double lambda)
{
  double dist_force_0 = 0;
  memset(coded_var_bins, 0, sizeof(*coded_var_bins) * MAX_NUM_ALF_CLASSES);

  for (int filt_idx = 0; filt_idx < num_filters; filt_idx++)
  {
    double cost_diff = (error_force_0_coeff_tab[filt_idx][0] + lambda * zero_bits_var_bin) - (error_force_0_coeff_tab[filt_idx][1] + lambda * bits_var_bin[filt_idx]);
    coded_var_bins[filt_idx] = cost_diff > 0 ? true : false;
    dist_force_0 += error_force_0_coeff_tab[filt_idx][coded_var_bins[filt_idx] ? 1 : 0];
  }

  return dist_force_0;
}

static double get_dist_force_0(const alf_aps *alf_param, channel_type channel, const int num_filters, double error_tab_force_0_coeff[MAX_NUM_ALF_CLASSES][2], bool* coded_var_bins, double lambda,
  int filter_coeff_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  int filter_clipp_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF]
)
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int bits_var_bin[MAX_NUM_ALF_CLASSES];

  for (int ind = 0; ind < num_filters; ++ind)
  {
    bits_var_bin[ind] = 0;
    for (int i = 0; i < num_coeff - 1; i++)
    {
      bits_var_bin[ind] += length_uvlc(abs(filter_coeff_set[ind][i]));
      if (abs(filter_coeff_set[ind][i]) != 0)
        bits_var_bin[ind] += 1;
    }
  }

  int zero_bits_var_bin = 0;
  for (int i = 0; i < num_coeff - 1; i++)
  {
    zero_bits_var_bin += length_uvlc(0);
  }
  if (alf_param->non_linear_flag[CHANNEL_TYPE_LUMA])
  {
    for (int ind = 0; ind < num_filters; ++ind)
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        if (!abs(filter_coeff_set[ind][i]))
        {
          filter_clipp_set[ind][i] = 0;
        }
      }
    }
  }
  double dist_force_0 = get_dist_coeff_force_0(coded_var_bins, error_tab_force_0_coeff, bits_var_bin, zero_bits_var_bin, num_filters, lambda);
  return dist_force_0;
}

static int get_cost_filter_coeff_force_0(const alf_aps *alf_param, channel_type channel, const int num_filters, bool* coded_var_bins,
  int p_diff_q_filter_coeff_int_pp[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  int filter_clipp_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF]
)
{
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int len = 0;

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
    else
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        len += length_uvlc(0); // alf_coeff_luma_delta[i][j]
      }
    }
  }

  if (alf_param->non_linear_flag[CHANNEL_TYPE_LUMA])
  {
    for (int ind = 0; ind < num_filters; ++ind)
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        if (!abs(p_diff_q_filter_coeff_int_pp[ind][i]))
        {
          filter_clipp_set[ind][i] = 0;
        }
        len += 2;
      }
    }
  }

  return len;
}

static int length_filter_coeffs(channel_type channel, const int num_filters, int filter_coeff[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF])
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

static int get_cost_filter_coeff(channel_type channel, const int num_filters, int p_diff_q_filter_coeff_int_pp[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF])
{
  return length_filter_coeffs(channel, num_filters, p_diff_q_filter_coeff_int_pp);  // alf_coeff_luma_delta[i][j];
}

static int get_cost_filter_clipp(channel_type channel, const int num_filters,
  int p_diff_q_filter_coeff_int_pp[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  int filter_clipp_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF])
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  for (int filter_idx = 0; filter_idx < num_filters; ++filter_idx)
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      if (!abs(p_diff_q_filter_coeff_int_pp[filter_idx][i]))
      {
        filter_clipp_set[filter_idx][i] = 0;
      }
    }
  }
  return (num_filters * (num_coeff - 1)) << 1;
}

static int get_non_filter_coeff_rate(alf_aps *aps)
{
  int len = 0                                        // alf_coefficients_delta_flag
    + 2                                        // slice_alf_chroma_idc                     u(2)
    + length_uvlc(aps->num_luma_filters - 1);  // alf_luma_num_filters_signalled_minus1   ue(v)

  if (aps->num_luma_filters > 1)
  {
    const int coeff_length = kvz_math_ceil_log2(aps->num_luma_filters); //#if JVET_O0491_HLS_CLEANUP
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
    {
      len += coeff_length;
    }
  }

  return len;
}

static double calculate_error_opt_filt(alf_covariance *cov, const int *clip)
{
  double c[MAX_NUM_ALF_LUMA_COEFF];
  return optimize_filter_gns_calc(cov, clip, c, cov->num_coeff);
}

static int get_chroma_coeff_rate(alf_aps* aps, int alt_idx)
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
  if (aps->non_linear_flag[CHANNEL_TYPE_CHROMA])
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

static double get_filtered_distortion(alf_covariance* cov, array_variables *arr_vars, const int num_classes, const int num_filters_minus1, const int num_coeff, const int bit_depth)
{
  double dist = 0;

  for (int class_idx = 0; class_idx < num_classes; class_idx++)
  {
    dist += calc_error_for_coeffs(&cov[class_idx], arr_vars->filter_clipp_set[class_idx], arr_vars->filter_coeff_set[class_idx], num_coeff, bit_depth);
  }

  return dist;
}

static double get_unfiltered_distortion_cov_classes(alf_covariance* cov, const int num_classes)
{
  double dist = 0;
  for (int class_idx = 0; class_idx < num_classes; class_idx++)
  {
    dist += cov[class_idx].pix_acc;
  }
  return dist;
}

static double get_unfiltered_distortion_cov_channel(alf_covariance* cov, channel_type channel)
{
  double dist = 0;
  if (channel == CHANNEL_TYPE_LUMA)
  {
    dist = get_unfiltered_distortion_cov_classes(cov, MAX_NUM_ALF_CLASSES);
  }
  else
  {
    dist = get_unfiltered_distortion_cov_classes(cov, 1);
  }
  return dist;
}

static void add_alf_cov(alf_covariance *dst, alf_covariance *src)
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

static void add_alf_cov_lhs_rhs(alf_covariance *dst, alf_covariance *lhs, alf_covariance *rhs)
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

static void get_frame_stat(alf_covariance* frame_cov, alf_covariance* ctb_cov, bool* ctb_enable_flags, uint8_t* ctb_alt_idx, const int num_classes, int alt_idx, const int32_t num_ctus)
{
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
          add_alf_cov(&frame_cov[is_luma ? class_idx : alt_idx], &ctb_cov[(ctu_idx * num_classes) + class_idx]);
        }
      }
    }
  }
}

static void get_frame_stats(alf_info_t *alf_info, channel_type channel, const int32_t num_ctus)
{
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? true : false;
  int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
  int num_alternatives = is_luma ? 1 : alf_info->alf_param_temp.num_alternatives_chroma;
  // When calling this function m_ctuEnableFlag shall be set to 0 for CTUs using alternative APS
  // Here we compute frame stats for building new alternative filters
  for (int alt_idx = 0; alt_idx < num_alternatives; ++alt_idx)
  {
    for (int i = 0; i < num_classes; i++)
    {
      is_luma ? reset_alf_covariance(&alf_info->alf_covariance_frame_luma[i], MAX_ALF_NUM_CLIPPING_VALUES) :
        reset_alf_covariance(&alf_info->alf_covariance_frame_chroma[alt_idx], MAX_ALF_NUM_CLIPPING_VALUES);
    }
    if (is_luma)
    {
      get_frame_stat(alf_info->alf_covariance_frame_luma, alf_info->alf_covariance_y, alf_info->ctu_enable_flag[COMPONENT_Y], NULL, num_classes, alt_idx, num_ctus);
    }
    else
    {
      get_frame_stat(alf_info->alf_covariance_frame_chroma, alf_info->alf_covariance_u, alf_info->ctu_enable_flag[COMPONENT_Cb], alf_info->ctu_alternative[COMPONENT_Cb], num_classes, alt_idx, num_ctus);
      get_frame_stat(alf_info->alf_covariance_frame_chroma, alf_info->alf_covariance_v, alf_info->ctu_enable_flag[COMPONENT_Cr], alf_info->ctu_alternative[COMPONENT_Cr], num_classes, alt_idx, num_ctus);
    }
  }
}

static void copy_cov(alf_covariance *dst, alf_covariance *src)
{
  dst->num_coeff = src->num_coeff;
  dst->num_bins = src->num_bins;
  memcpy(&dst->ee, &src->ee, sizeof(dst->ee));
  memcpy(&dst->y, &src->y, sizeof(dst->y));
  dst->pix_acc = src->pix_acc;
}

static void copy_alf_param(alf_aps *dst, alf_aps *src)
{
  memcpy(dst->enabled_flag, src->enabled_flag, sizeof(dst->enabled_flag));
  memcpy(dst->non_linear_flag, src->non_linear_flag, sizeof(dst->non_linear_flag));
  memcpy(dst->luma_coeff, src->luma_coeff, sizeof(dst->luma_coeff));
  memcpy(dst->luma_clipp, src->luma_clipp, sizeof(dst->luma_clipp));
  dst->num_alternatives_chroma = src->num_alternatives_chroma;
  memcpy(dst->chroma_coeff, src->chroma_coeff, sizeof(dst->chroma_coeff));
  memcpy(dst->chroma_clipp, src->chroma_clipp, sizeof(dst->chroma_clipp));
  memcpy(dst->filter_coeff_delta_idx, src->filter_coeff_delta_idx, sizeof(dst->filter_coeff_delta_idx));
  memcpy(dst->alf_luma_coeff_flag, src->alf_luma_coeff_flag, sizeof(dst->alf_luma_coeff_flag));
  dst->num_luma_filters = src->num_luma_filters;
  dst->alf_luma_coeff_delta_flag = src->alf_luma_coeff_delta_flag;
  memcpy(dst->new_filter_flag, src->new_filter_flag, sizeof(dst->new_filter_flag));
}

static void copy_cc_alf_param(cc_alf_filter_param *dst, cc_alf_filter_param *src)
{
  memcpy(dst->cc_alf_filter_enabled, src->cc_alf_filter_enabled, sizeof(dst->cc_alf_filter_enabled));
  memcpy(dst->cc_alf_filter_idx_enabled, src->cc_alf_filter_idx_enabled, sizeof(dst->cc_alf_filter_idx_enabled));
  memcpy(dst->cc_alf_filter_count, src->cc_alf_filter_count, sizeof(dst->cc_alf_filter_count));
  memcpy(dst->cc_alf_coeff, src->cc_alf_coeff, sizeof(dst->cc_alf_coeff));
  memcpy(dst->new_cc_alf_filter, src->new_cc_alf_filter, sizeof(dst->new_cc_alf_filter));
  dst->number_valid_components = src->number_valid_components;
}

static void copy_alf_param_w_channel(alf_aps* dst, alf_aps* src, channel_type channel)
{
  if (channel == CHANNEL_TYPE_LUMA)
  {
    copy_alf_param(dst, src);
  }
  else
  {
    dst->enabled_flag[COMPONENT_Cb] = src->enabled_flag[COMPONENT_Cb];
    dst->enabled_flag[COMPONENT_Cr] = src->enabled_flag[COMPONENT_Cr];
    dst->num_alternatives_chroma = src->num_alternatives_chroma;
    dst->non_linear_flag[CHANNEL_TYPE_CHROMA] = src->non_linear_flag[CHANNEL_TYPE_CHROMA];
    memcpy(dst->chroma_coeff, src->chroma_coeff, sizeof(dst->chroma_coeff));
    memcpy(dst->chroma_clipp, src->chroma_clipp, sizeof(dst->chroma_clipp));
  }
}

static void copy_aps(alf_aps *dst, alf_aps *src, bool cc_alf_enabled)
{
  dst->aps_id = src->aps_id;
  dst->temporal_id = src->temporal_id;
  dst->layer_id = src->layer_id;
  dst->aps_type = src->aps_type;
  copy_alf_param(dst, src);
  if (cc_alf_enabled) {
    copy_cc_alf_param(&dst->cc_alf_aps_param, &src->cc_alf_aps_param);
  }
}

static void copy_aps_to_map(param_set_map *dst, alf_aps *src, int8_t aps_id, bool cc_alf_enabled)
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
    if (cc_alf_enabled) {
      copy_cc_alf_param(&dst[aps_id + T_ALF_APS].parameter_set.cc_alf_aps_param, &src->cc_alf_aps_param);
    }
  }
}

static void init_alf_covariance(alf_covariance *alf, int num_coeffs) {
  alf->num_coeff = num_coeffs;
  alf->num_bins = MAX_ALF_NUM_CLIPPING_VALUES;
  alf->pix_acc = 0;
  memset(alf->y, 0, sizeof(alf->y));
  memset(alf->ee, 0, sizeof(alf->ee));
}

static void copy_pixels(kvz_pixel *src, int x_src_start, int y_src_start, int src_stride,
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

static void adjust_pixels(kvz_pixel *src, int x_start, int x_end, int y_start, int y_end, int stride, int pic_width, int pic_height)
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

static void adjust_pixels_chroma(kvz_pixel *src, int x_start, int x_end, int y_start, int y_end, int stride, int pic_width, int pic_height)
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
  //left bottom corner
  if (bottom_left) {
    const int y_px = y_end - 1;
    for (int x = -2; x < 0; x++) {
      src[(2 + y_px) * stride + x] =
        src[(1 + y_px) * stride + x] = src[stride * y_px];
    }
  }
  //right bottom corner
  if (bottom_right) {
    const int x_px = x_end - 1;
    const int y_px = y_end - 1;
    for (int x = x_end; x < x_end + 2; x++) {
      src[(2 + y_px) * stride + x] =
        src[(1 + y_px) * stride + x] = src[stride * y_px + x_px];
    }
  }
}

static void set_ctu_enable_flag(bool **flags, channel_type channel, uint8_t value, const int32_t num_ctus)
{
  if (channel == CHANNEL_TYPE_LUMA) {
    memset(flags[COMPONENT_Y], value, sizeof(bool) * num_ctus);
  }
  else {
    memset(flags[COMPONENT_Cr], value, sizeof(bool) * num_ctus);
    memset(flags[COMPONENT_Cb], value, sizeof(bool) * num_ctus);
  }
}

static void copy_ctu_enable_flag(bool **flags_dst, bool **flags_src, channel_type channel, const int32_t num_ctus)
{
  if (channel == CHANNEL_TYPE_LUMA) {
    memcpy(flags_dst[COMPONENT_Y], flags_src[COMPONENT_Y], sizeof(bool) * num_ctus);
  }
  else {
    memcpy(flags_dst[COMPONENT_Cr], flags_src[COMPONENT_Cr], sizeof(bool) * num_ctus);
    memcpy(flags_dst[COMPONENT_Cb], flags_src[COMPONENT_Cb], sizeof(bool) * num_ctus);
  }
}

//-------------------------cabac writer functions------------------------

static void alf_cabac_reset_bits(cabac_data_t * const data)
{
  data->low = 0;
  data->bits_left = 23;
  data->num_buffered_bytes = 0;
  data->buffered_byte = 0xff;
}

static void code_alf_ctu_enable_flag(encoder_state_t * const state,
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

    bool left_avail = ctu_rs_addr % frame_width_in_ctus ? 1 : 0;
    bool above_avail = ctu_rs_addr/ frame_width_in_ctus ? 1 : 0;

    int left_ctu_addr = left_avail ? ctu_rs_addr - 1 : -1;
    int above_ctu_addr = above_avail ? ctu_rs_addr - frame_width_in_ctus : -1;

    bool* ctb_alf_flag = state->tile->frame->alf_info->ctu_enable_flag[component_id];

    int ctx = 0;
    ctx += left_ctu_addr > -1 ? (ctb_alf_flag[left_ctu_addr] ? 1 : 0) : 0;
    ctx += above_ctu_addr > -1 ? (ctb_alf_flag[above_ctu_addr] ? 1 : 0) : 0;

    cabac->cur_ctx = &(cabac->ctx.alf_ctb_flag_model[component_id * 3 + ctx]);
    CABAC_BIN(cabac, ctb_alf_flag[ctu_rs_addr], "alf_ctb_flag");
  }
}

static void code_alf_ctu_enable_flags_component(encoder_state_t * const state,
  cabac_data_t * const cabac,
  alf_component_id component_id,
  alf_aps *aps)
{
  const int32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    code_alf_ctu_enable_flag(state, cabac, ctu_idx, component_id, aps);
  }
}

static void code_alf_ctu_enable_flags_channel(encoder_state_t * const state,
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

static void code_alf_ctu_filter_index(encoder_state_t * const state,
  cabac_data_t * const cabac,
  uint32_t ctu_rs_addr,
  bool alf_enable_luma)
{
  const encoder_control_t * const encoder = state->encoder_control;

  if (!encoder->cfg.alf_type || !alf_enable_luma)//(!cs.sps->getALFEnabledFlag()) || (!alfEnableLuma))
  {
    return;
  }

  if (!state->tile->frame->alf_info->ctu_enable_flag[COMPONENT_Y][ctu_rs_addr])
  {
    return;
  }

  const unsigned filter_set_idx = state->tile->frame->alf_info->alf_ctb_filter_index[ctu_rs_addr];
  unsigned num_aps = state->slice->tile_group_num_aps;
  unsigned num_available_filt_sets = num_aps + ALF_NUM_FIXED_FILTER_SETS;
  if (num_available_filt_sets > ALF_NUM_FIXED_FILTER_SETS)
  {
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
  }
  else
  {
    assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //Fixed set numavail < num_fixed
    kvz_cabac_encode_trunc_bin(cabac, filter_set_idx, ALF_NUM_FIXED_FILTER_SETS);
  }
}

static void code_alf_ctu_alternative_ctu(encoder_state_t * const state,
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
    bool* ctb_alf_flag = state->tile->frame->alf_info->ctu_enable_flag[comp_idx];

    if (ctb_alf_flag[ctu_rs_addr])
    {
      const int num_alts = alf_param_ref->num_alternatives_chroma;
      uint8_t* ctb_alf_alternative = state->tile->frame->alf_info->ctu_alternative[comp_idx];
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

static void code_alf_ctu_alternatives_component(encoder_state_t * const state,
  cabac_data_t * const cabac,
  alf_component_id comp_id,
  alf_aps* aps)
{
  if (comp_id == COMPONENT_Y)
    return;
  uint32_t num_ctus = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  bool* ctb_alf_flag = state->tile->frame->alf_info->ctu_enable_flag[comp_id];
  for (int ctu_idx = 0; ctu_idx < num_ctus; ctu_idx++)
  {
    if (ctb_alf_flag[ctu_idx])
    {
      code_alf_ctu_alternative_ctu(state, cabac, ctu_idx, comp_id, aps);
    }
  }
}

static void code_alf_ctu_alternatives_channel(encoder_state_t * const state,
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

static void code_cc_alf_filter_control_idc(encoder_state_t * const state,
  cabac_data_t * const cabac, uint8_t idc_val,
  const alf_component_id comp_id, const int ctu_idx,
  const uint8_t *filter_control_idc,
  const int filter_count)
{
  assert(!(idc_val > filter_count)); //Filter index is too large
  int width_in_lcu = state->tile->frame->width_in_lcu;

  bool left_avail = ctu_idx % width_in_lcu ? 1 : 0;
  bool above_avail = ctu_idx / width_in_lcu ? 1 : 0;
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
}

void kvz_encode_alf_bits(encoder_state_t * const state, const int ctu_idx)
{
  if (state->encoder_control->cfg.alf_type)
  {
    alf_info_t *alf_info = state->tile->frame->alf_info;
    cc_alf_filter_param *cc_filter_param = state->slice->cc_filter_param;
    bool **ctu_enable_flag = state->tile->frame->alf_info->ctu_enable_flag;
    for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
    {
      bool is_luma = comp_idx == COMPONENT_Y ? true : false;
      //Pitäisi poistaa//
      /*if (!is_luma)
      {
        state->slice->tile_group_alf_enabled_flag[comp_idx] = false;
      }*/
      //---------------//
      code_alf_ctu_enable_flag(state, &state->cabac, ctu_idx, comp_idx, NULL);
      if (is_luma)
      {
        if (ctu_enable_flag[comp_idx][ctu_idx])
        {
          //int num_aps = state->slice->tile_group_num_aps;
          //state->slice->tile_group_num_aps = 0;
          code_alf_ctu_filter_index(state, &state->cabac, ctu_idx, state->slice->tile_group_alf_enabled_flag[COMPONENT_Y]);
          //state->slice->tile_group_num_aps = num_aps;
        }
      }
      if (!is_luma)
      {
        bool* ctb_alf_flag = state->slice->tile_group_alf_enabled_flag[comp_idx] ? ctu_enable_flag[comp_idx] : NULL;
        if (ctb_alf_flag && ctb_alf_flag[ctu_idx])
        {
          code_alf_ctu_alternative_ctu(state, &state->cabac, ctu_idx, comp_idx, NULL);
        }
      }
    }

    if (state->encoder_control->cfg.alf_type == KVZ_ALF_FULL) {
      int num_components = state->encoder_control->chroma_format == KVZ_CSP_400 ? 1 : MAX_NUM_COMPONENT;
      for (int comp_idx = 1; comp_idx < num_components; comp_idx++) {
        if (cc_filter_param->cc_alf_filter_enabled[comp_idx - 1]) {
          const int filter_count = cc_filter_param->cc_alf_filter_count[comp_idx - 1];
          code_cc_alf_filter_control_idc(state, &state->cabac, alf_info->cc_alf_filter_control[comp_idx - 1][ctu_idx], comp_idx,
            ctu_idx, alf_info->cc_alf_filter_control[comp_idx - 1], filter_count);
        }
      }
    }
  }
}

static void encode_alf_aps_filter(encoder_state_t * const state,
  alf_aps* aps,
  const bool is_chroma,
  const int alt_idx)
{
  bitstream_t * const stream = &state->stream;
  const int num_coeff = is_chroma ? 7 : 13;
  const short* coeff = is_chroma ? aps->chroma_coeff[alt_idx] : aps->luma_coeff;
  const int16_t* clipp = is_chroma ? aps->chroma_clipp[alt_idx] : aps->luma_clipp;
  const int num_filters = is_chroma ? 1 : aps->num_luma_filters;

  // Filter coefficients
  for (int ind = 0; ind < num_filters; ++ind)
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      WRITE_UE(stream, abs(coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i]), is_chroma ? "alf_chroma_coeff_abs" : "alf_luma_coeff_abs"); //alf_coeff_chroma[i], alf_coeff_luma_delta[i][j]
      if (abs(coeff[ind* MAX_NUM_ALF_LUMA_COEFF + i]) != 0)
      {
        WRITE_U(stream, (coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] < 0) ? 1 : 0, 1, is_chroma ? "alf_chroma_coeff_sign" : "alf_luma_coeff_sign");
      }
    }
  }

  // Clipping values coding
  if (aps->non_linear_flag[is_chroma])
  {
    for (int ind = 0; ind < num_filters; ++ind)
    {
      for (int i = 0; i < num_coeff - 1; i++)
      {
        WRITE_U(stream, clipp[ind * MAX_NUM_ALF_LUMA_COEFF + i], 2, is_chroma ? "alf_chroma_clip_idx" : "alf_luma_clip_idx");
      }
    }
  }
}

static void encode_alf_aps_flags(encoder_state_t * const state,
  alf_aps* aps)
{
  bitstream_t * const stream = &state->stream;
  const bool cc_alf_enabled = state->encoder_control->cfg.alf_type == KVZ_ALF_FULL;

  WRITE_U(stream, aps->new_filter_flag[CHANNEL_TYPE_LUMA], 1, "alf_luma_new_filter");
  if (state->encoder_control->chroma_format != KVZ_CSP_400)
  {
    WRITE_U(stream, aps->new_filter_flag[CHANNEL_TYPE_CHROMA], 1, "alf_chroma_new_filter")
  }

  if (state->encoder_control->chroma_format != KVZ_CSP_400)
  {
    if (cc_alf_enabled) {
      WRITE_U(stream, aps->cc_alf_aps_param.new_cc_alf_filter[COMPONENT_Cb - 1], 1, "alf_cc_cb_filter_signal_flag");
      WRITE_U(stream, aps->cc_alf_aps_param.new_cc_alf_filter[COMPONENT_Cr - 1], 1, "alf_cc_cr_filter_signal_flag");
    }
    else {
      WRITE_U(stream, 0, 1, "alf_cc_cb_filter_signal_flag");
      WRITE_U(stream, 0, 1, "alf_cc_cr_filter_signal_flag");
    }
  }

  if (aps->new_filter_flag[CHANNEL_TYPE_LUMA])
  {
    WRITE_U(stream, aps->non_linear_flag[CHANNEL_TYPE_LUMA], 1, "alf_luma_clip");
    WRITE_UE(stream, aps->num_luma_filters - 1, "alf_luma_num_filters_signalled_minus1");
    if (aps->num_luma_filters > 1)
    {
      //const int length = ceilLog2(param.numLumaFilters);
      const int length = kvz_math_ceil_log2(aps->num_luma_filters);
      for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
      {
        WRITE_U(stream, aps->filter_coeff_delta_idx[i], length, "alf_luma_coeff_delta_idx");
      }
    }
    encode_alf_aps_filter(state, aps, false, 0);
  }
  if (aps->new_filter_flag[CHANNEL_TYPE_CHROMA])
  {
    WRITE_U(stream, aps->non_linear_flag[CHANNEL_TYPE_CHROMA], 1, "alf_nonlinear_enable_flag_chroma");
    if (MAX_NUM_ALF_ALTERNATIVES_CHROMA > 1)
    {
      //WRITE_UVLC(param.numAlternativesChroma - 1, "alf_chroma_num_alts_minus1");
      WRITE_UE(stream, aps->num_alternatives_chroma - 1, "alf_chroma_num_alts_minus1");
    }

    for (int alt_idx = 0; alt_idx < aps->num_alternatives_chroma; ++alt_idx)
    {
      encode_alf_aps_filter(state, aps, true, alt_idx);
    }
  }

  if (cc_alf_enabled) {
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
          int num_coeff = MAX_NUM_CC_ALF_CHROMA_COEFF; //CC_ALF_FILTER

          const short *coeff = aps->cc_alf_aps_param.cc_alf_coeff[cc_idx][filter_idx];
          // Filter coefficients
          for (int i = 0; i < num_coeff - 1; i++)
          {
            if (coeff[i] == 0)
            {
              WRITE_U(stream, 0, 3,
                cc_idx == 0 ? "alf_cc_cb_mapped_coeff_abs" : "alf_cc_cr_mapped_coeff_abs");
            }
            else
            {
              WRITE_U(stream, 1 + kvz_math_floor_log2(abs(coeff[i])), 3,
                cc_idx == 0 ? "alf_cc_cb_mapped_coeff_abs" : "alf_cc_cr_mapped_coeff_abs");
              WRITE_U(stream, coeff[i] < 0 ? 1 : 0, 1, cc_idx == 0 ? "alf_cc_cb_coeff_sign" : "alf_cc_cr_coeff_sign");
            }
          }
        }
      }
    }
  }
}


// ToDo: Fill in LMCS APS
static void encode_lmcs_aps(encoder_state_t* const state, lmcs_aps* aps)
{
  bitstream_t* const stream = &state->stream;
  //SliceReshapeInfo param = pcAPS->getReshaperAPSInfo();
  WRITE_UE(stream, 0/*param.reshaperModelMinBinIdx*/, "lmcs_min_bin_idx");
  WRITE_UE(stream, 16 - 1/*16 - 1 - param.reshaperModelMaxBinIdx*/, "lmcs_delta_max_bin_idx");

  WRITE_UE(stream, 7/*param.maxNbitsNeededDeltaCW - 1*/, "lmcs_delta_cw_prec_minus1");
  /*
  for (int i = param.reshaperModelMinBinIdx; i <= param.reshaperModelMaxBinIdx; i++)
  {
    int deltaCW = param.reshaperModelBinCWDelta[i];
    int signCW = (deltaCW < 0) ? 1 : 0;
    int absCW = (deltaCW < 0) ? (-deltaCW) : deltaCW;
    WRITE_CODE(absCW, param.maxNbitsNeededDeltaCW, "lmcs_delta_abs_cw[ i ]");
    if (absCW > 0)
    {
      WRITE_FLAG(signCW, "lmcs_delta_sign_cw_flag[ i ]");
    }
  }
  int deltaCRS = pcAPS->chromaPresentFlag ? param.chrResScalingOffset : 0;
  int signCRS = (deltaCRS < 0) ? 1 : 0;
  int absCRS = (deltaCRS < 0) ? (-deltaCRS) : deltaCRS;
  if (pcAPS->chromaPresentFlag)
  {
    WRITE_CODE(absCRS, 3, "lmcs_delta_abs_crs");
  }
  if (absCRS > 0)
  {
    WRITE_FLAG(signCRS, "lmcs_delta_sign_crs_flag");
  }
  */
}

static void encoder_state_write_adaptation_parameter_set(encoder_state_t * const state, alf_aps *aps)
{
#ifdef KVZ_DEBUG
  printf("=========== Adaptation Parameter Set  ===========\n");
#endif

  bitstream_t * const stream = &state->stream;

  WRITE_U(stream, (int)aps->aps_type, 3, "aps_params_type");
  WRITE_U(stream, aps->aps_id, 5, "adaptation_parameter_set_id");
  WRITE_U(stream, state->encoder_control->chroma_format != KVZ_CSP_400, 1, "aps_chroma_present_flag");

  if (aps->aps_type == T_ALF_APS)
  {
    encode_alf_aps_flags(state, aps);
  }
  else if (aps->aps_type == T_LMCS_APS)
  {
    //encode_lmcs_aps(state);
  }
  /*else if (aps->aps_type == T_SCALING_LIST_APS)
  {
    codeScalingListAps(pcAPS);
  }*/
  WRITE_U(stream, 0, 1, "aps_extension_flag"); //Implementation when this flag is equal to 1 should be added when it is needed. Currently in the spec we don't have case when this flag is equal to 1
  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encode_alf_aps(encoder_state_t * const state)
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

      if (write_aps)
      {
        kvz_nal_write(stream, NAL_UNIT_PREFIX_APS, 0, state->frame->first_nal);
        state->frame->first_nal = false;
        encoder_state_write_adaptation_parameter_set(state, &aps);

        aps_map[aps_id + T_ALF_APS].b_changed = false;
      }
    }
  }
}

void kvz_encode_alf_adaptive_parameter_set(encoder_state_t * const state)
{
  //send LMCS APS when LMCSModel is updated. It can be updated even current slice does not enable reshaper.
  //For example, in RA, update is on intra slice, but intra slice may not use reshaper
  //encode_alf_aps_lmcs(state);

  // only 1 SCALING LIST data for 1 picture
  //encode_alf_aps_scaling_list(state);

  encode_alf_aps(state);
}

//--------------------------------------------------------------------------

//-------------------------CC ALF encoding functions------------------------

static void filter_blk_cc_alf(encoder_state_t * const state,
  kvz_pixel *dst_buf, const kvz_pixel *rec_src,
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

  const kvz_pixel* src_buf = rec_src;
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

static void apply_cc_alf_filter(encoder_state_t * const state, alf_component_id comp_id, kvz_pixel *dst_buf,
  const kvz_pixel *rec_yuv_ext, const int luma_stride, uint8_t *filter_control,
  const short filter_set[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF],
  const int   selected_filter_idx,
  array_variables *arr_vars)
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

        {
          filter_blk_cc_alf(state, dst_buf, rec_yuv_ext, luma_stride, comp_id, filter_coeff, arr_vars->clp_rngs, alf_vb_luma_ctu_height,
            alf_vb_luma_pos, x_pos >> component_scale_x, y_pos >> component_scale_y,
            width >> component_scale_x, height >> component_scale_y);
        }
      }
      ctu_idx++;
    }
  }
}

static void setup_cc_alf_aps(encoder_state_t * const state,
  const int *cc_reuse_aps_id)
{
  cc_alf_filter_param *cc_filter_param = state->slice->cc_filter_param;
  if (cc_filter_param->cc_alf_filter_enabled[COMPONENT_Cb - 1])
  {
    int cc_alf_cb_aps_id = state->slice->tile_group_cc_alf_cb_aps_id;
    alf_aps *aps = &state->encoder_control->cfg.param_set_map[cc_alf_cb_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
    if (aps->aps_id >= 0 && aps->aps_id < ALF_CTB_MAX_NUM_APS)
    {
      //aps = m_apsMap->allocatePS((cc_alf_cb_aps_id << NUM_APS_TYPE_LEN) + ALF_APS);
      aps->temporal_id = 0; // cs.slice->getTLayer()
    }
    aps->cc_alf_aps_param.cc_alf_filter_enabled[COMPONENT_Cb - 1] = 1;
    aps->cc_alf_aps_param.cc_alf_filter_count[COMPONENT_Cb - 1] = cc_filter_param->cc_alf_filter_count[COMPONENT_Cb - 1];
    for (int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++)
    {
      aps->cc_alf_aps_param.cc_alf_filter_idx_enabled[COMPONENT_Cb - 1][filter_idx] =
        cc_filter_param->cc_alf_filter_idx_enabled[COMPONENT_Cb - 1][filter_idx];
      memcpy(aps->cc_alf_aps_param.cc_alf_coeff[COMPONENT_Cb - 1][filter_idx],
        cc_filter_param->cc_alf_coeff[COMPONENT_Cb - 1][filter_idx], sizeof(short) * MAX_NUM_CC_ALF_CHROMA_COEFF);
    }
    aps->aps_id = cc_alf_cb_aps_id;
    aps->aps_type = T_ALF_APS;
    if (cc_reuse_aps_id[COMPONENT_Cb - 1] < 0)
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
  if (cc_filter_param->cc_alf_filter_enabled[COMPONENT_Cr - 1])
  {
    int  cc_alf_cr_aps_id = state->slice->tile_group_cc_alf_cr_aps_id;
    alf_aps *aps = &state->encoder_control->cfg.param_set_map[cc_alf_cr_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
    if (aps->aps_id >= 0 && aps->aps_id < ALF_CTB_MAX_NUM_APS)
    {
      //aps = m_apsMap->allocatePS((cc_alf_cb_aps_id << NUM_APS_TYPE_LEN) + ALF_APS);
      aps->temporal_id = 0; // cs.slice->getTLayer()
    }
    aps->cc_alf_aps_param.cc_alf_filter_enabled[COMPONENT_Cr - 1] = 1;
    aps->cc_alf_aps_param.cc_alf_filter_count[COMPONENT_Cr - 1] = cc_filter_param->cc_alf_filter_count[COMPONENT_Cr - 1];
    for (int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++)
    {
      aps->cc_alf_aps_param.cc_alf_filter_idx_enabled[COMPONENT_Cr - 1][filter_idx] =
        cc_filter_param->cc_alf_filter_idx_enabled[COMPONENT_Cr - 1][filter_idx];
      memcpy(aps->cc_alf_aps_param.cc_alf_coeff[COMPONENT_Cr - 1][filter_idx],
        cc_filter_param->cc_alf_coeff[COMPONENT_Cr - 1][filter_idx], sizeof(short) * MAX_NUM_CC_ALF_CHROMA_COEFF);
    }
    aps->aps_id = cc_alf_cr_aps_id;
    aps->aps_type = T_ALF_APS;
    if (cc_reuse_aps_id[COMPONENT_Cr - 1] < 0)
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

static void round_filt_coeff_cc_alf(int16_t *filter_coeff_quant, double *filter_coeff, const int num_coeff, const int factor)
{
  for (int i = 0; i < num_coeff; i++)
  {
    int sign = filter_coeff[i] > 0 ? 1 : -1;
    double best_err = 128.0*128.0;
    int best_index = 0;
    for (int k = 0; k < CCALF_CANDS_COEFF_NR; k++)
    {
      double err = (filter_coeff[i] * sign * factor - cc_alf_small_tab[k]);
      err = err * err;
      if (err < best_err)
      {
        best_err = err;
        best_index = k;
      }
    }
    filter_coeff_quant[i] = cc_alf_small_tab[best_index] * sign;
  }
}

static int get_coeff_rate_cc_alf(short chroma_coeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF], bool filter_enabled[MAX_NUM_CC_ALF_FILTERS], uint8_t filter_count, alf_component_id comp_id)
{
  int bits = 0;

  if (filter_count > 0)
  {
    bits += length_uvlc(filter_count - 1);
    int signaled_filter_count = 0;
    for (int filterIdx = 0; filterIdx < MAX_NUM_CC_ALF_FILTERS; filterIdx++)
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

static void derive_cc_alf_filter_coeff(alf_covariance *alf_covariance_frame_cc_alf,
  short filter_coeff[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF],
  const uint8_t filter_idx)
{
  int forward_tab[CCALF_CANDS_COEFF_NR * 2 - 1] = { 0 };
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
    ky[k] = alf_covariance_frame_cc_alf[filter_idx].y[0][k];
    for (int l = 0; l < size; l++)
    {
      k_e[k][l] = alf_covariance_frame_cc_alf[filter_idx].ee[0][0][k][l];
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
  double err_ref = calc_error_for_cc_alf_coeffs(&alf_covariance_frame_cc_alf[filter_idx], filter_coeff_int, size, (7/*m_scaleBits*/ + 1));
  while (modified)
  {
    modified = 0;
    for (int i = 1; i > -2; i -= 2)
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
        if ((org_idx - delta < 0) || (org_idx - delta >= CCALF_CANDS_COEFF_NR * 2 - 1))
          continue;

        filter_coeff_int[k] = forward_tab[org_idx - delta];
        double error = calc_error_for_cc_alf_coeffs(&alf_covariance_frame_cc_alf[filter_idx], filter_coeff_int, size, (7/*m_scaleBits*/ + 1));
        if (error < err_min)
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

static void determine_control_idc_values(encoder_state_t *const state, const alf_component_id comp_id,
  const int ctu_width_c, const int ctu_height_c, const int pic_width_c,
  const int pic_height_c, double **unfiltered_distortion,
  uint64_t *training_distortion[MAX_NUM_CC_ALF_FILTERS],
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

  alf_info_t *alf_info = state->tile->frame->alf_info;
  cabac_data_t *cabac_estimator = &alf_info->cabac_estimator;
  cabac_data_t ctx_initial;
  cabac_data_t ctx_best;
  cabac_data_t ctx_start;
  memcpy(&ctx_initial, cabac_estimator, sizeof(ctx_initial));
  memcpy(&ctx_best, cabac_estimator, sizeof(ctx_best));
  cabac_estimator->only_count = 1;
  ctx_initial.only_count = 1;
  ctx_best.only_count = 1;
  ctx_start.only_count = 1;

  //enum kvz_chroma_format chroma_format = state->encoder_control->chroma_format;
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
      //const uint32_t threshold_s = MIN(pic_height_c - y_ctu, ctu_height_c) << (chroma_format != KVZ_CSP_420 ? 0 : 1);
      //const uint32_t number_of_chroma_samples = MIN(pic_height_c - y_ctu, ctu_height_c) * MIN(pic_height_c - x_ctu, ctu_width_c);
      //const uint32_t threshold_c = (number_of_chroma_samples >> 2);

      memcpy(cabac_estimator, &ctx_best, sizeof(*cabac_estimator));
      memcpy(&ctx_start, cabac_estimator, sizeof(ctx_start));

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

        memcpy(cabac_estimator, &ctx_start, sizeof(*cabac_estimator));
        alf_cabac_reset_bits(cabac_estimator);

        //const Position lumaPos = Position({ xCtu << getComponentScaleX(comp_id, cs.pcv->chrFormat),
        //  yCtu << getComponentScaleY(comp_id, cs.pcv->chrFormat) });
        code_cc_alf_filter_control_idc(state, cabac_estimator, filter_idc, comp_id, ctu_idx,
          filter_control, *cc_alf_filter_count);
        //rate = FRAC_BITS_SCALE * m_CABACEstimator->getEstFracBits();
        rate = (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3);
        cost = rate * lambda + ssd;

        bool limitation_exceeded = false;
        if (limit_cc_alf && filter_idx < MAX_NUM_CC_ALF_FILTERS)
        {
          assert(false); // should have returned from cc alf with limit_cc_alf == true
          //limitation_exceeded = limitation_exceeded || (luma_swing_greater_than_threshold_count[ctu_idx] >= threshold_s);
          //limitation_exceeded = limitation_exceeded || (chroma_sample_count_near_mid_point[ctu_idx] >= threshold_c);
        }
        if (cost < best_cost && !limitation_exceeded)
        {
          best_cost = cost;
          best_rate = rate;
          best_ssd = ssd;
          best_filter_idc = filter_idc;
          best_filter_idx = filter_idx;

          //ctx_best = SubCtx(Ctx::CcAlfFilterControlFlag, m_CABACEstimator->getCtx());
          memcpy(&ctx_best, cabac_estimator, sizeof(ctx_best));

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
    memcpy(cabac_estimator, &ctx_initial, sizeof(*cabac_estimator));
    //m_CABACEstimator->resetBits();
    alf_cabac_reset_bits(cabac_estimator);

    int ctu_idx = 0;
    for (int y = 0; y < pic_height_c; y += ctu_height_c)
    {
      for (int x = 0; x < pic_width_c; x += ctu_width_c)
      {
        const int filter_idx_plus1 = filter_control[ctu_idx];
        code_cc_alf_filter_control_idc(state, cabac_estimator, (filter_idx_plus1 == 0 ? 0
          : map_filter_idx_to_filter_idc[filter_idx_plus1 - 1]),
          comp_id, ctu_idx, filter_control, *cc_alf_filter_count);

        ctu_idx++;
      }
    }
    (*cur_total_rate) += (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3);
  }
#endif

  // restore for next iteration
  memcpy(cabac_estimator, &ctx_initial, sizeof(*cabac_estimator));
}

static void get_available_cc_alf_aps_ids(encoder_state_t *const state, alf_component_id compID,
  int *aps_ids_size, int *aps_ids)
{
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    param_set_map* param_set = &state->encoder_control->cfg.param_set_map[i + NUM_APS_TYPE_LEN + T_ALF_APS];
    if (param_set->b_changed && (param_set->parameter_set.aps_id >= 0 || param_set->parameter_set.aps_id < ALF_CTB_MAX_NUM_APS)) {
      copy_aps(&state->slice->apss[i], &param_set->parameter_set, true);
    }
  }

  int aps_id_checked = 0, cur_aps_id = state->tile->frame->alf_info->aps_id_start;
  if (cur_aps_id < ALF_CTB_MAX_NUM_APS)
  {
    while (aps_id_checked < ALF_CTB_MAX_NUM_APS &&
      !state->frame->is_irap &&
      (*aps_ids_size) < ALF_CTB_MAX_NUM_APS
      /*&& !cs.slice->getPendingRasInit()*/)
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

static void get_frame_stats_cc_alf(alf_covariance* alf_covariance_cc_alf,
  alf_covariance* alf_covariance_frame_cc_alf, int filter_idc, const int num_ctus_in_frame,
  uint8_t *training_cov_control)
{
  const int filter_idx = filter_idc - 1;

  // init Frame stats buffers
  reset_alf_covariance(&alf_covariance_frame_cc_alf[filter_idx], -1);

  for (int ctu_rs_addr = 0; ctu_rs_addr < num_ctus_in_frame; ctu_rs_addr++)
  {
    if (training_cov_control[ctu_rs_addr] == filter_idc)
    {
      add_alf_cov(&alf_covariance_frame_cc_alf[filter_idx],
        &alf_covariance_cc_alf[ctu_rs_addr]);
    }
  }
}

static void derive_cc_alf_filter(encoder_state_t * const state, alf_component_id comp_id,
  const kvz_picture *org_yuv, const kvz_picture *rec_dst_yuv,
  int *cc_reuse_aps_id)
{
  cc_alf_filter_param *cc_filter_param = state->slice->cc_filter_param;
  if (!state->slice->tile_group_alf_enabled_flag[COMPONENT_Y])
  {
    cc_filter_param->cc_alf_filter_enabled[comp_id - 1] = false;
    return;
  }

  bool limit_cc_alf = state->encoder_control->cfg.qp >= 37; // m_encCfg->getCCALFQpThreshold();
  if (limit_cc_alf) // && state->slice. cs.slice->getSliceQp() <= m_encCfg->getBaseQP() + 1)
  {
    cc_filter_param->cc_alf_filter_enabled[comp_id - 1] = false;
    return;
  }

  alf_info_t *alf_info = state->tile->frame->alf_info;
  alf_covariance *alf_covariance_cc_alf = state->tile->frame->alf_info->alf_covariance_cc_alf[comp_id - 1];
  alf_covariance *alf_covariance_frame_cc_alf = state->tile->frame->alf_info->alf_covariance_frame_cc_alf[comp_id - 1];
  uint8_t* training_cov_control = alf_info->training_cov_control;
  uint8_t* filter_control = alf_info->filter_control;
  uint8_t* best_filter_control = alf_info->best_filter_control;
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  uint8_t best_map_filter_idx_to_filter_idc[MAX_NUM_CC_ALF_FILTERS + 1];
  bool scale_x = (comp_id == COMPONENT_Y || chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool scale_y = (comp_id == COMPONENT_Y || chroma_fmt != KVZ_CSP_420) ? 0 : 1;
  const int ctu_width_c = LCU_WIDTH >> scale_x;
  const int ctu_height_c = LCU_WIDTH >> scale_y;
  const int pic_width_c = state->tile->frame->width >> scale_x;
  const int pic_height_c = state->tile->frame->height >> scale_y;
  //const int pic_stride_c = rec_dst_yuv->stride >> scale_x;
  //const int8_t bit_depth = state->encoder_control->bitdepth;
  const int max_training_iter_count = 15;
  int max_ctu_height_log2 = kvz_math_floor_log2(LCU_WIDTH);
  //int max_ctu_height_log2_chrma = kvz_math_floor_log2(LCU_WIDTH) >> scale_y;
  int max_ctu_width_log2 = kvz_math_floor_log2(LCU_WIDTH);
  //int max_ctu_width_log2_chrma = kvz_math_floor_log2(LCU_WIDTH) >> scale_x;
  int32_t ctus_in_width = state->tile->frame->width_in_lcu;
  const uint32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  short best_filter_coeff_set[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
  bool best_filter_idx_enabled[MAX_NUM_CC_ALF_FILTERS];
  uint8_t best_filter_count = 0;
  double lambda = state->frame->lambda;
  cabac_data_t *cabac_estimator = &alf_info->cabac_estimator;
  cabac_data_t ctx_start_cc_alf_filter_control_flag;

  //uint64_t* luma_swing_greater_than_threshold_count;
  //uint64_t* chroma_sample_count_near_mid_point;

  /*if (limit_cc_alf)
  {
    luma_swing_greater_than_threshold_count = malloc(num_ctus_in_pic * sizeof(*luma_swing_greater_than_threshold_count));

    count_luma_swing_greater_than_threshold(rec_dst_yuv->y, rec_dst_yuv->stride, rec_dst_yuv->height, rec_dst_yuv->width,
      max_ctu_width_log2, max_ctu_height_log2, luma_swing_greater_than_threshold_count,
      ctus_in_width, bit_depth);
  }
  if (limit_cc_alf)
  {
    chroma_sample_count_near_mid_point = malloc(num_ctus_in_pic * sizeof(*chroma_sample_count_near_mid_point));
    if (comp_id == COMPONENT_Cb)
    {
      count_chroma_sample_value_near_mid_point(rec_dst_yuv->u, pic_stride_c, pic_height_c, pic_width_c,
        max_ctu_width_log2_chrma, max_ctu_height_log2_chrma, chroma_sample_count_near_mid_point,
        ctus_in_width, bit_depth);

    }
    else if (comp_id == COMPONENT_Cr)
    {
      count_chroma_sample_value_near_mid_point(rec_dst_yuv->v, pic_stride_c, pic_height_c, pic_width_c,
        max_ctu_width_log2_chrma, max_ctu_height_log2_chrma, chroma_sample_count_near_mid_point,
        ctus_in_width, bit_depth);
    }
    else
    {
      assert(false); // Component ID not allowed.
    }
  }*/

  for (int filter_idx = 0; filter_idx <= MAX_NUM_CC_ALF_FILTERS; filter_idx++)
  {
    if (filter_idx < MAX_NUM_CC_ALF_FILTERS)
    {
      memset(best_filter_coeff_set[filter_idx], 0, sizeof(best_filter_coeff_set[filter_idx]));
      best_map_filter_idx_to_filter_idc[filter_idx] = filter_idx + 1;
    }
    else
    {
      best_map_filter_idx_to_filter_idc[filter_idx] = 0;
    }
  }

  memset(best_filter_control, 0, sizeof(uint8_t) * num_ctus_in_pic);
  int cc_alf_reuse_aps_id = -1;
  cc_reuse_aps_id[comp_id - 1] = -1;

  memcpy(&ctx_start_cc_alf_filter_control_flag, cabac_estimator, sizeof(ctx_start_cc_alf_filter_control_flag));
  ctx_start_cc_alf_filter_control_flag.only_count = 1;

  // compute cost of not filtering
  uint64_t unfiltered_distortion = 0;
  for (int ctb_idx = 0; ctb_idx < num_ctus_in_pic; ctb_idx++)
  {
    unfiltered_distortion += (uint64_t)alf_covariance_cc_alf[ctb_idx].pix_acc;
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

  for (int test_filter_idx = 0; test_filter_idx < (aps_ids_size + 1); test_filter_idx++)
  {
    bool referencing_existing_aps = (test_filter_idx < aps_ids_size) ? true : false;
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
      for (int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++)
      {
        cc_alf_filter_idx_enabled[filter_idx] = false;
        memset(cc_alf_filter_coeff[filter_idx], 0, sizeof(cc_alf_filter_coeff[filter_idx]));
      }
      if (referencing_existing_aps)
      {
        max_number_of_filters_being_tested =
          state->encoder_control->cfg.param_set_map[aps_ids[test_filter_idx] + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.cc_alf_aps_param.cc_alf_filter_count[comp_id - 1];
        cc_alf_filter_count = max_number_of_filters_being_tested;
        for (int filter_idx = 0; filter_idx < max_number_of_filters_being_tested; filter_idx++)
        {
          cc_alf_filter_idx_enabled[filter_idx] = true;
          memcpy(cc_alf_filter_coeff[filter_idx], cc_filter_param->cc_alf_coeff[comp_id - 1][filter_idx],
            sizeof(cc_alf_filter_coeff[filter_idx]));
        }
        memcpy(cc_alf_filter_coeff, state->encoder_control->cfg.param_set_map[aps_ids[test_filter_idx] + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.cc_alf_aps_param.cc_alf_coeff[comp_id - 1], sizeof(cc_alf_filter_coeff));
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
      assert(max_number_of_filters_being_tested != 0); //max_number_of_filters_being_tested should not be 0.
      const int column_size = (pic_width_c / max_number_of_filters_being_tested);
      for (int y = 0; y < pic_height_c; y += ctu_height_c)
      {
        for (int x = 0; x < pic_width_c; x += ctu_width_c)
        {
          training_cov_control[control_idx] = (x / column_size) + 1;
          control_idx++;
        }
      }

      // compute cost of filtering
      training_iter_count = 0;
      keep_training = true;
      improvement = false;
      prev_total_cost = MAX_DOUBLE;
      while (keep_training)
      {
        improvement = false;
        for (int filter_idx = 0; filter_idx < max_number_of_filters_being_tested; filter_idx++)
        {
          if (cc_alf_filter_idx_enabled[filter_idx])
          {
            if (!referencing_existing_aps)
            {
              get_frame_stats_cc_alf(alf_covariance_cc_alf, alf_covariance_frame_cc_alf, (filter_idx + 1), num_ctus_in_pic, training_cov_control);
              derive_cc_alf_filter_coeff(alf_covariance_frame_cc_alf, cc_alf_filter_coeff, filter_idx);
            }

            for (int y = 0; y < pic_height_c; y += (1 << log2_block_height))
            {
              for (int x = 0; x < pic_width_c; x += (1 << log2_block_width))
              {
                int ctu_idx = (y >> log2_block_height) * ctus_in_width + (x >> log2_block_width);
                alf_info->training_distortion[filter_idx][ctu_idx] =
                  (int)(alf_info->ctb_distortion_unfilter[comp_id][ctu_idx]
                    + calc_error_for_cc_alf_coeffs(&alf_covariance_cc_alf[ctu_idx],
                      cc_alf_filter_coeff[filter_idx], num_coeff, 7 + 1));
              }
            }
          }
        }

        memcpy(cabac_estimator, &ctx_start_cc_alf_filter_control_flag, sizeof(*cabac_estimator));

        cur_total_distortion = 0;
        cur_total_rate = 0;
        determine_control_idc_values(state, comp_id, ctu_width_c, ctu_height_c, pic_width_c, pic_height_c,
          alf_info->ctb_distortion_unfilter, alf_info->training_distortion,
          referencing_existing_aps, training_cov_control, filter_control, &cur_total_distortion,
          &cur_total_rate, cc_alf_filter_idx_enabled, map_filter_idx_to_filter_idc, &cc_alf_filter_count);

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
            memcpy(best_filter_control, filter_control, sizeof(uint8_t) * num_ctus_in_pic);
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
    memset(best_filter_control, 0, sizeof(uint8_t) * num_ctus_in_pic);
  }

  // save best coeff and control
  bool atleast_one_block_undergoes_fitlering = false;
  for (int controlIdx = 0; best_filter_count > 0 && controlIdx < num_ctus_in_pic; controlIdx++)
  {
    if (best_filter_control[controlIdx])
    {
      atleast_one_block_undergoes_fitlering = true;
      break;
    }
  }
  cc_filter_param->number_valid_components = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_COMPONENT;
  cc_filter_param->cc_alf_filter_enabled[comp_id - 1] = atleast_one_block_undergoes_fitlering;
  if (atleast_one_block_undergoes_fitlering)
  {
    // update the filter control indicators
    if (best_reuse_temporal_filter_coeff != 1)
    {
      short stored_best_filter_coeff_set[MAX_NUM_CC_ALF_FILTERS][MAX_NUM_CC_ALF_CHROMA_COEFF];
      for (int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++)
      {
        memcpy(stored_best_filter_coeff_set[filter_idx], best_filter_coeff_set[filter_idx], sizeof(best_filter_coeff_set[filter_idx]));
      }
      memcpy(filter_control, best_filter_control, sizeof(uint8_t) * num_ctus_in_pic);

      int filter_count = 0;
      for (int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++)
      {
        uint8_t cur_filter_idc = best_map_filter_idx_to_filter_idc[filter_idx];
        if (best_filter_idx_enabled[filter_idx])
        {
          for (int control_idx = 0; control_idx < num_ctus_in_pic; control_idx++)
          {
            if (filter_control[control_idx] == (filter_idx + 1))
            {
              best_filter_control[control_idx] = cur_filter_idc;
            }
          }
          memcpy(best_filter_coeff_set[cur_filter_idc - 1], stored_best_filter_coeff_set[filter_idx], sizeof(stored_best_filter_coeff_set[filter_idx]));
          filter_count++;
        }
        best_filter_idx_enabled[filter_idx] = (filter_idx < best_filter_count) ? true : false;
      }
      assert(filter_count == best_filter_count); //Number of filters enabled did not match the filter count
    }

    cc_filter_param->cc_alf_filter_count[comp_id - 1] = best_filter_count;
    // cleanup before copying
    memset(alf_info->cc_alf_filter_control[comp_id - 1], 0, sizeof(uint8_t) * num_ctus_in_pic);
    for (int filter_idx = 0; filter_idx < MAX_NUM_CC_ALF_FILTERS; filter_idx++)
    {
      memset(cc_filter_param->cc_alf_coeff[comp_id - 1][filter_idx], 0,
        sizeof(cc_filter_param->cc_alf_coeff[comp_id - 1][filter_idx]));
    }
    memset(cc_filter_param->cc_alf_filter_idx_enabled[comp_id - 1], false,
      sizeof(cc_filter_param->cc_alf_filter_idx_enabled[comp_id - 1]));
    for (int filter_idx = 0; filter_idx < best_filter_count; filter_idx++)
    {
      cc_filter_param->cc_alf_filter_idx_enabled[comp_id - 1][filter_idx] = best_filter_idx_enabled[filter_idx];
      memcpy(cc_filter_param->cc_alf_coeff[comp_id - 1][filter_idx], best_filter_coeff_set[filter_idx],
        sizeof(best_filter_coeff_set[filter_idx]));
    }
    memcpy(alf_info->cc_alf_filter_control[comp_id - 1], best_filter_control, sizeof(uint8_t) * num_ctus_in_pic);
    if (cc_alf_reuse_aps_id >= 0)
    {
      cc_reuse_aps_id[comp_id - 1] = cc_alf_reuse_aps_id;
      if (comp_id == COMPONENT_Cb)
      {
        state->slice->tile_group_cc_alf_cb_aps_id = cc_alf_reuse_aps_id;
      }
      else
      {
        state->slice->tile_group_cc_alf_cr_aps_id = cc_alf_reuse_aps_id;
      }
    }
  }

  /*if (luma_swing_greater_than_threshold_count)
  {
    FREE_POINTER(luma_swing_greater_than_threshold_count);
  }
  if (chroma_sample_count_near_mid_point)
  {
    FREE_POINTER(chroma_sample_count_near_mid_point);
  }*/

}

static void calc_covariance_cc_alf(int32_t e_local[MAX_NUM_CC_ALF_CHROMA_COEFF][1], const kvz_pixel *rec, const int stride, int vb_distance)
{
  const kvz_pixel *rec_y_m1 = rec - 1 * stride;
  const kvz_pixel *rec_y_0 = rec;
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

static void get_blk_stats_cc_alf(encoder_state_t * const state,
  alf_covariance *alf_covariance,
  const kvz_picture *org_yuv,
  const alf_component_id comp_id,
  const int x_pos, const int y_pos,
  const int width, const int height)
{
  alf_info_t *alf_info = state->tile->frame->alf_info;
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
  kvz_pixel *rec_y = &alf_info->alf_tmp_y[luma_rec_pos];
  kvz_pixel *rec_u = &alf_info->alf_tmp_u[chroma_rec_pos];
  kvz_pixel *rec_v = &alf_info->alf_tmp_v[chroma_rec_pos];

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

static void derive_stats_for_cc_alf_filtering(encoder_state_t * const state,
  const kvz_picture *org_yuv,
  const int comp_idx, const int mask_stride,
  const uint8_t filter_idc)
{
  alf_covariance **alf_covariance_cc_alf = state->tile->frame->alf_info->alf_covariance_cc_alf;
  alf_covariance *alf_covariance_frame_cc_alf = state->tile->frame->alf_info->alf_covariance_frame_cc_alf[comp_idx - 1];
  const int32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  const int filter_idx = filter_idc - 1;

  // init CTU stats buffers
  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    reset_alf_covariance(&alf_covariance_cc_alf[comp_idx - 1][(filter_idx * num_ctus_in_pic) + ctu_idx], -1);
  }

  // init Frame stats buffers
  reset_alf_covariance(&alf_covariance_frame_cc_alf[filter_idx], -1);

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
      get_blk_stats_cc_alf(state, &alf_covariance_cc_alf[comp_idx - 1][(filter_idx * num_ctus_in_pic) + ctu_rs_addr],
        org_yuv, comp_idx, x_pos, y_pos, width, height);
      add_alf_cov(&alf_covariance_frame_cc_alf[filter_idx],
        &alf_covariance_cc_alf[comp_idx - 1][(filter_idx * num_ctus_in_pic) + ctu_rs_addr]);
      ctu_rs_addr++;
    }
  }
}

static void count_luma_swing_greater_than_threshold(const kvz_pixel* luma,
  int luma_stride, int height, int width,
  int log2_block_width, int log2_block_height,
  uint64_t* luma_swing_greater_than_threshold_count,
  int luma_count_stride,
  int8_t input_bit_depth)
{
  const int threshold = (1 << (input_bit_depth - 2)) - 1;

  // 3x4 Diamond
  int x_support[] = { 0, -1, 0, 1, -1, 0, 1, 0 };
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

            if (p < min_val)
            {
              min_val = p;
            }
            if (p > max_val)
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

static void count_chroma_sample_value_near_mid_point(const kvz_pixel* chroma, int chroma_stride, int height, int width,
  int log2_block_width, int log2_block_height,
  uint64_t* chroma_sample_count_near_mid_point,
  int chroma_sample_count_near_mid_point_stride,
  int8_t input_bit_depth)
{
  const int mid_point = (1 << input_bit_depth) >> 1;
  const int threshold = 16;

  for (int y = 0; y < height; y += (1 << log2_block_height))
  {
    for (int x = 0; x < width; x += (1 << log2_block_width))
    {
      chroma_sample_count_near_mid_point[(y >> log2_block_height) * chroma_sample_count_near_mid_point_stride + (x >> log2_block_width)] = 0;

      for (int y_off = 0; y_off < (1 << log2_block_height); y_off++)
      {
        for (int x_off = 0; x_off < (1 << log2_block_width); x_off++)
        {
          if ((y + y_off) >= height || (x + x_off) >= width)
          {
            continue;
          }

          int distance_to_mid_point = abs(chroma[y_off * chroma_stride + x + x_off] - mid_point);
          if (distance_to_mid_point < threshold)
          {
            chroma_sample_count_near_mid_point[(y >> log2_block_height) * chroma_sample_count_near_mid_point_stride + (x >> log2_block_width)]++;
          }
        }
      }
    }
    chroma += (chroma_stride << log2_block_height);
  }
}

static void init_distortion_cc_alf(alf_covariance* alf_covariance_cc_alf[MAX_NUM_COMPONENT], double **ctb_distortion_unfilter, const int num_ctus)
{
  for (int comp = 1; comp < MAX_NUM_COMPONENT; comp++)
  {
    for (int ctb_idx = 0; ctb_idx < num_ctus; ctb_idx++)
    {
      ctb_distortion_unfilter[comp][ctb_idx] = alf_covariance_cc_alf[comp - 1][ctb_idx].pix_acc;
    }
  }
}

static void alf_reconstruct_coeff(encoder_state_t * const state,
  alf_aps *aps,
  channel_type channel,
  const bool is_rdo,
  const bool is_redo,
  array_variables *arr_vars)
{
  const int8_t bit_depth = state->encoder_control->bitdepth;
  int factor = is_rdo ? 0 : (1 << (bit_depth - 1));
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  alf_filter_type filter_type = is_luma ? ALF_FILTER_7X7 : ALF_FILTER_5X5;
  int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
  int num_coeff = filter_type == ALF_FILTER_5X5 ? 7 : 13;
  int num_coeff_minus1 = num_coeff - 1;
  const int num_alts = is_luma ? 1 : aps->num_alternatives_chroma;

  for (int alt_idx = 0; alt_idx < num_alts; ++alt_idx)
  {
    int num_filters = is_luma ? aps->num_luma_filters : 1;
    short* coeff = is_luma ? aps->luma_coeff : aps->chroma_coeff[alt_idx];
    int16_t* clipp = is_luma ? aps->luma_clipp : aps->chroma_clipp[alt_idx];

    for (int filter_idx = 0; filter_idx < num_filters; filter_idx++)
    {
      coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
    }

    if (!is_luma)
    {
      for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
      {
        arr_vars->chroma_coeff_final[alt_idx][coeff_idx] = coeff[coeff_idx];
        int clip_idx = aps->non_linear_flag[channel] ? clipp[coeff_idx] : 0;
        arr_vars->chroma_clipp_final[alt_idx][coeff_idx] = is_rdo ? clip_idx : arr_vars->alf_clipping_values[channel][clip_idx];
      }
      arr_vars->chroma_coeff_final[alt_idx][num_coeff_minus1] = factor;
      arr_vars->chroma_clipp_final[alt_idx][num_coeff_minus1] = is_rdo ? 0 : arr_vars->alf_clipping_values[channel][0];
      continue;
    }
    for (int class_idx = 0; class_idx < num_classes; class_idx++)
    {
      int filter_idx = aps->filter_coeff_delta_idx[class_idx];
      assert((filter_idx >= 0 && filter_idx <= aps->num_luma_filters)); // "Bad coeff delta idx in ALF"
      for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
      {
        arr_vars->coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] = coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx];
      }
      arr_vars->coeff_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
      arr_vars->clipp_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = is_rdo ? 0 : arr_vars->alf_clipping_values[channel][0];
      for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
      {
        int clip_idx = aps->non_linear_flag[channel] ? clipp[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] : 0;
        assert((clip_idx >= 0 && clip_idx < MAX_ALF_NUM_CLIPPING_VALUES)); // "Bad clip idx in ALF"
        arr_vars->clipp_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] = is_rdo ? clip_idx : arr_vars->alf_clipping_values[channel][clip_idx];
      }
      arr_vars->clipp_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] =
        is_rdo ? 0 :
        arr_vars->alf_clipping_values[channel][0];
    }
  }
}

static double alf_derive_ctb_alf_enable_flags(encoder_state_t * const state,
  channel_type channel,
  double *dist_unfilter,
  const int num_classes,
  const double chroma_weight,
  array_variables *arr_vars)
{
  alf_info_t *alf_info = state->tile->frame->alf_info;
  short* alf_ctb_filter_index = alf_info->alf_ctb_filter_index;
  cabac_data_t *cabac_estimator = &alf_info->cabac_estimator;
  cabac_data_t ctx_temp_start;
  cabac_data_t ctx_temp_best;
  cabac_data_t ctx_temp_alt_start;
  //cabac_data_t ctx_temp_alt_best;

  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;

  alf_aps *alf_param_temp = &alf_info->alf_param_temp;
  const kvz_pixel comp_id_first = is_luma ? COMPONENT_Y : COMPONENT_Cb;
  const kvz_pixel comp_id_last = is_luma ? COMPONENT_Y : COMPONENT_Cr;
  const int num_alts = is_luma ? 1 : alf_param_temp->num_alternatives_chroma;
  const int8_t bit_depth = state->encoder_control->bitdepth;
  const int32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  int num_coeff = is_luma ? 13 : 7;
  double cost = 0;
  double lambda = state->frame->lambda;
  *dist_unfilter = 0;

  if (is_luma) {
    alf_param_temp->enabled_flag[COMPONENT_Y] = 1;
  }
  else {
    alf_param_temp->enabled_flag[COMPONENT_Cb] = 1;
    alf_param_temp->enabled_flag[COMPONENT_Cr] = 1;
  }

  assert((chroma_weight <= 0.0) && (state->slice->start_in_rs == 0)); //"incompatible start CTU address, must be 0"

  alf_reconstruct_coeff(state, alf_param_temp, channel, true, is_luma, arr_vars);

  for (int alt_idx = 0; alt_idx < (is_luma ? 1 : MAX_NUM_ALF_ALTERNATIVES_CHROMA); alt_idx++)
  {
    for (int class_idx = 0; class_idx < (is_luma ? MAX_NUM_ALF_CLASSES : 1); class_idx++)
    {
      for (int i = 0; i < (is_luma ? MAX_NUM_ALF_LUMA_COEFF : MAX_NUM_ALF_CHROMA_COEFF); i++)
      {
        arr_vars->filter_coeff_set[is_luma ? class_idx : alt_idx][i] = is_luma ? arr_vars->coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + i] : arr_vars->chroma_coeff_final[alt_idx][i];
        arr_vars->filter_clipp_set[is_luma ? class_idx : alt_idx][i] = is_luma ? arr_vars->clipp_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + i] : arr_vars->chroma_clipp_final[alt_idx][i];
      }
    }
  }

  alf_covariance* alf_cov;
  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    for (int comp_id = comp_id_first; comp_id <= comp_id_last; comp_id++)
    {
      alf_cov = (is_luma ? alf_info->alf_covariance_y :
        (comp_id == COMPONENT_Cb ? alf_info->alf_covariance_u : alf_info->alf_covariance_v));

      const double ctu_lambda = chroma_weight > 0.0 ? (is_luma ? 0/*cs.picture->m_uEnerHpCtu[ctuIdx]*/ : 0/*cs.picture->m_uEnerHpCtu[ctuIdx]*/ / chroma_weight) : lambda;
      double dist_unfilter_ctu = get_unfiltered_distortion_cov_classes(&alf_cov[ctu_idx * num_classes], num_classes);

      //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
      memcpy(&ctx_temp_start, cabac_estimator, sizeof(ctx_temp_start));
      //m_CABACEstimator->resetBits();
      alf_cabac_reset_bits(cabac_estimator);
      cabac_estimator->only_count = 1;
      alf_info->ctu_enable_flag[comp_id][ctu_idx] = 1;
      code_alf_ctu_enable_flag(state, cabac_estimator, ctu_idx, comp_id, alf_param_temp);

      if (is_luma)
      {
        // Evaluate cost of signaling filter set index for convergence of filters enabled flag / filter derivation
        assert(alf_ctb_filter_index[ctu_idx] == ALF_NUM_FIXED_FILTER_SETS);
        assert(state->slice->tile_group_num_aps == 1);
        code_alf_ctu_filter_index(state, cabac_estimator, ctu_idx, alf_param_temp->enabled_flag[COMPONENT_Y]);
      }
      double cost_on = dist_unfilter_ctu + ctu_lambda * (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3);

      //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
      memcpy(&ctx_temp_best, cabac_estimator, sizeof(ctx_temp_best));

      if (is_luma)
      {
        cost_on += get_filtered_distortion(&alf_cov[ctu_idx * num_classes], arr_vars, num_classes, alf_param_temp->num_luma_filters - 1, num_coeff, bit_depth);
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
            memcpy(cabac_estimator, &ctx_temp_alt_start, sizeof(*cabac_estimator));
          }
          //m_CABACEstimator->resetBits();
          alf_cabac_reset_bits(cabac_estimator);
          cabac_estimator->only_count = 1;
          alf_info->ctu_alternative[comp_id][ctu_idx] = alt_idx;
          code_alf_ctu_alternative_ctu(state, cabac_estimator, ctu_idx, comp_id, alf_param_temp);
          double r_altCost = ctu_lambda * (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3); //frac_bits_scale * 0/*m_CABACEstimator->getEstFracBits()*/;

          double alt_dist = 0.;
          alt_dist += calc_error_for_coeffs(&alf_cov[ctu_idx * num_classes], arr_vars->filter_clipp_set[alt_idx], arr_vars->filter_coeff_set[alt_idx], num_coeff, bit_depth);

          double alt_cost = alt_dist + r_altCost;
          if (alt_cost < best_alt_cost)
          {
            best_alt_cost = alt_cost;
            best_alt_idx = alt_idx;
            //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
            memcpy(&ctx_temp_best, cabac_estimator, sizeof(ctx_temp_best));
          }
        }
        alf_info->ctu_alternative[comp_id][ctu_idx] = best_alt_idx;
        cost_on += best_alt_cost;
      }

      //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
      memcpy(cabac_estimator, &ctx_temp_start, sizeof(*cabac_estimator));
      //m_CABACEstimator->resetBits();
      alf_cabac_reset_bits(cabac_estimator);
      cabac_estimator->only_count = 1;
      alf_info->ctu_enable_flag[comp_id][ctu_idx] = 0;
      code_alf_ctu_enable_flag(state, cabac_estimator, ctu_idx, comp_id, alf_param_temp);
      double cost_off = dist_unfilter_ctu + ctu_lambda * (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3); //frac_bits_scale * 0;// m_CABACEstimator->getEstFracBits();

      if (cost_on < cost_off)
      {
        cost += cost_on;
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
        memcpy(cabac_estimator, &ctx_temp_best, sizeof(*cabac_estimator));
        alf_info->ctu_enable_flag[comp_id][ctu_idx] = 1;
      }
      else
      {
        cost += cost_off;
        alf_info->ctu_enable_flag[comp_id][ctu_idx] = 0;
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
      alf_param_temp->enabled_flag[compId] = false;
      for (int i = 0; i < num_ctus_in_pic; i++)
      {
        if (alf_info->ctu_enable_flag[compId][i])
        {
          alf_param_temp->enabled_flag[compId] = true;
          break;
        }
      }
    }
  }

  return cost;
}

static void alf_create_frame_buffer(encoder_state_t * const state, alf_info_t *alf_info)
{
  if (!alf_info->alf_fulldata_buf)
  {
    enum kvz_chroma_format chroma_format = state->encoder_control->chroma_format;
    const size_t simd_padding_width = 64;
    int width = state->tile->frame->width;
    int height = state->tile->frame->height;
    int stride = state->tile->frame->source->stride;
    unsigned int luma_size = (width + 8) * (height + 8);
    unsigned chroma_sizes[] = { 0, luma_size / 4, luma_size / 2, luma_size };
    unsigned chroma_size = chroma_sizes[chroma_format];

    alf_info->alf_fulldata_buf = MALLOC_SIMD_PADDED(kvz_pixel, (luma_size + 2 * chroma_size), simd_padding_width * 2);
    alf_info->alf_fulldata = &alf_info->alf_fulldata_buf[4 * (width + 8) + 4] + simd_padding_width / sizeof(kvz_pixel);
    alf_info->alf_tmp_y = &alf_info->alf_fulldata[0];

    if (chroma_format == KVZ_CSP_400) {
      alf_info->alf_tmp_u = NULL;
      alf_info->alf_tmp_v = NULL;
    }
    else {
      alf_info->alf_tmp_u = &alf_info->alf_fulldata[luma_size - (4 * (width + 8) + 4) + (2 * (stride / 2) + 2)];
      alf_info->alf_tmp_v = &alf_info->alf_fulldata[luma_size - (4 * (width + 8) + 4) + chroma_size + (2 * (stride / 2) + 2)];
    }
  }

}


void kvz_alf_create(videoframe_t *frame, enum kvz_chroma_format chroma_format)
{
  const int num_ctus_in_pic = frame->width_in_lcu * frame->height_in_lcu;
  const int pic_width = frame->width;
  const int pic_height = frame->height;
  const int luma_coeffs = 13;
  const int chroma_coeffs = 7;
  const int cc_alf_coeff = 8;
  int num_classes = 0;

  alf_info_t *alf_info = frame->alf_info;
  alf_info->aps_id_start = ALF_CTB_MAX_NUM_APS;

  alf_info->ctu_enable_flag[MAX_NUM_COMPONENT] = malloc(num_ctus_in_pic * MAX_NUM_COMPONENT * sizeof(*alf_info->ctu_enable_flag[MAX_NUM_COMPONENT]));
  memset(alf_info->ctu_enable_flag[MAX_NUM_COMPONENT], 0, num_ctus_in_pic * MAX_NUM_COMPONENT * sizeof(*alf_info->ctu_enable_flag[MAX_NUM_COMPONENT]));
  alf_info->ctu_enable_flag_tmp[MAX_NUM_COMPONENT] = malloc(num_ctus_in_pic * MAX_NUM_COMPONENT * sizeof(*alf_info->ctu_enable_flag_tmp[MAX_NUM_COMPONENT]));
  memset(alf_info->ctu_enable_flag_tmp[MAX_NUM_COMPONENT], 0, num_ctus_in_pic * MAX_NUM_COMPONENT * sizeof(*alf_info->ctu_enable_flag_tmp[MAX_NUM_COMPONENT]));

  alf_info->ctu_alternative[MAX_NUM_COMPONENT] = malloc(num_ctus_in_pic * (MAX_NUM_COMPONENT - 1) * sizeof(*alf_info->ctu_alternative[MAX_NUM_COMPONENT]));
  memset(alf_info->ctu_alternative[MAX_NUM_COMPONENT], 0, num_ctus_in_pic * (MAX_NUM_COMPONENT - 1) * sizeof(*alf_info->ctu_alternative[MAX_NUM_COMPONENT]));
  alf_info->ctu_alternative_tmp[MAX_NUM_COMPONENT] = malloc(num_ctus_in_pic * (MAX_NUM_COMPONENT - 1) * sizeof(*alf_info->ctu_alternative_tmp[MAX_NUM_COMPONENT]));
  memset(alf_info->ctu_alternative_tmp[MAX_NUM_COMPONENT], 0, num_ctus_in_pic * (MAX_NUM_COMPONENT - 1) * sizeof(*alf_info->ctu_alternative_tmp[MAX_NUM_COMPONENT]));

  alf_info->ctb_distortion_unfilter[MAX_NUM_COMPONENT] = malloc(num_ctus_in_pic * MAX_NUM_COMPONENT * sizeof(*alf_info->ctb_distortion_unfilter[MAX_NUM_COMPONENT]));
  memset(alf_info->ctb_distortion_unfilter[MAX_NUM_COMPONENT], 0, num_ctus_in_pic * MAX_NUM_COMPONENT * sizeof(*alf_info->ctb_distortion_unfilter[MAX_NUM_COMPONENT]));

  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    alf_info->ctu_enable_flag[comp_idx] = &alf_info->ctu_enable_flag[MAX_NUM_COMPONENT][comp_idx * num_ctus_in_pic];
    alf_info->ctu_enable_flag_tmp[comp_idx] = &alf_info->ctu_enable_flag_tmp[MAX_NUM_COMPONENT][comp_idx * num_ctus_in_pic];
    alf_info->ctb_distortion_unfilter[comp_idx] = &alf_info->ctb_distortion_unfilter[MAX_NUM_COMPONENT][comp_idx * num_ctus_in_pic];

    if (comp_idx == COMPONENT_Y)
    {
      alf_info->ctu_alternative[comp_idx] = NULL;
      alf_info->ctu_alternative_tmp[comp_idx] = NULL;
    }
    else
    {
      alf_info->ctu_alternative[comp_idx] = &alf_info->ctu_alternative[MAX_NUM_COMPONENT][(comp_idx - 1) * num_ctus_in_pic];
      alf_info->ctu_alternative_tmp[comp_idx] = &alf_info->ctu_alternative_tmp[MAX_NUM_COMPONENT][(comp_idx - 1) * num_ctus_in_pic];
    }
  }

  if (chroma_format != KVZ_CSP_400) {
    for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
    {
      num_classes += comp_idx ? 1 : MAX_NUM_ALF_CLASSES;
    }
  }
  else
  {
    num_classes = MAX_NUM_ALF_CLASSES;
  }

  const int num_covs = num_ctus_in_pic * num_classes;
  const int num_luma_covs = num_ctus_in_pic * MAX_NUM_ALF_CLASSES;
  alf_info->alf_covariance = malloc(num_covs * sizeof(alf_covariance));
  alf_info->alf_covariance_y = &alf_info->alf_covariance[0];

  for (int indx = 0; indx < num_luma_covs; indx++)
  {
    init_alf_covariance(&frame->alf_info->alf_covariance_y[indx], luma_coeffs);
  }
  for (int k = 0; k < MAX_NUM_ALF_CLASSES; k++)
  {
    init_alf_covariance(&alf_info->alf_covariance_frame_luma[k], luma_coeffs);
  }

  if (chroma_format != KVZ_CSP_400) {
    const int num_chroma_covs = num_ctus_in_pic;
    alf_info->alf_covariance_u = &alf_info->alf_covariance[num_luma_covs];
    alf_info->alf_covariance_v = &alf_info->alf_covariance[num_luma_covs + num_chroma_covs];
    for (int k = 0; k < num_chroma_covs; k++)
    {
      init_alf_covariance(&alf_info->alf_covariance_u[k], chroma_coeffs);
      init_alf_covariance(&alf_info->alf_covariance_v[k], chroma_coeffs);
    }

    for (int k = 0; k < MAX_NUM_ALF_ALTERNATIVES_CHROMA; k++)
    {
      init_alf_covariance(&alf_info->alf_covariance_frame_chroma[k], chroma_coeffs);
    }

    alf_info->alf_covariance_cc_alf[MAX_NUM_COMPONENT - 1] = malloc(num_ctus_in_pic * MAX_NUM_CC_ALF_FILTERS * (MAX_NUM_COMPONENT - 1) * sizeof(*alf_info->alf_covariance_cc_alf[MAX_NUM_COMPONENT - 1]));
    for (int comp_idx = 0; comp_idx < (MAX_NUM_COMPONENT - 1); comp_idx++)
    {
      alf_info->alf_covariance_cc_alf[comp_idx] = &alf_info->alf_covariance_cc_alf[MAX_NUM_COMPONENT - 1][comp_idx  * MAX_NUM_CC_ALF_FILTERS * num_ctus_in_pic];
    }
    for (int k = 0; k < num_ctus_in_pic * MAX_NUM_CC_ALF_FILTERS * (MAX_NUM_COMPONENT - 1); k++)
    {
      init_alf_covariance(&alf_info->alf_covariance_cc_alf[MAX_NUM_COMPONENT - 1][k], cc_alf_coeff);
    }
    for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT - 1; comp_idx++)
    {
      for (int k = 0; k < MAX_NUM_CC_ALF_FILTERS; k++)
      {
        init_alf_covariance(&alf_info->alf_covariance_frame_cc_alf[comp_idx][k], cc_alf_coeff);
      }
    }
  }

  for (int k = 0; k <= MAX_NUM_ALF_CLASSES + 1; k++)
  {
    init_alf_covariance(&alf_info->alf_covariance_merged[k], luma_coeffs);
  }

  alf_info->training_cov_control = malloc(num_ctus_in_pic * sizeof(*alf_info->training_cov_control));
  alf_info->training_distortion[MAX_NUM_CC_ALF_FILTERS] = malloc(num_ctus_in_pic * MAX_NUM_CC_ALF_FILTERS * sizeof(*alf_info->training_distortion[MAX_NUM_CC_ALF_FILTERS]));
  memset(alf_info->training_distortion[MAX_NUM_CC_ALF_FILTERS], 0, num_ctus_in_pic * MAX_NUM_CC_ALF_FILTERS * sizeof(*alf_info->training_distortion[MAX_NUM_CC_ALF_FILTERS]));
  for (int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++)
  {
    alf_info->training_distortion[i] = &alf_info->training_distortion[MAX_NUM_CC_ALF_FILTERS][num_ctus_in_pic * i];
  }

  alf_info->filter_control = malloc(num_ctus_in_pic * sizeof(*alf_info->filter_control));
  alf_info->best_filter_control = malloc(num_ctus_in_pic * sizeof(*alf_info->best_filter_control));

  // Classification
  alf_info->classifier = malloc(pic_height * sizeof(**alf_info->classifier));
  alf_info->classifier[0] = malloc(pic_height * pic_width * sizeof(*alf_info->classifier));

  for (int i = 1; i < pic_height; i++)
  {
    alf_info->classifier[i] = alf_info->classifier[0] + i * pic_width;
  }

  alf_info->cc_alf_filter_control[2] = malloc(2 * num_ctus_in_pic * sizeof(*alf_info->cc_alf_filter_control[2]));
  memset(alf_info->cc_alf_filter_control[2], 0, 2 * num_ctus_in_pic * sizeof(*alf_info->cc_alf_filter_control[2]));
  alf_info->cc_alf_filter_control[0] = &alf_info->cc_alf_filter_control[2][0];
  alf_info->cc_alf_filter_control[1] = &alf_info->cc_alf_filter_control[2][num_ctus_in_pic];

  alf_info->alf_ctb_filter_index = malloc(num_ctus_in_pic * sizeof(*alf_info->alf_ctb_filter_index));
  alf_info->alf_ctb_filter_set_index_tmp = malloc(num_ctus_in_pic * sizeof(*alf_info->alf_ctb_filter_set_index_tmp));

  alf_info->alf_fulldata_buf = NULL;
  alf_info->alf_fulldata = NULL;
  alf_info->alf_tmp_y = NULL;
  alf_info->alf_tmp_u = NULL;
  alf_info->alf_tmp_v = NULL;

}


void kvz_alf_destroy(videoframe_t * const frame)
{
  alf_info_t *alf_info = frame->alf_info;

  if (alf_info->alf_covariance_y)
  {
    alf_info->alf_covariance_y = NULL;
  }
  if (alf_info->alf_covariance_u)
  {
    alf_info->alf_covariance_u = NULL;
  }
  if (alf_info->alf_covariance_v)
  {
    alf_info->alf_covariance_v = NULL;
  }
  if (alf_info->alf_covariance)
  {
    FREE_POINTER(alf_info->alf_covariance);
  }

  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    if (alf_info->ctu_enable_flag[comp_idx])
    {
      alf_info->ctu_enable_flag[comp_idx] = NULL;
    }
    if (alf_info->ctu_enable_flag_tmp[comp_idx])
    {
      alf_info->ctu_enable_flag_tmp[comp_idx] = NULL;
    }
    if (alf_info->ctu_alternative[comp_idx])
    {
      alf_info->ctu_alternative[comp_idx] = NULL;
    }
    if (alf_info->ctu_alternative_tmp[comp_idx])
    {
      alf_info->ctu_alternative_tmp[comp_idx] = NULL;
    }
    if (alf_info->ctb_distortion_unfilter[comp_idx])
    {
      alf_info->ctb_distortion_unfilter[comp_idx] = NULL;
    }
    if (comp_idx > 0)
    {
      if (alf_info->alf_covariance_cc_alf[comp_idx - 1])
      {
        alf_info->alf_covariance_cc_alf[comp_idx - 1] = NULL;
      }
    }

  }
  if (alf_info->ctu_enable_flag[MAX_NUM_COMPONENT])
  {
    FREE_POINTER(alf_info->ctu_enable_flag[MAX_NUM_COMPONENT]);
  }
  if (alf_info->ctu_enable_flag_tmp[MAX_NUM_COMPONENT])
  {
    FREE_POINTER(alf_info->ctu_enable_flag_tmp[MAX_NUM_COMPONENT]);
  }
  if (alf_info->ctu_alternative[MAX_NUM_COMPONENT])
  {
    FREE_POINTER(alf_info->ctu_alternative[MAX_NUM_COMPONENT]);
  }
  if (alf_info->ctu_alternative_tmp[MAX_NUM_COMPONENT])
  {
    FREE_POINTER(alf_info->ctu_alternative_tmp[MAX_NUM_COMPONENT]);
  }
  if (alf_info->ctb_distortion_unfilter[MAX_NUM_COMPONENT])
  {
    FREE_POINTER(alf_info->ctb_distortion_unfilter[MAX_NUM_COMPONENT]);
  }
  if (alf_info->alf_covariance_cc_alf[MAX_NUM_COMPONENT - 1])
  {
    FREE_POINTER(alf_info->alf_covariance_cc_alf[MAX_NUM_COMPONENT - 1]);
  }

  if (alf_info->training_cov_control)
  {
    FREE_POINTER(alf_info->training_cov_control);
  }

  for (int i = 0; i < MAX_NUM_CC_ALF_FILTERS; i++)
  {
    if (alf_info->training_distortion[i])
    {
      alf_info->training_distortion[i] = NULL;
    }
  }
  if (alf_info->training_distortion[MAX_NUM_CC_ALF_FILTERS])
  {
    FREE_POINTER(alf_info->training_distortion[MAX_NUM_CC_ALF_FILTERS]);
  }

  if (alf_info->filter_control)
  {
    FREE_POINTER(alf_info->filter_control);
  }
  if (alf_info->best_filter_control)
  {
    FREE_POINTER(alf_info->best_filter_control);
  }

  if (alf_info->classifier)
  {
    FREE_POINTER(alf_info->classifier[0]);
    FREE_POINTER(alf_info->classifier);
  }

  if (alf_info->cc_alf_filter_control[0])
  {
    alf_info->cc_alf_filter_control[0] = NULL;
  }
  if (alf_info->cc_alf_filter_control[1])
  {
    alf_info->cc_alf_filter_control[1] = NULL;
  }
  if (alf_info->cc_alf_filter_control[2])
  {
    FREE_POINTER(alf_info->cc_alf_filter_control[2]);
  }

  if (alf_info->alf_ctb_filter_index)
  {
    FREE_POINTER(alf_info->alf_ctb_filter_index);
  }

  if (alf_info->alf_ctb_filter_set_index_tmp)
  {
    FREE_POINTER(alf_info->alf_ctb_filter_set_index_tmp);
  }

  if (alf_info->alf_tmp_y)
  {
    alf_info->alf_tmp_y = NULL;
  }
  if (alf_info->alf_tmp_u)
  {
    alf_info->alf_tmp_u = NULL;
  }
  if (alf_info->alf_tmp_v)
  {
    alf_info->alf_tmp_v = NULL;
  }
  if (alf_info->alf_fulldata)
  {
    alf_info->alf_fulldata = NULL;
  }
  if (alf_info->alf_fulldata_buf)
  {
    FREE_POINTER(alf_info->alf_fulldata_buf);
  }
}

static void alf_merge_classes(alf_aps *alf_aps,
  channel_type channel,
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
    copy_cov(&cov_merged[i], &cov[i]);
    cov_merged[i].num_bins = alf_aps->non_linear_flag[CHANNEL_TYPE_LUMA] ? MAX_ALF_NUM_CLIPPING_VALUES : 1;
  }

  // Try merging different covariance matrices

  // temporal AlfCovariance structure is allocated as the last element in covMerged array, the size of covMerged is MAX_NUM_ALF_CLASSES + 1
  alf_covariance* tmp_cov = &cov_merged[MAX_NUM_ALF_CLASSES];
  tmp_cov->num_bins = alf_aps->non_linear_flag[CHANNEL_TYPE_LUMA] ? MAX_ALF_NUM_CLIPPING_VALUES : 1;

  // init Clip
  for (int i = 0; i < num_classes; i++)
  {
    for (int val = 0; val < MAX_NUM_ALF_LUMA_COEFF; val++) {
      clip_merged[num_remaining - 1][i][val] = alf_aps->non_linear_flag[CHANNEL_TYPE_LUMA] ? MAX_ALF_NUM_CLIPPING_VALUES / 2 : 0;
    }
    if (alf_aps->non_linear_flag[CHANNEL_TYPE_LUMA])
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

            double error_merged = alf_aps->non_linear_flag[CHANNEL_TYPE_LUMA] ? optimize_filter_clip(tmp_cov, tmp_clip) : calculate_error_opt_filt(tmp_cov, tmp_clip);
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

static double alf_derive_coeff_quant(channel_type channel,
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

static double alf_derive_filter_coeffs(alf_aps *aps,
  channel_type channel,
  alf_covariance *cov,
  alf_covariance *covMerged,
  short* filter_indices,
  int num_filters,
  double error_tab_force_0_coeff[MAX_NUM_ALF_CLASSES][2],
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  int filter_coeff_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  int filter_clipp_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  const int bit_depth)
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;

  double error = 0.0;
  alf_covariance *tmp_cov = &covMerged[MAX_NUM_ALF_CLASSES];

  for (int filt_idx = 0; filt_idx < num_filters; filt_idx++)
  {
    reset_alf_covariance(tmp_cov, -1);
    bool found_clip = false;
    for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
    {
      if (filter_indices[class_idx] == filt_idx)
      {
        add_alf_cov(tmp_cov, &cov[class_idx]);
        if (!found_clip)
        {
          found_clip = true; // clip should be at the adress of shortest one
          memcpy(filter_clipp_set[filt_idx], clip_merged[num_filters - 1][class_idx], sizeof(int[MAX_NUM_ALF_LUMA_COEFF]));
        }
      }
    }

    // Find coeffcients
    assert(num_coeff == tmp_cov->num_coeff);
    error_tab_force_0_coeff[filt_idx][1] = tmp_cov->pix_acc + alf_derive_coeff_quant(channel, filter_clipp_set[filt_idx], filter_coeff_set[filt_idx], tmp_cov, bit_depth, false);
    error_tab_force_0_coeff[filt_idx][0] = tmp_cov->pix_acc;
    error += error_tab_force_0_coeff[filt_idx][1];
  }
  return error;
}

static int alf_derive_filter_coefficients_prediction_mode(const alf_aps *alf_param,
  channel_type channel,
  const int num_filters,
  int filter_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  int filter_clipp_set[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF])
{
  return (alf_param->non_linear_flag[CHANNEL_TYPE_LUMA] ? get_cost_filter_clipp(channel, num_filters, filter_set, filter_clipp_set) : 0) + get_cost_filter_coeff(channel, num_filters, filter_set);
}

static double alf_merge_filters_and_cost(encoder_state_t * const state,
  alf_aps *alf_aps,
  channel_type channel,
  int *ui_coeff_bits,
  alf_covariance *cov_frame,
  alf_covariance *cov_merged,
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF],
  array_variables *arr_vars)
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

  //clip_merged:iä ei tarvitse nollata ennen
  alf_merge_classes(alf_aps, channel, cov_frame, cov_merged, clip_merged, MAX_NUM_ALF_CLASSES, arr_vars->filter_indices);

  while (num_filters >= 1)
  {
    dist = alf_derive_filter_coeffs(alf_aps, channel, cov_frame, cov_merged, arr_vars->filter_indices[num_filters - 1], num_filters, error_force_0_coeff_tab, clip_merged, arr_vars->filter_coeff_set, arr_vars->filter_clipp_set, bit_depth);

    // filter coeffs are stored in m_filterCoeffSet
    dist_force0 = get_dist_force_0(alf_aps, channel, num_filters, error_force_0_coeff_tab, coded_var_bins, lambda, arr_vars->filter_coeff_set, arr_vars->filter_clipp_set);
    coeff_bits = alf_derive_filter_coefficients_prediction_mode(alf_aps, channel, num_filters, arr_vars->filter_coeff_set, arr_vars->filter_clipp_set);
    coeff_bits_force0 = get_cost_filter_coeff_force_0(alf_aps, channel, num_filters, coded_var_bins, arr_vars->filter_coeff_set, arr_vars->filter_clipp_set);

    cost = dist + lambda * coeff_bits;
    cost0 = dist_force0 + lambda * coeff_bits_force0;

    if (cost0 < cost)
    {
      cost = cost0;
    }

    if (cost <= cost_min)
    {
      cost_min = cost;
      num_filters_best = num_filters;
    }
    num_filters--;
  }

  dist = alf_derive_filter_coeffs(alf_aps, channel, cov_frame, cov_merged, arr_vars->filter_indices[num_filters_best - 1], num_filters_best, error_force_0_coeff_tab, clip_merged, arr_vars->filter_coeff_set, arr_vars->filter_clipp_set, bit_depth);

  coeff_bits = alf_derive_filter_coefficients_prediction_mode(alf_aps, channel, num_filters_best, arr_vars->filter_coeff_set, arr_vars->filter_clipp_set);
  dist_force0 = get_dist_force_0(alf_aps, channel, num_filters_best, error_force_0_coeff_tab, coded_var_bins, lambda, arr_vars->filter_coeff_set, arr_vars->filter_clipp_set);
  coeff_bits_force0 = get_cost_filter_coeff_force_0(alf_aps, channel, num_filters_best, coded_var_bins, arr_vars->filter_coeff_set, arr_vars->filter_clipp_set);

  cost = dist + lambda * coeff_bits;
  cost0 = dist_force0 + lambda * coeff_bits_force0;

  alf_aps->num_luma_filters = num_filters_best;
  double dist_return;
  if (cost <= cost0)
  {
    dist_return = dist;
    alf_aps->alf_luma_coeff_delta_flag = 0;
    *ui_coeff_bits = coeff_bits;
  }
  else
  {
    dist_return = dist_force0;
    alf_aps->alf_luma_coeff_delta_flag = 1;
    *ui_coeff_bits = coeff_bits_force0;
    memcpy(alf_aps->alf_luma_coeff_flag, coded_var_bins, sizeof(coded_var_bins));

    for (int var_ind = 0; var_ind < num_filters_best; var_ind++)
    {
      if (coded_var_bins[var_ind] == 0)
      {
        memset(arr_vars->filter_coeff_set[var_ind], 0, sizeof(int) * MAX_NUM_ALF_LUMA_COEFF);
        memset(arr_vars->filter_clipp_set[var_ind], 0, sizeof(int) * MAX_NUM_ALF_LUMA_COEFF);
      }
    }
  }

  for (int ind = 0; ind < alf_aps->num_luma_filters; ++ind)
  {
    for (int i = 0; i < num_coeff; i++)
    {
      alf_aps->luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = arr_vars->filter_coeff_set[ind][i];
      alf_aps->luma_clipp[ind * MAX_NUM_ALF_LUMA_COEFF + i] = arr_vars->filter_clipp_set[ind][i];
    }
  }

  memcpy(alf_aps->filter_coeff_delta_idx, arr_vars->filter_indices[num_filters_best - 1], sizeof(short) * MAX_NUM_ALF_CLASSES);
  *ui_coeff_bits += get_non_filter_coeff_rate(alf_aps);
  return dist_return;
}

static double alf_get_filter_coeff_and_cost(encoder_state_t * const state,
  channel_type channel,
  double dist_unfilter,
  int *ui_coeff_bits,
  bool b_re_collect_stat,
  bool only_filter_cost,
  array_variables *arr_vars)
{
  alf_info_t *alf_info = state->tile->frame->alf_info;
  alf_aps *alf_param_temp = &alf_info->alf_param_temp;
  cabac_data_t *cabac_estimator = &alf_info->cabac_estimator;
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  double lambda = state->frame->lambda;
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF];
  const int8_t bit_depth = state->encoder_control->bitdepth;
  const int32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  alf_covariance *alf_cov_frame = is_luma ? alf_info->alf_covariance_frame_luma : alf_info->alf_covariance_frame_chroma;

  //collect stat based on CTU decision
  if (b_re_collect_stat)
  {
    get_frame_stats(state->tile->frame->alf_info, channel, num_ctus_in_pic);
  }

  double dist = dist_unfilter;
  (*ui_coeff_bits) = 0;

  //get filter coeff
  if (is_luma)
  {
    const int fill_val = alf_param_temp->non_linear_flag[channel] ? MAX_ALF_NUM_CLIPPING_VALUES / 2 : 0;
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++) {
      for (int j = 0; j < MAX_NUM_ALF_CLASSES; j++) {
        for (int k = 0; k < MAX_NUM_ALF_LUMA_COEFF; k++) {
          clip_merged[i][j][k] = fill_val;
        }
      }
    }

    // Reset Merge Tmp Cov
    reset_alf_covariance(&alf_info->alf_covariance_merged[MAX_NUM_ALF_CLASSES], MAX_ALF_NUM_CLIPPING_VALUES);
    reset_alf_covariance(&alf_info->alf_covariance_merged[MAX_NUM_ALF_CLASSES + 1], MAX_ALF_NUM_CLIPPING_VALUES);
    //distortion
    dist += alf_merge_filters_and_cost(state, alf_param_temp, channel, ui_coeff_bits, alf_cov_frame, alf_info->alf_covariance_merged, clip_merged, arr_vars);
  }
  else
  {
    //distortion
    for (int alt_idx = 0; alt_idx < alf_param_temp->num_alternatives_chroma; ++alt_idx)
    {
      assert(num_coeff == alf_cov_frame[alt_idx].num_coeff);
      alf_aps best_slice_param;
      double best_cost = MAX_DOUBLE;
      double best_dist = MAX_DOUBLE;
      int best_coeff_bits = 0;
      const int non_linear_flag_max = state->encoder_control->cfg.alf_non_linear_chroma ? 2 : 1;

      for (int non_linear_flag = 0; non_linear_flag < non_linear_flag_max; non_linear_flag++)
      {
        int current_non_linear_flag = alf_param_temp->non_linear_flag[channel] ? 1 : 0;
        if (non_linear_flag != current_non_linear_flag)
        {
          continue;
        }

        int fill_val = non_linear_flag ? MAX_ALF_NUM_CLIPPING_VALUES / 2 : 0;
        for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++) {
          arr_vars->filter_clipp_set[alt_idx][i] = fill_val;
        }

        double dist = alf_cov_frame[alt_idx].pix_acc + alf_derive_coeff_quant(channel, arr_vars->filter_clipp_set[alt_idx], arr_vars->filter_coeff_set[alt_idx], &alf_cov_frame[alt_idx], bit_depth, non_linear_flag);
        for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
        {
          alf_param_temp->chroma_coeff[alt_idx][i] = arr_vars->filter_coeff_set[alt_idx][i];
          alf_param_temp->chroma_clipp[alt_idx][i] = arr_vars->filter_clipp_set[alt_idx][i];
        }
        int coeff_bits = get_chroma_coeff_rate(alf_param_temp, alt_idx);
        double cost = dist + lambda * coeff_bits;
        if (cost < best_cost)
        {
          best_cost = cost;
          best_dist = dist;
          best_coeff_bits = coeff_bits;
          copy_alf_param(&best_slice_param, alf_param_temp);
        }
      }
      *ui_coeff_bits += best_coeff_bits;
      dist += best_dist;
      copy_alf_param(alf_param_temp, &best_slice_param);
    }
    (*ui_coeff_bits) += length_uvlc(alf_param_temp->num_alternatives_chroma - 1);
    (*ui_coeff_bits)++;
  }

  if (only_filter_cost)
  {
    return dist + lambda * (*ui_coeff_bits);
  }
  double rate = (*ui_coeff_bits);
  //m_CABACEstimator->resetBits();
  alf_cabac_reset_bits(cabac_estimator);
  //m_CABACEstimator->codeAlfCtuEnableFlags(cs, channel, &m_alfParamTemp);
  code_alf_ctu_enable_flags_channel(state, cabac_estimator, channel, alf_param_temp);

  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    if (is_luma)
    {
      // Evaluate cost of signaling filter set index for convergence of filters enabled flag / filter derivation
      assert(alf_info->alf_ctb_filter_index[ctu_idx] == ALF_NUM_FIXED_FILTER_SETS);
      assert(state->slice->tile_group_num_aps == 1);
      //m_CABACEstimator->codeAlfCtuFilterIndex(cs, ctu_idx, &m_alfParamTemp.enabledFlag[COMPONENT_Y]);
      code_alf_ctu_filter_index(state, cabac_estimator, ctu_idx, alf_param_temp->enabled_flag[COMPONENT_Y]);
    }
  }
  //m_CABACEstimator->codeAlfCtuAlternatives(cs, channel, &m_alfParamTemp);
  code_alf_ctu_alternatives_channel(state, cabac_estimator, channel, alf_param_temp);
  rate += (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3); //frac_bits_scale * 0;/*(double)m_CABACEstimator->getEstFracBits();*/
  return dist + lambda * rate;
}

static void alf_encoder(encoder_state_t * const state,
  alf_aps *aps,
  channel_type channel,
  const double lambda_chroma_weight, // = 0.0
  array_variables *arr_vars)
{
  alf_info_t *alf_info = state->tile->frame->alf_info;
  bool **ctu_enable_flags = alf_info->ctu_enable_flag;
  bool **ctu_enable_flags_tmp = alf_info->ctu_enable_flag_tmp;
  uint8_t **ctu_alternatives = alf_info->ctu_alternative;
  uint8_t **ctu_alternatives_tmp = alf_info->ctu_alternative_tmp;
  alf_aps *alf_param_temp = &alf_info->alf_param_temp;

  cabac_data_t ctx_start;
  cabac_data_t *cabac_estimator = &alf_info->cabac_estimator;
  memcpy(&ctx_start, cabac_estimator, sizeof(ctx_start));
  //TempCtx        ctxBest(m_CtxCache);
  cabac_data_t ctx_best;
  memcpy(&ctx_best, &ctx_start, sizeof(ctx_best));

  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  alf_covariance *alf_cov_frame = is_luma ? alf_info->alf_covariance_frame_luma : alf_info->alf_covariance_frame_chroma;
  kvz_config cfg = state->encoder_control->cfg;

  double cost_min = MAX_DOUBLE;
  double lambda = state->frame->lambda;
  unsigned *bits_new_filter = arr_vars->bits_new_filter;
  bits_new_filter[channel] = 0;
  const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
  const int32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  int ui_coeff_bits = 0;

  //m_alfSliceParamTemp = alfSliceParam;
  copy_alf_param(alf_param_temp, aps);

  //1. get unfiltered distortion
  if (!is_luma)
  {
    alf_param_temp->num_alternatives_chroma = 1;
  }
  double cost = get_unfiltered_distortion_cov_channel(alf_cov_frame, channel);
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
    set_ctu_enable_flag(ctu_enable_flags_tmp, channel, 0, num_ctus_in_pic);
    if (!is_luma)
    {
      for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++) {
        ctu_alternatives_tmp[COMPONENT_Cb][ctu_idx] = 0;
        ctu_alternatives_tmp[COMPONENT_Cr][ctu_idx] = 0;
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
        alf_param_temp->num_alternatives_chroma = num_alternatives;
      }
      //2. all CTUs are on
      if (is_luma)
      {
        alf_param_temp->enabled_flag[COMPONENT_Y] = 1;
      }
      else
      {
        alf_param_temp->enabled_flag[COMPONENT_Cb] = 1;
        alf_param_temp->enabled_flag[COMPONENT_Cr] = 1;
      }
      alf_param_temp->non_linear_flag[channel] = non_linear_flag;

      //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
      memcpy(cabac_estimator, &ctx_start, sizeof(*cabac_estimator));
      //setCtuEnableFlag(m_ctuEnableFlag, channel, 1);
      set_ctu_enable_flag(ctu_enable_flags, channel, 1, num_ctus_in_pic);
      // all alternatives are on
      if (!is_luma)
      {
        init_ctu_alternative_chroma(alf_param_temp, ctu_alternatives, num_ctus_in_pic);
      }
      cost = alf_get_filter_coeff_and_cost(state, channel, 0, &ui_coeff_bits, true, false, arr_vars);

      if (cost < cost_min)
      {
        bits_new_filter[channel] = ui_coeff_bits;
        cost_min = cost;
        copy_alf_param_w_channel(aps, alf_param_temp, channel);
        //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
        memcpy(&ctx_best, cabac_estimator, sizeof(ctx_best));
        //setCtuEnableFlag(m_ctuEnableFlagTmp, channel, 1);
        set_ctu_enable_flag(ctu_enable_flags_tmp, channel, 1, num_ctus_in_pic);
        if (!is_luma)
        {
          memcpy(ctu_alternatives_tmp[COMPONENT_Cb], ctu_alternatives[COMPONENT_Cb], sizeof(uint8_t) * num_ctus_in_pic);
          memcpy(ctu_alternatives_tmp[COMPONENT_Cr], ctu_alternatives[COMPONENT_Cr], sizeof(uint8_t) * num_ctus_in_pic);
        }
      }

      //3. CTU decision
      double dist_unfilter = 0;
      double prev_it_cost = MAX_DOUBLE;
      const int iter_num = is_luma ? (2 * 4 + 1) : (2 * (2 + alf_param_temp->num_alternatives_chroma - 1) + 1);

      for (int iter = 0; iter < iter_num; iter++)
      {
        if ((iter & 0x01) == 0)
        {
          //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
          memcpy(cabac_estimator, &ctx_start, sizeof(*cabac_estimator));
          cost = lambda * ui_coeff_bits;
          cost += alf_derive_ctb_alf_enable_flags(state, channel, &dist_unfilter, num_classes, lambda_chroma_weight, arr_vars);
          if (cost < cost_min)
          {
            bits_new_filter[channel] = ui_coeff_bits;
            cost_min = cost;
            //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
            memcpy(&ctx_best, cabac_estimator, sizeof(ctx_best));
            //copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, channel);
            copy_ctu_enable_flag(ctu_enable_flags_tmp, ctu_enable_flags, channel, num_ctus_in_pic);
            if (!is_luma)
            {
              for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++) {
                ctu_alternatives_tmp[COMPONENT_Cb][ctu_idx] = ctu_alternatives[COMPONENT_Cb][ctu_idx];
                ctu_alternatives_tmp[COMPONENT_Cr][ctu_idx] = ctu_alternatives[COMPONENT_Cr][ctu_idx];
              }
            }
            copy_alf_param_w_channel(aps, alf_param_temp, channel);
          }
          else if (cost >= prev_it_cost)
          {
            // High probability that we have converged or we are diverging
            break;
          }
          prev_it_cost = cost;
        }
        else
        {
          // unfiltered distortion is added due to some CTBs may not use filter
          // no need to reset CABAC here, since uiCoeffBits is not affected
          /*cost = */alf_get_filter_coeff_and_cost(state, channel, dist_unfilter, &ui_coeff_bits, true, false, arr_vars);
        }
      }//for iter
    // Decrease number of alternatives and reset ctu params and filters
    }
  }//for non_linea_flag

  memcpy(cabac_estimator, &ctx_best, sizeof(*cabac_estimator));
  if (!is_luma)
  {
    memcpy(ctu_alternatives[COMPONENT_Cb], ctu_alternatives_tmp[COMPONENT_Cb], sizeof(uint8_t) * num_ctus_in_pic);
    memcpy(ctu_alternatives[COMPONENT_Cr], ctu_alternatives_tmp[COMPONENT_Cr], sizeof(uint8_t) * num_ctus_in_pic);
  }
  copy_ctu_enable_flag(ctu_enable_flags, ctu_enable_flags_tmp, channel, num_ctus_in_pic);
}

static void alf_get_avai_aps_ids_luma(encoder_state_t * const state,
  int *new_aps_id,
  int *aps_ids,
  int *size_of_aps_ids,
  short alf_clipping_values[MAX_NUM_CHANNEL_TYPE][MAX_ALF_NUM_CLIPPING_VALUES])
{
  //alf_aps *apss = state->slice->apss;
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    param_set_map* param_set = &state->encoder_control->cfg.param_set_map[i + NUM_APS_TYPE_LEN + T_ALF_APS];
    if (param_set->b_changed && (param_set->parameter_set.aps_id >= 0 || param_set->parameter_set.aps_id < ALF_CTB_MAX_NUM_APS)) {
      copy_aps(&state->slice->apss[i], &param_set->parameter_set, false);
    }
  }

  //std::vector<int> result;
  int aps_id_checked = 0, cur_aps_id = state->tile->frame->alf_info->aps_id_start;
  if (cur_aps_id < ALF_CTB_MAX_NUM_APS)
  {
    while ((aps_id_checked < ALF_CTB_MAX_NUM_APS) && !state->frame->is_irap && *size_of_aps_ids < ALF_CTB_MAX_NUM_APS /*&& /*!cs.slice->getPendingRasInit()*/)
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
  *new_aps_id = state->tile->frame->alf_info->aps_id_start - 1;
  if (*new_aps_id < 0)
  {
    *new_aps_id = (int)ALF_CTB_MAX_NUM_APS - 1;
  }
  assert(*new_aps_id < (int)ALF_CTB_MAX_NUM_APS); //Wrong APS index assignment in getAvaiApsIdsLuma
}

static void alf_calc_covariance(int16_t e_local[MAX_NUM_ALF_LUMA_COEFF][MAX_ALF_NUM_CLIPPING_VALUES],
  const kvz_pixel *rec,
  const int stride,
  const channel_type channel,
  const int transpose_idx,
  int vb_distance,
  short alf_clipping_values[MAX_NUM_CHANNEL_TYPE][MAX_ALF_NUM_CLIPPING_VALUES])
{
  static const int alf_pattern_5[13] = {
              0,
          1,  2,  3,
      4,  5,  6,  5,  4,
          3,  2,  1,
              0
  };

  static const int alf_pattern_7[25] = {
                0,
            1,  2,  3,
        4,  5,  6,  7,  8,
    9, 10, 11, 12, 11, 10, 9,
        8,  7,  6,  5,  4,
            3,  2,  1,
                0
  };

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
  const short* clip = alf_clipping_values[channel];
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

static void alf_get_blk_stats(encoder_state_t * const state,
  channel_type channel,
  alf_covariance *alf_covariance,
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
  int vb_pos,
  short alf_clipping_values[MAX_NUM_CHANNEL_TYPE][MAX_ALF_NUM_CLIPPING_VALUES])
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
      alf_calc_covariance(e_local, rec + j, rec_stride, channel, transpose_idx, vb_distance, alf_clipping_values);
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
                alf_covariance[class_idx].ee[b0][b1][k][l] += weight * (e_local[k][b0] * (double)e_local[l][b1]);
              }
              else
              {
                alf_covariance[class_idx].ee[b0][b1][k][l] += e_local[k][b0] * (double)e_local[l][b1];
              }
            }
          }
        }
        for (int b = 0; b < num_bins; b++)
        {
          if (0/*m_alfWSSD*/)
          {
            alf_covariance[class_idx].y[b][k] += weight * (e_local[k][b] * (double)y_local);
          }
          else
          {
            alf_covariance[class_idx].y[b][k] += e_local[k][b] * (double)y_local;
          }
        }
      }
      if (0/*m_alfWSSD*/)
      {
        alf_covariance[class_idx].pix_acc += weight * (y_local * (double)y_local);
      }
      else
      {
        alf_covariance[class_idx].pix_acc += y_local * (double)y_local;
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
            alf_covariance[class_idx].ee[b0][b1][k][l] = alf_covariance[class_idx].ee[b1][b0][l][k];
          }
        }
      }
    }
  }
}

static void alf_derive_stats_for_filtering(encoder_state_t * const state,
  short alf_clipping_values[MAX_NUM_CHANNEL_TYPE][MAX_ALF_NUM_CLIPPING_VALUES])
{
  alf_info_t *alf_info = state->tile->frame->alf_info;
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  bool chroma_scale_x = (chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool chroma_scale_y = (chroma_fmt != KVZ_CSP_420) ? 0 : 1;

  const int32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  const int alf_vb_luma_ctu_height = LCU_WIDTH;
  const int alf_vb_chma_ctu_height = (LCU_WIDTH >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0));
  const int alf_vb_luma_pos = LCU_WIDTH - ALF_VB_POS_ABOVE_CTUROW_LUMA;
  const int alf_vb_chma_pos = (LCU_WIDTH >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0)) - ALF_VB_POS_ABOVE_CTUROW_CHMA;
  int32_t pic_width = state->tile->frame->width;
  int32_t pic_height = state->tile->frame->height;
  int ctu_rs_addr = 0;

  const int number_of_components = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_COMPONENT;

  // init CTU stats buffers
  {
    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
    {
      for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
      {
        reset_alf_covariance(&state->tile->frame->alf_info->alf_covariance_y[(ctu_idx * MAX_NUM_ALF_CLASSES) + class_idx], MAX_ALF_NUM_CLIPPING_VALUES);
      }
    }

    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
    {
      reset_alf_covariance(&state->tile->frame->alf_info->alf_covariance_u[ctu_idx], MAX_ALF_NUM_CLIPPING_VALUES);
      reset_alf_covariance(&state->tile->frame->alf_info->alf_covariance_v[ctu_idx], MAX_ALF_NUM_CLIPPING_VALUES);
    }
  }

  // init Frame stats buffers
  const int number_of_channels = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_CHANNEL_TYPE;
  for (int channel_idx = 0; channel_idx < number_of_channels; channel_idx++)
  {
    const channel_type channel_id = channel_idx;
    const int num_classes = channel_id == CHANNEL_TYPE_LUMA ? MAX_NUM_ALF_CLASSES : 1;
    for (int class_idx = 0; class_idx < num_classes; class_idx++)
    {
      reset_alf_covariance(&alf_info->alf_covariance_frame_luma[class_idx], MAX_ALF_NUM_CLIPPING_VALUES);
    }
    reset_alf_covariance(&alf_info->alf_covariance_frame_chroma[0], MAX_ALF_NUM_CLIPPING_VALUES);
  }

  alf_covariance* alf_cov;
  alf_covariance* alf_cov_frame;
  for (int y_pos = 0; y_pos < pic_height; y_pos += LCU_WIDTH)
  {
    for (int x_pos = 0; x_pos < pic_width; x_pos += LCU_WIDTH)
    {
      const int width = (x_pos + LCU_WIDTH > pic_width) ? (pic_width - x_pos) : LCU_WIDTH;
      const int height = (y_pos + LCU_WIDTH > pic_height) ? (pic_height - y_pos) : LCU_WIDTH;
      for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
      {
        alf_cov = comp_idx == COMPONENT_Y ? state->tile->frame->alf_info->alf_covariance_y :
          comp_idx == COMPONENT_Cb ? state->tile->frame->alf_info->alf_covariance_u :
          comp_idx == COMPONENT_Cr ? state->tile->frame->alf_info->alf_covariance_v : NULL;

        if (alf_cov == NULL) {
          assert(0);
        }

        const bool is_luma = comp_idx == COMPONENT_Y ? 1 : 0;
        channel_type ch_type = is_luma ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;
        alf_cov_frame = is_luma ? alf_info->alf_covariance_frame_luma : alf_info->alf_covariance_frame_chroma;

        int blk_w = is_luma ? width : width >> chroma_scale_x;
        int blk_h = is_luma ? height : height >> chroma_scale_y;
        int pos_x = is_luma ? x_pos : x_pos >> chroma_scale_x;
        int pos_y = is_luma ? y_pos : y_pos >> chroma_scale_y;

        int32_t org_stride = is_luma ? state->tile->frame->source->stride : state->tile->frame->source->stride >> chroma_scale_x;
        int32_t rec_stride = is_luma ? state->tile->frame->rec->stride : state->tile->frame->rec->stride >> chroma_scale_x;

        kvz_pixel *org = comp_idx ? (comp_idx - 1 ? &state->tile->frame->source->v[pos_x + pos_y * org_stride] : &state->tile->frame->source->u[pos_x + pos_y * org_stride]) : &state->tile->frame->source->y[pos_x + pos_y * org_stride];
        kvz_pixel *rec = comp_idx ? (comp_idx - 1 ? &state->tile->frame->rec->v[pos_x + pos_y * rec_stride] : &state->tile->frame->rec->u[pos_x + pos_y * rec_stride]) : &state->tile->frame->rec->y[pos_x + pos_y * rec_stride];

        const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
        const int cov_index = ctu_rs_addr * num_classes;
        alf_get_blk_stats(state, ch_type,
          &alf_cov[cov_index],
          comp_idx ? NULL : alf_info->classifier,
          org, org_stride, rec, rec_stride, pos_x, pos_y, pos_x, pos_y, blk_w, blk_h,
          (is_luma ? alf_vb_luma_ctu_height : alf_vb_chma_ctu_height),
          (is_luma) ? alf_vb_luma_pos : alf_vb_chma_pos,
          alf_clipping_values
        );

        for (int class_idx = 0; class_idx < num_classes; class_idx++)
        {
          add_alf_cov(&alf_cov_frame[is_luma ? class_idx : 0],
            &alf_cov[cov_index + class_idx]
          );
        }
      }
      ctu_rs_addr++;
    }
  }
}

static void alf_reconstruct_coeff_aps(encoder_state_t * const state, bool luma, bool chroma, bool is_rdo,
  array_variables *arr_vars)
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

      alf_reconstruct_coeff(state, &alf_param_tmp, CHANNEL_TYPE_LUMA, is_rdo, true, arr_vars);
      memcpy(arr_vars->coeff_aps_luma[i], arr_vars->coeff_final, sizeof(arr_vars->coeff_final));
      memcpy(arr_vars->clipp_aps_luma[i], arr_vars->clipp_final, sizeof(arr_vars->clipp_final));
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
    alf_reconstruct_coeff(state, &alf_param_tmp, CHANNEL_TYPE_CHROMA, is_rdo, true, arr_vars);
  }
}

static void alf_encoder_ctb(encoder_state_t * const state,
  alf_aps *aps,
  const double lambda_chroma_weight,
  array_variables *arr_vars)
{
  alf_info_t *alf_info = state->tile->frame->alf_info;
  bool **ctu_enable_flag = alf_info->ctu_enable_flag;
  bool **ctu_enable_flag_tmp = alf_info->ctu_enable_flag_tmp;
  uint8_t **ctu_alternatives = alf_info->ctu_alternative;
  uint8_t **ctu_alternatives_tmp = alf_info->ctu_alternative_tmp;
  double **ctb_distortions_unfilter = alf_info->ctb_distortion_unfilter;
  unsigned *arr_bits_new_filter = arr_vars->bits_new_filter;
  short *alf_ctb_filter_index = alf_info->alf_ctb_filter_index;
  short *alf_ctb_filter_set_index_tmp = alf_info->alf_ctb_filter_set_index_tmp;
  alf_aps *alf_param_temp = &alf_info->alf_param_temp;

  cabac_data_t ctx_start;
  cabac_data_t *cabac_estimator = &alf_info->cabac_estimator;
  memcpy(&ctx_start, cabac_estimator, sizeof(ctx_start));
  cabac_data_t ctx_best;
  memcpy(&ctx_best, &ctx_start, sizeof(ctx_best));
  cabac_data_t ctx_temp_start;
  cabac_data_t ctx_temp_best;
  memcpy(&ctx_temp_best, &ctx_start, sizeof(ctx_temp_best));
  cabac_data_t ctx_temp_alt_start;
  //cabac_data_t ctx_temp_alt_best;

  int best_aps_ids[ALF_CTB_MAX_NUM_APS] = { -1, -1, -1, -1, -1, -1, -1, -1 };
  int filter_tmp[MAX_NUM_ALF_LUMA_COEFF] = { 0 };
  int g_clip_tmp[MAX_NUM_ALF_LUMA_COEFF] = { 0 };
  int size_of_best_aps_ids = 0;
  int clip_default[MAX_NUM_ALF_LUMA_COEFF] = { 0 };
  const int8_t bit_depth = state->encoder_control->bitdepth;
  double lambda = state->frame->lambda;
  int size_of_aps_ids = 0;
  const int32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  alf_aps alf_aps_temp_nl;
  alf_covariance *alf_cov_chroma;
  int cov_indx = 0;

  //AlfSliceParam  alfSliceParamNewFiltersBest = alfSliceParamNewFilters;
  alf_aps alf_aps_new_filters_best;
  copy_alf_param(&alf_aps_new_filters_best, aps);
  alf_aps* apss = state->slice->apss;

  bool has_new_filters[2] = { aps->enabled_flag[COMPONENT_Y] , aps->enabled_flag[COMPONENT_Cb] || aps->enabled_flag[COMPONENT_Cr] };

  //initDistortion();
  for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
  {
    ctb_distortions_unfilter[COMPONENT_Y][ctu_idx] = get_unfiltered_distortion_cov_classes(&alf_info->alf_covariance_y[ctu_idx], MAX_NUM_ALF_CLASSES);
    ctb_distortions_unfilter[COMPONENT_Cb][ctu_idx] = get_unfiltered_distortion_cov_classes(&alf_info->alf_covariance_u[ctu_idx], 1);
    ctb_distortions_unfilter[COMPONENT_Cr][ctu_idx] = get_unfiltered_distortion_cov_classes(&alf_info->alf_covariance_v[ctu_idx], 1);
  }

  //luma
  copy_alf_param(alf_param_temp, aps);
  memset(ctu_enable_flag[COMPONENT_Y], 1, sizeof(bool) * num_ctus_in_pic);
  get_frame_stats(alf_info, CHANNEL_TYPE_LUMA, num_ctus_in_pic);
  memset(ctu_enable_flag[COMPONENT_Y], 0, sizeof(bool) * num_ctus_in_pic);

  double cost_off = get_unfiltered_distortion_cov_channel(alf_info->alf_covariance_frame_luma, CHANNEL_TYPE_LUMA);

  int new_aps_id;
  int aps_ids[ALF_CTB_MAX_NUM_APS];
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    aps_ids[i] = -1;
  }

  alf_get_avai_aps_ids_luma(state, &new_aps_id, aps_ids, &size_of_aps_ids, arr_vars->alf_clipping_values);

  double cost_min = MAX_DOUBLE;
  alf_reconstruct_coeff_aps(state, true, false, true, arr_vars);

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
        bits_new_filter = arr_bits_new_filter[CHANNEL_TYPE_LUMA];
        alf_reconstruct_coeff(state, aps, CHANNEL_TYPE_LUMA, true, true, arr_vars);
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
        copy_alf_param(alf_param_temp, aps);
        alf_param_temp->enabled_flag[CHANNEL_TYPE_LUMA] = true;
        double cur_cost = 3 * lambda;

        if (iter > 0)  //re-derive new filter-set
        {
          double d_dist_org_new_filter = 0;
          int blocks_using_new_filter = 0;
          for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
          {
            if (ctu_enable_flag[COMPONENT_Y][ctu_idx] && alf_ctb_filter_index[ctu_idx] != ALF_NUM_FIXED_FILTER_SETS)
            {
              ctu_enable_flag[COMPONENT_Y][ctu_idx] = 0;
            }
            else if (ctu_enable_flag[COMPONENT_Y][ctu_idx] && alf_ctb_filter_index[ctu_idx] == ALF_NUM_FIXED_FILTER_SETS)
            {
              blocks_using_new_filter++;
              d_dist_org_new_filter += ctb_distortions_unfilter[COMPONENT_Y][ctu_idx];
              cov_indx = ctu_idx * MAX_NUM_ALF_CLASSES;
              for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
              {
                short* p_coeff = arr_vars->coeff_final;
                int16_t* p_clipp = arr_vars->clipp_final;
                for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
                {
                  filter_tmp[i] = p_coeff[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                  g_clip_tmp[i] = p_clipp[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                }
                d_dist_org_new_filter += calc_error_for_coeffs(&alf_info->alf_covariance_y[cov_indx + class_idx], g_clip_tmp, filter_tmp, MAX_NUM_ALF_LUMA_COEFF, bit_depth);
              }
            }
          } //for ctb
          if (blocks_using_new_filter > 0 && blocks_using_new_filter < num_ctus_in_pic)
          {
            int bit_nl[2] = { 0, 0 };
            double err_nl[2] = { 0.0, 0.0 };
            alf_param_temp->non_linear_flag[CHANNEL_TYPE_LUMA] = 1;
            if (state->encoder_control->cfg.alf_non_linear_luma)
            {
              err_nl[1] = alf_get_filter_coeff_and_cost(state, CHANNEL_TYPE_LUMA, 0, &bit_nl[1], true, true, arr_vars);
              copy_alf_param(&alf_aps_temp_nl, alf_param_temp);
            }
            else
            {
              err_nl[1] = MAX_DOUBLE;
            }
            alf_param_temp->non_linear_flag[CHANNEL_TYPE_LUMA] = 0;
            err_nl[0] = alf_get_filter_coeff_and_cost(state, CHANNEL_TYPE_LUMA, 0, &bit_nl[0], true, true, arr_vars);

            int bits_new_filter_temp_luma = bit_nl[0];
            double err = err_nl[0];
            if (err_nl[1] < err_nl[0])
            {
              err = err_nl[1];
              bits_new_filter_temp_luma = bit_nl[1];
              copy_alf_param(alf_param_temp, &alf_aps_temp_nl);
            }
            if (d_dist_org_new_filter + lambda * arr_bits_new_filter[CHANNEL_TYPE_LUMA] < err) //re-derived filter is not good, skip
            {
              continue;
            }
            alf_reconstruct_coeff(state, alf_param_temp, CHANNEL_TYPE_LUMA, true, true, arr_vars);
            bits_new_filter = bits_new_filter_temp_luma;
          }
          else //no blocks using new filter, skip
          {
            continue;
          }
        }

        //m_CABACEstimator->getCtx() = ctxStart;
        memcpy(cabac_estimator, &ctx_start, sizeof(*cabac_estimator));
        for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
        {
          double dist_unfilter_ctb = ctb_distortions_unfilter[COMPONENT_Y][ctu_idx];
          //ctb on
          ctu_enable_flag[COMPONENT_Y][ctu_idx] = 1;
          double cost_on = MAX_DOUBLE;
          //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
          memcpy(&ctx_temp_start, cabac_estimator, sizeof(ctx_temp_start));
          ctx_temp_start.only_count = 1;
          int i_best_filter_set_idx = 0;
          int first_filter_set_idx = 0;
          if (!state->encoder_control->cfg.alf_allow_predefined_filters)
          {
            first_filter_set_idx = ALF_NUM_FIXED_FILTER_SETS;
          }
          for (int filter_set_idx = first_filter_set_idx; filter_set_idx < num_filter_set; filter_set_idx++)
          {
            //rate
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
            memcpy(cabac_estimator, &ctx_temp_start, sizeof(*cabac_estimator));
            //m_CABACEstimator->resetBits();
            alf_cabac_reset_bits(cabac_estimator);
            //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, COMPONENT_Y, &m_alfSliceParamTemp);
            code_alf_ctu_enable_flag(state, cabac_estimator, ctu_idx, COMPONENT_Y, alf_param_temp);
            alf_ctb_filter_index[ctu_idx] = filter_set_idx;
            code_alf_ctu_filter_index(state, cabac_estimator, ctu_idx, alf_param_temp->enabled_flag[COMPONENT_Y]);

            double rate_on = (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3); //frac_bits_scale * 0; /*(double)m_CABACEstimator->getEstFracBits()*/ ;
            //distortion
            double dist = dist_unfilter_ctb;
            cov_indx = ctu_idx * MAX_NUM_ALF_CLASSES;
            for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
            {
              if (filter_set_idx < ALF_NUM_FIXED_FILTER_SETS)
              {
                int filter_idx = g_class_to_filter_mapping[filter_set_idx][class_idx];
                dist += calc_error_for_coeffs(&alf_info->alf_covariance_y[cov_indx + class_idx], clip_default, g_fixed_filter_set_coeff[filter_idx], MAX_NUM_ALF_LUMA_COEFF, bit_depth);
              }
              else
              {
                short *p_coeff;
                int16_t *p_clipp;
                if (use_new_filter && filter_set_idx == ALF_NUM_FIXED_FILTER_SETS)
                {
                  p_coeff = arr_vars->coeff_final;
                  p_clipp = arr_vars->clipp_final;
                }
                else if (use_new_filter)
                {
                  p_coeff = arr_vars->coeff_aps_luma[filter_set_idx - 1 - ALF_NUM_FIXED_FILTER_SETS];
                  p_clipp = arr_vars->clipp_aps_luma[filter_set_idx - 1 - ALF_NUM_FIXED_FILTER_SETS];
                }
                else
                {
                  p_coeff = arr_vars->coeff_aps_luma[filter_set_idx - ALF_NUM_FIXED_FILTER_SETS];
                  p_clipp = arr_vars->clipp_aps_luma[filter_set_idx - ALF_NUM_FIXED_FILTER_SETS];
                }
                for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
                {
                  filter_tmp[i] = p_coeff[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                  g_clip_tmp[i] = p_clipp[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                }
                dist += calc_error_for_coeffs(&alf_info->alf_covariance_y[cov_indx + class_idx], g_clip_tmp, filter_tmp, MAX_NUM_ALF_LUMA_COEFF, bit_depth);
              }
            }
            //cost
            double cost_on_tmp = dist + lambda * rate_on;
            if (cost_on_tmp < cost_on)
            {
              //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
              memcpy(&ctx_temp_best, cabac_estimator, sizeof(ctx_temp_best));
              ctx_temp_best.only_count = 1;
              cost_on = cost_on_tmp;
              i_best_filter_set_idx = filter_set_idx;
            }
          }
          //ctb off
          ctu_enable_flag[COMPONENT_Y][ctu_idx] = 0;
          //rate
          //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
          memcpy(cabac_estimator, &ctx_temp_start, sizeof(*cabac_estimator));
          //m_CABACEstimator->resetBits();
          alf_cabac_reset_bits(cabac_estimator);
          //m_CABACEstimator->codeAlfCtuEnableFlag(cs, ctbIdx, COMPONENT_Y, &m_alfSliceParamTemp);
          code_alf_ctu_enable_flag(state, cabac_estimator, ctu_idx, COMPONENT_Y, alf_param_temp);

          //cost
          double cost_off = dist_unfilter_ctb + lambda * (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3);// frac_bits_scale * 0; /* (double)m_CABACEstimator->getEstFracBits()*/ ;
          if (cost_on < cost_off)
          {
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
            memcpy(cabac_estimator, &ctx_temp_best, sizeof(*cabac_estimator));
            ctu_enable_flag[COMPONENT_Y][ctu_idx] = 1;
            alf_ctb_filter_index[ctu_idx] = i_best_filter_set_idx;
            cur_cost += cost_on;
          }
          else
          {
            ctu_enable_flag[COMPONENT_Y][ctu_idx] = 0;
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
          copy_alf_param(&alf_aps_new_filters_best, alf_param_temp);

          //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
          memcpy(&ctx_best, cabac_estimator, sizeof(ctx_best));
          //copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, CHANNEL_TYPE_LUMA);
          memcpy(ctu_enable_flag_tmp[COMPONENT_Y], ctu_enable_flag[COMPONENT_Y], sizeof(bool) * num_ctus_in_pic);
          for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
          {
            alf_ctb_filter_set_index_tmp[ctu_idx] = alf_ctb_filter_index[ctu_idx];
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
      memset(ctu_enable_flag[i], 0, sizeof(bool) * num_ctus_in_pic);
    }
    return;
  }
  else
  {
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
    copy_ctu_enable_flag(ctu_enable_flag, ctu_enable_flag_tmp, CHANNEL_TYPE_LUMA, num_ctus_in_pic);
    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
    {
      alf_ctb_filter_index[ctu_idx] = alf_ctb_filter_set_index_tmp[ctu_idx];
    }

    if (alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA])
    {
      //APS* newAPS = m_apsMap->getPS((new_aps_id << NUM_APS_TYPE_LEN) + T_ALF_APS);
      alf_aps* new_aps = &state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
      if (new_aps->aps_id < 0 || new_aps->aps_id >= ALF_CTB_MAX_NUM_APS) // new_aps == NULL
      {
        //newAPS = m_apsMap->allocatePS(new_aps_id);
        assert(new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS < ALF_CTB_MAX_NUM_APS); //Invalid PS id
        bool found = false;
        for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++) {
          if (state->encoder_control->cfg.param_set_map[i].parameter_set.aps_id == new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS) {
            found = true;
          }
        }
        if (!found) {
          state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
          //state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN+ T_ALF_APS].p_nalu_data = 0;
          //state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set = malloc(sizeof(alf_aps));
          state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.aps_id = new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS;
        }
        //copy_alf_param(new_aps, &state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set);
        new_aps->aps_id = new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS;
        new_aps->aps_type = T_ALF_APS;
      }
      copy_alf_param(new_aps, &alf_aps_new_filters_best);
      new_aps->temporal_id = state->slice->id;
      new_aps->new_filter_flag[CHANNEL_TYPE_CHROMA] = false;
      state->encoder_control->cfg.param_set_map[new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
      alf_info->aps_id_start = new_aps_id;
    }

    int8_t *aps_ids = state->slice->tile_group_luma_aps_id;
    for (int i = 0; i < state->slice->tile_group_num_aps; i++)
    {
      copy_aps(&apss[aps_ids[i]], &state->encoder_control->cfg.param_set_map[aps_ids[i] + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set, false);
    }
  }

  //chroma
  if (state->encoder_control->chroma_format != KVZ_CSP_400)
  {
    copy_alf_param(alf_param_temp, &alf_aps_new_filters_best);
    if (alf_param_temp->num_alternatives_chroma < 1)
    {
      alf_param_temp->num_alternatives_chroma = 1;
    }
    for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
    {
      ctu_alternatives[COMPONENT_Cb][ctu_idx] = 0;
      ctu_alternatives[COMPONENT_Cr][ctu_idx] = 0;
    }
    set_ctu_enable_flag(ctu_enable_flag, CHANNEL_TYPE_CHROMA, 1, num_ctus_in_pic);
    get_frame_stats(alf_info, CHANNEL_TYPE_CHROMA, num_ctus_in_pic);
    cost_off = get_unfiltered_distortion_cov_channel(alf_info->alf_covariance_frame_chroma, CHANNEL_TYPE_CHROMA);
    cost_min = MAX_DOUBLE;
    //m_CABACEstimator->getCtx() = AlfCtx(ctxBest);
    memcpy(cabac_estimator, &ctx_best, sizeof(*cabac_estimator));
    //ctxStart = AlfCtx(m_CABACEstimator->getCtx());
    memcpy(&ctx_start, cabac_estimator, sizeof(ctx_start));
    ctx_start.only_count = 1;
    int new_aps_id_chroma = -1;
    if (alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA] && (alf_aps_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_aps_new_filters_best.enabled_flag[COMPONENT_Cr]))
    {
      new_aps_id_chroma = new_aps_id;
    }
    else if (alf_aps_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_aps_new_filters_best.enabled_flag[COMPONENT_Cr])
    {
      int cur_id = alf_info->aps_id_start;
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
      const bool reuse_existing_aps = cur_aps_id != new_aps_id_chroma;

      if ((/*(cs.slice->getPendingRasInit() ||*/ state->frame->is_irap) && reuse_existing_aps)
      {
        continue;
      }

      alf_aps* cur_aps = &state->encoder_control->cfg.param_set_map[cur_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
      double cur_cost = lambda * 3;

      if (!reuse_existing_aps)
      {
        copy_alf_param(alf_param_temp, aps);
        cur_cost += lambda * arr_bits_new_filter[CHANNEL_TYPE_CHROMA];
      }
      else if (cur_aps && cur_aps->temporal_id <= state->slice->id
        && cur_aps->layer_id != 0
        && cur_aps->new_filter_flag[CHANNEL_TYPE_CHROMA])
      {
        //g_alf_slice_aps_temp = cur_aps;
        copy_alf_param(alf_param_temp, cur_aps);
      }
      else
      {
        continue;
      }
      alf_reconstruct_coeff(state, alf_param_temp, CHANNEL_TYPE_CHROMA, true, true, arr_vars);
      //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
      memcpy(cabac_estimator, &ctx_start, sizeof(*cabac_estimator));
      for (int comp_id = 1; comp_id < MAX_NUM_COMPONENT; comp_id++)
      {
        alf_cov_chroma = comp_id == COMPONENT_Cb ? alf_info->alf_covariance_u : alf_info->alf_covariance_v;
        alf_param_temp->enabled_flag[comp_id] = true;
        for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
        {
          double dist_unfilter_ctu = ctb_distortions_unfilter[comp_id][ctu_idx];
          //cost on
          ctu_enable_flag[comp_id][ctu_idx] = 1;
          memcpy(&ctx_temp_start, cabac_estimator, sizeof(ctx_temp_start));
          ctx_temp_start.only_count = 1;
          //rate
          //memcpy(&cabac_estimator, &ctx_temp_start, sizeof(cabac_estimator));
          alf_cabac_reset_bits(cabac_estimator);
          //ctb flag
          code_alf_ctu_enable_flag(state, cabac_estimator, ctu_idx, comp_id, alf_param_temp);
          double rate_on = (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3); //frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;
          const double ctu_lambda = lambda_chroma_weight > 0.0 ? 0/*cs.picture->m_uEnerHpCtu[ctbIdx]*/ / lambda_chroma_weight : lambda;
          double dist = MAX_DOUBLE;
          int num_alts = alf_param_temp->num_alternatives_chroma;
          //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
          memcpy(&ctx_temp_best, cabac_estimator, sizeof(ctx_temp_best));
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
              memcpy(cabac_estimator, &ctx_temp_alt_start, sizeof(*cabac_estimator));
            }
            //m_CABACEstimator->resetBits();
            alf_cabac_reset_bits(cabac_estimator);
            ctu_alternatives[comp_id][ctu_idx] = alt_idx;
            //m_CABACEstimator->codeAlfCtuAlternative(cs, ctbIdx, compId, &m_alfParamTemp);
            code_alf_ctu_alternative_ctu(state, cabac_estimator, ctu_idx, comp_id, alf_param_temp);
            double alt_rate = (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3); //frac_bits_scale * 0/*m_CABACEstimator->getEstFracBits()*/;
            double r_alt_cost = ctu_lambda * alt_rate;

            //distortion
            for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
            {
              filter_tmp[i] = arr_vars->chroma_coeff_final[alt_idx][i];
              g_clip_tmp[i] = arr_vars->chroma_clipp_final[alt_idx][i];
            }
            double alt_dist = calc_error_for_coeffs(&alf_cov_chroma[ctu_idx], g_clip_tmp, filter_tmp, MAX_NUM_ALF_CHROMA_COEFF, bit_depth);
            double alt_cost = alt_dist + r_alt_cost;
            if (alt_cost < best_alt_cost)
            {
              best_alt_cost = alt_cost;
              best_alt_idx = alt_idx;
              best_alt_rate = alt_rate;
              //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
              memcpy(&ctx_temp_best, cabac_estimator, sizeof(ctx_temp_best));
              ctx_temp_best.only_count = 1;
              dist = alt_dist;
            }
          }
          ctu_alternatives[comp_id][ctu_idx] = best_alt_idx;
          rate_on += best_alt_rate;
          dist += dist_unfilter_ctu;
          //cost
          double cost_on = dist + ctu_lambda * rate_on;
          //cost off
          ctu_enable_flag[comp_id][ctu_idx] = 0;
          //rate
          memcpy(cabac_estimator, &ctx_temp_start, sizeof(*cabac_estimator));
          alf_cabac_reset_bits(cabac_estimator);
          code_alf_ctu_enable_flag(state, cabac_estimator, ctu_idx, comp_id, alf_param_temp);
          //cost
          double cost_off = dist_unfilter_ctu + lambda * (23 - cabac_estimator->bits_left) + (cabac_estimator->num_buffered_bytes << 3); //frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;
          if (cost_on < cost_off)
          {
            //m_CABACEstimator->getCtx() = AlfCtx(ctxTempBest);
            memcpy(cabac_estimator, &ctx_temp_best, sizeof(*cabac_estimator));
            ctu_enable_flag[comp_id][ctu_idx] = 1;
            cur_cost += cost_on;
          }
          else
          {
            ctu_enable_flag[comp_id][ctu_idx] = 0;
            cur_cost += cost_off;
          }
        }//ctb_idx
      }//comp_id
      //chroma idc
      //setEnableFlag(m_alfSliceParamTemp, CHANNEL_TYPE_CHROMA, m_ctuEnableFlag);
      for (int comp_id = COMPONENT_Cb; comp_id <= COMPONENT_Cr; comp_id++)
      {
        alf_param_temp->enabled_flag[comp_id] = false;
        for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
        {
          if (ctu_enable_flag[comp_id][ctu_idx])
          {
            alf_param_temp->enabled_flag[comp_id] = true;
            break;
          }
        }
      }

      if (cur_cost < cost_min)
      {
        cost_min = cur_cost;
        state->slice->tile_group_chroma_aps_id = cur_aps_id;
        state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = alf_param_temp->enabled_flag[COMPONENT_Cb];
        state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = alf_param_temp->enabled_flag[COMPONENT_Cr];
        copy_ctu_enable_flag(ctu_enable_flag_tmp, ctu_enable_flag, CHANNEL_TYPE_CHROMA, num_ctus_in_pic);

        for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
        {
          ctu_alternatives_tmp[COMPONENT_Cb][ctu_idx] = ctu_alternatives[COMPONENT_Cb][ctu_idx];
          ctu_alternatives_tmp[COMPONENT_Cr][ctu_idx] = ctu_alternatives[COMPONENT_Cr][ctu_idx];
        }
      }
    }

    if (new_aps_id_chroma >= 0)
    {
      state->slice->tile_group_cc_alf_cb_aps_id = new_aps_id_chroma;
      state->slice->tile_group_cc_alf_cr_aps_id = new_aps_id_chroma;
    }
    if (cost_off < cost_min)
    {
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = false;
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = false;
      set_ctu_enable_flag(ctu_enable_flag, CHANNEL_TYPE_CHROMA, 0, num_ctus_in_pic);
    }
    else
    {
      copy_ctu_enable_flag(ctu_enable_flag, ctu_enable_flag_tmp, CHANNEL_TYPE_CHROMA, num_ctus_in_pic);
      for (int ctu_idx = 0; ctu_idx < num_ctus_in_pic; ctu_idx++)
      {
        ctu_alternatives[COMPONENT_Cb][ctu_idx] = ctu_alternatives_tmp[COMPONENT_Cb][ctu_idx];
        ctu_alternatives[COMPONENT_Cr][ctu_idx] = ctu_alternatives_tmp[COMPONENT_Cr][ctu_idx];
      }
      if (state->slice->tile_group_chroma_aps_id == new_aps_id_chroma)  //new filter
      {
        //APS* newAPS = m_apsMap->getPS(new_aps_id_chroma);
        alf_aps* new_aps = &state->encoder_control->cfg.param_set_map[new_aps_id_chroma + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set;
        if (new_aps->aps_id < 0 || new_aps->aps_id >= ALF_CTB_MAX_NUM_APS) //new_aps == NULL
        {
          //newAPS = m_apsMap->allocatePS(new_aps_id);
          assert(new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS < ALF_CTB_MAX_NUM_APS); //Invalid PS id
          bool found = false;
          for (int i = 0; i < (sizeof(state->encoder_control->cfg.param_set_map) / sizeof(state->encoder_control->cfg.param_set_map[0])); i++) {
            if (state->encoder_control->cfg.param_set_map[i].parameter_set.aps_id == new_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS) {
              found = true;
            }
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

        new_aps->num_alternatives_chroma = aps->num_alternatives_chroma;
        new_aps->non_linear_flag[CHANNEL_TYPE_CHROMA] = aps->non_linear_flag[CHANNEL_TYPE_CHROMA];
        new_aps->temporal_id = state->slice->id;
        for (int alt_idx = 0; alt_idx < MAX_NUM_ALF_ALTERNATIVES_CHROMA; ++alt_idx)
        {
          for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
          {
            new_aps->chroma_coeff[alt_idx][i] = aps->chroma_coeff[alt_idx][i];
            new_aps->chroma_clipp[alt_idx][i] = aps->chroma_clipp[alt_idx][i];
          }
        }
        state->encoder_control->cfg.param_set_map[new_aps_id_chroma + NUM_APS_TYPE_LEN + T_ALF_APS].b_changed = true;
        alf_info->aps_id_start = new_aps_id_chroma;
      }
      apss[state->slice->tile_group_chroma_aps_id].aps_id = state->encoder_control->cfg.param_set_map[state->slice->tile_group_chroma_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.aps_id;
      apss[state->slice->tile_group_chroma_aps_id].aps_type = state->encoder_control->cfg.param_set_map[state->slice->tile_group_chroma_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set.aps_type;
      copy_alf_param(&apss[state->slice->tile_group_chroma_aps_id], &state->encoder_control->cfg.param_set_map[state->slice->tile_group_chroma_aps_id + NUM_APS_TYPE_LEN + T_ALF_APS].parameter_set);
    }
  }
}

static void alf_filter_block(encoder_state_t * const state,
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
      p_class = state->tile->frame->alf_info->classifier[blk_dst_y + i] + blk_dst_x;
    }

    for (int j = 0; j < end_width - start_width; j += cls_size_x)
    {
      if (!chroma)
      {
        alf_classifier cl = p_class[j];
        transpose_idx = cl.transpose_idx;
        coef = filter_set + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
        clip = fClipSet + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
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

          if (!(is_near_vb_above || is_near_vb_below))
          {
            sum = (sum + offset) >> shift;
          }
          else
          {
            sum = (sum + (1 << ((shift + 3) - 1))) >> (shift + 3);
          }
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

static void alf_reconstruct(encoder_state_t * const state,
  array_variables *arr_vars)
{
  if (!state->slice->tile_group_alf_enabled_flag[COMPONENT_Y])
  {
    return;
  }

  alf_reconstruct_coeff_aps(state, true, state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] || state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr], false, arr_vars);

  alf_info_t *alf_info = state->tile->frame->alf_info;
  bool **ctu_enable_flags = alf_info->ctu_enable_flag;
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
  memcpy(&alf_info->alf_tmp_y[index_luma], &state->tile->frame->rec->y[index_luma],
    sizeof(kvz_pixel) * luma_stride * (luma_height + MAX_ALF_PADDING_SIZE * 2));
  memcpy(&alf_info->alf_tmp_u[index_chroma], &state->tile->frame->rec->u[index_chroma],
    sizeof(kvz_pixel) * chroma_stride * (chroma_height + chroma_padding * 2));
  memcpy(&alf_info->alf_tmp_v[index_chroma], &state->tile->frame->rec->v[index_chroma],
    sizeof(kvz_pixel) * chroma_stride * (chroma_height + chroma_padding * 2));

  for (int y_pos = 0; y_pos < luma_height; y_pos += max_cu_height)
  {
    for (int x_pos = 0; x_pos < luma_width; x_pos += max_cu_width)
    {

      const int width = (x_pos + max_cu_width > luma_width) ? (luma_width - x_pos) : max_cu_width;
      const int height = (y_pos + max_cu_height > luma_height) ? (luma_height - y_pos) : max_cu_height;

      bool ctu_enable_flag = ctu_enable_flags[COMPONENT_Y][ctu_idx];
      for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
      {
        ctu_enable_flag |= ctu_enable_flags[comp_idx][ctu_idx] > 0;
      }

      {
        if (ctu_enable_flags[COMPONENT_Y][ctu_idx])
        {
          short filter_set_index = alf_info->alf_ctb_filter_index[ctu_idx];
          short *coeff;
          int16_t *clip;
          if (filter_set_index >= ALF_NUM_FIXED_FILTER_SETS)
          {
            coeff = arr_vars->coeff_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
            clip = arr_vars->clipp_aps_luma[filter_set_index - ALF_NUM_FIXED_FILTER_SETS];
          }
          else
          {
            coeff = arr_vars->fixed_filter_set_coeff_dec[filter_set_index];
            clip = arr_vars->clip_default;
          }
          alf_filter_block(state,
            alf_info->alf_tmp_y, state->tile->frame->rec->y,
            luma_stride, luma_stride,
            coeff, clip, arr_vars->clp_rngs.comp[COMPONENT_Y], COMPONENT_Y,
            width, height, x_pos, y_pos, x_pos, y_pos,
            alf_vb_luma_pos, alf_vb_luma_ctu_height);
        }
        for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
        {
          alf_component_id comp_id = comp_idx;

          if (ctu_enable_flags[comp_idx][ctu_idx])
          {
            kvz_pixel *dst_pixels = comp_id - 1 ? state->tile->frame->rec->v : state->tile->frame->rec->u;
            const kvz_pixel *src_pixels = comp_id - 1 ? alf_info->alf_tmp_v : alf_info->alf_tmp_u;

            const int alt_num = alf_info->ctu_alternative[comp_id][ctu_idx];
            alf_filter_block(state,
              src_pixels, dst_pixels,
              chroma_stride, chroma_stride,
              arr_vars->chroma_coeff_final[alt_num], arr_vars->chroma_clipp_final[alt_num], arr_vars->clp_rngs.comp[comp_idx], comp_idx,
              width >> chroma_scale_x, height >> chroma_scale_y,
              x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
              x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
              alf_vb_chma_pos, alf_vb_chma_ctu_height);
          }
        }
      }
      ctu_idx++;
    }
  }
}

static void alf_derive_classification_blk(encoder_state_t * const state,
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
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;

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

    int *p_y_ver = laplacian[ALF_VER][i];
    int *p_y_hor = laplacian[ALF_HOR][i];
    int *p_y_dig_0 = laplacian[ALF_DIAG0][i];
    int *p_y_dig_1 = laplacian[ALF_DIAG1][i];

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
      p_y_dig_0[j] = abs(y0 - p_y_down[-1] - p_y_up[1]) + abs(y_up1 - p_y[0] - p_y_up2[2]);
      p_y_dig_1[j] = abs(y0 - p_y_up[-1] - p_y_down[1]) + abs(y_up1 - p_y_up2[0] - p_y[2]);

      if (j > 4 && (j - 6) % 4 == 0)
      {
        int j_m_6 = j - 6;
        int j_m_4 = j - 4;
        int j_m_2 = j - 2;

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
        activity = alf_clip3_int(0, max_activity, (temp_act * 96) >> shift);
      }
      else
      {
        activity = alf_clip3_int(0, max_activity, (temp_act * 64) >> shift);
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
      if ((uint32_t)d1 * (uint32_t)hv0 > (uint32_t)hv1 * (uint32_t)d0)
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

static void alf_derive_classification(encoder_state_t * const state,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  const int blk_dst_x,
  const int blk_dst_y)
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

  adjust_pixels(state->tile->frame->rec->y, x_pos, pic_width, y_pos, pic_height, state->tile->frame->rec->stride,
    pic_width, pic_height);
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

      alf_derive_classification_blk(state, state->encoder_control->cfg.input_bitdepth + 4, n_height, n_width, j, i,
        j - x_pos + blk_dst_x, i - y_pos + blk_dst_y,
        alf_vb_luma_ctu_height,
        alf_vb_luma_pos);
    }
  }
}

void kvz_alf_enc_process(encoder_state_t *const state)
{
  alf_info_t *alf_info = state->tile->frame->alf_info;
  alf_create_frame_buffer(state, alf_info);

  if (1 /*!layerIdx*/ && (false/*cs.slice->getPendingRasInit()*/ || (state->frame->pictype == KVZ_NAL_IDR_W_RADL || state->frame->pictype == KVZ_NAL_IDR_N_LP)))
  {
    for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++) {
      reset_aps(&state->slice->apss[i], state->encoder_control->cfg.alf_type == KVZ_ALF_FULL);
      if (state->encoder_control->cfg.param_set_map[i + T_ALF_APS].b_changed)
      {
        alf_aps* alf_aps = &state->encoder_control->cfg.param_set_map[i + T_ALF_APS].parameter_set;
        state->encoder_control->cfg.param_set_map[i + T_ALF_APS].b_changed = false;
        reset_aps(alf_aps, state->encoder_control->cfg.alf_type == KVZ_ALF_FULL);
      }
    }
    alf_info->aps_id_start = ALF_CTB_MAX_NUM_APS;
  }

  alf_aps alf_param;
  reset_alf_param(&alf_param);
  cc_alf_filter_param *cc_filter_param = state->slice->cc_filter_param;


  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  bool chroma_scale_x = (chroma_fmt == KVZ_CSP_444) ? 0 : 1;
  bool chroma_scale_y = (chroma_fmt != KVZ_CSP_420) ? 0 : 1;
  int8_t kvz_bit_depth = state->encoder_control->bitdepth;
  const int32_t num_ctus_in_pic = state->tile->frame->width_in_lcu * state->tile->frame->height_in_lcu;
  const int8_t input_bitdepth = state->encoder_control->bitdepth;
  double lambda_chroma_weight = 0.0;

  cabac_data_t ctx_start;
  cabac_data_t ctx_start_cc_alf;
  cabac_data_t *cabac_estimator = &alf_info->cabac_estimator;
  memcpy(cabac_estimator, &state->cabac, sizeof(*cabac_estimator));
  memcpy(&ctx_start, &state->cabac, sizeof(ctx_start));
  memcpy(&ctx_start_cc_alf, cabac_estimator, sizeof(ctx_start_cc_alf));
  cabac_estimator->only_count = 1;
  ctx_start.only_count = 1;
  ctx_start_cc_alf.only_count = 1;

  // derive classification
  const int luma_height = state->tile->frame->height;
  const int luma_width = state->tile->frame->width;

  static array_variables arr_vars;
  static bool init_values = false;

  if (!init_values)
  {
    assert(MAX_ALF_NUM_CLIPPING_VALUES > 0); //"g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] must be at least one"
    arr_vars.alf_clipping_values[CHANNEL_TYPE_LUMA][0] = 1 << input_bitdepth;
    int shift_luma = input_bitdepth - 8;
    for (int i = 1; i < MAX_ALF_NUM_CLIPPING_VALUES; ++i)
    {
      arr_vars.alf_clipping_values[CHANNEL_TYPE_LUMA][i] = 1 << (7 - 2 * i + shift_luma);
    }

    assert(MAX_ALF_NUM_CLIPPING_VALUES > 0); //"g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] must be at least one"
    arr_vars.alf_clipping_values[CHANNEL_TYPE_CHROMA][0] = 1 << input_bitdepth;
    int shift_chroma = input_bitdepth - 8;
    for (int i = 1; i < MAX_ALF_NUM_CLIPPING_VALUES; ++i)
    {
      arr_vars.alf_clipping_values[CHANNEL_TYPE_CHROMA][i] = 1 << (7 - 2 * i + shift_chroma);
    }

    for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF * MAX_NUM_ALF_CLASSES; i++)
    {
      arr_vars.clip_default[i] = arr_vars.alf_clipping_values[CHANNEL_TYPE_LUMA][0];
    }

    for (int filter_set_index = 0; filter_set_index < ALF_NUM_FIXED_FILTER_SETS; filter_set_index++)
    {
      for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
      {
        int fixed_filter_idx = g_class_to_filter_mapping[filter_set_index][class_idx];
        for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF - 1; i++)
        {
          arr_vars.fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + i] = g_fixed_filter_set_coeff[fixed_filter_idx][i];
        }
        arr_vars.fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + MAX_NUM_ALF_LUMA_COEFF - 1] = (1 << (input_bitdepth - 1));
      }
    }

    //Default clp_rng
    arr_vars.clp_rngs.comp[COMPONENT_Y].min = arr_vars.clp_rngs.comp[COMPONENT_Cb].min = arr_vars.clp_rngs.comp[COMPONENT_Cr].min = 0;
    arr_vars.clp_rngs.comp[COMPONENT_Y].max = (1 << kvz_bit_depth) - 1;
    arr_vars.clp_rngs.comp[COMPONENT_Y].bd = kvz_bit_depth;
    arr_vars.clp_rngs.comp[COMPONENT_Y].n = 0;
    arr_vars.clp_rngs.comp[COMPONENT_Cb].max = arr_vars.clp_rngs.comp[COMPONENT_Cr].max = (1 << kvz_bit_depth) - 1;
    arr_vars.clp_rngs.comp[COMPONENT_Cb].bd = arr_vars.clp_rngs.comp[COMPONENT_Cr].bd = kvz_bit_depth;
    arr_vars.clp_rngs.comp[COMPONENT_Cb].n = arr_vars.clp_rngs.comp[COMPONENT_Cr].n = 0;
    arr_vars.clp_rngs.used = arr_vars.clp_rngs.chroma = false;

    init_values = true;
  }

  for (int y_pos = 0; y_pos < luma_height; y_pos += LCU_WIDTH)
  {
    for (int x_pos = 0; x_pos < luma_width; x_pos += LCU_WIDTH)
    {
      const int width = (x_pos + LCU_WIDTH > luma_width) ? (luma_width - x_pos) : LCU_WIDTH;
      const int height = (y_pos + LCU_WIDTH > luma_height) ? (luma_height - y_pos) : LCU_WIDTH;
      {
        alf_derive_classification(state, width, height, x_pos, y_pos, x_pos, y_pos);
      }
    }
  }


  // get CTB stats for filtering 
  alf_derive_stats_for_filtering(state, arr_vars.alf_clipping_values);

  for (int ctb_iIdx = 0; ctb_iIdx < num_ctus_in_pic; ctb_iIdx++)
  {
    alf_info->alf_ctb_filter_index[ctb_iIdx] = ALF_NUM_FIXED_FILTER_SETS;
  }

  // consider using new filter (only)
  alf_param.new_filter_flag[CHANNEL_TYPE_LUMA] = true;
  alf_param.new_filter_flag[CHANNEL_TYPE_CHROMA] = true;
  state->slice->tile_group_num_aps = 1; // Only new filter for RD cost optimization

  // derive filter (luma)
  alf_encoder(state,
    &alf_param, CHANNEL_TYPE_LUMA,
    lambda_chroma_weight,
    &arr_vars
  );

  // derive filter (chroma)
  if (state->encoder_control->chroma_format != KVZ_CSP_400) {
    alf_encoder(state,
      &alf_param, CHANNEL_TYPE_CHROMA,
      lambda_chroma_weight,
      &arr_vars
    );
  }
  // let alfEncoderCtb decide now
  alf_param.new_filter_flag[CHANNEL_TYPE_LUMA] = false;
  alf_param.new_filter_flag[CHANNEL_TYPE_CHROMA] = false;
  state->slice->tile_group_num_aps = 0;

  //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
  memcpy(cabac_estimator, &ctx_start, sizeof(*cabac_estimator));
  alf_encoder_ctb(state, &alf_param, lambda_chroma_weight, &arr_vars);

  //for (int s = 0; s < state.; s++) //numSliceSegments
  {
    if (state->encoder_control->cfg.lossless)
    {
      for (uint32_t ctb_idx = 0; ctb_idx < num_ctus_in_pic; ctb_idx++) //pcPic->slices[s]->getNumCtuInSlice()
      {
        //uint32_t ctuRsAddr = pcPic->slices[s]->getCtuAddrInSlice(ctuIdx);
        state->tile->frame->alf_info->ctu_enable_flag[COMPONENT_Y][ctb_idx] = 0;
        state->tile->frame->alf_info->ctu_enable_flag[COMPONENT_Cb][ctb_idx] = 0;
        state->tile->frame->alf_info->ctu_enable_flag[COMPONENT_Cr][ctb_idx] = 0;
      }
    }
  }

  alf_reconstruct(state, &arr_vars);

  if (state->encoder_control->cfg.alf_type != KVZ_ALF_FULL)
  {
    return;
  }

  // Do not transmit CC ALF if it is unchanged
  if (state->slice->tile_group_alf_enabled_flag[COMPONENT_Y])
  {
    for (int32_t luma_alf_aps_id = 0; luma_alf_aps_id < state->slice->tile_group_num_aps; luma_alf_aps_id++)
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

  const kvz_picture *org_yuv = state->tile->frame->source;
  const kvz_picture *rec_yuv = state->tile->frame->rec;

  const int luma_stride = state->tile->frame->rec->stride;
  const int chroma_stride = luma_stride >> chroma_scale_x;
  const int chroma_height = luma_height >> chroma_scale_y;
  const int chroma_padding = MAX_ALF_PADDING_SIZE >> chroma_scale_x;

  const int index_chroma = -(chroma_stride * chroma_padding + chroma_padding);

  //Copy reconstructed samples to a buffer.
  memcpy(&alf_info->alf_tmp_u[index_chroma], &state->tile->frame->rec->u[index_chroma],
    sizeof(kvz_pixel) * chroma_stride * (chroma_height + chroma_padding * 2));
  memcpy(&alf_info->alf_tmp_v[index_chroma], &state->tile->frame->rec->v[index_chroma],
    sizeof(kvz_pixel) * chroma_stride * (chroma_height + chroma_padding * 2));

  adjust_pixels_chroma(alf_info->alf_tmp_u,
    0,
    rec_yuv->width >> chroma_scale_x,
    0,
    rec_yuv->height >> chroma_scale_y,
    rec_yuv->stride >> chroma_scale_x,
    rec_yuv->width >> chroma_scale_x,
    rec_yuv->height >> chroma_scale_y);
  adjust_pixels_chroma(alf_info->alf_tmp_v,
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
  init_distortion_cc_alf(alf_info->alf_covariance_cc_alf, alf_info->ctb_distortion_unfilter, num_ctus_in_pic);

  memcpy(cabac_estimator, &ctx_start_cc_alf, sizeof(*cabac_estimator));
  derive_cc_alf_filter(state, COMPONENT_Cb, org_yuv, rec_yuv, arr_vars.cc_reuse_aps_id);
  memcpy(cabac_estimator, &ctx_start_cc_alf, sizeof(*cabac_estimator));
  derive_cc_alf_filter(state, COMPONENT_Cr, org_yuv, rec_yuv, arr_vars.cc_reuse_aps_id);

  setup_cc_alf_aps(state, arr_vars.cc_reuse_aps_id);

  for (alf_component_id comp_idx = 1; comp_idx < (state->encoder_control->chroma_format == KVZ_CSP_400 ? 1 : MAX_NUM_COMPONENT); comp_idx++)
  {
    if (cc_filter_param->cc_alf_filter_enabled[comp_idx - 1])
    {
      kvz_pixel* rec_uv = comp_idx == COMPONENT_Cb ? rec_yuv->u : rec_yuv->v;
      const int luma_stride = rec_yuv->stride;
      apply_cc_alf_filter(state, comp_idx, rec_uv, alf_info->alf_tmp_y, luma_stride, alf_info->cc_alf_filter_control[comp_idx - 1],
        cc_filter_param->cc_alf_coeff[comp_idx - 1], -1, &arr_vars);
    }
  }
}