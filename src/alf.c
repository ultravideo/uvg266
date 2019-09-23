

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

bool is_crossed_by_virtual_boundaries(const int x_pos, const int y_pos, const int width, const int height, bool* clip_top, bool* clip_bottom, bool* clip_left, bool* clip_right, 
                                      int* num_hor_vir_bndry, int* num_ver_vir_bndry, int hor_vir_bndry_pos[], int ver_vir_bndry_pos[], encoder_state_t *const state)
{
  *clip_top = false; *clip_bottom = false; *clip_left = false; *clip_right = false;
  *num_hor_vir_bndry = 0; *num_ver_vir_bndry = 0;
  if (state->encoder_control->cfg.loop_filter_across_virtual_boundaries_disabled_flag)
  {
    for (int i = 0; i < state->slice->num_hor_virtual_boundaries; i++)
    {
      if (state->slice->virtual_boundaries_pos_y[i] == y_pos)
      {
        *clip_top = true;
      }
      else if (state->slice->virtual_boundaries_pos_y[i] == y_pos + height)
      {
        *clip_bottom = true;
      }
      else if (y_pos < state->slice->virtual_boundaries_pos_y[i] && state->slice->virtual_boundaries_pos_y[i] < y_pos + height)
      {
        hor_vir_bndry_pos[*num_hor_vir_bndry++] = state->slice->virtual_boundaries_pos_y[i];
      }
    }
    for (int i = 0; i < state->slice->num_ver_virtual_boundaries; i++)
    {
      if (state->slice->virtual_boundaries_pos_x[i] == x_pos)
      {
        *clip_left = true;
      }
      else if (state->slice->virtual_boundaries_pos_x[i] == x_pos + width)
      {
        *clip_right = true;
      }
      else if (x_pos < state->slice->virtual_boundaries_pos_x[i] && state->slice->virtual_boundaries_pos_x[i] < x_pos + width)
      {
        ver_vir_bndry_pos[*num_ver_vir_bndry++] = state->slice->virtual_boundaries_pos_x[i];
      }
    }
  }
  return *num_hor_vir_bndry > 0 || *num_ver_vir_bndry > 0 || *clip_top || *clip_bottom || *clip_left || *clip_right;
}

int clip_alf(const int clip, const short ref, const short val0, const short val1)
{
  return alf_clip3(-clip, +clip, val0 - ref) + alf_clip3(-clip, +clip, val1 - ref);
}

int alf_clip_pixel(const int a, const clp_rng clp_rng)
{
  return MIN(MAX(clp_rng.min, a), clp_rng.max);
}

int alf_clip3(const int minVal, const int maxVal, const int a)
{
  return MIN(MAX(minVal, a), maxVal);
}

void get_clip_max(alf_covariance *cov, int *clip_max)
{
  const int num_coeff = cov->num_coeff;
  for (int k = 0; k < num_coeff - 1; ++k)
  {
    clip_max[k] = 0;

    bool inc = true;
    while (inc && clip_max[k] + 1 < cov->num_bins && cov->y[clip_max[k] + 1][k] == cov->y[clip_max[k]][k])
    {
      for (int l = 0; inc && l < num_coeff; ++l)
        if (cov->ee[clip_max[k]][0][k][l] != cov->ee[clip_max[k] + 1][0][k][l])
        {
          inc = false;
        }
      if (inc)
      {
        ++clip_max[k];
      }
    }
  }
  clip_max[num_coeff - 1] = 0;
}

void reduce_clip_cost(alf_covariance *cov, int *clip)
{
  for (int k = 0; k < cov->num_coeff - 1; ++k)
  {
    bool dec = true;
    while (dec && clip[k] > 0 && cov->y[clip[k] - 1][k] == cov->y[clip[k]][k])
    {
      for (int l = 0; dec && l < cov->num_coeff; ++l)
        if (cov->ee[clip[k]][clip[l]][k][l] != cov->ee[clip[k] - 1][clip[l]][k][l])
        {
          dec = false;
        }
      if (dec)
      {
        --clip[k];
      }
    }
  }
}

void set_ey_from_clip(alf_covariance *cov,const int* clip, double ee[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double y[MAX_NUM_ALF_LUMA_COEFF], int size)
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

double optimize_filter(alf_covariance *cov, int* clip, double *f, bool optimize_clip)
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
  static double inv_diag[MAX_NUM_ALF_LUMA_COEFF];  /* Vector of the inverse of diagonal entries of outMatr */

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

int gns_solve_by_chol(double lhs[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF], double rhs[MAX_NUM_ALF_LUMA_COEFF], double *x, int num_eq)
{
  static double aux[MAX_NUM_ALF_LUMA_COEFF];     /* Auxiliary vector */
  static double u[MAX_NUM_ALF_LUMA_COEFF][MAX_NUM_ALF_LUMA_COEFF];    /* Upper triangular Cholesky factor of lhs */
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

double calc_error_for_coeffs(alf_covariance *cov, const int *clip, const int *coeff, const int num_coeff, const int bit_depth)
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

/*#if !JVET_O0216_ALF_COEFF_EG3 || !JVET_O0064_SIMP_ALF_CLIP_CODING
int get_golomb_k_min(channel_type channel, const int num_filters, int k_min_tab[MAX_NUM_ALF_LUMA_COEFF], int bits_coeff_scan[m_MAX_SCAN_VAL][m_MAX_EXP_GOLOMB])
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
}*/

int length_golomb(int coeff_val, int k, bool signed_coeff)
{
  int num_bins = 0;
  unsigned int symbol = abs(coeff_val);
  while (symbol >= (unsigned int)(1 << k))
  {
    num_bins++;
    symbol -= 1 << k;
    k++;
  }
  num_bins += (k + 1);
  if (signed_coeff && coeff_val != 0)
  {
    num_bins++;
  }
  return num_bins;
}

double get_dist_coeff_force_0(bool* coded_var_bins, double error_force_0_coeff_tab[MAX_NUM_ALF_CLASSES][2], int* bits_var_bin, const int num_filters)
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

double get_dist_force_0(channel_type channel, const int num_filters, double error_tab_force_0_coeff[MAX_NUM_ALF_CLASSES][2], bool* coded_var_bins)
{
  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  const int *golomb_idx = channel == CHANNEL_TYPE_LUMA ? alf_golomb_idx_7 : alf_golomb_idx_5;

  static int bits_var_bin[MAX_NUM_ALF_CLASSES];

  for (int ind = 0; ind < num_filters; ++ind)
  {
    bits_var_bin[ind] = 0;
    for (int i = 0; i < num_coeff - 1; i++)
    {
      bits_var_bin[ind] += length_golomb(abs(g_filter_coeff_set[ind][i]), 3, true);
    }
  }

  if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
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
        len += length_golomb(abs(p_diff_q_filter_coeff_int_pp[ind][i]), 3, true); // alf_coeff_luma_delta[i][j]
      }
    }
  }

  if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
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
  int num_luma_filters = aps->num_luma_filters;
  //int fixed_filter_pattern = aps->fixed_filter_pattern;
  //int fixed_filter_set_index = aps->fixed_filter_set_index;

  int len = 1   // alf_coefficients_delta_flag
          + length_truncated_unary(0, 3)    // chroma_idc = 0, it is signalled when ALF is enabled for luma
          + get_tb_length(num_luma_filters - 1, MAX_NUM_ALF_CLASSES);   //numLumaFilters

  if (num_luma_filters > 1)
  {
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
    {
      len += get_tb_length((int)filter_coeff_delta_idx[i], num_luma_filters);  //filter_coeff_delta[i]
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
      bit_cnt += length_golomb(abs(filter_coeff[ind][i]), 3, true);
    }
  }
  return bit_cnt;
}

double calculate_error(alf_covariance *cov, const int *clip, const double *coeff)
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

int get_coeff_rate(alf_aps* alf_aps, bool is_chroma)
{
  int iBits = 0;
  assert(is_chroma);
  int num_coeff = is_chroma ? 7 : 13;
  // Filter coefficients
  for (int i = 0; i < num_coeff - 1; i++)
  {
    iBits += length_golomb(alf_aps->chroma_coeff[i], 3, true);  // alf_coeff_chroma[i], alf_coeff_luma_delta[i][j]
  }

  if (g_alf_aps_temp.non_linear_flag[is_chroma])
  {
    for (int i = 0; i < num_coeff - 1; i++)
    {
      if (!abs(alf_aps->chroma_coeff[i]))
      {
        alf_aps->chroma_clipp[i] = 0;
      }
    }
    iBits += ((num_coeff - 1) << 1);
  }
  return iBits;
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

double get_filtered_distortion(alf_covariance* cov, const int num_classes, const int num_filters_minus1, const int num_coeff)
{
  double dist = 0;

  for (int class_idx = 0; class_idx < num_classes; class_idx++)
  {
    dist += calc_error_for_coeffs(&cov[class_idx], g_filter_clipp_set[class_idx], g_filter_coeff_set[class_idx], num_coeff, ALF_NUM_BITS);
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
    dist = get_unfiltered_distortion_cov_classes(cov, 1) + length_truncated_unary(0, 3) * g_lambda[COMPONENT_Cb];
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

void get_frame_stats(channel_type channel, int i_shape_idx)
{
  int num_classes = channel == CHANNEL_TYPE_LUMA ? MAX_NUM_ALF_CLASSES : 1;
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
}

void get_frame_stat(alf_covariance* frame_cov, alf_covariance** ctb_cov, uint8_t* ctb_enable_flags, const int num_classes)
{
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
}

void copy_alf_param(alf_aps *dst, alf_aps *src)
{
  memcpy(dst->enabled_flag, src->enabled_flag, sizeof(dst->enabled_flag));
  memcpy(dst->non_linear_flag, src->non_linear_flag, sizeof(dst->non_linear_flag));
  memcpy(dst->luma_coeff, src->luma_coeff, sizeof(dst->luma_coeff));
  memcpy(dst->luma_clipp, src->luma_clipp, sizeof(dst->luma_clipp));
  memcpy(dst->chroma_coeff, src->chroma_coeff, sizeof(dst->chroma_coeff));
  memcpy(dst->chroma_clipp, src->chroma_clipp, sizeof(dst->chroma_clipp));
  memcpy(dst->filter_coeff_delta_idx, src->filter_coeff_delta_idx, sizeof(dst->filter_coeff_delta_idx));
  memcpy(dst->alf_luma_coeff_flag, src->alf_luma_coeff_flag, sizeof(dst->alf_luma_coeff_flag));
  dst->num_luma_filters = src->num_luma_filters;
  dst->alf_luma_coeff_delta_flag = src->alf_luma_coeff_delta_flag;
  //dst->alf_luma_coeff_delta_prediction_flag = src->alf_luma_coeff_delta_prediction_flag;
  dst->t_layer = src->t_layer;
  memcpy(dst->new_filter_flag, src->new_filter_flag, sizeof(dst->new_filter_flag));
  //dst->fixed_filter_pattern = src->fixed_filter_pattern;
  //memcpy(dst->fixed_filter_idx, src->fixed_filter_idx, sizeof(dst->fixed_filter_idx));
  //dst->fixed_filter_set_index = src->fixed_filter_set_index;
}

void reset_alf_param(alf_aps *src)
{
  memset(src->enabled_flag, false, sizeof(src->enabled_flag));
  memset(src->non_linear_flag, false, sizeof(src->non_linear_flag));
  memset(src->luma_coeff, 0, sizeof(src->luma_coeff));
  memset(src->luma_clipp, 0, sizeof(src->luma_clipp));
  memset(src->chroma_coeff, 0, sizeof(src->chroma_coeff));
  memset(src->chroma_clipp, 0, sizeof(src->chroma_clipp));
  memset(src->filter_coeff_delta_idx, 0, sizeof(src->filter_coeff_delta_idx));
  memset(src->alf_luma_coeff_flag, true, sizeof(src->alf_luma_coeff_flag));
  src->num_luma_filters = 1;
  src->alf_luma_coeff_delta_flag = false;
  src->t_layer = 0;
  memset(src->new_filter_flag, 0, sizeof(src->new_filter_flag));
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

void reset_alf_covariance(alf_covariance *alf, int num_bins) {
  if (num_bins > 0) { alf->num_bins = num_bins; }
  alf->pix_acc = 0;
  memset(alf->y, 0, sizeof(alf->y));
  memset(alf->ee, 0, sizeof(alf->ee));
}

//-------------------------------------------------------------------

//-------------------------encoding functions------------------------

void kvz_alf_enc_process(encoder_state_t *const state,
  const lcu_order_element_t *const lcu)
{
  //kerran muualla 
  if (false/*cs.slice->getPendingRasInit()*/ || true/*cs.slice->isIDRorBLA()*/)
  {
    //memset(state->slice->apss, 0, sizeof(alf_aps) * ALF_CTB_MAX_NUM_APS);
    for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++) {
      state->slice->apss[i].aps_id = 0;
      state->slice->apss[i].aps_type = 0;
      reset_alf_param(&state->slice->apss[i]);
      state->slice->apss[i].num_luma_filters = 0;

    }
    g_aps_id_start = ALF_CTB_MAX_NUM_APS;
    //state->slice->aps = malloc(sizeof(alf_aps*) * MAX_NUM_APS);
    
    //FREE_POINTER(state->slice->param_set_map);

    /* T‰ss‰ ei pit‰isi koskaan tapahtua mit‰‰n, koska yl‰puolella olisi FREE_POINTER
    for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
    {
      alf_aps* alf_aps = &state->slice->param_set_map[i + T_ALF_APS].parameter_set;
      state->slice->param_set_map[i + T_ALF_APS].b_changed = false;
      if (alf_aps)
      {
        reset_alf_param(alf_aps);
        alf_aps = NULL;
      }
    }*/
  }

  //luotava vain kerran?
  alf_aps alf_param;
  reset_alf_param(&alf_param);

  //const TempCtx  ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));

  //muualla kerran
  // set CTU ALF enable flags, it was already reset before ALF process
  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++) {
      g_ctu_enable_flag[comp_idx][ctu_idx] = 0;
    }
  }
  int width = state->tile->frame->width;
  int height = state->tile->frame->height;

  //turha
  //reset_alf_param(&alf_param);

  //int shift_luma = 2 * 0;// DISTORTION_PRECISION_ADJUSTMENT(m_inputBitDepth[CHANNEL_TYPE_LUMA]);
  //int shift_chroma = 2 * 0;// DISTORTION_PRECISION_ADJUSTMENT(m_inputBitDepth[CHANNEL_TYPE_CHROMA]);
  g_lambda[COMPONENT_Y] = state->lambda;// 981.62883931057581 // lambdas[COMPONENT_Y] * double(1 << shiftLuma);
  g_lambda[COMPONENT_Cb] = state->lambda;// 1236.7748378312663 // lambdas[COMPONENT_Cb] * double(1 << shiftChroma);
  g_lambda[COMPONENT_Cr] = state->lambda;// 1236.7748378312663 //lambdas[COMPONENT_Cr] * double(1 << shiftChroma);

  //PelUnitBuf orgYuv = cs.getOrgBuf();

  //m_tempBuf.copyFrom(cs.getRecoBuf());
  //PelUnitBuf recYuv = m_tempBuf.getBuf(cs.area);
  //recYuv.extendBorderPel(MAX_ALF_FILTER_LENGTH >> 1);

  // derive classification
  //const kvz_pixel rec_luma = state->tile->frame->rec->y;
  //const PreCalcValues& pcv = *cs.pcv;

  //turhat
  const int lumaHeight = state->tile->frame->height;
  const int lumaWidth = state->tile->frame->width;
  const int maxCUHeight = LCU_WIDTH;
  const int maxCUWidth = LCU_WIDTH;

  bool clip_top = false, clip_bottom = false, clip_left = false, clip_right = false;
  int num_hor_vir_bndry = 0, num_ver_vir_bndry = 0;
  int hor_vir_bndry_pos[] = { 0, 0, 0 };
  int ver_vir_bndry_pos[] = { 0, 0, 0 };
  
  //turha
  for (int y_pos = 0; y_pos < lumaHeight; y_pos += maxCUHeight)
  {
    for (int x_pos = 0; x_pos < lumaWidth; x_pos += maxCUWidth)
    {
      //lcu->size n‰iden tilalle
      const int width = (x_pos + maxCUWidth > lumaWidth) ? (lumaWidth - x_pos) : maxCUWidth;
      const int height = (y_pos + maxCUHeight > lumaHeight) ? (lumaHeight - y_pos) : maxCUHeight;

      //T‰t‰ algoritmia pit‰‰ viel‰ viilata
      if (is_crossed_by_virtual_boundaries(x_pos, y_pos, width, height, &clip_top, &clip_bottom, &clip_left, &clip_right, &num_hor_vir_bndry, &num_ver_vir_bndry, hor_vir_bndry_pos, ver_vir_bndry_pos, state))
      {
        int y_start = y_pos;
        for (int i = 0; i <= num_hor_vir_bndry; i++)
        {
          const int y_end = i == num_hor_vir_bndry ? y_pos + height : hor_vir_bndry_pos[i];
          const int h = y_end - y_start;
          const bool clip_t = (i == 0 && clip_top) || (i > 0) || (y_start == 0);
          const bool clip_b = (i == num_hor_vir_bndry && clip_bottom) || (i < num_hor_vir_bndry) || (y_end == lumaHeight);

          int x_start = x_pos;
          for (int j = 0; j <= num_ver_vir_bndry; j++)
          {
            const int x_end = j == num_ver_vir_bndry ? x_pos + width : ver_vir_bndry_pos[j];
            const int w = x_end - x_start;
            const bool clip_l = (j == 0 && clip_left) || (j > 0) || (x_start == 0);
            const bool clip_r = (j == num_ver_vir_bndry && clip_right) || (j < num_ver_vir_bndry) || (x_end == lumaWidth);

            //const int w_buf = w + (clip_l ? 0 : MAX_ALF_PADDING_SIZE) + (clip_r ? 0 : MAX_ALF_PADDING_SIZE);
            //const int h_buf = h + (clip_t ? 0 : MAX_ALF_PADDING_SIZE) + (clip_b ? 0 : MAX_ALF_PADDING_SIZE);
            //PelUnitBuf buf = m_tempBuf2.subBuf(UnitArea(cs.area.chromaFormat, Area(0, 0, w_buf, h_buf)));
            //buf.copyFrom(recYuv.subBuf(UnitArea(cs.area.chromaFormat, Area(x_start - (clip_l ? 0 : MAX_ALF_PADDING_SIZE), y_start - (clip_t ? 0 : MAX_ALF_PADDING_SIZE), w_buf, h_buf))));
            //buf.extendBorderPel(MAX_ALF_PADDING_SIZE);
            //buf = buf.subBuf(UnitArea(cs.area.chromaFormat, Area(clip_l ? 0 : MAX_ALF_PADDING_SIZE, clip_t ? 0 : MAX_ALF_PADDING_SIZE, w, h)));

            //const Area blkSrc(0, 0, w, h);
            //const Area blkDst(xStart, yStart, w, h);
            kvz_alf_derive_classification(state, lcu, w, h, x_start, y_start, x_start, y_start);
            //Area blkPCM(xStart, yStart, w, h);
            kvz_alf_reset_pcm_blk_class_info(state, lcu, w, h, x_start, y_start);

            x_start = x_end;
          }

          y_start = y_end;
        }
      }
      else
      {
        //Area blk(x_pos, y_pos, width, height);
        kvz_alf_derive_classification(state, lcu, width, height, x_pos, y_pos, x_pos, y_pos);
        //Area blkPCM(x_pos, y_pos, width, height);
        kvz_alf_reset_pcm_blk_class_info(state, lcu, width, height, x_pos, y_pos);
      }

    } //
  } //

  // get CTB stats for filtering 
  kvz_alf_derive_stats_for_filtering(state, lcu); //checked

  // derive filter (luma)
  kvz_alf_encoder(state, lcu, &alf_param, CHANNEL_TYPE_LUMA); //ulkopuolelle

  // derive filter (chroma)
  kvz_alf_encoder(state, lcu, &alf_param, CHANNEL_TYPE_CHROMA); //ulkopuolelle

  //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
  kvz_alf_encoder_ctb(state, lcu, &alf_param);

  kvz_alf_reconstructor(state, lcu);
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

  if (is_luma) {
    g_alf_aps_temp.enabled_flag[COMPONENT_Y] = 1;
  }
  else {
    g_alf_aps_temp.enabled_flag[COMPONENT_Cb] = 1;
    g_alf_aps_temp.enabled_flag[COMPONENT_Cr] = 1;
  }

  kvz_alf_reconstruct_coeff(state, &g_alf_aps_temp, channel, true, is_luma);

  for (int class_idx = 0; class_idx < (is_luma ? MAX_NUM_ALF_CLASSES : 1); class_idx++)
  {
    for (int i = 0; i < (is_luma ? MAX_NUM_ALF_LUMA_COEFF : MAX_NUM_ALF_CHROMA_COEFF); i++)
    {
      g_filter_coeff_set[class_idx][i] = is_luma ? g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + i] : g_chroma_coeff_final[i];
      g_filter_clipp_set[class_idx][i] = is_luma ? g_clipp_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + i] : g_chroma_clipp_final[i];
    }
  }

  for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++)
  {
    for (int comp_id = comp_id_first; comp_id <= comp_id_last; comp_id++)
    {
      double dist_unfilter_ctu = get_unfiltered_distortion_cov_classes(g_alf_covariance[comp_id][i_shape_idx][ctu_idx], num_classes);

      //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
      //m_CABACEstimator->resetBits();
      g_ctu_enable_flag[comp_id][ctu_idx] = 1;
      code_alf_ctu_enable_flag(state, lcu, ctu_idx, comp_id, &g_alf_aps_temp);
      double cost_on = dist_unfilter_ctu + get_filtered_distortion(g_alf_covariance[comp_id][i_shape_idx][ctu_idx], num_classes, g_alf_aps_temp.num_luma_filters - 1, num_coeff);

      const double ctu_lambda = g_lambda[comp_id];

      cost_on += ctu_lambda * frac_bits_scale * 0;// m_CABACEstimator->getEstFracBits();

      //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());

      //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
      //m_CABACEstimator->resetBits();
      g_ctu_enable_flag[comp_id][ctu_idx] = 0;
      code_alf_ctu_enable_flag(state, lcu, ctu_idx, comp_id, &g_alf_aps_temp);
      double cost_off = dist_unfilter_ctu + ctu_lambda * frac_bits_scale * 0;// m_CABACEstimator->getEstFracBits();

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
    const alf_component_id compIDFirst = COMPONENT_Cb;
    const alf_component_id compIDLast = COMPONENT_Cr;
    for (int compId = compIDFirst; compId <= compIDLast; compId++)
    {
      g_alf_aps_temp.enabled_flag[compId] = false;
      for (int i = 0; i < g_num_ctus_in_pic; i++)
      {
        if (g_ctu_enable_flag[compId][i])
        {
          g_alf_aps_temp.enabled_flag[compId] = true;
          break;
        }
      }
    }
    const int alf_chroma_idc = g_alf_aps_temp.enabled_flag[COMPONENT_Cb] * 2 + g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
    cost += length_truncated_unary(alf_chroma_idc, 3) * g_lambda[channel];
  }

  return cost;
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
    //m_alfCovarianceFrame[ch_type] = new AlfCovariance*[m_filterShapes[ch_type].size()];
    g_alf_covariance_frame[ch_type] = malloc(sizeof(**g_alf_covariance_frame[ch_type]));
    for (int i = 0; i != 1/*m_filterShapes[ch_type].size()*/; i++)
    {
      g_alf_covariance_frame[ch_type][i] = malloc(num_classes * sizeof(alf_covariance));
      for (int k = 0; k < num_classes; k++)
      {
        g_alf_covariance_frame[ch_type][i][k].num_coeff = num_coeffs;
        g_alf_covariance_frame[ch_type][i][k].num_bins = g_max_alf_num_clipping_values;
        memset(g_alf_covariance_frame[ch_type][i][k].y, 0, sizeof(g_alf_covariance_frame[ch_type][i][k].y));
        memset(g_alf_covariance_frame[ch_type][i][k].ee, 0, sizeof(g_alf_covariance_frame[ch_type][i][k].ee));
      }
    }
  }

  for (int comp_idx = 0; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
  {
    g_ctu_enable_flag[comp_idx] = malloc(g_num_ctus_in_pic * sizeof(*g_ctu_enable_flag[comp_idx]));
    g_ctu_enable_flag_tmp[comp_idx] = malloc(g_num_ctus_in_pic * sizeof(*g_ctu_enable_flag_tmp[comp_idx]));

    int num_classes = comp_idx ? 1 : MAX_NUM_ALF_CLASSES;
    int num_coeffs = comp_idx ? 7 : 13;

    g_alf_covariance[comp_idx] = malloc(sizeof(***g_alf_covariance[comp_idx]));

    for (int i = 0; i != 1/*m_filterShapes[ch_type].size()*/; i++)
    {
      g_alf_covariance[comp_idx][i] = malloc(g_num_ctus_in_pic * sizeof(**g_alf_covariance[comp_idx][i]));
      for (int j = 0; j < g_num_ctus_in_pic; j++)
      {
        g_alf_covariance[comp_idx][i][j] = malloc(num_classes * sizeof(alf_covariance));
        for (int k = 0; k < num_classes; k++)
        {
          g_alf_covariance[comp_idx][i][j][k].num_coeff = num_coeffs;
          g_alf_covariance[comp_idx][i][j][k].num_bins = g_max_alf_num_clipping_values;
          memset(g_alf_covariance[comp_idx][i][j][k].y, 0, sizeof(g_alf_covariance[comp_idx][i][j][k].y));
          memset(g_alf_covariance[comp_idx][i][j][k].ee, 0, sizeof(g_alf_covariance[comp_idx][i][j][k].ee));
        }
      }
    }
  }

  for (int i = 0; i != 1/*m_filterShapes[COMPONENT_Y].size()*/; i++)
  {
    for (int j = 0; j <= MAX_NUM_ALF_CLASSES + 1; j++)
    {
      g_alf_covariance_merged[i][j].num_coeff = 13;
      g_alf_covariance_merged[i][j].num_bins = g_max_alf_num_clipping_values;
      memset(g_alf_covariance_merged[i][j].y, 0, sizeof(g_alf_covariance_merged[i][j].y));
      memset(g_alf_covariance_merged[i][j].ee, 0, sizeof(g_alf_covariance_merged[i][j].ee));
    }
  }
  g_filter_coeff_set = malloc(MAX_NUM_ALF_CLASSES * sizeof(int*));
  g_filter_clipp_set = malloc(MAX_NUM_ALF_CLASSES * sizeof(int*));
  g_diff_filter_coeff = malloc(MAX_NUM_ALF_CLASSES * sizeof(int*));

  for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
  {
    g_filter_coeff_set[i] = malloc(MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
    g_filter_clipp_set[i] = malloc(MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
    g_diff_filter_coeff[i] = malloc(MAX_NUM_ALF_LUMA_COEFF * sizeof(int));
  }
  g_aps_id_start = ALF_CTB_MAX_NUM_APS;

  //g_ctb_distortion_fixed_filter = malloc(g_num_ctus_in_pic * sizeof(double));
  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    g_ctb_distortion_unfilter[comp] = malloc(g_num_ctus_in_pic * sizeof(double));
  }

  g_alf_ctb_filter_index = malloc(g_num_ctus_in_pic * sizeof(*g_alf_ctb_filter_index));
  g_alf_ctb_filter_set_index_tmp = malloc(g_num_ctus_in_pic * sizeof(*g_alf_ctb_filter_set_index_tmp));
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
      for (int i = 0; i != 1/*m_filterShapes[ch_type].size()*/; i++)
      {
        /*for (int k = 0; k < num_classes; k++)
        {
          
          g_alf_covariance_frame[channel_idx][i][k].pix_acc = 0;
          memset(g_alf_covariance_frame[channel_idx][i][k].y, 0, sizeof(g_alf_covariance_frame[channel_idx][i][k].y));
          memset(g_alf_covariance_frame[channel_idx][i][k].ee, 0, sizeof(g_alf_covariance_frame[channel_idx][i][k].ee));
        }*/
        FREE_POINTER(g_alf_covariance_frame[channel_idx][i]);
      }
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

    if (g_alf_covariance[comp_idx])
    {
      channel_type chType = comp_idx ? CHANNEL_TYPE_CHROMA : CHANNEL_TYPE_LUMA;
      int numClasses = comp_idx ? 1 : MAX_NUM_ALF_CLASSES;

      for (int i = 0; i != 1/*m_filterShapes[ch_type].size()*/; i++)
      {
        for (int j = 0; j < g_num_ctus_in_pic; j++)
        {
          /*
          for (int k = 0; k < num_classes; k++)
          {
            g_alf_covariance[comp_idx][i][j][k].pix_acc = 0;
            memset(g_alf_covariance[comp_idx][i][j][k].y, 0, sizeof(g_alf_covariance[comp_idx][i][j][k].y));
            memset(g_alf_covariance[comp_idx][i][j][k].ee, 0, sizeof(g_alf_covariance[comp_idx][i][j][k].ee));
          }*/
          FREE_POINTER(g_alf_covariance[comp_idx][i][j]);
        }
        FREE_POINTER(g_alf_covariance[comp_idx][i]);
      }
      FREE_POINTER(g_alf_covariance[comp_idx]);
    }
  }

  /*
  for (int i = 0; i != 1/*m_filterShapes[COMPONENT_Y].size()*//*; i++)
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
  }*/

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

  if (g_diff_filter_coeff)
  {
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++)
    {
      FREE_POINTER(g_diff_filter_coeff[i]);
    }
    FREE_POINTER(g_diff_filter_coeff);
  }

  /*if (g_ctb_distortion_fixed_filter != NULL) {
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
  channel_type channel)
{
  //const TempCtx  ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));
  //TempCtx        ctxBest(m_CtxCache);

  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  kvz_config cfg = state->encoder_control->cfg;

  double cost_min = MAX_DOUBLE;
  g_bits_new_filter[channel] = 0;
  const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
  int ui_coeff_bits = 0;

  for (int i_shape_idx = 0; i_shape_idx < 1/*alfFilterShape.size()*/; i_shape_idx++)
  {
    //m_alfSliceParamTemp = alfSliceParam;
    copy_alf_param(&g_alf_aps_temp, aps);

    //1. get unfiltered distortion
    double cost = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[channel][i_shape_idx], channel);
    cost /= 1.001; // slight preference for unfiltered choice

    if (cost < cost_min)
    {
      cost_min = cost;
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
        memset(g_ctu_enable_flag_tmp[COMPONENT_Y], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
      }
      else {
        memset(g_ctu_enable_flag_tmp[COMPONENT_Cb], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
        memset(g_ctu_enable_flag_tmp[COMPONENT_Cr], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
      }
    }

    const int non_linear_flag_max = (is_luma ? cfg.alf_non_linear_luma : cfg.alf_non_linear_chroma) ? 2 : 1;

    for (int non_linear_flag = 0; non_linear_flag < non_linear_flag_max; non_linear_flag++)
    {
      //2. all CTUs are on
      //setEnableFlag(m_alfSliceParamTemp, channel, true);
      if (is_luma) {
        g_alf_aps_temp.enabled_flag[COMPONENT_Y] = 1;
      }
      else {
        g_alf_aps_temp.enabled_flag[COMPONENT_Cb] = 1;
        g_alf_aps_temp.enabled_flag[COMPONENT_Cr] = 1;
      }
      g_alf_aps_temp.non_linear_flag[channel] = non_linear_flag;

      //m_CABACEstimator->getCtx() = AlfCtx(ctxStart);
      //setCtuEnableFlag(m_ctuEnableFlag, channel, 1);
      if (is_luma) {
        memset(g_ctu_enable_flag[COMPONENT_Y], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
      }
      else {
        memset(g_ctu_enable_flag[COMPONENT_Cb], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
        memset(g_ctu_enable_flag[COMPONENT_Cr], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
      }
      cost = kvz_alf_get_filter_coeff_and_cost(state, lcu, channel, 0, &ui_coeff_bits, i_shape_idx, non_linear_flag != 0, false);

      if (cost < cost_min)
      {
        g_bits_new_filter[channel] = ui_coeff_bits;
        cost_min = cost;
        //copyAlfSliceParam(alfSliceParam, m_alfSliceParamTemp, channel);
        if (is_luma)
        {
          copy_alf_param(aps, &g_alf_aps_temp);
        }
        else
        {
          aps->non_linear_flag[channel] = g_alf_aps_temp.non_linear_flag[channel];
          aps->enabled_flag[COMPONENT_Cb] = g_alf_aps_temp.enabled_flag[COMPONENT_Cb];
          aps->enabled_flag[COMPONENT_Cr] = g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
          memcpy(aps->chroma_coeff, g_alf_aps_temp.chroma_coeff, sizeof(aps->chroma_coeff));
          memcpy(aps->chroma_clipp, g_alf_aps_temp.chroma_clipp, sizeof(aps->chroma_clipp));
        }
        //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
        //setCtuEnableFlag(m_ctuEnableFlagTmp, channel, 1);
        if (is_luma) {
          memset(g_ctu_enable_flag_tmp[COMPONENT_Y], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
        }
        else {
          memset(g_ctu_enable_flag_tmp[COMPONENT_Cb], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
          memset(g_ctu_enable_flag_tmp[COMPONENT_Cr], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
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
              memcpy(g_ctu_enable_flag_tmp[COMPONENT_Y], g_ctu_enable_flag[COMPONENT_Y], sizeof(uint8_t) * g_num_ctus_in_pic);
            }
            else {
              memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cr], g_ctu_enable_flag[COMPONENT_Cr], sizeof(uint8_t) * g_num_ctus_in_pic);
              memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cb], g_ctu_enable_flag[COMPONENT_Cb], sizeof(uint8_t) * g_num_ctus_in_pic);
            }

            //copyAlfSliceParam(alfSliceParam, m_alfSliceParamTemp, channel);
            if (is_luma)
            {
              copy_alf_param(aps, &g_alf_aps_temp);
            }
            else
            {
              aps->non_linear_flag[channel] = g_alf_aps_temp.non_linear_flag[channel];
              aps->enabled_flag[COMPONENT_Cb] = g_alf_aps_temp.enabled_flag[COMPONENT_Cb];
              aps->enabled_flag[COMPONENT_Cr] = g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
              memcpy(aps->chroma_coeff, g_alf_aps_temp.chroma_coeff, sizeof(aps->chroma_coeff));
              memcpy(aps->chroma_clipp, g_alf_aps_temp.chroma_clipp, sizeof(aps->chroma_clipp));
            }
          }
        }
        else
        {
          cost = kvz_alf_get_filter_coeff_and_cost(state, lcu, channel, dist_unfilter, &ui_coeff_bits, i_shape_idx, true, false);
        }
      }//for iter
      // Decrease number of alternatives and reset ctu params and filters
    }//for non_linea_flag
  }//for shape_idx
  //m_CABACEstimator->getCtx() = AlfCtx(ctxBest);
}

void kvz_alf_get_avai_aps_ids_luma(encoder_state_t *const state, 
  int *newApsId, 
  int aps_ids[ALF_CTB_MAX_NUM_APS],
  int *size_of_aps_ids)
{
  param_set_map *aps_set = state->slice->param_set_map;
  for (int i = 0; i < ALF_CTB_MAX_NUM_APS; i++)
  {
    state->slice->apss[i] = aps_set[i + T_ALF_APS].parameter_set;
  }

  //std::vector<int> result;
  int aps_id_checked = 0, cur_aps_id = g_aps_id_start;
  if (cur_aps_id < ALF_CTB_MAX_NUM_APS)
  {
    while (aps_id_checked < ALF_CTB_MAX_NUM_APS && /*!cs.slice->isIntra() &&*/ *size_of_aps_ids < ALF_CTB_MAX_NUM_APS /*&& /*!cs.slice->getPendingRasInit() &&*/ /*!cs.slice->isIDRorBLA()*/)
    {
      alf_aps *cur_aps = &state->slice->apss[cur_aps_id];

      if (cur_aps && cur_aps->t_layer/*cur_aps->getTemporalId()*/ <= state->slice->id/*cs.slice->getTLayer()*/ && cur_aps->new_filter_flag[CHANNEL_TYPE_LUMA])
      {
        //result.push_back(cur_aps_id);
        aps_ids[*size_of_aps_ids] = cur_aps_id;
        *size_of_aps_ids++;
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

  *newApsId = g_aps_id_start - 1;
  if (*newApsId < 0)
  {
    *newApsId = (int)ALF_CTB_MAX_NUM_APS - 1;
  }
  assert(*newApsId < (int)MAX_NUM_APS); //Wrong APS index assignment in getAvaiApsIdsLuma
}

void kvz_alf_derive_stats_for_filtering(encoder_state_t *const state,
  const lcu_order_element_t *const lcu)
{
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;

  //turhat
  int32_t pic_width = state->tile->frame->rec->width;
  int32_t pic_height = state->tile->frame->rec->height;
  int ctu_rs_addr = 0;

  const int number_of_components = (chroma_fmt == KVZ_CSP_400) ? 1 : MAX_NUM_COMPONENT;

  // init CTU stats buffers
  for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
  {
    bool is_luma = comp_idx == 0 ? 1 : 0;
    const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;

    for (int shape = 0; shape != 1 /*m_filterShapes[toChannelType(comp_id)].size()*/; shape++)
    {
      for (int class_idx = 0; class_idx < num_classes; class_idx++)
      {
        //t‰m‰ luup pois
        for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++)
        {

          reset_alf_covariance(&g_alf_covariance[comp_idx][shape][ctu_idx][class_idx],
            g_alf_num_clipping_values[comp_idx == COMPONENT_Y ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA]);
        
        } //
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
    for (int shape = 0; shape != 1/*m_filterShapes[channel_idx].size()*/; shape++)
    {
      for (int class_idx = 0; class_idx < num_classes; class_idx++)
      {
        reset_alf_covariance(&g_alf_covariance_frame[channel_idx][shape][class_idx], g_alf_num_clipping_values[channel_id]);
      }
    }
  }

  bool clip_top = false, clip_bottom = false, clip_left = false, clip_right = false;
  int num_hor_vir_bndry = 0, num_ver_vir_bndry = 0;
  int hor_vir_bndry_pos[] = { 0, 0, 0 };
  int ver_vir_bndry_pos[] = { 0, 0, 0 };

  int max_cu_width = LCU_WIDTH;
  int max_cu_height = LCU_WIDTH;

  //turhat
  for (int y_pos = 0; y_pos < pic_height; y_pos += max_cu_height)
  {
    for (int x_pos = 0; x_pos < pic_width; x_pos += max_cu_width)
    {//

      //lcu->size
      int width = (x_pos + max_cu_width > pic_width) ? (pic_width - x_pos) : max_cu_width;
      int height = (y_pos + max_cu_height > pic_height) ? (pic_height - y_pos) : max_cu_height;

      if (is_crossed_by_virtual_boundaries(x_pos, y_pos, width, height, &clip_top, &clip_bottom, &clip_left, &clip_right, &num_hor_vir_bndry, &num_ver_vir_bndry, hor_vir_bndry_pos, ver_vir_bndry_pos, state))
      {
        int y_start = y_pos;
        for (int i = 0; i <= num_hor_vir_bndry; i++)
        {
          const int y_end = i == num_hor_vir_bndry ? y_pos + height : hor_vir_bndry_pos[i];
          const int h = y_end - y_start;
          const bool clip_t = (i == 0 && clip_top) || (i > 0) || (y_start == 0);
          const bool clip_b = (i == num_hor_vir_bndry && clip_bottom) || (i < num_hor_vir_bndry) || (y_end == pic_height);

          int x_start = x_pos;
          for (int j = 0; j <= num_ver_vir_bndry; j++)
          {
            const int x_end = j == num_ver_vir_bndry ? x_pos + width : ver_vir_bndry_pos[j];
            const int w = x_end - x_start;
            const bool clip_l = (j == 0 && clip_left) || (j > 0) || (x_start == 0);
            const bool clip_r = (j == num_ver_vir_bndry && clip_right) || (j < num_ver_vir_bndry) || (x_end == pic_width);

            const int w_buf = w + (clip_l ? 0 : MAX_ALF_PADDING_SIZE) + (clip_r ? 0 : MAX_ALF_PADDING_SIZE);
            const int h_buf = h + (clip_t ? 0 : MAX_ALF_PADDING_SIZE) + (clip_b ? 0 : MAX_ALF_PADDING_SIZE);
            //PelUnitBuf recBuf = m_tempBuf2.subBuf(UnitArea(cs.area.chromaFormat, Area(0, 0, w_buf, h_buf)));
            //recBuf.copyFrom(recYuv.subBuf(UnitArea(cs.area.chromaFormat, Area(x_start - (clip_l ? 0 : MAX_ALF_PADDING_SIZE), y_start - (clip_t ? 0 : MAX_ALF_PADDING_SIZE), w_buf, h_buf))));
            //recBuf.extendBorderPel(MAX_ALF_PADDING_SIZE);
            //recBuf = recBuf.subBuf(UnitArea(cs.area.chromaFormat, Area(clip_l ? 0 : MAX_ALF_PADDING_SIZE, clip_t ? 0 : MAX_ALF_PADDING_SIZE, w, h)));

            //const UnitArea area(m_chromaFormat, Area(0, 0, w, h));
            //const UnitArea areaDst(m_chromaFormat, Area(x_start, y_start, w, h));
            for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
            {
              const bool is_luma = comp_idx == COMPONENT_Y ? 1 : 0;

              const int chroma_scale_x = ((is_luma) || (chroma_fmt == KVZ_CSP_444)) ? 0 : 1;
              const int chroma_scale_y = ((is_luma) || (chroma_fmt != KVZ_CSP_420)) ? 0 : 1;
              channel_type ch_type = is_luma ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;

              int blk_w = is_luma ? w : w >> chroma_scale_x;
              int blk_h = is_luma ? h : h >> chroma_scale_y;
              int pos_x = is_luma ? x_start : x_start >> chroma_scale_x;
              int pos_y = is_luma ? y_start : y_start >> chroma_scale_y;

              int32_t org_stride = is_luma ? state->tile->frame->source->stride : state->tile->frame->source->stride >> chroma_scale_x;
              int32_t rec_stride = is_luma ? state->tile->frame->rec->stride    : state->tile->frame->rec->stride >> chroma_scale_x;

              kvz_pixel *org = comp_idx ? (comp_idx - 1 ? &state->tile->frame->source->v[pos_x + pos_y * org_stride] : &state->tile->frame->source->u[pos_x + pos_y * org_stride]) : &state->tile->frame->source->y[pos_x + pos_y * org_stride];
              kvz_pixel *rec = comp_idx ? (comp_idx - 1 ? &state->tile->frame->rec->v[pos_x + pos_y * org_stride]    : &state->tile->frame->rec->u[pos_x + pos_y * org_stride])    : &state->tile->frame->rec->y[pos_x + pos_y * org_stride];

              for (int shape = 0; shape !=1/*m_filterShapes[ch_type].size()*/; shape++)
              {
                kvz_alf_get_blk_stats(state, lcu, ch_type, &g_alf_covariance[comp_idx][shape][ctu_rs_addr], comp_idx ? NULL : classifier,
                  org, org_stride, rec, rec_stride, pos_x, pos_y, pos_x, pos_y, blk_w, blk_h, 
                  ((comp_idx == 0) ? g_alf_vb_luma_ctu_height : g_alf_vb_chma_ctu_height), 
                  ((y_pos + max_cu_height >= pic_height) ? pic_height : ((comp_idx == 0) ? g_alf_vb_luma_pos : g_alf_vb_chma_pos))
                );
              }
            }

            x_start = x_end;
          }

          y_start = y_end;
        }

        for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
        {
          const bool is_luma = comp_idx == COMPONENT_Y ? 1 : 0;
          channel_type ch_type = is_luma ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;

          for (int shape = 0; shape != 1/*m_filterShapes[chType].size()*/; shape++)
          {
            const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;

            for (int class_idx = 0; class_idx < num_classes; class_idx++)
            {
              add_alf_cov(&g_alf_covariance_frame[ch_type][shape][class_idx], &g_alf_covariance[comp_idx][shape][ctu_rs_addr][class_idx]);
            }
          }
        }
      }
      else
      {
        for (int comp_idx = 0; comp_idx < number_of_components; comp_idx++)
        {
          const bool is_luma = comp_idx == COMPONENT_Y ? 1 : 0;

          const int chroma_scale_x = ((is_luma) || (chroma_fmt == KVZ_CSP_444)) ? 0 : 1;
          const int chroma_scale_y = ((is_luma) || (chroma_fmt != KVZ_CSP_420)) ? 0 : 1;
          channel_type ch_type = is_luma ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;

          int blk_w = is_luma ? width : width >> chroma_scale_x;
          int blk_h = is_luma ? height : height >> chroma_scale_y;
          int pos_x = is_luma ? x_pos : x_pos >> chroma_scale_x;
          int pos_y = is_luma ? y_pos : y_pos >> chroma_scale_y;

          int32_t org_stride = is_luma ? state->tile->frame->source->stride : state->tile->frame->source->stride >> chroma_scale_x;
          int32_t rec_stride = is_luma ? state->tile->frame->rec->stride    : state->tile->frame->rec->stride >> chroma_scale_x;

          kvz_pixel *org = comp_idx ? (comp_idx - 1 ? &state->tile->frame->source->v[pos_x + pos_y * org_stride] : &state->tile->frame->source->u[pos_x + pos_y * org_stride]) : &state->tile->frame->source->y[pos_x + pos_y * org_stride];
          kvz_pixel *rec = comp_idx ? (comp_idx - 1 ? &state->tile->frame->rec->v[pos_x + pos_y * org_stride]    : &state->tile->frame->rec->u[pos_x + pos_y * org_stride])    : &state->tile->frame->rec->y[pos_x + pos_y * org_stride];

          for (int shape = 0; shape != 1/*m_filterShapes[ch_type].size()*/; shape++)
          {
            g_alf_covariance[comp_idx][shape][ctu_rs_addr][0];
            kvz_alf_get_blk_stats(state, lcu, ch_type, &g_alf_covariance[comp_idx][shape][ctu_rs_addr], comp_idx ? NULL : classifier, org, org_stride, rec, rec_stride, pos_x, pos_y, pos_x, pos_y, blk_w, blk_h
              , (is_luma ? g_alf_vb_luma_ctu_height : g_alf_vb_chma_ctu_height)
              , ((y_pos + max_cu_height >= pic_height) ? pic_height : ((is_luma) ? g_alf_vb_luma_pos : g_alf_vb_chma_pos)));

            const int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;

            for (int class_idx = 0; class_idx < num_classes; class_idx++)
            {
              add_alf_cov(&g_alf_covariance_frame[ch_type][shape][class_idx], &g_alf_covariance[comp_idx][shape][ctu_rs_addr][class_idx]);
            }
          }
        }
      }
      //turhat
      ctu_rs_addr++;
    }
  }//
}

void kvz_alf_get_blk_stats(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  alf_covariance **alf_covariance,
  alf_classifier **classifier,
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
  static int e_local[MAX_NUM_ALF_LUMA_COEFF][MAX_ALF_NUM_CLIPPING_VALUES];

  const int num_bins = g_alf_num_clipping_values[channel];

  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int transpose_idx = 0;
  int class_idx = 0;

  for (int i = 0; i < height; i++)
  {
    int vb_distance = ((y_dst + i) % vb_ctu_height) - vb_pos;
    for (int j = 0; j < width; j++)
    {
      if (classifier && classifier[y_dst + i][x_dst + j].class_idx == ALF_UNUSED_CLASS_IDX && classifier[y_dst + i][x_dst + j].transpose_idx == ALF_UNUSED_TRANSPOSE_IDX)
      {
        continue;
      }
      memset(e_local, 0, sizeof(e_local));
      if (classifier)
      {
        alf_classifier* cl = &classifier[y_dst + i][x_dst + j];
        transpose_idx = cl->transpose_idx;
        class_idx = cl->class_idx;
      }

      double weight = 1.0;
      if (0/*m_alfWSSD*/)
      {
        weight = g_luma_level_to_weight_plut[org[j]];
      }
      int y_local = org[j] - rec[j];
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
                (*alf_covariance)[class_idx].ee[b0][b1][k][l] += weight * (double)(e_local[k][b0] * e_local[l][b1]);
              }
              else
              {
                (*alf_covariance)[class_idx].ee[b0][b1][k][l] += e_local[k][b0] * e_local[l][b1];
              }
            }
          }
        }
        for (int b = 0; b < num_bins; b++)
        {
          if (0/*m_alfWSSD*/)
          {
            (*alf_covariance)[class_idx].y[b][k] += weight * (double)(e_local[k][b] * y_local);
          }
          else
          {
            (*alf_covariance)[class_idx].y[b][k] += e_local[k][b] * y_local;
          }
        }
      }
      if (0/*m_alfWSSD*/)
      {
        (*alf_covariance)[class_idx].pix_acc += weight * (double)(y_local * y_local);
      }
      else
      {
        (*alf_covariance)[class_idx].pix_acc += y_local * y_local;
      }
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

void kvz_alf_calc_covariance(int e_local[MAX_NUM_ALF_LUMA_COEFF][MAX_ALF_NUM_CLIPPING_VALUES],
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
  const int num_bins = g_alf_num_clipping_values[channel];

  int k = 0;

  const short curr = rec[0];

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

double kvz_alf_get_filter_coeff_and_cost(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  double dist_unfilter,
  int *ui_coeff_bits,
  int i_shape_idx,
  bool b_re_collect_stat,
  bool only_filter_cost)
{
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;

  //collect stat based on CTU decision
  if (b_re_collect_stat)
  {
    get_frame_stats(channel, i_shape_idx);
  }

  double dist = dist_unfilter;
  *ui_coeff_bits = 0;
  int ui_slice_flag = 0;
  
  //AlfFilterShape& alfFilterShape = m_alfSliceParamTemp.filterShapes[channel][iShapeIdx];
  //get filter coeff
  if (is_luma)
  {
    //std::fill_n(m_alfClipMerged[iShapeIdx][0][0], MAX_NUM_ALF_LUMA_COEFF*MAX_NUM_ALF_CLASSES*MAX_NUM_ALF_CLASSES, m_alfParamTemp.non_linear_flag[channel] ? AlfNumClippingValues[CHANNEL_TYPE_LUMA] / 2 : 0);
    const int fill_val = g_alf_aps_temp.non_linear_flag[channel] ? g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] / 2 : 0;
    for (int i = 0; i < MAX_NUM_ALF_CLASSES; i++) {
      for (int j = 0; j < MAX_NUM_ALF_CLASSES; j++) {
        for (int k = 0; k < MAX_NUM_ALF_LUMA_COEFF; k++) {
          g_alf_clip_merged[i_shape_idx][i][j][k] = fill_val;
        }
      }
    }
    // Reset Merge Tmp Cov
    reset_alf_covariance(&g_alf_covariance_merged[i_shape_idx][MAX_NUM_ALF_CLASSES], g_alf_num_clipping_values[channel]);
    reset_alf_covariance(&g_alf_covariance_merged[i_shape_idx][MAX_NUM_ALF_CLASSES + 1], g_alf_num_clipping_values[channel]);
    //distortion
   dist += kvz_alf_merge_filters_and_cost(state, &g_alf_aps_temp, lcu, channel, ui_coeff_bits, g_alf_covariance_frame[channel][i_shape_idx], g_alf_covariance_merged[i_shape_idx], g_alf_clip_merged[i_shape_idx]);
  }
  else
  {
    //distortion
    assert(num_coeff == g_alf_covariance_frame[channel][i_shape_idx][0].num_coeff);
    //std::fill_n(m_filterClippSet[0], MAX_NUM_ALF_CHROMA_COEFF, m_alfParamTemp.non_linear_flag[channel] ? AlfNumClippingValues[CHANNEL_TYPE_CHROMA] / 2 : 0);
    const int fill_val = g_alf_aps_temp.non_linear_flag[channel] ? g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] / 2 : 0;
    for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++) {
      g_filter_clipp_set[0][i] = fill_val;
    }

    dist += g_alf_covariance_frame[channel][i_shape_idx][0].pix_acc + kvz_alf_derive_coeff_quant(channel, g_filter_clipp_set[0], g_filter_coeff_set[0], &g_alf_covariance_frame[channel][i_shape_idx][0], ALF_NUM_BITS, g_alf_aps_temp.non_linear_flag[channel]);

    //setEnableFlag( m_alfSliceParamTemp, channel, m_ctuEnableFlag );
    const int alf_chroma_idc = g_alf_aps_temp.enabled_flag[COMPONENT_Cb] * 2 + g_alf_aps_temp.enabled_flag[COMPONENT_Cr];

    for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
    {
      g_alf_aps_temp.chroma_coeff[i] = g_filter_coeff_set[0][i];
      g_alf_aps_temp.chroma_clipp[i] = g_filter_clipp_set[0][i];
    }
    ui_coeff_bits += get_coeff_rate(&g_alf_aps_temp, true);
    ui_slice_flag = length_truncated_unary(alf_chroma_idc, 3);
  }

  if (only_filter_cost)
  {
    return dist + g_lambda[channel] * *ui_coeff_bits;
  }

  double rate = *ui_coeff_bits + ui_slice_flag;
  //m_CABACEstimator->resetBits();
  //m_CABACEstimator->codeAlfCtuEnableFlags(cs, channel, &m_alfParamTemp);
  if (is_luma)
  {
    if (g_alf_aps_temp.enabled_flag[COMPONENT_Y])
      code_alf_ctu_enable_flags(state, lcu, COMPONENT_Y, &g_alf_aps_temp);
  }
  else
  {
    if (g_alf_aps_temp.enabled_flag[COMPONENT_Cb])
      code_alf_ctu_enable_flags(state, lcu, COMPONENT_Cb, &g_alf_aps_temp);
    if (g_alf_aps_temp.enabled_flag[COMPONENT_Cr])
      code_alf_ctu_enable_flags(state, lcu, COMPONENT_Cr, &g_alf_aps_temp);
  }

  rate += frac_bits_scale * 0;/*(double)m_CABACEstimator->getEstFracBits();*/ // <-- Not supported
  return dist + g_lambda[channel] * rate;
}

//#if JVET_O0669_REMOVE_ALF_COEFF_PRED
int kvz_alf_derive_filter_coefficients_prediction_mode(channel_type channel,
  int **filter_set,
  int** filter_coeff_diff,
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
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;

  static int tmp_clip[MAX_NUM_ALF_LUMA_COEFF];
  static int best_merge_clip[MAX_NUM_ALF_LUMA_COEFF];
  static double err[MAX_NUM_ALF_CLASSES];
  static double best_merge_err;
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
    cov_merged[i].num_coeff = cov[i].num_coeff;
    memcpy(cov_merged[i].ee, cov[i].ee, sizeof(cov_merged[i].ee));
    memcpy(cov_merged[i].y, cov[i].y, sizeof(cov_merged[i].y));
    cov_merged[i].pix_acc = cov[i].pix_acc;
    cov_merged[i].num_bins = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? g_alf_num_clipping_values[COMPONENT_Y] : 1;
  }

  // Try merging different covariance matrices

  // temporal AlfCovariance structure is allocated as the last element in covMerged array, the size of covMerged is MAX_NUM_ALF_CLASSES + 1
  alf_covariance* tmp_cov = &cov_merged[MAX_NUM_ALF_CLASSES];
  tmp_cov->num_bins = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? g_alf_num_clipping_values[COMPONENT_Y] : 1;

  // init Clip
  for (int i = 0; i < num_classes; i++)
  {
    //std::fill_n(clipMerged[numRemaining - 1][i], MAX_NUM_ALF_LUMA_COEFF, m_alfParamTemp.non_linear_flag[CHANNEL_TYPE_LUMA] ? AlfNumClippingValues[CHANNEL_TYPE_LUMA] / 2 : 0);
    for (int val = 0; val < MAX_NUM_ALF_LUMA_COEFF; val++) {
      clip_merged[num_remaining - 1][i][val] = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] / 2 : 0;
    }
    
    if (g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA])
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

            //tmp_cov.add(cov_merged[i], cov_merged[j]);
            for (int b0 = 0; b0 < cov_merged[i].num_bins; b0++)
            {
              for (int b1 = 0; b1 < cov_merged[i].num_bins; b1++)
              {
                for (int j = 0; j < cov_merged[i].num_coeff; j++)
                {
                  for (int i = 0; i < cov_merged[i].num_coeff; i++)
                  {
                    tmp_cov->ee[b0][b1][j][i] = cov_merged[i].ee[b0][b1][j][i] + cov_merged[j].ee[b0][b1][j][i];
                  }
                }
              }
            }
            for (int b = 0; b < cov_merged[i].num_bins; b++)
            {
              for (int j = 0; j < cov_merged[i].num_coeff; j++)
              {
                tmp_cov->y[b][j] = cov_merged[i].y[b][j] + cov_merged[j].y[b][j];
              }
            }
            tmp_cov->pix_acc = cov_merged[i].pix_acc + cov_merged[j].pix_acc;

            for (int l = 0; l < MAX_NUM_ALF_LUMA_COEFF; ++l)
            {
              tmp_clip[l] = (clip_merged[num_remaining - 1][i][l] + clip_merged[num_remaining - 1][j][l] + 1) >> 1;
            }

            double error_merged = g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] ? optimize_filter_clip(tmp_cov, tmp_clip) : calculate_error_opt_filt(tmp_cov, tmp_clip);
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

    //cov_merged[best_to_merge_idx1] += cov_merged[best_to_merge_idx2];
    for (int b0 = 0; b0 < cov_merged[best_to_merge_idx1].num_bins; b0++)
    {
      for (int b1 = 0; b1 < cov_merged[best_to_merge_idx1].num_bins; b1++)
      {
        for (int j = 0; j < cov_merged[best_to_merge_idx1].num_coeff; j++)
        {
          for (int i = 0; i < cov_merged[best_to_merge_idx1].num_coeff; i++)
          {
            cov_merged[best_to_merge_idx1].ee[b0][b1][j][i] += cov_merged[best_to_merge_idx2].ee[b0][b1][j][i];
          }
        }
      }
    }
    for (int b = 0; b < cov_merged[best_to_merge_idx1].num_bins; b++)
    {
      for (int j = 0; j < cov_merged[best_to_merge_idx1].num_coeff; j++)
      {
        cov_merged[best_to_merge_idx1].y[b][j] += cov_merged[best_to_merge_idx2].y[b][j];
      }
    }
    cov_merged[best_to_merge_idx1].pix_acc += cov_merged[best_to_merge_idx2].pix_acc;

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
  alf_aps *alf_aps,
  const lcu_order_element_t *const lcu,
  channel_type channel,
  int *ui_coeff_bits,
  alf_covariance *cov_frame,
  alf_covariance *cov_merged, 
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF])
{
  const int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int num_filters_best = 0;
  int num_filters = MAX_NUM_ALF_CLASSES;
  static bool coded_var_bins[MAX_NUM_ALF_CLASSES];
  static double error_force_0_coeff_tab[MAX_NUM_ALF_CLASSES][2];

  double cost, cost0, dist, dist_force0, cost_min = MAX_DOUBLE;
  int coeff_bits, coeff_bits_force0;

  kvz_alf_merge_classes(channel, cov_frame, cov_merged, clip_merged, MAX_NUM_ALF_CLASSES, g_filter_indices);

  while (num_filters >= 1)
  {
    dist = kvz_alf_derive_filter_coeffs(alf_aps, channel, cov_frame, cov_merged, g_filter_indices[num_filters-1], num_filters, error_force_0_coeff_tab, clip_merged);

    // filter coeffs are stored in m_filterCoeffSet
    dist_force0 = get_dist_force_0(channel, num_filters, error_force_0_coeff_tab, coded_var_bins);
    coeff_bits = kvz_alf_derive_filter_coefficients_prediction_mode(channel, g_filter_coeff_set, g_diff_filter_coeff, num_filters);
    coeff_bits_force0 = get_cost_filter_coeff_force_0(channel, g_filter_coeff_set, num_filters, coded_var_bins);

    cost = dist + g_lambda[COMPONENT_Y] * coeff_bits;
    cost0 = dist_force0 + g_lambda[COMPONENT_Y] * coeff_bits_force0;

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

  dist = kvz_alf_derive_filter_coeffs(alf_aps, channel, cov_frame, cov_merged, g_filter_indices[num_filters_best - 1], num_filters_best, error_force_0_coeff_tab, clip_merged);

  coeff_bits = kvz_alf_derive_filter_coefficients_prediction_mode(channel, g_filter_coeff_set, g_diff_filter_coeff, num_filters_best);
  dist_force0 = get_dist_force_0(channel, num_filters_best, error_force_0_coeff_tab, coded_var_bins);
  coeff_bits_force0 = get_cost_filter_coeff_force_0(channel, g_filter_coeff_set, num_filters_best, coded_var_bins);

  cost = dist + g_lambda[COMPONENT_Y] * coeff_bits;
  cost0 = dist_force0 + g_lambda[COMPONENT_Y] * coeff_bits_force0;

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
      alf_aps->luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_filter_coeff_set[ind][i];
      alf_aps->luma_clipp[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_filter_clipp_set[ind][i];
      /* #if !JVET_O0669_REMOVE_ALF_COEFF_PRED
      if (alf_aps->alf_luma_coeff_delta_prediction_flag)
      {
        alf_aps->luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_diff_filter_coeff[ind][i];
      }
      else
      {
        alf_aps->luma_coeff[ind * MAX_NUM_ALF_LUMA_COEFF + i] = g_filter_coeff_set[ind][i];
      }*/
    }
  }

  memcpy(alf_aps->filter_coeff_delta_idx, g_filter_indices[num_filters_best - 1], sizeof(short) * MAX_NUM_ALF_CLASSES);
  ui_coeff_bits += get_non_filter_coeff_rate(alf_aps);
  return dist_return;
}

double kvz_alf_derive_filter_coeffs(alf_aps *aps,
  channel_type channel,
  alf_covariance *cov,
  alf_covariance *covMerged,
  short* filter_indices,
  int num_filters,
  double error_tab_force_0_coeff[MAX_NUM_ALF_CLASSES][2],
  int clip_merged[MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_CLASSES][MAX_NUM_ALF_LUMA_COEFF])
{
  // #if !JVET_O0669_REMOVE_ALF_COEFF_PRED
  //int *fixed_filter_pattern = &aps->fixed_filter_pattern;
  //int *fixed_filter_idx = aps->fixed_filter_idx;
  //int *fixed_filter_set_index = &aps->fixed_filter_set_index;

  int num_coeff = channel == CHANNEL_TYPE_LUMA ? 13 : 7;
  int *weights = channel == CHANNEL_TYPE_LUMA ? alf_weights_7 : alf_weights_5;

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
    //tmp_cov.reset();
    if (tmp_cov->num_bins > 0)
      tmp_cov->num_bins = -1;
    tmp_cov->pix_acc = 0;
    memset(tmp_cov->y, 0, sizeof(*tmp_cov->y));
    memset(tmp_cov->ee, 0, sizeof(*tmp_cov->ee));

    bool found_clip = false;
    for( int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++ )
    {
      if( filter_indices[class_idx] == filt_idx )
      {
        //tmp_cov += cov[class_idx];
        for (int b0 = 0; b0 < tmp_cov->num_bins; b0++)
        {
          for (int b1 = 0; b1 < tmp_cov->num_bins; b1++)
          {
            for (int j = 0; j < tmp_cov->num_coeff; j++)
            {
              for (int i = 0; i < tmp_cov->num_coeff; i++)
              {
                tmp_cov->ee[b0][b1][j][i] += cov[class_idx].ee[b0][b1][j][i];
              }
            }
          }
        }
        for (int b = 0; b < tmp_cov->num_bins; b++)
        {
          for (int j = 0; j < tmp_cov->num_coeff; j++)
          {
            tmp_cov->y[b][j] += cov[class_idx].y[b][j];
          }
        }
        tmp_cov->pix_acc += cov[class_idx].pix_acc;

         /* #if !JVET_O0669_REMOVE_ALF_COEFF_PRED  //adjust stat
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
    error_tab_force_0_coeff[filt_idx][1] = tmp_cov->pix_acc + kvz_alf_derive_coeff_quant(channel, g_filter_clipp_set[filt_idx], g_filter_coeff_set[filt_idx], tmp_cov, ALF_NUM_BITS, false);
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
  int *weights = is_luma ? alf_weights_7 : alf_weights_5;
  const int factor = 1 << (ALF_NUM_BITS - 1);
  const int max_value = factor - 1;
  const int min_value = -factor + 1;

  static double filter_coeff[MAX_NUM_ALF_LUMA_COEFF];

  optimize_filter(cov, filter_clipp, filter_coeff, optimize_clip);

  //roundFiltCoeff(filterCoeffQuant, filter_coeff, num_coeff, factor);
  for (int i = 0; i < num_coeff; i++)
  {
    int sign = filter_coeff[i] > 0 ? 1 : -1;
    filter_coeff_quant[i] = (int)(filter_coeff[i] * sign * factor + 0.5) * sign;
  }

  for (int i = 0; i < num_coeff; i++)
  {
    filter_coeff_quant[i] += MIN(max_value, MAX(min_value, filter_coeff_quant[i]));
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
          continue;

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

void kvz_alf_encoder_ctb(encoder_state_t *const state,
  const lcu_order_element_t *const lcu,
  alf_aps *aps)
{
  /*
  TempCtx        ctxStart(m_CtxCache, AlfCtx(m_CABACEstimator->getCtx()));
  TempCtx        ctxBest(m_CtxCache);
  TempCtx        ctxTempStart(m_CtxCache);
  TempCtx        ctxTempBest(m_CtxCache);
  */

  //AlfSliceParam  alfSliceParamNewFiltersBest = alfSliceParamNewFilters;
  alf_aps alf_aps_new_filters_best;
  copy_alf_param(&alf_aps_new_filters_best, aps);
  alf_aps* apss = state->slice->apss;

  alf_component_id component_id = state->tile->frame->alf_info->component_id;
  bool is_luma = component_id == COMPONENT_Y ? true : false;

  bool has_new_filters[2] = { aps->enabled_flag[COMPONENT_Y] , aps->enabled_flag[COMPONENT_Cb] || aps->enabled_flag[COMPONENT_Cr] };
  //initDistortion();
  for (int comp = 0; comp < MAX_NUM_COMPONENT; comp++)
  {
    for (int ctb_idx = 0; ctb_idx < g_num_ctus_in_pic; ctb_idx++)
    {
      g_ctb_distortion_unfilter[comp][ctb_idx] = get_unfiltered_distortion_cov_classes(g_alf_covariance[comp][0][ctb_idx], comp == 0 ? MAX_NUM_ALF_CLASSES : 1);
    }
  }

  //luma
  memset(g_ctu_enable_flag[COMPONENT_Y], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
  get_frame_stats(CHANNEL_TYPE_LUMA, 0);
  memset(g_ctu_enable_flag[COMPONENT_Y], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
  double cost_off = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[CHANNEL_TYPE_LUMA][0], CHANNEL_TYPE_LUMA);

  int new_aps_id;
  int aps_ids[ALF_CTB_MAX_NUM_APS];
  int size_of_aps_ids = 0;

  kvz_alf_get_avai_aps_ids_luma(state, &new_aps_id, &aps_ids[ALF_CTB_MAX_NUM_APS], &size_of_aps_ids);

  int best_aps_ids[8] = { -1, -1, -1, -1, -1, -1, -1. -1 };
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
      //memset(state->slice->tile_group_num_aps, num_temporal_aps + use_new_filter, sizeof(int));
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
        double cur_cost = 3 * g_lambda[CHANNEL_TYPE_LUMA];

        if (iter > 0)  //re-derive new filter-set
        {
          double d_dist_org_new_filter = 0;
          int blocks_using_new_filter = 0;
          for (int ctb_idx = 0; ctb_idx < g_num_ctus_in_pic; ctb_idx++)
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
                short* p_clipp = g_clipp_final;
                for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF; i++)
                {
                  g_filter_tmp[i] = p_coeff[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                  g_clip_tmp[i] = p_clipp[class_idx * MAX_NUM_ALF_LUMA_COEFF + i];
                }
                d_dist_org_new_filter += calc_error_for_coeffs(&g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx], g_clip_tmp, g_filter_tmp, MAX_NUM_ALF_LUMA_COEFF, ALF_NUM_BITS);
              }
            }
          }
          if (blocks_using_new_filter > 0 && blocks_using_new_filter < g_num_ctus_in_pic)
          {
            int bit_nl[2] = { 0, 0 };
            double err_nl[2] = { 0.0, 0.0 };
            err_nl[1] = MAX_DOUBLE;
            g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] = 1;
            if (state->encoder_control->cfg.alf_non_linear_luma)
            {
              err_nl[1] = kvz_alf_get_filter_coeff_and_cost(state, lcu, CHANNEL_TYPE_LUMA, 0, &bit_nl[1], 0, true, true);
              copy_alf_param(&g_alf_aps_temp_nl, &g_alf_aps_temp);
            }
            else
            {
              err_nl[1] = MAX_DOUBLE;
            }
            g_alf_aps_temp.non_linear_flag[CHANNEL_TYPE_LUMA] = 0;

            //errNL[0] = getFilterCoeffAndCost(cs, 0, CHANNEL_TYPE_LUMA, true, 0, bitNL[0], true);
            err_nl[0] = kvz_alf_get_filter_coeff_and_cost(state, lcu, CHANNEL_TYPE_LUMA, 0, &bit_nl[0], 0, true, true);

            int bitsNewFilterTempLuma = bit_nl[0];
            int bits_new_filter_temp_luma = bit_nl[0];
            double err = err_nl[0];
            if (err_nl[1]  < err_nl[0])
            {
              err = err_nl[1];
              bits_new_filter_temp_luma = bit_nl[1];
              copy_alf_param(&g_alf_aps_temp, &g_alf_aps_temp_nl);
            }
            if (d_dist_org_new_filter + g_lambda[CHANNEL_TYPE_LUMA] * g_bits_new_filter[CHANNEL_TYPE_LUMA] < err) //re-derived filter is not good, skip
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
        for (int ctb_idx = 0; ctb_idx < g_num_ctus_in_pic; ctb_idx++)
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
            code_alf_ctu_enable_flag(state, lcu, ctb_idx, COMPONENT_Y, &g_alf_aps_temp);
            g_alf_ctb_filter_index[ctb_idx] = filter_set_idx;
            code_alf_ctu_filter_index(state, lcu, ctb_idx, g_alf_aps_temp.enabled_flag[COMPONENT_Y]);

            double rate_on = 4.0255737304687500/*FracBitsScale *(double)m_CABACEstimator->getEstFracBits()*/;
            //distortion
            double dist = dist_unfilter_ctb;
            for (int class_idx = 0; class_idx < MAX_NUM_ALF_CLASSES; class_idx++)
            {
              if (filter_set_idx < ALF_NUM_FIXED_FILTER_SETS)
              {
                int filter_idx = g_class_to_filter_mapping[filter_set_idx][class_idx];
                dist += calc_error_for_coeffs(&g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx], g_clip_default_enc, g_fixed_filter_set_coeff[filter_idx], MAX_NUM_ALF_LUMA_COEFF, ALF_NUM_BITS);
              }
              else
              {
                short *p_coeff;
                short *p_clipp;
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
                dist += calc_error_for_coeffs(&g_alf_covariance[COMPONENT_Y][0][ctb_idx][class_idx], g_clip_tmp, g_filter_tmp, MAX_NUM_ALF_LUMA_COEFF, ALF_NUM_BITS);
              }
            }
            //cost
            double cost_on_tmp = dist + g_lambda[COMPONENT_Y] * rate_on;
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
          code_alf_ctu_enable_flag(state, lcu, ctb_idx, COMPONENT_Y, &g_alf_aps_temp);

          //cost
          double cost_off = dist_unfilter_ctb + (1243132.9728229337 - 1237410.0000000000)/*m_lambda[COMPONENT_Y] * FracBitsScale*(double)m_CABACEstimator->getEstFracBits()*/;
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

        int tmp_bits = bits_new_filter + 3 * (num_filter_set - ALF_NUM_FIXED_FILTER_SETS);
        cur_cost += tmp_bits * g_lambda[COMPONENT_Y];
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
          copy_alf_param(&alf_aps_new_filters_best, &g_alf_aps_temp);

          //ctxBest = AlfCtx(m_CABACEstimator->getCtx());
          //copyCtuEnableFlag(m_ctuEnableFlagTmp, m_ctuEnableFlag, CHANNEL_TYPE_LUMA);
          memcpy(g_ctu_enable_flag_tmp[COMPONENT_Y], g_ctu_enable_flag[COMPONENT_Y], sizeof(uint8_t) * g_num_ctus_in_pic);

          for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++)
          {
            g_alf_ctb_filter_set_index_tmp[ctu_idx] = g_alf_ctb_filter_index[ctu_idx];
          }
          alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA] = use_new_filter;
        }
      }//for (int iter = 0; iter < numIter; iter++)
    }// for (int numTemporalAps = 0; numTemporalAps < apsIds.size(); numTemporalAps++)
  }//for (int useNewFilter = 0; useNewFilter <= 1; useNewFilter++)

  if (cost_off <= cost_min)
  {
    memset(state->slice->tile_group_alf_enabled_flag, 0, sizeof(state->slice->tile_group_alf_enabled_flag));
    state->slice->tile_group_num_aps = 0;
    for (int i = 0; i < MAX_NUM_COMPONENT; i++) {
      memset(g_ctu_enable_flag[i], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
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
    int size_of_best_aps_ids = 0;
    for (int i = 0; i < 8; i++) {
      if (best_aps_ids[i] != -1) {
        size_of_best_aps_ids++;
      }
    }
    state->slice->tile_group_num_aps = size_of_best_aps_ids;
    //cs.slice->setAPSs(bestApsIds);
    for (int i = 0; i < size_of_best_aps_ids; i++)
    {
      state->slice->tile_group_luma_aps_id[i] = best_aps_ids[i];
    }

    //copyCtuEnableFlag(m_ctuEnableFlag, m_ctuEnableFlagTmp, CHANNEL_TYPE_LUMA);
    memcpy(g_ctu_enable_flag[COMPONENT_Y], g_ctu_enable_flag_tmp[COMPONENT_Y], sizeof(uint8_t) * g_num_ctus_in_pic);
    for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++)
    {
      g_alf_ctb_filter_index[ctu_idx] = g_alf_ctb_filter_set_index_tmp[ctu_idx];
    }

    if (alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA])
    {
      //APS* newAPS = m_apsMap->getPS((newApsId << NUM_APS_TYPE_LEN) + T_ALF_APS);
      alf_aps* new_aps = &state->slice->param_set_map[new_aps_id + T_ALF_APS].parameter_set;
      if (new_aps == NULL)
      {
        //newAPS = m_apsMap->allocatePS(newApsId);
        assert(new_aps_id  + T_ALF_APS < MAX_NUM_APS); //Invalid PS id
        bool found = false;
        for (int i = 0; i < (sizeof(state->slice->param_set_map) / sizeof(state->slice->param_set_map[0])); i++) {
          if (state->slice->param_set_map[i].parameter_set.aps_id == new_aps_id  + T_ALF_APS) {
            found = true;
          }
        }
        if (!found) {
          state->slice->param_set_map[new_aps_id + T_ALF_APS].b_changed = true;
          //state->slice->param_set_map[new_aps_id + T_ALF_APS].p_nalu_data = 0;
          //state->slice->param_set_map[new_aps_id + T_ALF_APS].parameter_set = malloc(sizeof(alf_aps));
          state->slice->param_set_map[new_aps_id + T_ALF_APS].parameter_set.aps_id = new_aps_id + T_ALF_APS;
        }
        copy_alf_param(new_aps, &state->slice->param_set_map[new_aps_id + T_ALF_APS].parameter_set);
        new_aps->aps_id = new_aps_id;
        new_aps->aps_type = T_ALF_APS;
      }
      copy_alf_param(new_aps, &alf_aps_new_filters_best);
      new_aps->new_filter_flag[CHANNEL_TYPE_CHROMA] = false;
      state->slice->param_set_map[new_aps_id + T_ALF_APS].b_changed = true;   
      g_aps_id_start = new_aps_id;
    }

    //std::vector<int> apsIds = cs.slice->getTileGroupApsIdLuma();
    int* aps_ids = state->slice->tile_group_luma_aps_id;
    for (int i = 0; i < state->slice->tile_group_num_aps; i++)
    {
      //apss[apsIds[i]] = m_apsMap->getPS((apsIds[i] << NUM_APS_TYPE_LEN) + T_ALF_APS);
      state->slice->apss[aps_ids[i]].aps_id = state->slice->param_set_map[aps_ids[new_aps_id + T_ALF_APS]].parameter_set.aps_id;
      state->slice->apss[aps_ids[i]].aps_type = state->slice->param_set_map[aps_ids[new_aps_id + T_ALF_APS]].parameter_set.aps_type;
      copy_alf_param(&state->slice->apss[aps_ids[i]], &state->slice->param_set_map[aps_ids[new_aps_id + T_ALF_APS]].parameter_set);
    }
  }

  //chroma
  memset(g_ctu_enable_flag[COMPONENT_Cb], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
  memset(g_ctu_enable_flag[COMPONENT_Cr], 1, sizeof(uint8_t) * g_num_ctus_in_pic);
  get_frame_stats(CHANNEL_TYPE_CHROMA, 0);
  cost_off = get_unfiltered_distortion_cov_channel(g_alf_covariance_frame[CHANNEL_TYPE_CHROMA][0], CHANNEL_TYPE_CHROMA) + g_lambda[CHANNEL_TYPE_CHROMA]/*981.62883931057581*/ * 1.0;
  cost_min = MAX_DOUBLE;
  //m_CABACEstimator->getCtx() = AlfCtx(ctxBest);
  //ctxStart = AlfCtx(m_CABACEstimator->getCtx());
  int new_aps_id_chroma = -1;
  if (alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA] && (alf_aps_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_aps_new_filters_best.enabled_flag[COMPONENT_Cr]))
  {
    new_aps_id_chroma = new_aps_id;
  }
  else if (alf_aps_new_filters_best.enabled_flag[COMPONENT_Cb] || alf_aps_new_filters_best.enabled_flag[COMPONENT_Cr])
  {
    int cur_id = g_aps_id_start;
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

  for (int cur_aps_id = 0; cur_aps_id < ALF_CTB_MAX_NUM_APS; cur_aps_id++)
  {
    if (1/*(cs.slice->getPendingRasInit() || cs.slice->isIDRorBLA() || cs.slice->isIntra())*/ && cur_aps_id != new_aps_id_chroma)
    {
      continue;
    }
    //APS* cur_aps = m_apsMap->getPS(cur_aps_id);
    alf_aps* cur_aps = &state->slice->param_set_map[cur_aps_id + T_ALF_APS].parameter_set;
    double cur_cost = g_lambda[CHANNEL_TYPE_CHROMA]/*981.62883931057581*/ * 3;

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
    for (int comp_id = 1; comp_id < MAX_NUM_COMPONENT; comp_id++)
    {
      g_alf_aps_temp.enabled_flag[comp_id] = true;
      for (int ctb_idx = 0; ctb_idx < g_num_ctus_in_pic; ctb_idx++)
      {
        double dist_unfilter_ctu = g_ctb_distortion_unfilter[comp_id][ctb_idx];
        //cost on
        g_ctu_enable_flag[comp_id][ctb_idx] = 1;
        //ctxTempStart = AlfCtx(m_CABACEstimator->getCtx());
        //rate
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
        //m_CABACEstimator->resetBits();
        //ctb flag
        code_alf_ctu_enable_flag(state, lcu, ctb_idx, comp_id, &g_alf_aps_temp);
        double rate_on = frac_bits_scale*(double)838/*m_CABACEstimator->getEstFracBits()*/;
        //distortion
        for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
        {
          g_filter_tmp[i] = g_chroma_coeff_final[i];
          g_clip_tmp[i] = g_chroma_clipp_final[i];
        }
        double dist = dist_unfilter_ctu + calc_error_for_coeffs(&g_alf_covariance[comp_id][0][ctb_idx][0], g_clip_tmp, g_filter_tmp, MAX_NUM_ALF_CHROMA_COEFF, ALF_NUM_BITS);
        double cost_on = dist + g_lambda[comp_id] * rate_on;
        //ctxTempBest = AlfCtx(m_CABACEstimator->getCtx());
        //cost off
        g_ctu_enable_flag[comp_id][ctb_idx] = 0;
        //rate
        //m_CABACEstimator->getCtx() = AlfCtx(ctxTempStart);
        //m_CABACEstimator->resetBits();
        code_alf_ctu_enable_flag(state, lcu, ctb_idx, comp_id, &g_alf_aps_temp);
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
    for (int comp_id = COMPONENT_Cb; comp_id <= COMPONENT_Cr; comp_id++)
    {
      g_alf_aps_temp.enabled_flag[comp_id] = false;
      for (int i = 0; i < g_num_ctus_in_pic; i++)
      {
        if (g_ctu_enable_flag[comp_id][i])
        {
          g_alf_aps_temp.enabled_flag[comp_id] = true;
          break;
        }
      }
    }

    const int alf_chroma_idc = g_alf_aps_temp.enabled_flag[COMPONENT_Cb] * 2 + g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
    cur_cost += length_truncated_unary(alf_chroma_idc, 3) * g_lambda[CHANNEL_TYPE_CHROMA];

    if (cur_cost < cost_min)
    {
      cost_min = cur_cost;
      state->slice->tile_group_chroma_aps_id = cur_aps_id;
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = g_alf_aps_temp.enabled_flag[COMPONENT_Cb];
      state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = g_alf_aps_temp.enabled_flag[COMPONENT_Cr];
      memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cb], g_ctu_enable_flag[COMPONENT_Cb], sizeof(uint8_t) * g_num_ctus_in_pic);
      memcpy(g_ctu_enable_flag_tmp[COMPONENT_Cr], g_ctu_enable_flag[COMPONENT_Cr], sizeof(uint8_t) * g_num_ctus_in_pic);
    }
  }
  if (cost_off < cost_min)
  {
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] = false;
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr] = false;
    memset(g_ctu_enable_flag[COMPONENT_Cb], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
    memset(g_ctu_enable_flag[COMPONENT_Cr], 0, sizeof(uint8_t) * g_num_ctus_in_pic);
  }
  else
  {
    memcpy(g_ctu_enable_flag[COMPONENT_Cb], g_ctu_enable_flag_tmp[COMPONENT_Cb], sizeof(uint8_t) * g_num_ctus_in_pic);
    memcpy(g_ctu_enable_flag[COMPONENT_Cr], g_ctu_enable_flag_tmp[COMPONENT_Cr], sizeof(uint8_t) * g_num_ctus_in_pic);
    if (state->slice->tile_group_chroma_aps_id == new_aps_id_chroma)  //new filter
    {
      //APS* newAPS = m_apsMap->getPS(new_aps_id_chroma);
      alf_aps* new_aps = &state->slice->param_set_map[new_aps_id_chroma + T_ALF_APS].parameter_set;
      if (new_aps == NULL)
      {
        //newAPS = m_apsMap->allocatePS(newApsId);
        assert(new_aps_id + T_ALF_APS < MAX_NUM_APS); //Invalid PS id
        bool found = false;
        for (int i = 0; i < (sizeof(state->slice->param_set_map) / sizeof(state->slice->param_set_map[0])); i++) {
          if (state->slice->param_set_map[i].parameter_set.aps_id == new_aps_id + T_ALF_APS) {
            found = true;
          }
        }
        if (!found) {
          state->slice->param_set_map[new_aps_id + T_ALF_APS].b_changed = true;
          //state->slice->param_set_map[new_aps_id + 1 + T_ALF_APS].p_nalu_data = 0;
          //state->slice->param_set_map[new_aps_id + 1 + T_ALF_APS].parameter_set = malloc(sizeof(alf_aps));
          state->slice->param_set_map[new_aps_id + T_ALF_APS].parameter_set.aps_id = new_aps_id + T_ALF_APS;
        }
        copy_alf_param(new_aps, &state->slice->param_set_map[new_aps_id + T_ALF_APS].parameter_set);
        new_aps->aps_id = new_aps_id;
        new_aps->aps_type = T_ALF_APS;

        //newAPS->getAlfAPSParam().reset();
        memset(new_aps->enabled_flag, false, sizeof(new_aps->enabled_flag));
        memset(new_aps->non_linear_flag, false, sizeof(new_aps->non_linear_flag));
        memset(new_aps->luma_coeff, 0, sizeof(new_aps->luma_coeff));
        memset(new_aps->luma_clipp, 0, sizeof(new_aps->luma_coeff));
        memset(new_aps->chroma_coeff, 0, sizeof(new_aps->chroma_coeff));
        memset(new_aps->chroma_clipp, 0, sizeof(new_aps->chroma_coeff));
        memset(new_aps->filter_coeff_delta_idx, 0, sizeof(new_aps->filter_coeff_delta_idx));
        memset(new_aps->alf_luma_coeff_flag, true, sizeof(new_aps->alf_luma_coeff_flag));
        new_aps->num_luma_filters = 1;
        new_aps->alf_luma_coeff_delta_flag = false;
        //new_aps->alf_luma_coeff_delta_prediction_flag = false;
        new_aps->t_layer = 0;
        memset(new_aps->new_filter_flag, 0, sizeof(new_aps->new_filter_flag));
        //new_aps->fixed_filter_pattern = 0;
        //memset(new_aps->fixed_filter_idx, 0, sizeof(new_aps->fixed_filter_idx));
        //new_aps->fixed_filter_set_index = 0;;
      }
      new_aps->new_filter_flag[CHANNEL_TYPE_CHROMA] = true;
      if (!alf_aps_new_filters_best.new_filter_flag[CHANNEL_TYPE_LUMA]) {
        new_aps->new_filter_flag[CHANNEL_TYPE_LUMA];
      }
      new_aps->non_linear_flag[CHANNEL_TYPE_CHROMA] = aps->non_linear_flag[CHANNEL_TYPE_CHROMA];
      new_aps->t_layer = state->slice->id;
      for (int i = 0; i < MAX_NUM_ALF_CHROMA_COEFF; i++)
      {
        new_aps->chroma_coeff[i] = aps->chroma_coeff[i];
        new_aps->chroma_clipp[i] = aps->chroma_clipp[i];
      }
      state->slice->param_set_map[new_aps_id_chroma + T_ALF_APS].b_changed = true;
      g_aps_id_start = new_aps_id_chroma;
    }
    apss[state->slice->tile_group_chroma_aps_id].aps_id = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id + T_ALF_APS].parameter_set.aps_id;
    apss[state->slice->tile_group_chroma_aps_id].aps_type = state->slice->param_set_map[state->slice->tile_group_chroma_aps_id + T_ALF_APS].parameter_set.aps_type;
    copy_alf_param(&apss[state->slice->tile_group_chroma_aps_id], &state->slice->param_set_map[state->slice->tile_group_chroma_aps_id + T_ALF_APS].parameter_set);
  }
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
  const int max_cu_width = LCU_WIDTH;
  const int max_cu_height = LCU_WIDTH;

  int ctu_idx = 0;
  bool clip_top = false, clip_bottom = false, clip_left = false, clip_right = false;
  int num_hor_vir_bndry = 0, num_ver_vir_bndry = 0;
  int hor_vir_bndry_pos[] = { 0, 0, 0 };
  int ver_vir_bndry_pos[] = { 0, 0, 0 };

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
      if (ctu_enable_flag && is_crossed_by_virtual_boundaries(x_pos, y_pos, width, height, &clip_top, &clip_bottom, &clip_left, &clip_right, &num_hor_vir_bndry, &num_ver_vir_bndry, hor_vir_bndry_pos, ver_vir_bndry_pos, state))
      {
        int y_start = y_pos;
        for (int i = 0; i <= num_hor_vir_bndry; i++)
        {
          const int y_end = i == num_hor_vir_bndry ? y_pos + height : hor_vir_bndry_pos[i];
          const int h = y_end - y_start;
          const bool clip_t = (i == 0 && clip_top) || (i > 0) || (y_start == 0);
          const bool clip_b = (i == num_hor_vir_bndry && clip_bottom) || (i < num_hor_vir_bndry) || (y_end == luma_height);

          int x_start = x_pos;
          for (int j = 0; j <= num_ver_vir_bndry; j++)
          {
            const int x_end = j == num_ver_vir_bndry ? x_pos + width : ver_vir_bndry_pos[j];
            const int w = x_end - x_start;
            const bool clip_l = (j == 0 && clip_left) || (j > 0) || (x_start == 0);
            const bool clip_r = (j == num_ver_vir_bndry && clip_right) || (j < num_ver_vir_bndry) || (x_end == luma_width);

            const int w_buf = w + (clip_l ? 0 : MAX_ALF_PADDING_SIZE) + (clip_r ? 0 : MAX_ALF_PADDING_SIZE);
            const int h_buf = h + (clip_t ? 0 : MAX_ALF_PADDING_SIZE) + (clip_b ? 0 : MAX_ALF_PADDING_SIZE);
            //PelUnitBuf buf = m_tempBuf2.subBuf(UnitArea(cs.area.chromaFormat, Area(0, 0, w_buf, h_buf)));
            //buf.copyFrom(recExtBuf.subBuf(UnitArea(cs.area.chromaFormat, Area(x_start - (clip_l ? 0 : MAX_ALF_PADDING_SIZE), y_start - (clip_t ? 0 : MAX_ALF_PADDING_SIZE), w_buf, h_buf))));
            //buf.extendBorderPel(MAX_ALF_PADDING_SIZE);
            //buf = buf.subBuf(UnitArea(cs.area.chromaFormat, Area(clip_l ? 0 : MAX_ALF_PADDING_SIZE, clip_t ? 0 : MAX_ALF_PADDING_SIZE, w, h)));

            if (g_ctu_enable_flag[COMPONENT_Y][ctu_idx])
            {
              //const Area blkSrc(0, 0, w, h);
              //const Area blkDst(x_start, y_start, w, h);
              short filter_set_index = alf_ctu_filter_index[ctu_idx];
              short *coeff;
              short *clip;
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
              kvz_alf_filter_block(state, lcu, coeff, clip, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y,
                w, h, x_start, y_start, x_start, y_start,
                ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_luma_pos),
                g_alf_vb_luma_ctu_height);
            }

            for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
            {
              alf_component_id comp_id = comp_idx;
              const int chroma_scale_x = ((comp_id == COMPONENT_Y) || (chroma_fmt == KVZ_CSP_444)) ? 0 : 1;
              const int chroma_scale_y = ((comp_id == COMPONENT_Y) || (chroma_fmt != KVZ_CSP_420)) ? 0 : 1;

              if (g_ctu_enable_flag[comp_idx][ctu_idx])
              {
                //const Area blkSrc(0, 0, w >> chromaScaleX, h >> chromaScaleY);
                //const Area blkDst(x_start >> chromaScaleX, y_start >> chromaScaleY, w >> chromaScaleX, h >> chromaScaleY);
                kvz_alf_filter_block(state, lcu, g_chroma_coeff_final, g_chroma_clipp_final, g_clp_rngs.comp[comp_idx], comp_id,
                  w >> chroma_scale_x, h >> chroma_scale_y,
                  x_start >> chroma_scale_x, y_start >> chroma_scale_y,
                  x_start >> chroma_scale_x, y_start >> chroma_scale_y,
                  ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
                  g_alf_vb_chma_ctu_height);
              }
            }

            x_start = x_end;
          }

          y_start = y_end;
        }
      }
      else
      {
        if (g_ctu_enable_flag[COMPONENT_Y][ctu_idx])
        {
          //Area blk(x_pos, y_pos, width, height);
          short filter_set_index = g_alf_ctb_filter_index[ctu_idx];
          short *coeff;
          short *clip;
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
          kvz_alf_filter_block(state, lcu, coeff, clip, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y,
            width, height, x_pos, y_pos, x_pos, y_pos,
            ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_luma_pos),
            g_alf_vb_luma_ctu_height);
        }
        for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
        {
          alf_component_id comp_id = comp_idx;
          const int chroma_scale_x = ((comp_id == COMPONENT_Y) || (chroma_fmt == KVZ_CSP_444)) ? 0 : 1;
          const int chroma_scale_y = ((comp_id == COMPONENT_Y) || (chroma_fmt != KVZ_CSP_420)) ? 0 : 1;

          if (g_ctu_enable_flag[comp_idx][ctu_idx])
          {
            //Area blk(x_pos >> chroma_scale_x, y_pos >> chroma_scale_y, width >> chroma_scale_x, height >> chroma_scale_y);
            //m_filter5x5Blk(m_classifier, recYuv, tmpYuv, blk, comp_id, m_chromaCoeffFinal, clp_rngs.comp[comp_idx], cs);
            kvz_alf_filter_block(state, lcu, g_chroma_coeff_final, g_chroma_clipp_final, g_clp_rngs.comp[comp_idx], comp_idx,
              width >> chroma_scale_x, height >> chroma_scale_y,
              x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
              x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
              ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
              g_alf_vb_chma_ctu_height);
          }
        }
      }
      ctu_idx++;
    }
  }
}

//----------------------------------------------------------------------

//-------------------------cabac writer functions------------------------

void code_alf_ctu_enable_flags(encoder_state_t * const state,
  const lcu_order_element_t *const lcu,
  int8_t component_id,
  alf_aps *aps)
{
  for (int ctu_idx = 0; ctu_idx < g_num_ctus_in_pic; ctu_idx++)
  {
    code_alf_ctu_enable_flag(state, lcu, ctu_idx, component_id, aps);
  }
}

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

    int ctx = 0;
    ctx += left_ctu_addr > -1 ? (ctb_alf_flag[left_ctu_addr] ? 1 : 0) : 0;
    ctx += above_ctu_addr > -1 ? (ctb_alf_flag[above_ctu_addr] ? 1 : 0) : 0;

    //m_BinEncoder.encodeBin(ctbAlfFlag[ctuRsAddr], Ctx::ctbAlfFlag(comp_idx * 3 + ctx));
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
  }
  else
  {
    assert(filter_set_idx < ALF_NUM_FIXED_FILTER_SETS); //Fixed set numavail < num_fixed
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

  kvz_alf_reconstruct_coeff_aps(state,
    true,
    state->slice->tile_group_alf_enabled_flag[COMPONENT_Cb] || state->slice->tile_group_alf_enabled_flag[COMPONENT_Cr],
    false);

  //PelUnitBuf recYuv = cs.getRecoBuf();
  kvz_picture *rec_yuv = state->tile->frame->rec;
  //m_tempBuf.copyFrom(recYuv);
  //PelUnitBuf tmpYuv = m_tempBuf.getBuf(cs.area);
  //tmpYuv.extendBorderPel(MAX_ALF_FILTER_LENGTH >> 1);

  int luma_height = state->tile->frame->height;
  int luma_width = state->tile->frame->width;

  int ctu_idx = 0;
  int max_cu_width = LCU_WIDTH;
  int max_cu_height = LCU_WIDTH;

  bool clip_top = false, clip_bottom = false, clip_left = false, clip_right = false;
  int num_hor_vir_bndry = 0, num_ver_vir_bndry = 0;
  int hor_vir_bndry_pos[] = { 0, 0, 0 };
  int ver_vir_bndry_pos[] = { 0, 0, 0 };

  for (int y_pos = 0; y_pos < luma_height; y_pos += max_cu_width)
  {
    for (int x_pos = 0; x_pos < luma_width; x_pos += max_cu_width)
    {
      const int width = (x_pos + max_cu_width > luma_width) ? (luma_width - x_pos) : max_cu_width;
      const int height = (y_pos + max_cu_width > luma_height) ? (luma_height - y_pos) : max_cu_width;

      bool ctuEnableFlag = g_ctu_enable_flag[COMPONENT_Y][ctu_idx];
      for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
      {
        ctuEnableFlag |= g_ctu_enable_flag[comp_idx][ctu_idx] > 0;
      }

      if(ctuEnableFlag && is_crossed_by_virtual_boundaries(x_pos, y_pos, width, height, &clip_top, &clip_bottom, &clip_left, &clip_right, &num_hor_vir_bndry, &num_ver_vir_bndry, hor_vir_bndry_pos, ver_vir_bndry_pos, state))
      {
        int y_start = y_pos;
        for (int i = 0; i <= num_hor_vir_bndry; i++)
        {
          const int y_end = i == num_hor_vir_bndry ? y_pos + height : hor_vir_bndry_pos[i];
          const int h = y_end - y_start;
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
              kvz_alf_derive_classification(state, lcu, w, h, x_start, y_start, x_start, y_start);
              //const Area blkPCM(x_start, y_start, w, h);
              //resetPCMBlkClassInfo(cs, m_classifier, buf.get(COMPONENT_Y), blkPCM);
              kvz_alf_reset_pcm_blk_class_info(state, lcu, w, h, x_start, y_start);
              short filter_set_index = g_alf_ctb_filter_index[ctu_idx];
              short *coeff;
              short *clip;
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
              kvz_alf_filter_block(state, lcu, coeff, clip, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y, 
                                    w, h, x_start, y_start, x_start, y_start, 
                                  ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_luma_pos),
                                  g_alf_vb_luma_ctu_height);
            }

            for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
            {
              alf_component_id comp_id = comp_idx;
              const int chroma_scale_x = ((comp_id == COMPONENT_Y) || (chroma_fmt == KVZ_CSP_444)) ? 0 : 1;
              const int chroma_scale_y = ((comp_id == COMPONENT_Y) || (chroma_fmt != KVZ_CSP_420)) ? 0 : 1;

              if (g_ctu_enable_flag[comp_idx][ctu_idx])
              {
                //const Area blkSrc(0, 0, w >> chromaScaleX, h >> chromaScaleY);
                //const Area blkDst(x_start >> chromaScaleX, y_start >> chromaScaleY, w >> chromaScaleX, h >> chromaScaleY);
                kvz_alf_filter_block(state, lcu, g_chroma_coeff_final, g_chroma_clipp_final, g_clp_rngs.comp[comp_idx], 
                  comp_id, w >> chroma_scale_x, h >> chroma_scale_y, 
                  x_start >> chroma_scale_x, y_start >> chroma_scale_y, 
                  x_start >> chroma_scale_x, y_start >> chroma_scale_y, 
                  ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
                  g_alf_vb_chma_ctu_height);
              }
            }

            x_start = x_end;
          }

          y_start = y_end;
        }
      }
      else
      {
        //const UnitArea area(cs.area.chromaFormat, Area(x_pos, y_pos, width, height));
        if (g_ctu_enable_flag[COMPONENT_Y][ctu_idx])
        {
          //Area blk(x_pos, y_pos, width, height);
          //deriveClassification(m_classifier, tmpYuv.get(COMPONENT_Y), blk, blk);
          kvz_alf_derive_classification(state, lcu, width, height, x_pos, y_pos, x_pos, y_pos);
          //Area blkPCM(x_pos, y_pos, width, height);
          //resetPCMBlkClassInfo(cs, m_classifier, tmpYuv.get(COMPONENT_Y), blkPCM);
          kvz_alf_reset_pcm_blk_class_info(state, lcu, width, height, x_pos, y_pos);
          short filter_set_index = g_alf_ctb_filter_index[ctu_idx];
          short *coeff;
          short *clip;
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
          kvz_alf_filter_block(state, lcu, coeff, clip, g_clp_rngs.comp[COMPONENT_Y], COMPONENT_Y,
                               width, height, x_pos, y_pos, x_pos, y_pos,
                              ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_luma_pos), 
                               g_alf_vb_luma_ctu_height);
        }
        for (int comp_idx = 1; comp_idx < MAX_NUM_COMPONENT; comp_idx++)
        {
          alf_component_id comp_id = comp_idx;
          const int chroma_scale_x = ((comp_id == COMPONENT_Y) || (chroma_fmt == KVZ_CSP_444)) ? 0 : 1;
          const int chroma_scale_y = ((comp_id == COMPONENT_Y) || (chroma_fmt != KVZ_CSP_420)) ? 0 : 1;

          if (g_ctu_enable_flag[comp_idx][ctu_idx])
          {
            //Area blk(x_pos >> chroma_scale_x, y_pos >> chroma_scale_y, width >> chroma_scale_x, height >> chroma_scale_y);
            //m_filter5x5Blk(m_classifier, recYuv, tmpYuv, blk, comp_id, m_chromaCoeffFinal, clp_rngs.comp[comp_idx], cs);
            kvz_alf_filter_block(state, lcu, g_chroma_coeff_final, g_chroma_clipp_final, g_clp_rngs.comp[comp_idx], comp_idx,
                                 width >> chroma_scale_x, height >> chroma_scale_y,
                                 x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
                                 x_pos >> chroma_scale_x, y_pos >> chroma_scale_y,
                                ((y_pos + max_cu_height >= luma_height) ? luma_height : g_alf_vb_chma_pos),
                                 g_alf_vb_chma_ctu_height);
          }
        }
      }
      ctu_idx++;
    }
  }
}

void kvz_alf_reconstruct_coeff_aps(encoder_state_t *const state, bool luma, bool chroma, bool is_rdo)
{
  //luma
  alf_aps* apss = state->slice->apss;
  //AlfSliceParam alfSliceParamTmp;
  alf_aps* alf_param_tmp;
  //APS* cur_aps;
  alf_aps* cur_aps;

  if (luma)
  {
    for (int i = 0; i < state->slice->tile_group_num_aps /* 1,  cs.slice->getTileGroupNumAps()*/; i++) {
      int aps_idx = state->slice->tile_group_luma_aps_id[i];
      cur_aps = &apss[aps_idx];

      assert(cur_aps != NULL); // "invalid APS"
      alf_param_tmp = cur_aps;

      kvz_alf_reconstruct_coeff(state, alf_param_tmp, CHANNEL_TYPE_LUMA, is_rdo, true);
      memcpy(g_coeff_aps_luma[i], g_coeff_final, sizeof(g_coeff_final));
      memcpy(g_clipp_aps_luma[i], g_clipp_final, sizeof(g_clipp_final));
    }
  }

  //chroma
  if (chroma)
  {
    int aps_idx_chroma = state->slice->tile_group_chroma_aps_id;
    cur_aps = &apss[aps_idx_chroma];
    alf_param_tmp = cur_aps;
    kvz_alf_reconstruct_coeff(state, alf_param_tmp, CHANNEL_TYPE_CHROMA, is_rdo, true);
  }
}

//void reconstructCoeff(AlfSliceParam& alfSliceParam, ChannelType channel, const bool isRdo, const bool isRedo)
void kvz_alf_reconstruct_coeff(encoder_state_t *const state,
  alf_aps *aps,
  channel_type channel,
  const bool is_rdo,
  const bool is_redo)
{
  int num_luma_filters = aps->num_luma_filters;
  //bool *alf_luma_coeff_delta_prediction_flag = &aps->alf_luma_coeff_delta_prediction_flag;
  short *luma_coeff = aps->luma_coeff;
  short *luma_clipp = aps->luma_clipp;
  short *chroma_coeff = aps->chroma_coeff;
  short *chroma_clipp = aps->chroma_clipp;
  short *filter_coeff_delta_idx = aps->filter_coeff_delta_idx;
  //int *fixed_filter_set_index = &aps->fixed_filter_set_index;

  int factor = is_rdo ? 0 : (1 << (ALF_NUM_BITS - 1));
  bool is_luma = channel == CHANNEL_TYPE_LUMA ? 1 : 0;

  alf_filter_type filter_type = is_luma ? ALF_FILTER_7X7 : ALF_FILTER_5X5;
  int num_classes = is_luma ? MAX_NUM_ALF_CLASSES : 1;
  int num_coeff = filter_type == ALF_FILTER_5X5 ? 7 : 13;
  int num_coeff_minus1 = num_coeff - 1;
  int num_filters = is_luma ? num_luma_filters : 1;
  short* coeff = is_luma ? luma_coeff : chroma_coeff;
  short* clipp = is_luma ? luma_clipp : chroma_clipp;

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
  }*/

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
    }*/
    for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
    {
      g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] = coeff[filter_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx];
      /*#if !JVET_O0669_REMOVE_ALF_COEFF_PRED
      //fixed filter
      if (fixed_filter_idx >= 0)
      {
        g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + coeff_idx] += g_fixed_filter_set_coeff[fixed_filter_idx][coeff_idx];
      }*/
    }
    g_coeff_final[class_idx * MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = factor;
    g_clipp_final[class_idx* MAX_NUM_ALF_LUMA_COEFF + num_coeff_minus1] = is_rdo ? 0 : g_alf_clipping_values[channel][0];
    for (int coeff_idx = 0; coeff_idx < num_coeff_minus1; ++coeff_idx)
    {
      int clip_idx = aps->non_linear_flag[channel] ? (clipp + filter_idx * MAX_NUM_ALF_LUMA_COEFF)[coeff_idx] : 0;
      (g_clipp_final + class_idx * MAX_NUM_ALF_LUMA_COEFF)[coeff_idx] = is_rdo ? clip_idx : g_alf_clipping_values[channel][clip_idx];
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

void kvz_alf_create(encoder_state_t const *state,
  const lcu_order_element_t const *lcu) //,
  //const int format)
{
  //memcpy(m_inputBitDepth, inputBitDepth, sizeof(m_inputBitDepth));

  int ***laplacian = state->tile->frame->alf_info->laplacian;
  alf_classifier ***classifier = &state->tile->frame->alf_info->classifier;
  bool *created = &state->tile->frame->alf_info->created;
  enum kvz_chroma_format chroma_fmt = state->encoder_control->chroma_format;

  const int pic_width = state->tile->frame->width;
  const int pic_height = state->tile->frame->height;
  const int max_cu_width = LCU_WIDTH; //128
  const int max_cu_height = LCU_WIDTH; //128

  g_num_ctus_in_width = (pic_width / max_cu_width) + ((pic_width % max_cu_width) ? 1 : 0);
  g_num_ctus_in_height = (pic_height / max_cu_height) + ((pic_height % max_cu_height) ? 1 : 0);
  g_num_ctus_in_pic = g_num_ctus_in_width * g_num_ctus_in_height;

  g_alf_vb_luma_pos = max_cu_height - ALF_VB_POS_ABOVE_CTUROW_LUMA;
  g_alf_vb_chma_pos = (max_cu_height >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0)) - ALF_VB_POS_ABOVE_CTUROW_CHMA;

  g_alf_vb_luma_ctu_height = max_cu_height;
  g_alf_vb_chma_ctu_height = (max_cu_height >> ((chroma_fmt == KVZ_CSP_420) ? 1 : 0));

  assert(g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] > 0); //"g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] must be at least one"
  for (int i = 0; i < g_alf_num_clipping_values[CHANNEL_TYPE_LUMA]; ++i)
  {
    g_alf_clipping_values[CHANNEL_TYPE_LUMA][i] =
      (kvz_pixel)round(pow(2., g_input_bit_depth[CHANNEL_TYPE_LUMA] * 
      (g_alf_num_clipping_values[CHANNEL_TYPE_LUMA] - i) / g_alf_num_clipping_values[CHANNEL_TYPE_LUMA]));
  }

  assert(g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] > 0); //"g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] must be at least one"
  g_alf_clipping_values[CHANNEL_TYPE_CHROMA][0] = 1 << g_input_bit_depth[CHANNEL_TYPE_CHROMA];
  for (int i = 1; i < g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA]; ++i)
  {
    g_alf_clipping_values[CHANNEL_TYPE_CHROMA][i] =
      (kvz_pixel)round(pow(2., g_input_bit_depth[CHANNEL_TYPE_CHROMA] - 8 + 8. * 
      (g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] - i - 1) / (g_alf_num_clipping_values[CHANNEL_TYPE_CHROMA] - 1)));
  }

  if (*created)
  {
    return;
  }

  //m_tempBuf.destroy();
  //m_tempBuf.create(format, Area(0, 0, pic_width, pic_height), max_cu_width, MAX_ALF_FILTER_LENGTH >> 1, 0, false);
  //m_tempBuf2.destroy();
  //m_tempBuf2.create(format, Area(0, 0, maxCUWidth + (MAX_ALF_PADDING_SIZE << 1), maxCUHeight + (MAX_ALF_PADDING_SIZE << 1)), maxCUWidth, MAX_ALF_PADDING_SIZE, 0, false);

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
      for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF - 1; i++)
      {
        g_fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + i] = g_fixed_filter_set_coeff[fixed_filter_idx][i];
      }
      g_fixed_filter_set_coeff_dec[filter_set_index][class_idx * MAX_NUM_ALF_LUMA_COEFF + MAX_NUM_ALF_LUMA_COEFF - 1] = (1 << (ALF_NUM_BITS - 1));
    }
  }

  for (int i = 0; i < MAX_NUM_ALF_LUMA_COEFF * MAX_NUM_ALF_CLASSES; i++)
  {
    g_clip_default[i] = g_alf_clipping_values[CHANNEL_TYPE_LUMA][0];
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

  /*
  m_tempBuf.destroy();
  m_tempBuf2.destroy();
  */

  alf->created = false;
}

void kvz_alf_derive_classification(encoder_state_t *const state,
  const lcu_order_element_t * const lcu,
  const int width,
  const int height,
  int x_pos,
  int y_pos,
  const int blk_dst_x,
  const int blk_dst_y)//,
  //alf_classifier** classifier)
{
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;
  int32_t pic_height = state->tile->frame->rec->height;

  int max_height = y_pos + height;
  int max_width = x_pos + width;

  for (int i = y_pos; i < max_height; i += CLASSIFICATION_BLK_SIZE)
  {
    int n_height = MIN(i + CLASSIFICATION_BLK_SIZE, max_height) - i;

    for (int j = x_pos; j < max_width; j += CLASSIFICATION_BLK_SIZE)
    {
      int n_width = MIN(j + CLASSIFICATION_BLK_SIZE, max_width) - j;
      //kvz_alf_derive_classification_blk(state, lcu, classifier, alf_input_bit_depth[CHANNEL_TYPE_LUMA] + 4, n_height, n_width, j, i);
      kvz_alf_derive_classification_blk(state, lcu, classifier, g_input_bit_depth[CHANNEL_TYPE_LUMA] + 4, n_height, n_width, j, i, 
        j - x_pos + blk_dst_x, i - y_pos + blk_dst_y, g_alf_vb_luma_ctu_height, 
        ((i - y_pos + blk_dst_y + n_height >= pic_height) ? pic_height : g_alf_vb_luma_pos));
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
  const int blk_pos_x,
  const int blk_pos_y,
  const int blk_dst_x,
  const int blk_dst_y,
  const int vb_ctu_height,
  int vb_pos)
{
  videoframe_t* const frame = state->tile->frame;
  int ***laplacian = state->tile->frame->alf_info->laplacian;

  //CHECK((vb_ctu_height & (vb_ctu_height - 1)) != 0, "vb_ctu_height must be a power of 2");

  static const int th[16] = { 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
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
  const short* filter_set,
  const short *fClipSet,
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
  videoframe_t* const frame = state->tile->frame;
  alf_filter_type const filter_type = component_id == COMPONENT_Y ? ALF_FILTER_7X7 : ALF_FILTER_5X5;
  const bool chroma = component_id == COMPONENT_Y ? 0 : 1;
  alf_classifier **classifier = state->tile->frame->alf_info->classifier;

  //CHECK((vb_ctu_height & (vb_ctu_height - 1)) != 0, "vb_ctu_height must be a power of 2");

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
  kvz_pixel *dst = state->tile->frame->rec->y;// +start_height * dst_stride;

  const kvz_pixel *p_img_y_pad_0, *p_img_y_pad_1, *p_img_y_pad_2, *p_img_y_pad_3, *p_img_y_pad_4, *p_img_y_pad_5, *p_img_y_pad_6;
  const kvz_pixel *p_img_0, *p_img_1, *p_img_2, *p_img_3, *p_img_4, *p_img_5, *p_img_6;

  const short *coef = filter_set;
  const short *clip = fClipSet;

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
  int filter_coeff[MAX_NUM_ALF_LUMA_COEFF];
  //std::array<int, MAX_NUM_ALF_LUMA_COEFF> filterClipp;
  int filter_clipp[MAX_NUM_ALF_LUMA_COEFF];

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
      p_class = classifier[blk_dst_y + i] + blk_dst_x;
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
        clip = fClipSet + cl.class_idx * MAX_NUM_ALF_LUMA_COEFF;
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
            /*
            Position pos(j + blkDst.x + blkX, i + blkDst.y + blkY);
#if JVET_O0090_ALF_CHROMA_FILTER_ALTERNATIVES_CTB && !JVET_O0050_LOCAL_DUAL_TREE
            const CodingUnit* cu = isDualTree ? cs.getCU(pos, CH_C) : cs.getCU(recalcPosition(nChromaFormat, CH_C, CH_L, pos), CH_L);
#else
            CodingUnit* cu = isDualTree ? cs.getCU(pos, CH_C) : cs.getCU(recalcPosition(nChromaFormat, CH_C, CH_L, pos), CH_L);
#endif
#if JVET_O0050_LOCAL_DUAL_TREE
            cu = cu->isSepTree() ? cs.getCU( pos, CH_C ) : cu;
#endif
            *flags++ = cu->ipcm ? 1 : 0;
            */
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

        const int yVb = (blk_dst_y + i + ii) & (vb_ctu_height - 1);
        if (yVb < vb_pos && (yVb >= vb_pos - (chroma ? 2 : 4)))   // above
        {
          p_img_1 = (yVb == vb_pos - 1) ? p_img_0 : p_img_1;
          p_img_3 = (yVb >= vb_pos - 2) ? p_img_1 : p_img_3;
          p_img_5 = (yVb >= vb_pos - 3) ? p_img_3 : p_img_5;

          p_img_2 = (yVb == vb_pos - 1) ? p_img_0 : p_img_2;
          p_img_4 = (yVb >= vb_pos - 2) ? p_img_2 : p_img_4;
          p_img_6 = (yVb >= vb_pos - 3) ? p_img_4 : p_img_6;
        }
        else if (yVb >= vb_pos && (yVb <= vb_pos + (chroma ? 1 : 3)))   // bottom
        {
          p_img_2 = (yVb == vb_pos) ? p_img_0 : p_img_2;
          p_img_4 = (yVb <= vb_pos + 1) ? p_img_2 : p_img_4;
          p_img_6 = (yVb <= vb_pos + 2) ? p_img_4 : p_img_6;

          p_img_1 = (yVb == vb_pos) ? p_img_0 : p_img_1;
          p_img_3 = (yVb <= vb_pos + 1) ? p_img_1 : p_img_3;
          p_img_5 = (yVb <= vb_pos + 2) ? p_img_3 : p_img_5;
        }


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