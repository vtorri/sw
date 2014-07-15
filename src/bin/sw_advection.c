/*
 * SW - a Spline wavelet library
 * Copyright (c) 2013-2014, Vincent Torri
 * All rights reserved
 *
 * This software is governed by the CeCILL-C license under French law
 * and abiding by the rules of distribution of free software. You can
 * use, modify and/or redistribute the software under the terms of the
 * CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
 * URL: "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided
 * only with a limited warranty and the software's author, the holder of
 * the economic rights, and the successive licensors have only limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading, using, modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean that it is complicated to manipulate, and that also
 * therefore means that it is reserved for developers and experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards
 * their requirements in conditions enabling the security of their
 * systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <math.h>

#include <sw.h>

static double
mod(double x, double per)
{
  if (x >= per)
    {
      while (x >= per)
        x -= per;
    }
  else
    {
      if (x < 0.0)
        {
          while (x < 0.0)
            x += per;
        }
    }

  return x;
}

/* x translate par multiple de per pour etre >= inf et < inf + per */
static double
mod2(double inf, double per, double x)
{
  if (x < inf)
    {
      while (x < inf)
	x += per;
    }

  if (x >= inf + per)
    {
      while (x >= inf + per)
	x -= per;
    }

  return x;
}

static double
gauss(double lambda, double x)
{
  return exp(-lambda * x * x);
}

static void
gaussian(double *f, double lambda, double x0, int32_t scale)
{
  int32_t size;
  int32_t i;

  size = 1 << scale;

  for (i = 0; i < size; i++)
    {
      double x;

      x = (double)i / (double)size;
      f[i] = gauss(lambda, x - x0);
    }
}

void gaussian_period(double *f, double lambda, double x0, int scale)
{
  int size;
  int i;

  size = 1 << scale;

  for (i = 0; i < size; i++)
    {
      double x1;
      double x2;

      x1 = (double)i / (double)size;
      x2 = mod2(-0.5 + x0, 1, x1);
      f[i] = gauss(lambda, x2 - x0);
    }
}

static void
file_save(const double *t, int size, const char *filename)
{
  FILE *f;
  double *iter;
  int i;

  printf("size 2 %s : %d\n", filename, size);

  f = fopen(filename, "wb");
  if (!f) return;

  iter = (double *)t;
  for (i = 0; i < size; i++, iter++)
    {
      fprintf (f, "%f %f\n", (double)i / size, *iter);
      printf ("%f %f\n", (double)i / size, *iter);
    }
  fprintf (f, "\n");
  fclose (f);
}

static void
_advection(const sw_mra_t *mra, double delta_t, double *ps0, double *ps1)
{
  const sw_scale_fct_t *sf;
  const sw_scale_fct_dual_t *sfd;
  const double *weights;
  int32_t degree;
  int32_t k;
  int32_t l;
  int32_t n;
  int32_t i_k;
  int32_t i_l;

  double c = 1.0;

  sf = sw_mra_scale_fct_get(mra);
  sfd = sw_mra_scale_fct_dual_get(mra);
  weights = sw_scale_fct_dual_weights_get(sfd, &degree);

  printf("sf(0.0) = %f\n", sw_scale_fct_value_get(sf, 0.0));
  printf("sf(0.5) = %f\n", sw_scale_fct_value_get(sf, 0.5));
  printf("sf(1.0) = %f\n", sw_scale_fct_value_get(sf, 1.0));

  printf("weights :\n");
  for (k = 0; k < degree; k++)
  {
      printf("[%d] %E\n", k, weights[k]);
  }

  degree--;

  printf(" supp : %d %d\n", sw_mra_size_inf_x_get(mra), sw_mra_size_sup_x_get(mra));
  i_k = 0;
  for (k = sw_mra_size_inf_x_get(mra); k <= sw_mra_size_sup_x_get(mra); k++)
    {
      ps1[i_k] = 0.0;
      i_l = 0;
      for (l = sw_mra_size_inf_x_get(mra); l <= sw_mra_size_sup_x_get(mra); l++)
        {
          double v = 0.0;
          for (n = -(degree >> 1); n <= degree - (degree >> 1); n++)
            {
                printf("[%d] [%d] [%d] : %E   %E\n", k, l, n,
                       k + n - l - c * delta_t * sw_mra_size_x_get(mra),
                       sw_scale_fct_value_get(sf, k + n - l - c * delta_t * sw_mra_size_x_get(mra)));
                v += weights[n + (degree >> 1)] * sw_scale_fct_value_get(sf, k + n - l - c * delta_t * sw_mra_size_x_get(mra));
                //v += weights[n + (degree >> 1)] * sw_scale_fct_value_get(sf, mod(k + n - l - c * delta_t * sw_mra_size_x_get(mra), sw_mra_size_x_get(mra)));
            }
          ps1[i_k] += v * ps0[i_l];
          i_l++;
        }
      i_k++;
    }
}

static void
sw_advection(int32_t order, int32_t order_dual, int32_t scale, int32_t degree)
{
  sw_mra_t *mra;
  double *f0;
  double *f1;
  double *ps0;
  double *ps1;
  int32_t size;

  int i;

  double lambda;
  double delta_t;

  lambda = 100;

  size = 1 << scale;
  mra = sw_mra_new(order, order_dual, 2, scale, SW_WEIGHTS_TYPE_LAGRANGE, degree);
  printf("mra : %p\n", mra);

  printf("size 1 : %d\n", size);

  f0 = sw_new(size);
  f1 = sw_new(size);
  ps0 = sw_new(size);
  ps1 = sw_new(size);

  gaussian(f0, lambda, 0.5, scale);
  file_save(f0, size, "sw_advection_0.dat");

  gaussian_period(f1, lambda, 0.8, scale);
  file_save(f1, size, "sw_advection_per_0_8.dat");

  delta_t = 0.5;

  sw_mra_proj_x_forward(mra, f0, ps0);
  file_save(ps0, size, "sw_ps_0.dat");

  printf("avant advection\n");
  _advection(mra, delta_t, ps0, ps1);
  file_save(ps1, size, "sw_ps_0_5.dat");
  printf("apres advection\n");

  sw_mra_proj_x_backward(mra, ps1, f1);
  file_save(f1, size, "sw_advection_0_5.dat");

  sw_free(ps1);
  sw_free(ps0);
  sw_free(f1);
  sw_free(f0);

  sw_mra_del(mra);
}

int main()
{
    sw_advection(1, 1, 3, 9);
    return 0;
}
