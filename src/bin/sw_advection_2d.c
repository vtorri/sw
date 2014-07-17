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
gaussian_2d(double *f, double lambda_x, double lambda_y, double x0, double y0, int32_t scale)
{
    int32_t size;
    int32_t i;
    int32_t j;

    size = 1 << scale;

    for (i = 0; i < size; i++)
    {
        double x;

        x = (double)i / (double)size;

        for (j = 0; j < size; j++)
        {
            double y;

            y = (double)j / (double)size;

            f[i * size + j] = gauss(lambda_x, x - x0) * gauss(lambda_y, y - y0);
        }
    }
}

static void
gaussian_2d_period(double *f, double lambda_x, double lambda_y, double x0, double y0, int32_t scale)
{
    int32_t size;
    int32_t i;
    int32_t j;

    size = 1 << scale;

    for (i = 0; i < size; i++)
    {
        double x1;
        double x2;

        x1 = (double)i / (double)size;
        x2 = mod2(-0.5 + x0, 1, x1);

        for (j = 0; j < size; j++)
        {
            double y1;
            double y2;

            y1 = (double)j / (double)size;
            y2 = mod2(-0.5 + y0, 1, y1);

            f[i * size + j] = gauss(lambda_x, x2 - x0) * gauss(lambda_y, y2 - y0);
        }
    }
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

  degree--;

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
file_save(const double *t, int size, const char *filename)
{
    FILE *f;
    //double *iter;
    int i;
    int j;

    f = fopen(filename, "wb");
    if (!f) return;

    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            fprintf (f, "%f \n", t[i * size + j]);
        }
        fprintf (f, "\n");
    }
  fclose (f);
}

static void
sw_advection_2d(int32_t order, int32_t order_dual, int32_t scale, int32_t degree)
{
    sw_mra_t *mra;
    double *f0;
    double *f1;
    double *f2;
    double *ps0;
    double *ps1;
    double *tmp;
    double x0;
    double y0;
    double lambda_x;
    double lambda_y;
    double delta_t;
    int32_t size;
    int i;
    int j;

    x0 = 0.5;
    y0 = 0.5;
    lambda_x = 90;
    lambda_y = 90;
    delta_t = 0.1;

    size = 1 << scale;

    f0 = sw_new(size * size);
    f1 = sw_new(size * size);
    f2 = sw_new(size * size);
    ps0 = sw_new(size);
    ps1 = sw_new(size);
    tmp = sw_new(size);

    gaussian_2d(f0, lambda_x, lambda_y, x0, y0, scale);

    file_save(f0, size, "sw_advection_2d_0.dat");

    mra = sw_mra_lagrange_new(order, order_dual, 3, scale, degree);

    /* advection en x */
    for (j = 0; j < size; j++)
    {
        for (i = 0; i < size; i++)
        {
            tmp[i] = f0[i * size + j];
        }
        sw_mra_proj_x_forward(mra, tmp, ps0);
        _advection(mra, delta_t, ps0, ps1);
        sw_mra_proj_x_backward(mra, ps1, tmp);
        for (i = 0; i < size; i++)
        {
            f1[i * size + j] = tmp[i];
        }
    }

    /* advection en y */
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            tmp[j] = f1[i * size + j];
        }
        sw_mra_proj_x_forward(mra, tmp, ps0);
        _advection(mra, delta_t, ps0, ps1);
        sw_mra_proj_x_backward(mra, ps1, tmp);
        for (j = 0; j < size; j++)
        {
            f2[i * size + j] = tmp[j];
        }
    }

    file_save(f2, size, "sw_advection_2d_0_1.dat");

    sw_free(tmp);
    sw_free(ps1);
    sw_free(ps0);
    sw_free(f2);
    sw_free(f1);
    sw_free(f0);

    sw_mra_del(mra);
}

int main()
{
    sw_advection_2d(2, 2, 8, 9);
    return 0;
}
