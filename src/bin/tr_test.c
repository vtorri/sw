#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <math.h>

#include <sw.h>
#include <tr.h>

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
double mod2(double inf, double per, double x)
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

double gauss(double lambda, double x)
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

  f = fopen(filename, "wb");
  if (!f) return;

  iter = (double *)t;
  for (i = 0; i < size; i++, iter++)
    {
      fprintf (f, "%f %f\n", (double)i / size, *iter);
    }
  fprintf (f, "\n");
  fclose (f);
}

static void
test_rk()
{
  double y;
  int i;

  tr_field_data_set(0.0, 1.0, 1.0, 0.0, 1);

  y = 1.0;
  for (i = 0; i < 100; i++)
  y = tr_rk4(tr_field, y, 0.01);

  printf("%f %f %E\n", y, exp(-1.0), y - exp(-1.0));
}

static void
advection(const mra_t *mra, double delta_t, double *ps0, double *ps1)
{
  const scale_fct_t *sf;
  const scale_fct_dual_t *sfd;
  const double *weights;
  int32_t degree;
  int32_t k;
  int32_t l;
  int32_t n;
  int32_t i_k;
  int32_t i_l;

  sf = mra_scale_fct_get(mra);
  sfd = mra_scale_fct_dual_get(mra);
  weights = scale_fct_dual_weights_get(sfd, &degree);
  degree--;

  i_k = 0;
  for (k = mra_size_inf_x_get (mra); k <= mra_size_sup_x_get (mra); k++)
    {
      ps1[i_k] = 0.0;
      i_l = 0;
      for (l = mra_size_inf_x_get (mra); l <= mra_size_sup_x_get (mra); l++)
        {
          double v = 0.0;
          for (n = -(degree >> 1); n <= degree - (degree >> 1); n++)
            {
              //v += weights[n + (degree >> 1)] * scale_fct_value_get(sf, mod(k + n - l - c * delta_t * mra_size_x_get(mra), mra_size_x_get (mra)));
              //v += weights[n + (degree >> 1)] * scale_fct_value_get(sf, k + n - l - c * delta_t * mra_size_x_get(mra));
              double x;

	      //x = (k + n - c * delta_t * mra_size_x_get(mra)) / (double)mra_size_x_get(mra);

	      x = tr_rk4(tr_field, (double)(k+n) / (double)mra_size_x_get(mra), delta_t);

              v += weights[n + (degree >> 1)] * scale_fct_value_get(sf, mra_size_x_get(mra) * mod(x, 1) - l);
              //v += weights[n + (degree >> 1)] * scale_fct_value_get(sf, mra_size_x_get(mra) * x - l);
            }
          ps1[i_k] += v * ps0[i_l];
          i_l++;
        }
      i_k++;
    }
}

static void
test_advection(int32_t order, int32_t order_dual, int32_t scale, int32_t degree)
{
  mra_t *mra;
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
  mra = mra_new (order, order_dual, 3, scale, SW_WEIGHTS_TYPE_LAGRANGE, degree);

  f0 = sw_new(size);
  f1 = sw_new(size);
  ps0 = sw_new(size);
  ps1 = sw_new(size);

  gaussian (f0, lambda, 0.5, scale);
  mra_proj_x_forward(mra, f0, ps0);

  /* On fait l'advection ici */
  delta_t = 0.1;

  advection(mra, delta_t, ps0, ps1);
  ps0 = ps1;
  advection(mra, delta_t, ps0, ps1);
  ps0 = ps1;
  /* advection(mra, delta_t, ps0, ps1); */
  /* ps0 = ps1; */
  /* advection(mra, delta_t, ps0, ps1); */
  /* ps0 = ps1; */
  /* advection(mra, delta_t, ps0, ps1); */
  /* for (i = 0; i < 1; i++) */
  /*   { */
  /*     printf(" * %d\n", i); */
  /*     advection(mra, delta_t, ps0, ps1); */
  /*     ps0 = ps1; */
  /*   } */

  mra_proj_x_backward(mra, ps1, f1);

  file_save(f1, size, "f1.dat");

  /* centree en 0.9 */
  gaussian_period(f0, lambda, 0.55275, scale);
  file_save(f0, size, "f2.dat");

  printf ("adv err : %e\n", sw_error(f0, f1, size));

  sw_free(ps1);
  sw_free(ps0);
  sw_free(f1);
  sw_free(f0);

  mra_del(mra);
}

int main()
{
  test_rk();
  test_advection(6, 2, 10, 9);

  return 0;
}
