#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "sw.h"


/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

/**
 * @cond SW_LOCAL
 */

/**
 * @endcond SW_LOCAL
 */


/******************************************************************************/
/*                                                                            */
/*                                   GLOBAL                                   */
/*                                                                            */
/******************************************************************************/


/******************************************************************************/
/*                                                                            */
/*                                    API                                     */
/*                                                                            */
/******************************************************************************/

mra_t *
mra_new (int32_t order,
         int32_t order_dual,
         int32_t scale_coarse,
         int32_t scale_fine,
         weights_type_t type,
         ...)
{
    mra_t  *mra;
    va_list va;

    if (scale_coarse >= scale_fine)
        return NULL;

    mra = (mra_t *)malloc (sizeof (mra_t));
    if (!mra)
        return NULL;

    mra->scale_fct = scale_fct_new (order);
    if (!mra->scale_fct)
        goto no_scale_fct;

    scale_fct_type_set(mra->scale_fct, type);

    mra->scale_fct_dual = scale_fct_dual_new (order, order_dual);
    if (!mra->scale_fct_dual)
        goto no_scale_fct_dual;

    va_start (va, type);

    switch (type)
    {
     case SW_WEIGHTS_TYPE_LAGRANGE:
     {
         int32_t degree;

         degree = va_arg (va, int32_t);
         va_end (va);
         scale_fct_dual_type_set(mra->scale_fct_dual, type, degree);
         break;
     }
     case SW_WEIGHTS_TYPE_SWELDENS:
     {
         int32_t r;
         int32_t s;

         r = va_arg(va, int32_t);
         s = va_arg(va, int32_t);
         va_end (va);
         scale_fct_dual_type_set(mra->scale_fct_dual, type, r, s);
         break;
     }
     default:
         va_end (va);
         goto no_scale_fct_dual;
    }

    mra->wavelet = wavelet_new (mra->scale_fct, mra->scale_fct_dual);
    if (!mra->wavelet)
        goto no_wavelet;

    mra->wavelet_dual = wavelet_dual_new (mra->scale_fct, mra->scale_fct_dual);
    if (!mra->wavelet_dual)
        goto no_wavelet_dual;

    mra->scale_coarse = scale_coarse;
    mra->scale_fine = scale_fine;
    mra->lambda_fine_x_inf = 0;
    mra->lambda_fine_x_sup = (1 << scale_fine) - 1;
    mra->size_x = 1 << scale_fine;
    mra->lambda_fine_v_inf = -(1 << scale_fine);
    mra->lambda_fine_v_sup = (1 << scale_fine) - 1;
    mra->size_v = 1 << (scale_fine + 1);

    return mra;

  no_wavelet_dual:
    wavelet_del (mra->wavelet);
  no_wavelet:
    scale_fct_dual_del (mra->scale_fct_dual);
  no_scale_fct_dual:
    scale_fct_del (mra->scale_fct);
  no_scale_fct:
    free (mra);

    return NULL;
}

void
mra_del (mra_t *mra)
{
    if (!mra)
        return;

    wavelet_dual_del (mra->wavelet_dual);
    wavelet_del (mra->wavelet);
    printf("del mra->scale_fct_dual\n");
    scale_fct_dual_del (mra->scale_fct_dual);
    scale_fct_del (mra->scale_fct);
    free (mra);
}

const scale_fct_t *
mra_scale_fct_get (const mra_t *mra)
{
    if (!mra)
        return NULL;

    return mra->scale_fct;
}

const scale_fct_dual_t *
mra_scale_fct_dual_get (const mra_t *mra)
{
    if (!mra)
        return NULL;

    return mra->scale_fct_dual;
}

/* Projections 1D */

void
mra_proj_x_forward (const mra_t  *mra,
                    const double *function,
                    double       *ps)
{
    if (!mra)
        return;

    scale_fct_dual_proj_periodic_forward (mra->scale_fct_dual,
                                          mra->scale_fine,
                                          function,
                                          ps);
}

void
mra_proj_x_backward (const mra_t  *mra,
                     const double *ps,
                     double       *function)
{
    if (!mra)
        return;

    scale_fct_proj_periodic_backward(mra->scale_fct,
                                     mra->scale_fine,
                                     ps,
                                     function);
}

void
mra_proj_v_forward (const mra_t  *mra,
                    const double *function,
                    double       *ps)
{
    if (!mra)
        return;

    scale_fct_dual_proj_dirichlet_forward(mra->scale_fct_dual,
                                          mra->scale_fine,
                                          function,
                                          ps);
}

void
mra_proj_v_backward (const mra_t  *mra,
                     const double *ps,
                     double       *function)
{
    if (!mra)
        return;

    scale_fct_proj_dirichlet_backward(mra->scale_fct,
                                      mra->scale_fine,
                                      ps,
                                      function);
}

/* Projections 2D */

void
mra_proj_2d_x_forward (const mra_t  *mra,
                       const double *function,
                       double       *ps,
                       double       *tmp1,
                       double       *tmp2)
{
    int32_t i;
    int32_t j;

    if (!mra)
        return;

    for (j = 0; j < mra->size_v; j++)
    {
        for (i = 0; i < mra->size_x; i++)
        {
            tmp1[i] = function[i * mra->size_v + j];
        }
    
        scale_fct_dual_proj_periodic_forward (mra->scale_fct_dual,
                                              mra->scale_fine,
                                              tmp1,
                                              tmp2);

        for (i = 0; i < mra->size_x; i++)
        {
            ps[i * mra->size_v + j] = tmp2[i];
        }
    }
}

void
mra_proj_2d_x_backward (const mra_t  *mra,
                        const double *ps,
                        double       *function,
                        double       *tmp1,
                        double       *tmp2)
{
    int32_t i;
    int32_t j;

    if (!mra)
        return;

    for (j = 0; j < mra->size_v; j++)
    {
        for (i = 0; i < mra->size_x; i++)
        {
            tmp1[i] = ps[i * mra->size_v + j];
        }
    
        scale_fct_proj_periodic_backward(mra->scale_fct,
                                         mra->scale_fine,
                                         tmp1,
                                         tmp2);

        for (i = 0; i < mra->size_x; i++)
        {
            function[i * mra->size_v + j] = tmp2[i];
        }
    }
}

void
mra_proj_2d_v_forward (const mra_t  *mra,
                       const double *function,
                       double       *ps)
{
    const double *iter_f;
    double       *iter_ps;
    int32_t       i;

    if (!mra)
        return;

    iter_f = function;
    iter_ps = ps;
    for (i = 0; i < mra->size_x; i++, iter_f += mra->size_v, iter_ps += mra->size_v)
    {
        scale_fct_dual_proj_dirichlet_forward(mra->scale_fct_dual,
                                              mra->scale_fine,
                                              iter_f,
                                              iter_ps);
    }
}

void
mra_proj_2d_v_backward (const mra_t  *mra,
                        const double *ps,
                        double       *function)
{
    double       *iter_f;
    const double *iter_ps;
    int32_t       i;

    if (!mra)
        return;

    iter_f = function;
    iter_ps = ps;
    for (i = 0; i < mra->size_x; i++, iter_f += mra->size_v, iter_ps += mra->size_v)
    {
        scale_fct_proj_dirichlet_backward (mra->scale_fct,
                                           mra->scale_fine,
                                           iter_ps,
                                           iter_f);
    }
}

void
mra_proj_2d_forward (const mra_t  *mra,
                     const double *function,
                     double       *ps,
                     double       *tmp1,
                     double       *tmp2)
{
    const double *iter_f;
    double       *iter_ps;
    int32_t       i;
    int32_t       j;

    if (!mra)
        return;

    for (j = 0; j < mra->size_v; j++)
    {
        for (i = 0; i < mra->size_x; i++)
        {
            tmp1[i] = function[i * mra->size_v + j];
        }
    
        scale_fct_dual_proj_periodic_forward (mra->scale_fct_dual,
                                              mra->scale_fine,
                                              tmp1,
                                              tmp2);

        for (i = 0; i < mra->size_x; i++)
        {
            ps[i * mra->size_v + j] = tmp2[i];
        }
    }

    iter_f = function;
    iter_ps = ps;
    for (i = 0; i < mra->size_x; i++, iter_f += mra->size_v, iter_ps += mra->size_v)
    {
        scale_fct_dual_proj_dirichlet_forward(mra->scale_fct_dual,
                                              mra->scale_fine,
                                              iter_f,
                                              iter_ps);
    }
}

void
mra_proj_2d_backward (const mra_t  *mra,
                      const double *ps,
                      double       *function,
                      double       *tmp1,
                      double       *tmp2)
{
    double       *iter_f;
    const double *iter_ps;
    int32_t       i;
    int32_t       j;

    if (!mra)
        return;

    for (j = 0; j < mra->size_v; j++)
    {
        for (i = 0; i < mra->size_x; i++)
        {
            tmp1[i] = ps[i * mra->size_v + j];
        }
    
        scale_fct_proj_periodic_backward(mra->scale_fct,
                                         mra->scale_fine,
                                         tmp1,
                                         tmp2);

        for (i = 0; i < mra->size_x; i++)
        {
            function[i * mra->size_v + j] = tmp2[i];
        }
    }

    iter_f = function;
    iter_ps = ps;
    for (i = 0; i < mra->size_x; i++, iter_f += mra->size_v, iter_ps += mra->size_v)
    {
        scale_fct_proj_dirichlet_backward (mra->scale_fct,
                                           mra->scale_fine,
                                           iter_ps,
                                           iter_f);
    }
}

/* double * */
/* mra_scale_fct_ps_get (const mra_t *mra, double coef) */
/* { */
/*   const double *weights; */
/*   double *ps; */
/*   double *coord_v; */
/*   int32_t ps_offset; */
/*   int32_t N1; */
/*   int32_t degree; */
/*   int32_t j; */
/*   int32_t k; */
/*   int32_t l; */
/*   int32_t n; */

/*   weights = scale_fct_dual_weights_get (mra->scale_fct_dual, &degree); */
/*   degree--; */

/*   ps = (double *)malloc (sizeof (double) * mra->size_v * mra->size_x * mra->size_x); */
/*   if (!ps) */
/*     return NULL; */

/*   coord_v = (double *)malloc (sizeof (double) * mra->size_v); */
/*   if (!coord_v) */
/*     { */
/*       free (ps); */
/*       return NULL; */
/*     } */

/*   for (j = 0; j < mra->size_v; j++) */
/*     coord_v[j] = (double)(-mra->size_v + 2 * j) / (double)mra->size_v; */

/*   ps_offset = 0; */
/*   N1 = scale_fct_base_N1_get ((scale_fct_base_t const *)mra->scale_fct); */

/*   for (j = 0; j < mra->size_v; j++) { */
/*     double  offset; */

/*     offset = coef * coord_v[j] * (double)mra->size_x; */
/*     for (k = mra->lambda_fine_x_inf; k <= mra->lambda_fine_x_sup; k++) { */
/*       for (l = mra->lambda_fine_x_inf; l <= mra->lambda_fine_x_sup; l++) { */
/*         double val; */

/*         val = 0.0; */
/*         for (n = -(degree >> 1); n <= degree - (degree >> 1); n++) { */
/*           double  x; */

/*           x = offset + k + n - l; */
/*           if (x < N1) */
/*             { */
/*               do { */
/*                 x += mra->size_x; */
/*               } while (x < N1); */
/*             } */
/*           if (x >= mra->lambda_fine_x_sup + 1 + N1) */
/*             { */
/*               do { */
/*                 x -= mra->size_x; */
/*               } while (x >= mra->lambda_fine_x_sup + 1 + N1); */
/*             } */
/*           val += weights[n + (degree >> 1)] * scale_fct_value_get (mra->scale_fct, x); */
/*         } */
/*         ps[ps_offset] = val; */
/*         ps_offset++; */
/*       } */
/*     } */
/*   } */

/*   free (coord_v); */

/*   return ps; */
/* } */

void
mra_advection_v (const mra_t  *mra,
                 double        coef,
                 const double *ps_in,
                 double       *ps_out)
{
    const double *weights;
    double        offset;
    int32_t       N1;
    int32_t       N2;
    int32_t       N1d;
    int32_t       N2d;
    int32_t       degree;
    int32_t       k;
    int32_t       l;
    int32_t       n;

    weights = scale_fct_dual_weights_get (mra->scale_fct_dual, &degree);
    degree--;

    N1 = scale_fct_base_N1_get ((const scale_fct_base_t *)mra->scale_fct);
    N2 = scale_fct_base_N2_get ((const scale_fct_base_t *)mra->scale_fct);
    N1d = scale_fct_base_N1_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    N2d = scale_fct_base_N2_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    offset = coef * (double)mra->size_x;

    for (k = (-mra->size_v / 2); k < (mra->size_v / 2); k++)
    {
        double val_ps;

        val_ps = 0.0;

        for (l = -N2 + N1d + k - 1; l <= N2d - N1 + k + 1; l++)
        {
            if ((l >= (-mra->size_v / 2)) && (l < (mra->size_v / 2)))
            {
                double val;

                val = 0.0;
                for (n = -(degree >> 1); n <= degree - (degree >> 1); n++)
                {
                    double tmp;

                    tmp = k - l + n + offset;
                    val += weights[n + (degree >> 1)] * scale_fct_value_get (mra->scale_fct, tmp);
                }
                val_ps += val * ps_in[l + (mra->size_v / 2)];
            }
        }
        ps_out[k + mra->size_v / 2] = val_ps;
    }
}

void
mra_fwt_x_forward (const mra_t *mra,
                   const double *proj_coef_fine,
                   double       *proj_coef_coarse,
                   double       *wavelet_coef)
{
    double *tmp;
    double *tmp_j_plus_1;
    double *tmp_j;
    const double *sf_filter_dual;
    const double *w_filter_dual;
    int32_t N1d;
    int32_t N2d;
    int32_t M1d;
    int32_t M2d;
    int32_t index;
    int32_t j;
    int32_t k;
    int32_t l;

    tmp = (double *)malloc (sizeof (double) * mra->size_x);
    if (!tmp)
        return;

    tmp_j_plus_1 = (double *)proj_coef_fine;
    tmp_j = tmp;

    N1d = scale_fct_base_N1_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    N2d = scale_fct_base_N2_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    sf_filter_dual = scale_fct_base_filter_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    M1d = wavelet_base_N1_get ((const wavelet_base_t *)mra->wavelet_dual);
    M2d = wavelet_base_N2_get ((const wavelet_base_t *)mra->wavelet_dual);
    w_filter_dual = wavelet_base_filter_get ((const wavelet_base_t *)mra->wavelet_dual);

    {
        printf ("\n\nfilter:\n");
        for (j = N1d; j <= N2d; j++)
            printf (" # %f\n", sf_filter_dual[j - N1d]);
        printf ("\n");
    }

    for (j = (mra->scale_fine - 1), index = 0; j >= mra->scale_coarse; index += (1 << j), j--)
    {
        printf ("index : %d %d\n", j, index);
        for (k = 0; k < (1 << j); k++)
        {
            double val;

            /* projection */
            val = 0;
            for (l = (-N2d + 2 * k); l <= (-N1d + 2 * k); l++)
            {
                if (l < 0 || l >= (1 << (j + 1)))
                    printf ("zut : %d %d %d %d %d\n", j, k , 1<<(j+1), l, l & ((1 << (j + 1)) - 1));
                val += sf_filter_dual[l - 2 * k + N2d] * tmp_j_plus_1[l & ((1 << (j + 1)) - 1)];
            }
            tmp_j[k] = val;

            /* wavelet */
            val = 0;
            for (l = (-M2d + 2 * k); l <= (-M1d + 2 * k); l++)
            {
                val += w_filter_dual[l - 2 * k + M2d] * tmp_j_plus_1[l & ((1 << (j + 1)) - 1)];
            }
            wavelet_coef[index + k] = val;
        }
        tmp_j_plus_1 = tmp_j;
    }
    memcpy (proj_coef_coarse, tmp_j, sizeof (double) * (1 << mra->scale_coarse));
    free (tmp);
}

void
mra_fwt_x_backward (const mra_t  *mra,
                    const double *proj_coef_coarse,
                    const double *wavelet_coef,
                    double       *proj_coef_fine)
{
    double *tmp;
    double *tmp_j_plus_1;
    double *tmp_j;
    const double *sf_filter_dual;
    const double *w_filter_dual;
    int32_t N1d;
    int32_t N2d;
    int32_t M1d;
    int32_t M2d;
    int32_t index;
    int32_t j;
    int32_t k;
    int32_t l;

    tmp = (double *)malloc (sizeof (double) * mra->size_x);
    if (!tmp)
        return;

    tmp_j_plus_1 = (double *)proj_coef_fine;
    tmp_j = tmp;
    memcpy(tmp_j, proj_coef_coarse, sizeof (double) * (1 << mra->scale_coarse));

    N1d = scale_fct_base_N1_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    N2d = scale_fct_base_N2_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    sf_filter_dual = scale_fct_base_filter_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    M1d = wavelet_base_N1_get ((const wavelet_base_t *)mra->wavelet_dual);
    M2d = wavelet_base_N2_get ((const wavelet_base_t *)mra->wavelet_dual);
    w_filter_dual = wavelet_base_filter_get ((const wavelet_base_t *)mra->wavelet_dual);

    for (j = mra->scale_coarse, index = ((1 << mra->scale_fine) - (1 << (mra->scale_coarse + 1))); j < mra->scale_fine; j++, index -= (1 << j))
    {
        printf (" ifwt : %d %d\n", j, index);
        for (k = 0; k < (1 << (j + 1)); k++)
        {
            double val;

            val = 0.0;

            for (l = (k - N2d) / 2 - 1; l <= (k - N1d) / 2 + 1; l++)
                if ((k - 2 * l >= N1d) && (k - 2 * l <= N2d))
                    val += sf_filter_dual[k - 2 * l - N1d] * tmp_j[l & ((1 << (j)) - 1)];

            for (l = (k - M2d) / 2 - 1; l <= (k - M1d) / 2 + 1; l++)
                if ((k - 2 * l >= M1d) && (k - 2 * l <= M2d))
                    val += w_filter_dual[k - 2 * l - M1d] * tmp_j[l & ((1 << (j)) - 1)];

/*       for (l = 0; l < (1 << j); l++) */
/*         if ((k - 2 * l >= N1d) && (k - 2 * l <= N2d)) */
/*           val += sf_filter_dual[k - 2 * l - N1d] * tmp_j[l]; */

/*       for (l = 0; l < (1 << j); l++) */
/*         if ((k - 2 * l >= M1d) && (k - 2 * l <= M2d)) */
/*           val += sf_filter_dual[k - 2 * l - M1d] * wavelet_coef[index + l]; */

            tmp_j_plus_1[k] = val;
        }
        tmp_j = tmp_j_plus_1;
    }
    memcpy(proj_coef_fine, tmp_j_plus_1, sizeof (double) * (1 << mra->scale_fine));
    free(tmp);
}

void
mra_fwt_v_forward (const mra_t *mra,
                   const double *proj_coef_fine,
                   double       *proj_coef_coarse,
                   double       *wavelet_coef)
{
    double *tmp;
    double *tmp_j_plus_1;
    double *tmp_j;
    const double *sf_filter_dual;
    const double *w_filter_dual;
    int32_t N1d;
    int32_t N2d;
    int32_t M1d;
    int32_t M2d;
    int32_t index;
    int32_t j;
    int32_t k;
    int32_t l;

    tmp = (double *)malloc (sizeof (double) * mra->size_x);
    if (!tmp)
        return;

    tmp_j_plus_1 = (double *)proj_coef_fine;
    tmp_j = tmp;

    N1d = scale_fct_base_N1_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    N2d = scale_fct_base_N2_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    sf_filter_dual = scale_fct_base_filter_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    M1d = wavelet_base_N1_get ((const wavelet_base_t *)mra->wavelet_dual);
    M2d = wavelet_base_N2_get ((const wavelet_base_t *)mra->wavelet_dual);
    w_filter_dual = wavelet_base_filter_get ((const wavelet_base_t *)mra->wavelet_dual);

    {
        printf ("\n\nfilter:\n");
        for (j = N1d; j <= N2d; j++)
            printf (" # %f\n", sf_filter_dual[j - N1d]);
        printf ("\n");
    }

    for (j = mra->scale_fine, index = 0; j >= mra->scale_coarse + 1; index += (1 << j), j--)
    {
        printf ("index : %d %d\n", j, index);
        for (k = 0; k < (1 << j); k++)
        {
            double val;

            /* projection */
            val = 0;
            for (l = (-N2d + 2 * k); l <= (-N1d + 2 * k); l++)
            {
                if ((l >= 0) && (l < (1 << (j + 1))))
                    val += sf_filter_dual[l - 2 * k + N2d] * tmp_j_plus_1[l];
            }
            tmp_j[k] = val;

            /* wavelet */
            val = 0;
            for (l = (-M2d + 2 * k); l <= (-M1d + 2 * k); l++)
            {
                if ((l >= 0) && (l < (1 << (j + 1))))
                    val += w_filter_dual[l - 2 * k + M2d] * tmp_j_plus_1[l];
            }
            wavelet_coef[index + k] = val;
        }
        tmp_j_plus_1 = tmp_j;
    }
    memcpy (proj_coef_coarse, tmp_j, sizeof (double) * (1 << mra->scale_coarse));
    free (tmp);
}

void
mra_fwt_v_backward (const mra_t  *mra,
                    const double *proj_coef_coarse,
                    const double *wavelet_coef,
                    double       *proj_coef_fine)
{
    double *tmp;
    double *tmp_j_plus_1;
    double *tmp_j;
    const double *sf_filter_dual;
    const double *w_filter_dual;
    int32_t N1d;
    int32_t N2d;
    int32_t M1d;
    int32_t M2d;
    int32_t index;
    int32_t j;
    int32_t k;
    int32_t l;

    tmp = (double *)malloc (sizeof (double) * mra->size_x);
    if (!tmp)
        return;

    tmp_j_plus_1 = (double *)proj_coef_fine;
    tmp_j = tmp;
    memcpy(tmp_j, proj_coef_coarse, sizeof (double) * (1 << mra->scale_coarse));

    N1d = scale_fct_base_N1_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    N2d = scale_fct_base_N2_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    sf_filter_dual = scale_fct_base_filter_get ((const scale_fct_base_t *)mra->scale_fct_dual);
    M1d = wavelet_base_N1_get ((const wavelet_base_t *)mra->wavelet_dual);
    M2d = wavelet_base_N2_get ((const wavelet_base_t *)mra->wavelet_dual);
    w_filter_dual = wavelet_base_filter_get ((const wavelet_base_t *)mra->wavelet_dual);

    for (j = mra->scale_coarse, index = ((1 << mra->scale_fine) - (1 << (mra->scale_coarse + 1))); j < mra->scale_fine; j++, index -= (1 << j))
    {
        printf (" ifwt : %d %d\n", j, index);
        for (k = 0; k < (1 << (j + 1)); k++)
        {
            double val;

            val = 0.0;

            for (l = (k - N2d) / 2 - 1; l <= (k - N1d) / 2 + 1; l++)
                if (((k - 2 * l >= N1d) && (k - 2 * l <= N2d)) &&
                    ((l >= 0) && (l < (1 << j))))
                    val += sf_filter_dual[k - 2 * l - N1d] * tmp_j[l];

            for (l = (k - M2d) / 2 - 1; l <= (k - M1d) / 2 + 1; l++)
                if (((k - 2 * l >= M1d) && (k - 2 * l <= M2d)) &&
                    ((l >= 0) && (l < (1 << j))))
                    val += w_filter_dual[k - 2 * l - M1d] * tmp_j[l];

            tmp_j_plus_1[k] = val;
        }
        tmp_j = tmp_j_plus_1;
    }
    memcpy(proj_coef_fine, tmp_j_plus_1, sizeof (double) * (1 << mra->scale_fine));
    free(tmp);
}
