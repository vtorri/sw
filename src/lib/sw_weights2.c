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

#include <stdlib.h>
#include <stdio.h>

#include "sw_polynomial.h"
#include "sw_system.h"
#include "sw.h"


/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

/**
 * @cond SW_LOCAL
 */

struct sw_weights_s
{
    double *weights;
    int32_t size;
};

/* Lagrange */

/* compute \int x^i\phi for i = 0 to order (inclusive) */
static rational_t *
_sw_moments_lagrange_get(const sw_scale_fct_base_t *sfb, int32_t order)
{
    const rational_t *filter_rat;
    rational_t       *moments;
    int32_t           N1;
    int32_t           N2;
    int32_t           m;
    int32_t           n;
    int32_t           k;
    int32_t           j;

    moments = (rational_t *)malloc((order + 1) * sizeof(rational_t));
    if (!moments)
        return NULL;

    N1 = sw_scale_fct_base_N1_get(sfb);
    N2 = sw_scale_fct_base_N2_get(sfb);
    filter_rat = sw_scale_fct_base_filter_rat_get(sfb);

    moments[0] = rat_new(1, 1, 0);
    for (m = 1; m <= order; m++)
    {
        rational_t num;
        rational_t den;

        num = rat_new(0, 1, 0);
        for (n = 0; n < m; n++)
        {
            rational_t tmp;
            rational_t coef;
            rational_t val;
            rational_t rpuiss;

            tmp = rat_new(sw_binomial(m, n), 1, 0);
            coef = rat_mul(&tmp, &moments[n]);
            val = rat_new(0, 1, 0);

            for (k = N1; k <= N2; k++)
            {
                int64_t npuiss;

                for (j = 0, npuiss = 1; j < m - n; j++, npuiss *= k) { }
                rpuiss = rat_new(npuiss, 1, 0);
                 tmp = rat_mul(&rpuiss, &filter_rat[k - N1]);
                val = rat_add(&val, &tmp);
            }
            val = rat_mul(&val, &coef);
            num = rat_add(&num, &val);
        }

        den = rat_new((1 << m) - 1, 1, 0);
        moments[m] = rat_div(&num, &den);
    }

    return moments;
}

/* Sweldens */

static rational_t *
_sw_coefs_sweldens_get(int32_t r, rational_t *lambda, int32_t *size)
{
    rational_t *res;
    rational_t *w1;
    rational_t *w2;
    rational_t *w3;
    rational_t  r1;
    rational_t  r2;
    int         i;

    if (r <= 0)
        return NULL;

    res = (rational_t *)malloc((sizeof(rational_t) * (r + 1) * (r + 2)) / 2);
    if (!res)
        return NULL;

    *size = ((r + 1) * (r + 2)) / 2;

    w1 = res;

    /* i = 0 */
    w1[0] = rat_new(1, 1, 0);

    /* i = 1 */
    w1[1] = *lambda;
    w1[2] = rat_new(1, 1, 0);

    /* i = 2 */

    /* 2 * lambda * lambda - 3 */
    r1 = rat_new(2, 1, 0);
    r1 = rat_mul(&r1, lambda);
    r1 = rat_mul(&r1, lambda);
    r2 = rat_new(-3, 1, 0);
    r1 = rat_add(&r1, &r2);
    w1[3] = r1;
    /* 4 * lambda */
    r1 = rat_new(4, 1, 0);
    r1 = rat_mul(&r1, lambda);
    w1[4] = r1;
    w1[5] = rat_new(1, 1, 0);

    w1 = res + 3;
    w2 = res + 1;
    w3 = res;

    for (i = 3; i <= r; i++)
    {
        int         j;

        w1 += i;
        w2 += i - 1;
        w3 += i - 2;

        r1 = rat_new(2, 1, 0);
        r1 = rat_mul(&r1, lambda);
        r2 = rat_new(-4, 1, 0);
        r1 = rat_add(&r1, &r2);
        r1 = rat_mul(&r1, w2);
        r1 = rat_add(&r1, w2 + 1);
        w1[0] = r1;
        r1 = rat_new(2, 1, 0);
        r1 = rat_mul(&r1, w2);
        r2 = rat_new(2, 1, 0);
        r2 = rat_mul(&r2, lambda);
        r2 = rat_mul(&r2, w2 + 1);
        r1 = rat_add(&r1, &r2);
        r1 = rat_add(&r1, w2 + 2);
        r2 = rat_new(-4, 1, 0);
        r2 = rat_mul(&r2, w3 + 1);
        r1 = rat_add(&r1, &r2);
        w1[1] = r1;

        for (j = 2; j < (i - 1); j++)
        {
            r1 = w2[j - 1];
            r1 = rat_add(&r1, w2 + j + 1);
            r2 = rat_new(2, 1, 0);
            r2 = rat_mul(&r2, lambda);
            r2 = rat_mul(&r2, w2 + j);
            r1 = rat_add(&r1, &r2);
            r2 = rat_new(-4, 1, 0);
            r2 = rat_mul(&r2, w3 + j);
            r1 = rat_add(&r1, &r2);
            w1[j] = r1;
        }
        r1 = w2[i - 2];
        r2 = rat_new(2, 1, 0);
        r2 = rat_mul(&r2, lambda);
        r2 = rat_mul(&r2, w2 + i - 1);
        r1 = rat_add(&r1, &r2);
        w1[i - 1] = r1;
        w1[i] = rat_new(1, 1, 0);
    }

    return res;
}

static rational_t *
_sw_moments_sweldens_get(const sw_scale_fct_base_t *sfb, int r)
{
    const rational_t *filter_rat;
    rational_t       *moments;
    int32_t           N1;
    int32_t           N2;
    int32_t           L;
    int32_t           i;
    int32_t           j;
    int32_t           k;
    int32_t           idx;

    moments = (rational_t *)malloc((r + 1) * sizeof(rational_t));
    if (!moments)
        return NULL;

    N1 = sw_scale_fct_base_N1_get(sfb);
    N2 = sw_scale_fct_base_N2_get(sfb);
    L = N2 - N1;
    filter_rat = sw_scale_fct_base_filter_rat_get(sfb);

    moments[0] = rat_new(1, 1, 0);
    idx = 0;
    for (i = 1; i <= r; i++)
    {
        rational_t num;
        rational_t den;

        idx += i;
        num = rat_new(0, 1, 0);

        for (k = N1; k <= N2; k++)
        {
            rational_t *coefs;
            rational_t  lambda;
            rational_t  tmp;
            rational_t  val;
            int         size;

            tmp = rat_new(1, 1, 0);
            lambda = rat_new(2 * (k - N1), L, 1);
            lambda = rat_sub(&lambda, &tmp);

            coefs = _sw_coefs_sweldens_get(r, &lambda, &size);
            if (!coefs)
            {
                free(moments);
                return NULL;
            }

            val = rat_new(0, 1, 0);
            for (j = 0; j < i; j++)
            {
                tmp = rat_mul(coefs + idx + j, moments + j);
                val = rat_add(&val, &tmp);
            }
            val = rat_mul(&val, &filter_rat[k - N1]);
            num = rat_add(&num, &val);
            free(coefs);
        }

        den = rat_new((1 << i) - 1, 1, 0);
        moments[i] = rat_div(&num, &den);
    }

    return moments;
}

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

sw_weights_t *
sw_weights_lagrange_new(const sw_scale_fct_base_t *sfb, int32_t degree)
{
    sw_weights_t *weights;
    rational_t   *moments;
    int64_t      *x;
    int32_t       i;

    if (!sfb || (degree <= 0))
        return NULL;

    weights = (sw_weights_t *)malloc(sizeof (sw_weights_t));
    if (!weights)
        return NULL;

    weights->weights = (double *)malloc((degree + 1) * sizeof(double));
    if (!weights->weights)
        goto free_weights;
    weights->size = degree + 1;

    moments = _sw_moments_lagrange_get(sfb, degree);
    if (!moments)
        goto free_weights_weights;

    x = (int64_t *)malloc((degree + 1) * sizeof (int64_t));
    if (!x)
        goto free_moments;

    for (i = 0; i <= degree; i++)
        x[i] = i - (degree >> 1);

    for (i = 0; i <= degree; i++)
    {
        rpol_t    *lagrange;
        rational_t r;
        int32_t    k;

        lagrange = rpol_lagrange_base_get (x, degree + 1, i);
        r = rat_new (0, 1, 0);
        for (k = 0; k <= degree; k++)
        {
            rational_t coef;

            coef = rpol_coef_get (lagrange, k);
            coef = rat_mul (&moments[k], &coef);
            r = rat_add (&r, &coef);
        }
        weights->weights[i] = rat_double_get (&r);
        rpol_delete (lagrange);
    }

    free (x);
    free (moments);

    return weights;

  free_moments:
    free(moments);
  free_weights_weights:
    free(weights->weights);
  free_weights:
    free(weights);

    return NULL;
}

sw_weights_t *
sw_weights_sweldens_new(const sw_scale_fct_base_t *sfb, int32_t r, int32_t s)
{
    sw_weights_t *weights;
    system_t     *system;
    rpol_t      **cheb;
    rational_t   *moments;
    rational_t   *matrix;
    rational_t   *iter;
    rational_t   *abscisses; /* x_k^* */
    rational_t   *sol;
    rational_t    tau_s;     /* Tau^* */
    int32_t       N1;
    int32_t       N2;
    int32_t       L;
    int32_t       i;
    int32_t       k;

    if (!sfb || (r < 0) || (s < 0))
        return NULL;

    N1 = sw_scale_fct_base_N1_get(sfb);
    N2 = sw_scale_fct_base_N2_get(sfb);
    L = N2 - N1;

    /* tau_s = (dr_s - 1) / 2 */
    tau_s = rat_new((r - 1) - L * (1 << (s - 1)), L * (1 << s), 1);

    weights = (sw_weights_t *)malloc(sizeof (sw_weights_t));
    if (!weights)
        return NULL;

    weights->weights = (double *)malloc(r * sizeof(double));
    if (!weights->weights)
        goto free_weights;
    weights->size = r;

    moments = _sw_moments_sweldens_get(sfb, r - 1);
    if (!moments)
        goto free_weights_weights;

    cheb = sw_rpol_chebychev_pols_get(r - 1);
    if (!cheb)
        goto free_moments;

    abscisses = (rational_t *)malloc(r * sizeof(rational_t));
    if (!abscisses)
        goto free_cheb;

    for (k = 0; k < r; k++)
    {
        /* d_k^* = 2.k/(L.2^s) */
        abscisses[k] = rat_new(k, L * (1 << (s - 1)), 1);
        abscisses[k] = rat_sub(abscisses + k, &tau_s);
    }

    matrix = (rational_t *)malloc(r * r * sizeof(rational_t));
    if (!matrix)
        goto free_abscisses;

    iter = matrix;
    for (i = 0; i < r; i++)
    {
        for (k = 0; k < r; k++, iter++)
        {
            *iter = rpol_value_get(cheb[i], abscisses + k);
        }
    }

    system = system_new(matrix, moments, r);
    if (!system)
        goto free_matrix;

    system_solve(system);
    sol = (rational_t *)system_solution_get(system);

    for (i = 0; i < r; i++)
        weights->weights[i] = rat_double_get(sol + i);

    system_delete(system);
    free(abscisses);
    for (i = 0; i < r; i++)
        rpol_delete(cheb[i]);
    free(cheb);

    return weights;

  free_matrix:
    free(matrix);
  free_abscisses:
    free(abscisses);
  free_cheb:
    for (i = 0; i < r; i++)
        rpol_delete(cheb[i]);
    free(cheb);
  free_moments:
    free(moments);
  free_weights_weights:
    free(weights->weights);
  free_weights:
    free(weights);

    return NULL;
}

void
sw_weights_del (sw_weights_t *w)
{
    if (!w)
        return;

    if (w->weights)
        free (w->weights);
    free (w);
}

const double *
sw_weights_get(const sw_weights_t *w,
               int32_t            *weights_size)
{
    if (!w || !weights_size)
        return NULL;
    *weights_size = w->size;
    return w->weights;
}
