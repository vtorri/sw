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

#include "sw_weights.h"
#include "sw_sweldens.h"
#include "sw_polynomial.h"


/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

/**
 * @cond SW_LOCAL
 */

struct weights
{
    weights_type_t          type;
    const scale_fct_base_t *sfb;

    union
    {
        struct
        {
            int32_t             degree;
        }lagrange;

        struct
        {
            int32_t             order;
            int32_t             scale;
            rational_t          tau;
        }sweldens;
    } data;

    double                 *weights;
    double                 *absciss;
    int32_t                 size;
};


/*
 * Util
 */

static int64_t
binomial (uint32_t n, uint32_t p)
{
    int64_t num;
    int64_t den;
    int64_t m;
    int64_t i;

    num = 1;
    den = 1;
    m = (p < (n >> 1)) ? p : n - p;
    for (i = 0; i < m; i++)
    {
        num *= (n - i);
        den *= (i + 1);
    }

    return num / den;
}

/*
 * Weight methods
 */

static rational_t *
_weights_lagrange_moments_get (const weights_t *w)
{
    const rational_t *filter_rat;
    rational_t *moments;
    int32_t     N1;
    int32_t     N2;
    int32_t     degree;
    int32_t     m;
    int32_t     n;
    int32_t     k;
    int32_t     j;

    degree = w->data.lagrange.degree;

    moments = (rational_t *)malloc (sizeof (rational_t) * (degree + 1));
    if (!moments)
        return NULL;

    N1 = scale_fct_base_N1_get (w->sfb);
    N2 = scale_fct_base_N2_get (w->sfb);
    filter_rat = scale_fct_base_filter_rat_get (w->sfb);

    moments[0] = rat_new (1, 1, 0);
    for (m = 1; m <= degree; m++)
    {
        rational_t num;
        rational_t den;

        num = rat_new (0, 1, 0);
        for (n = 0; n < m; n++)
        {
            rational_t tmp;
            rational_t coef;
            rational_t val;
            rational_t rpuiss;

            tmp = rat_new (binomial (m, n), 1, 0);
            coef = rat_mul (&tmp, &moments[n]);
            val = rat_new (0, 1, 0);

            for (k = N1; k <= N2; k++)
            {
                int64_t npuiss;

                for (j = 0, npuiss = 1; j < m - n; j++, npuiss *= k) { }
                rpuiss = rat_new (npuiss, 1, 0);
                tmp = rat_mul (&rpuiss, &filter_rat[k - N1]);
                val = rat_add (&val, &tmp);
            }
            val = rat_mul (&val, &coef);
            num = rat_add (&num, &val);
        }

        den = rat_new ((1 << m) - 1, 1, 0);
        moments[m] = rat_div (&num, &den);
    }

    return moments;
}

static double *
_weights_lagrange_weights_get (const weights_t *w)
{
    rational_t *moments;
    double     *weights;
    int64_t    *x;
    int32_t     degree;
    int32_t     i;

    moments = _weights_lagrange_moments_get (w);
    if (!moments)
        return NULL;

    degree = w->data.lagrange.degree;

    x = (int64_t *)malloc (sizeof (int64_t) * (degree + 1));
    if (!x)
    {
        free (moments);
        return NULL;
    }

    printf("size lag : %ld\n", (sizeof (double) * (degree + 1)));
    weights = (double *)malloc (sizeof (double) * (degree + 1));
    if (!weights)
    {
        free (moments);
        free (x);
        return NULL;
    }

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
        weights[i] = rat_double_get (&r);
        rpol_delete (lagrange);
    }

    free (moments);
    free (x);

    return weights;
}

static double *
_weights_sweldens_weights_get (const weights_t *w)
{
    rational_t  tau;
    rational_t *rat_weights;
    double     *weights;
    sweldens_t *sweldens;
    int32_t     order;
    int32_t     scale;
    int32_t     i;

    order = w->data.sweldens.order;
    scale = w->data.sweldens.scale;
    tau   = w->data.sweldens.tau;

    sweldens = sweldens_new (w->sfb, order, scale, &tau);

    rat_weights = sweldens_weights_get (sweldens);

    weights = (double *)malloc (sizeof (double) * (order));
    if (!weights)
    {
        sweldens_delete (sweldens);
        return NULL;
    }

    for (i = 0; i < order; i++)
        weights[i] = rat_double_get (&rat_weights[i]);

    sweldens_delete (sweldens);

    return weights;
}

/**
 * @endcond SW_LOCAL
 */


/******************************************************************************/
/*                                                                            */
/*                                   GLOBAL                                   */
/*                                                                            */
/******************************************************************************/

#define __UNUSED__ __attribute__ ((unused))

/******************************************************************************/
/*                                                                            */
/*                                    API                                     */
/*                                                                            */
/******************************************************************************/

weights_t *
weights_new (weights_type_t          type,
	     const scale_fct_base_t *sfb)
{
    weights_t *weights;

    if (!sfb)
        return NULL;

    weights = (weights_t *)malloc (sizeof (weights_t));
    if (!weights)
        return NULL;

    weights->type = type;
    weights->sfb = sfb;
    weights->weights = NULL;
    weights->absciss = NULL;

    switch (type)
    {
     case SW_WEIGHTS_TYPE_LAGRANGE:
         weights->size = weights->data.lagrange.degree + 1;
	 printf("new weights->weights\n");
         weights->weights = _weights_lagrange_weights_get(weights);
         break;
     case SW_WEIGHTS_TYPE_SWELDENS:
         weights->size = weights->data.sweldens.order;
         weights->weights = _weights_sweldens_weights_get(weights);
         break;
    }

    return weights;
}

void
weights_del (weights_t *w)
{
    if (w->weights)
    {
        printf("del weights->weights\n");
        free (w->weights);
    }
    free (w);
}

void
weights_lagrange_data_set (weights_t *w,
			   int32_t    degree)
{
    if (!w)
        return;

    if (w->type != SW_WEIGHTS_TYPE_LAGRANGE)
    {
        fprintf (stderr, "[sw] [weights] ERROR: wrong type\n");
        return;
    }

    w->data.lagrange.degree = degree;
    w->size = w->data.lagrange.degree + 1;
    w->weights = _weights_lagrange_weights_get(w);
}

void
weights_sweldens_data_set (weights_t        *w,
			   int32_t           order,
			   int32_t           scale,
			   const rational_t *tau)
{
    if (!w)
        return;

    if (w->type != SW_WEIGHTS_TYPE_SWELDENS)
    {
        fprintf (stderr, "[sw] [weights] ERROR: wrong type\n");
        return;
    }

    w->data.sweldens.order = order;
    w->data.sweldens.scale = scale;
    w->data.sweldens.tau = *tau;
}

const double *
weights_get (const weights_t *w,
	     int32_t         *weights_size)
{
    if (!w || !weights_size)
        return NULL;
    *weights_size = w->size;
    return w->weights;
}

double *
weights_absciss_get (const weights_t *w,
		     int32_t          start __UNUSED__,
		     int32_t          end __UNUSED__)
{
    double  *absciss;
    int32_t  i;

    absciss = (double *)malloc (sizeof (double) * w->size);
    if (!absciss)
        return NULL;

    switch (w->type)
    {
     case SW_WEIGHTS_TYPE_LAGRANGE:
         return absciss;
     case SW_WEIGHTS_TYPE_SWELDENS:
         for (i = 0; i < w->size; i++)
         {
             absciss[i] = (double)(i << w->data.sweldens.scale) - rat_double_get (&w->data.sweldens.tau);
         }
         return absciss;
     default:
         free (absciss);
         return NULL;
    }
}
