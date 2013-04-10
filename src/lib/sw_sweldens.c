#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sw_system.h"
#include "sw.h"

/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

struct sweldens
{
    scale_fct_base_t const *sfb;
    rational_t             *weights;
    rational_t             *absciss;
    int32_t                 order;
    int32_t                 scale;
};


static int32_t
_lookup0 (sweldens_t const *s, int32_t i, int32_t j)
{
    /* int32_t r; */

    /* r = s->order; */

    return i * (i - 1) / 2 + j;
}

static int32_t
_lookup1 (sweldens_t const *s, int32_t k, int32_t i, int32_t j)
{
    int32_t N1;
    int32_t r;

    N1 = scale_fct_base_N1_get (s->sfb);
    r = s->order;

    return (k - N1) * (r + 1) * (r + 2) / 2 + i * (i + 1) / 2 + j;
}

static int32_t
_lookup2 (int32_t m, int32_t i, int32_t j)
{
    return m * (m + 1) * (m + 2) / 6 + i * (m - i + 1) + i * (i + 1) / 2 + j;
}

static rational_t *
_coefs_get (sweldens_t const *s)
{
    rational_t *coefs;
    int32_t     N1;
    int32_t     N2;
    int32_t     r;
    int32_t     k;

    N1 = scale_fct_base_N1_get (s->sfb);
    N2 = scale_fct_base_N2_get (s->sfb);
    r = s->order;

    coefs = (rational_t *)malloc (sizeof (rational_t) * (N2 - N1 + 1) * (r + 1) * (r + 2) / 2);
    if (!coefs)
        return NULL;

    for (k = N1; k <= N2; k++)
    {
        rational_t lambda;
        rational_t un;
        rational_t deux;
        rational_t trois;
        rational_t quatre;
        rational_t tmp1;
        rational_t tmp2;
        int32_t    i;

        un     = rat_new (1, 1, 0);
        deux   = rat_new (2, 1, 0);
        trois  = rat_new (3, 1, 0);
        quatre = rat_new (4, 1, 0);

        lambda = rat_new (2 * k - N1 - N2, N2 - N1, 1);

        /* i = 0 */
        coefs[_lookup1 (s, k, 0, 0)] = un;

        /* i = 1 */
        coefs[_lookup1 (s, k, 1, 0)] = lambda;
        coefs[_lookup1 (s, k, 1, 1)] = un;

        /* i = 2 */
        tmp1 = rat_mul (&lambda, &lambda);
        tmp1 = rat_mul (&deux, &tmp1);
        tmp1 = rat_sub (&tmp1, &trois);
        coefs[_lookup1 (s, k, 2, 0)] = tmp1;
        tmp1 = rat_mul (&quatre, &lambda);
        coefs[_lookup1 (s, k, 2, 1)] = tmp1;
        coefs[_lookup1 (s, k, 2, 2)] = un;

        for (i = 3; i <= r; i++)
        {
            int32_t    j;

            /* j = 0 */
            tmp1 = rat_mul (&deux, &lambda);
            tmp1 = rat_sub (&tmp1, &quatre);
            tmp1 = rat_mul (&tmp1, &coefs[_lookup1 (s, k, i - 1, 0)]);
            tmp1 = rat_add (&tmp1, &coefs[_lookup1 (s, k, i - 1, 1)]);
            coefs[_lookup1 (s, k, i, 0)] = tmp1;

            /* j = 1 */
            tmp1 = rat_mul (&deux, &coefs[_lookup1 (s, k, i - 1, 0)]);
            tmp1 = rat_add (&tmp1, &coefs[_lookup1 (s, k, i - 1, 2)]);
            tmp2 = rat_mul (&deux, &lambda);
            tmp2 = rat_mul (&tmp2, &coefs[_lookup1 (s, k, i - 1, 1)]);
            tmp1 = rat_add (&tmp1, &tmp2);
            tmp2 = rat_mul (&quatre, &coefs[_lookup1 (s, k, i - 2, 1)]);
            tmp1 = rat_sub (&tmp1, &tmp2);
            coefs[_lookup1 (s, k, i, 1)] = tmp1;

            for (j = 2; j <= i - 2; j++)
            {
                tmp1 = rat_mul (&deux, &lambda);
                tmp1 = rat_mul (&tmp1, &coefs[_lookup1 (s, k, i - 1, j)]);
                tmp2 = rat_mul (&quatre, &coefs[_lookup1 (s, k, i - 2, j)]);
                tmp1 = rat_sub (&tmp1, &tmp2);
                tmp1 = rat_add (&tmp1, &coefs[_lookup1 (s, k, i - 1, j + 1)]);
                tmp1 = rat_add (&tmp1, &coefs[_lookup1 (s, k, i - 1, j - 1)]);
                coefs[_lookup1 (s, k, i, j)] = tmp1;
            }

            /* j = i - 1 */
            tmp1 = rat_mul (&deux, &lambda);
            tmp1 = rat_mul (&tmp1, &coefs[_lookup1 (s, k, i - 1, i - 1)]);
            tmp1 = rat_add (&tmp1, &coefs[_lookup1 (s, k, i - 1, i - 2)]);
            coefs[_lookup1 (s, k, i, i - 1)] = tmp1;

            /* j = i */
            coefs[_lookup1 (s, k, i, i)] = un;
        }
    }

    {
        for (k = N1; k <= N2; k++)
        {
            int32_t i;
            for (i = 0; i <= r; i++)
            {
                int32_t j;
                for (j = 0; j <= i; j++)
                {
                    printf ("k=%d, i=%d, j=%d : ", k, i, j);
                    rat_disp (&coefs[_lookup1 (s, k, i, j)], 1);
                }
            }
        }
    }

    return coefs;
}

static rational_t *
_weights_get (sweldens_t const *s)
{
    rational_t       *coefs;
    rational_t       *weights;
    rational_t const *filter;
    int32_t           N1;
    int32_t           N2;
    int32_t           r;
    int32_t           i;

    N1 = scale_fct_base_N1_get (s->sfb);
    N2 = scale_fct_base_N2_get (s->sfb);
    filter = scale_fct_base_filter_rat_get (s->sfb);
    r = s->order;

    coefs = _coefs_get (s);
    if (!coefs)
        return NULL;

    weights = (rational_t *)malloc (sizeof (rational_t) * (r + 1) * (r + 2) / 2);
    if (!weights)
    {
        free (coefs);
        return NULL;
    }

    weights[0] = rat_new (1, 1, 0);

    for (i = 1; i <= r; i++)
    {
        int32_t    j;

        for (j = 0; j < i; j++)
        {
            rational_t res;
            int32_t    k;

            res = rat_new (0, 1, 0);
            for (k = N1; k <= N2; k++)
            {
                rational_t tmp;

                tmp = rat_mul (&coefs[_lookup1 (s, k, i, j)], &filter[k - N1]);
                res = rat_add (&res, &tmp);
            }
            weights[_lookup0 (s, i, j)] = res;
        }
    }

    free (coefs);

    {
        int32_t ii;
        for (ii = 0; ii <= r; ii++)
        {
            int32_t j;
            for (j = 0; j < ii; j++)
            {
                printf ("i=%d, j=%d : ", ii, j);
                rat_disp (&weights[_lookup0 (s, ii, j)], 1);
            }
        }
    }

    return weights;
}

static rational_t *
_moments_get (sweldens_t const *s)
{
    rational_t *weights;
    rational_t *moments;
    int32_t     i;

    weights = _weights_get (s);
    if (!weights)
        return NULL;

    moments = (rational_t *)malloc (sizeof (rational_t) * (s->order + 1));
    if (!moments)
    {
        free (weights);
        return NULL;
    }

    moments[0] = rat_new (1, 1, 0);

    for (i = 1; i <= s->order; i++)
    {
        rational_t  puiss;
        rational_t  res;
        int32_t     j;

        puiss = rat_new (1, (1 << i) - 1, 0);
        res = rat_new (0, 1, 0);

        for (j = 0; j < i; j++)
        {
            rational_t  tmp;

            tmp = rat_mul (&moments[j], &weights[_lookup0 (s, i, j)]);
            printf ("j : %d : ", j);
            rat_disp (&moments[j], 0);
            printf ("  ");
            rat_disp (&weights[_lookup0 (s, i, j)], 1);
            res = rat_add (&res, &tmp);
        }
        res = rat_mul (&res, &puiss);
        moments[i] = res;
    }

    free (weights);

    {
        int32_t ii;

        for (ii = 0; ii <= s->order; ii++)
        {
            printf ("moment %d : ", ii);
            rat_disp (&moments[ii], 1);
        }
    }

    return moments;
}

static rational_t *
_coefs_pol_get (sweldens_t const *s)
{
    rational_t *coefs;
    int32_t     N1;
    int32_t     N2;
    int32_t     r;
    int32_t     size;
    int32_t     m;

    N1 = scale_fct_base_N1_get (s->sfb);
    N2 = scale_fct_base_N2_get (s->sfb);
    r = s->order;
    size = (1 + r) * (12 + r * (10 + 2 * r)) / 12;

    coefs = (rational_t *)calloc (size, sizeof (rational_t));
    if (!coefs)
    {
        return NULL;
    }

    coefs[0] = rat_new (1, 1, 0);

    for (m = 1; m < size; m++)
        coefs[m] = rat_new (0, 1, 0);

    for (m = 1; m <= r; m++)
    {
        rational_t dm;
        rational_t un;
        int32_t    i;

        un = rat_new (1, 1, 0);
        dm = rat_new (2 * (m - 1), (N2 - N1) * r, 1);
        dm = rat_sub (&dm, &un);

        for (i = 0; i <= m; i++)
        {
            int32_t j;

            for (j = 0; j <= (m - i); j++)
            {

                if (i >= 1)
                {
                    printf ("a\n");
                    coefs[_lookup2(m, i, j)] = rat_add (&coefs[_lookup2(m, i, j)], &coefs[_lookup2(m - 1, i - 1, j)]);
                }
                if ((i <= (m - 2)) &&
                    (j <= (m - i - 2)))
                {
                    printf ("b\n");
                    coefs[_lookup2(m, i, j)] = rat_add (&coefs[_lookup2(m, i, j)], &coefs[_lookup2(m - 1, i + 1, j)]);
                }
                if ((i <= (m - 1)) &&
                    (j >= 1))
                {
                    printf ("c\n");
                    coefs[_lookup2(m, i, j)] = rat_add (&coefs[_lookup2(m, i, j)], &coefs[_lookup2(m - 1, i, j - 1)]);
                }
                if (j <= (m - i - 2))
                {
                    printf ("d\n");
                    coefs[_lookup2(m, i, j)] = rat_add (&coefs[_lookup2(m, i, j)], &coefs[_lookup2(m - 1, i, j + 1)]);
                }
                if ((i <= (m - 1)) &&
                    (j <= (m - i - 1)))
                {
                    rational_t tmp;

                    printf ("e\n");
                    tmp = rat_new (2, 1, 0);
                    tmp = rat_mul (&tmp, &dm);
                    tmp = rat_mul (&tmp, &coefs[_lookup2(m - 1, i, j)]);
                    coefs[_lookup2(m, i, j)] = rat_add (&coefs[_lookup2(m, i, j)], &tmp);
                }
                if (i == 1)
                {
                    printf ("f\n");
                    coefs[_lookup2(m, i, j)] = rat_add (&coefs[_lookup2(m, i, j)], &coefs[_lookup2(m - 1, 0, j)]);
                }
                if (j == 1)
                {
                    printf ("g : %d\n", _lookup2(m - 1, i, 0));
                    rat_disp (&coefs[_lookup2(m - 1, i, 0)], 1);
                    coefs[_lookup2(m, i, j)] = rat_add (&coefs[_lookup2(m, i, j)], &coefs[_lookup2(m - 1, i, 0)]);
                }
                printf ("fin : %d %d %d\n", m, i, j);
            }
        }
    }

    {
        int32_t mm;

        for (mm = 0; mm <= r; mm++)
        {
            int32_t    i;

            for (i = 0; i <= mm; i++)
            {
                int32_t j;

                for (j = 0; j <= (mm - i); j++)
                {
                    printf ("coefs pol %d %d %d : ", mm, i, j);
                    rat_disp (&coefs[_lookup2(mm, i, j)], 1);
                }
            }
        }
    }

    return coefs;
}

static rational_t
_sweldens_chebyshev_value_get (int32_t n, rational_t const *x)
{
    rational_t un;
    rational_t deux;
    rational_t tmp1;
    rational_t tmp2;
    rational_t res;
    int32_t    i;

    un = rat_new (1, 1, 0);

    if (n == 0)
        return un;

    if (n == 1)
        return *x;

    deux = rat_new (2, 1, 0);
    tmp1 = *x;
    tmp2 = un;
    for (i = 2; i <= n; i++)
    {
        res = rat_mul (&deux, &tmp1);
        res = rat_mul (&res, x);
        res = rat_sub (&res, &tmp2);
        tmp2 = tmp1;
        tmp1 = res;
    }

    return res;
}

static int8_t
_sweldens_weights_get (sweldens_t const *s,
                       rational_t const *tau)
{
    rational_t *moments;
    rational_t *matrix;
    system_t   *system;
    int32_t     N1;
    int32_t     N2;
    int32_t     L;
    int32_t     r;
    int32_t     scale;
    int32_t     i;
    int32_t     k;

    N1 = scale_fct_base_N1_get (s->sfb);
    N2 = scale_fct_base_N2_get (s->sfb);
    L = N2 - N1;

    r = s->order;
    scale = s->scale;

    if (!s)
        return 0;

    {
        rational_t dr;

        dr = rat_new (0, 1, 0);
        if (rat_greater (tau, &dr))
            return 0;

        dr = rat_new (2 * (r - 1) - (L << scale), L << scale, 1);
        if (rat_lesser (tau, &dr))
            return 0;
    }

    moments = _moments_get (s);
    if (!moments)
        return 0;

    matrix = (rational_t *)malloc (sizeof (rational_t) * r * r);
    if (!matrix)
    {
        free (moments);
        return 0;
    }

    for (k = 0; k < r; k++)
    {
        rational_t dk;

        dk = rat_new (2 * k - (L << scale), L << scale, 1);
        dk = rat_sub (&dk, tau);
        s->absciss[k] = dk;
        for (i = 0; i < r; i++)
        {
            matrix[k * r + i] = _sweldens_chebyshev_value_get (i, &dk);
        }
    }

    system = system_new (matrix, moments, r);
    system_solve (system);
    memcpy (s->weights, system_solution_get (system), sizeof (rational_t) * r);
    system_delete (system);

    return 1;
}


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


sweldens_t *
sweldens_new (scale_fct_base_t const *sfb,
              int32_t                 order,
              int32_t                 scale,
              rational_t const       *tau)
{
    sweldens_t *s;

    if (!sfb || (order <= 0))
        return NULL;

    s = (sweldens_t *)malloc (sizeof (sweldens_t));
    if (!s)
        return NULL;

    s->sfb = sfb;
    s->order = order;
    s->scale = scale;

    s->weights = (rational_t *)malloc (sizeof (rational_t) * order);
    if (!s->weights)
    {
        free (s);
        return NULL;
    }

    s->absciss = (rational_t *)malloc (sizeof (rational_t) * order);
    if (!s->absciss)
    {
        free (s->weights);
        free (s);
        return NULL;
    }

    _sweldens_weights_get (s, tau);

    _coefs_pol_get (s);

    return s;
}

void
sweldens_delete (sweldens_t *s)
{
    if (!s)
        return;

    free (s);
}

int32_t
sweldens_order_get (sweldens_t const *s)
{
    return s->order;
}

rational_t *
sweldens_absciss_get (sweldens_t const *s)
{
    if (!s)
        return NULL;

    return s->absciss;
}

rational_t *
sweldens_weights_get (sweldens_t const *s)
{
    if (!s)
        return NULL;

    return s->weights;
}
