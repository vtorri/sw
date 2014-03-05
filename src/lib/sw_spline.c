#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

struct spline_s
{
    int32_t     order;
    int64_t     x1;
    int64_t     x2;
    rational_t *rat_coefs;
    double     *coefs;
    rational_t *rat_coefs_primitive;
    double     *coefs_primitive;
};

static rational_t *
_sw_spline_primitive_rat_get(sw_spline_t *spline)
{
    rational_t *primitive;
    rational_t  tmp1;
    rational_t  tmp2;
    rational_t  tmp3;
    rational_t  val;
    int32_t     i;
    int32_t     k;
    int32_t     l;

    primitive = (rational_t *)malloc(spline->order * (spline->order + 1) * sizeof(rational_t));
    if (!primitive)
        return NULL;

    for (k = 0; k < spline->order; k++)
        primitive[k * (spline->order + 1)] = rat_new(0, 1, 0);

    for (k = 0; k < spline->order; k++)
        for (i = 1; i <= spline->order; i++)
        {
            tmp1 = rat_new(i, 1, 0);
            val = rat_div(&spline->rat_coefs[k * spline->order + i - 1], &tmp1);
            primitive[k * (spline->order + 1) + i] = val;
        }

    val = rat_new(spline->x1, 1, 0);
    tmp1 = rat_new(0, 1, 0);
    for (i = 1; i <= spline->order; i++)
    {
        tmp2 = primitive[i];
        for (l = 1; l < i; l++)
            tmp2 = rat_mul(&tmp2, &val);
        tmp1 = rat_sub(&tmp1, &tmp2);
    }
    primitive[0] = tmp1;

    for (k = 1; k < spline->order; k++)
    {
        val = rat_new(spline->x1 + k, 1, 0);
        tmp1 = rat_new(0, 1, 0);
        for (i = 0; i <= spline->order; i++)
        {
            tmp2 = primitive[(k - 1) * (spline->order + 1) + i];
            for (l = 1; l < i; l++)
                tmp2 = rat_mul(&tmp2, &val);
            tmp1 = rat_add(&tmp1, &tmp2);
        }
        tmp3 = rat_new(0, 1, 0);
        for (i = 1; i <= spline->order; i++)
        {
            tmp2 = primitive[k * (spline->order + 1) + i];
            for (l = 1; l < i; l++)
                tmp2 = rat_mul(&tmp2, &val);
            tmp3 = rat_add(&tmp3, &tmp2);
        }
        primitive[k * (spline->order + 1)] = rat_sub(&tmp1, &tmp3);
    }

    return primitive;
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

sw_spline_t *
sw_spline_new (int32_t order)
{
    rational_t        integral;
    sw_spline_t      *spline;
    system_t         *s;
    rpol_t           *pol;
    rational_t       *matrix;
    rational_t       *vector;
    rational_t       *tmp;
    const rational_t *sol;
    int               offset_smooth;
    int               size;
    int               i;
    int               k;
    int               l;

    if (order <= 0)
        return NULL;

    tmp = (rational_t *)malloc(sizeof(rational_t) * order);
    if (!tmp)
        return NULL;
    i = order - 1;
    do
    {
        tmp[i] = rat_new(1, 1, 0);
    } while (--i >= 0);

    pol = rpol_new(tmp, order - 1);
    if (!pol)
    {
        free(tmp);
        return NULL;
    }

    size = order * order;
    matrix = (rational_t *)calloc(size * size, sizeof(rational_t));
    if (!matrix)
    {
        rpol_delete(pol);
        return NULL;
    }

    vector = (rational_t *)calloc(size, sizeof(rational_t));
    if (!vector)
    {
        rpol_delete(pol);
        free(matrix);
        return NULL;
    }

    /* built of the matrix */
    i = size * size - 1;
    do
    {
        matrix[i] = rat_new(0, 1, 0);
    } while (--i >= 0);
    *matrix = rat_new(1, 1, 0);
    offset_smooth = 0;
    for (l = 0; l < order - 1; l++)
    {
        rational_t tmp1;
        rational_t tmp2;
        int        offset_pt = 0;

        for (i = 0; i < order; i++)
        {
            tmp1 = rat_new(-(order / 2), 1, 0);
            matrix[(1 + offset_smooth) * size + i] = rpol_monom_get(pol, i, l, &tmp1);
        }
        for (k = 2; k <= order; k++)
        {
            for (i = 0; i < order; i++)
            {
                tmp1 = rat_new(-(order / 2) + k - 1, 1, 0);
                matrix[(k + offset_smooth) * size + i + offset_pt] = rpol_monom_get(pol, i, l, &tmp1);
                tmp1 = rat_new(-(order / 2) + k - 1, 1, 0);
                tmp2 = rpol_monom_get(pol, i, l, &tmp1);
                tmp1 = rat_new(0, 1, 0);
                tmp2 = rat_sub(&tmp1, &tmp2);
                matrix[(k + offset_smooth) * size + i + offset_pt + order] = tmp2;
            }
            offset_pt += order;
        }
        for (i = 1; i <= order; i++)
        {
            tmp1 = rat_new(order - (order / 2), 1, 0);
            matrix[(offset_smooth + order + 1) * size + i + offset_pt - 1] = rpol_monom_get(pol, i - 1, l, &tmp1);
        }
        offset_smooth += order + 1;
        free(pol);
        tmp[l] = rat_new(0, 1, 0);
        for (k = l + 1; k < order;  k++)
        {
            tmp1 = rat_new(k - l, 1, 0);
            tmp[k] = rat_mul(&tmp[k], &tmp1);
        }
        pol = rpol_new(tmp, order - 1);
    }
    rpol_delete(pol);


    /* built of the vector */
    i = size - 1;
    do
    {
        vector[i] = rat_new(0, 1, 0);
    } while (--i >= 0);
    vector[0] = rat_new(1, 1, 0);

    s = system_new(matrix, vector, size);
    if (!s)
    {
        rpol_delete(pol);
        free(vector);
        free(matrix);
        return NULL;
    }

    system_solve(s);
    sol = system_solution_get(s);
    if (!sol)
    {
        system_delete(s);
        return NULL;
    }

    spline = (sw_spline_t *)malloc(sizeof(sw_spline_t));
    if (!spline)
    {
        system_delete(s);
        return NULL;
    }

    spline->order = order;
    spline->x1 = -(order >> 1);
    spline->x2 = order - (order >> 1);

    spline->rat_coefs = (rational_t *)malloc(sizeof(rational_t) * size);
    if (!spline->rat_coefs)
    {
        free(spline);
        system_delete(s);
        return NULL;
    }
    memcpy(spline->rat_coefs, sol, size * sizeof(rational_t));

    spline->coefs = (double *)malloc(sizeof(double) * size);
    if (!spline->coefs)
    {
        free(spline->rat_coefs);
        free(spline);
        system_delete(s);
        return NULL;
    }

    /* Normalisation */
    integral = rat_new(0, 1, 0);
    for (i = 0; i < order; i++)
    {
        rational_t tmp_int;

        pol = rpol_new(sol + i * order, order - 1);
        tmp_int = rpol_integrate(pol, i - (order >> 1), i + 1 - (order >> 1));
        free(pol);
        integral = rat_add(&integral, &tmp_int);
    }

    /* getting real values */

    for (i = 0; i < size; i++)
    {
        rational_t val;
        val = rat_div(&sol[i], &integral);
        spline->rat_coefs[i] = val;
        spline->coefs[i] = rat_double_get(&val);
    }

    spline->rat_coefs_primitive = _sw_spline_primitive_rat_get(spline);
    if (!spline->rat_coefs_primitive)
    {
        free(spline->coefs);
        free(spline->rat_coefs);
        free(spline);
        system_delete(s);
        return NULL;
    }

    /* getting real values */

    spline->coefs_primitive = (double *)malloc(spline->order * (spline->order + 1) * sizeof(double));
    if (!spline->coefs)
    {
        free(spline->rat_coefs_primitive);
        free(spline->coefs);
        free(spline->rat_coefs);
        free(spline);
        system_delete(s);
        return NULL;
    }

    for (i = 0; i < spline->order * (spline->order + 1); i++)
    {
        spline->coefs_primitive[i] = rat_double_get(&spline->rat_coefs_primitive[i]);
    }

    system_delete(s);

    return spline;
}

void
sw_spline_del(sw_spline_t *spline)
{
    if (!spline)
        return;

    free(spline->coefs_primitive);
    free(spline->rat_coefs_primitive);
    free(spline->coefs);
    free(spline->rat_coefs);
    free(spline);
}

int32_t
sw_spline_order_get(const sw_spline_t *spline)
{
    if (!spline)
    {
        fprintf(stderr, "[spl] spline undefined");
        return 0;
    }

    return spline->order;
}

int64_t
sw_spline_x1_get(const sw_spline_t *spline)
{
    if (!spline)
    {
        fprintf(stderr, "[spl] spline undefined");
        return 0;
    }

    return spline->x1;
}

int64_t
sw_spline_x2_get(const sw_spline_t *spline)
{
    if (!spline)
    {
        fprintf(stderr, "[spl] spline undefined");
        return 0;
    }

    return spline->x2;
}

const rational_t *
sw_spline_rat_coef_get(const sw_spline_t *spline)
{
    if (!spline)
        return NULL;

    return spline->rat_coefs;
}

const double *
sw_spline_coef_get(const sw_spline_t *spline)
{
    if (!spline)
        return NULL;

    return spline->coefs;
}

rational_t
sw_spline_value_rat_get(const sw_spline_t *spline,
                        const rational_t  *val)
{
    rational_t tmp;
    rational_t res;
    rational_t rx1;
    rational_t rx2;
    rpol_t    *pol;
    int64_t    k;

    rx1 = rat_new(spline->x1, 1, 0);
    rx2 = rat_new(spline->x2, 1, 0);
    if (rat_lesser_or_equal(val, &rx1) ||
        rat_greater_or_equal(val, &rx2))
        return rat_new(0, 1, 0);

    tmp = rat_sub(val, &rx1);
    k = (int64_t)floor(rat_double_get(&tmp));
    pol = rpol_new(spline->rat_coefs + k * spline->order, spline->order - 1);
    res = rpol_value_get(pol, val);
    free(pol);

    return res;
}

double
sw_spline_value_get(const sw_spline_t *spline,
                    double             val)
{
    const double *coefs;
    const double *tmp;
    double        res;
    int64_t       k;
    int32_t       i;

    if ((val <= (double)spline->x1) ||
        (val >= (double)spline->x2))
        return 0.0;

    i = spline->order - 1;
    if (i == 0)
        return 1.0;

    k = (int64_t)floor(val - (double)spline->x1);

    coefs = spline->coefs + k * spline->order;
    tmp = coefs + i;
    res = *tmp;
    do
    {
        res *= val;
        res += *--tmp;
    } while (--i > 0);

    return res;
}

double
sw_spline_integral_value_get(const sw_spline_t *spline,
                             double x1,
                             double x2)
{
    const double *coefs;
    const double *tmp;
    double        res;
    double        epsilon;
    double        val1;
    double        val2;
    int64_t       k1;
    int64_t       k2;
    int32_t       i;
    int32_t       k;

    res = 0.0;
    epsilon = 1.0;

    if (x2 < x1)
    {
        double t;

        t = x1;
        x1 = x2;
        x2 = t;
        epsilon = -1.0;
    }

    if (x1 < (double)spline->x1)
        x1 = (double)spline->x1;

    if (x2 > (double)spline->x2)
        x2 = (double)spline->x2;

    k1 = (int64_t)floor(x1);
    k2 = (int64_t)floor(x2);
    if (k2 == spline->x2)
        k2--;

    /* x1 and x2 are in the same interval */
    if (k1 == k2)
    {
        coefs = spline->coefs_primitive + (k1 - spline->x1) * (spline->order + 1);
        tmp = coefs + spline->order;
        val1 = *tmp;
        val2 = *tmp;

        i = spline->order;
        do
        {
            val1 *= (double)x1;
            val2 *= (double)x2;
            val1 += *--tmp;
            val2 += *tmp;
        } while (--i > 0);

        return epsilon * (val2 - val1);
    }

    /* x1 and x2 are not in the same interval */
    /* |---x1---|-------|---x2---|  */

    coefs = spline->coefs_primitive + (k1 - spline->x1) * (spline->order + 1);
    tmp = coefs + spline->order;
    val1 = *tmp;
    val2 = *tmp;
    i = spline->order;
    do
    {
        val1 *= (double)x1;
        val2 *= (double)(k1 + 1);
        val1 += *--tmp;
        val2 += *tmp;
    } while (--i > 0);
    res += val2 - val1;

    k1++;
    for (k = k1; k < k2; k++)
    {
        coefs += spline->order + 1;
        tmp = coefs + spline->order;
        val1 = *tmp;
        val2 = *tmp;
        i = spline->order;
        do
        {
            val1 *= (double)k;
            val2 *= (double)(k + 1);
            val1 += *--tmp;
            val2 += *tmp;
        } while (--i > 0);
        res += val2 - val1;
    }

    coefs += spline->order + 1;
    tmp = coefs + spline->order;
    val1 = *tmp;
    val2 = *tmp;
    i = spline->order;
    do
    {
        val1 *= (double)k2;
        val2 *= (double)x2;
        val1 += *--tmp;
        val2 += *tmp;
    } while (--i > 0);
    res += val2 - val1;

    return epsilon * res;
}

rational_t
sw_spline_integral_value_rat_get(const sw_spline_t *spline,
                                 int x1,
                                 int x2)
{
    const rational_t *coefs;
    const rational_t *tmp;
    rational_t        res;
    rational_t        epsilon;
    rational_t        val1;
    rational_t        val2;
    int32_t           i;
    int32_t           k;

    if (x1 == x2)
        return rat_new(0, 1, 0);

    res = rat_new(0, 1, 0);
    epsilon = rat_new(1, 1, 0);

    if (x2 < x1)
    {
        int t;

        t = x1;
        x1 = x2;
        x2 = t;
        epsilon = rat_new(-1, 1, 0);
    }

    if (x1 < spline->x1)
        x1 = spline->x1;

    if (x2 > spline->x2)
        x2 = spline->x2;

    coefs = spline->rat_coefs_primitive + (x1 - spline->x1) * (spline->order + 1);
    tmp = coefs + spline->order;
    for (k = x1; k < x2; k++)
    {
        val1 = *tmp;
        val2 = *tmp;
        i = spline->order;
        do
        {
            rational_t rk;

            rk = rat_new(k, 1, 0);
            val1 = rat_mul(&val1, &rk);
            rk = rat_new(k + 1, 1, 0);
            val2 = rat_mul(&val2, &rk);
            val1 = rat_add(&val1, --tmp);
            val2 = rat_add(&val2, tmp);
        } while (--i > 0);
        val2 = rat_sub(&val2, &val1);
        res = rat_add(&res, &val2);

        coefs += spline->order + 1;
        tmp = coefs + spline->order;
    }

    res = rat_mul(&res, &epsilon);

    return res;
}
