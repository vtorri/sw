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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "sw_polynomial.h"


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


rpol_t *
rpol_new(const rational_t *coefs,
         int32_t           degree)
{
    rpol_t *pol;

    if (!coefs)
        return NULL;

    if (degree < 0)
        degree = -1;

    pol = (rpol_t *)malloc(sizeof(rpol_t));
    if (!pol)
        return NULL;

    pol->coefs = coefs;
    pol->degree = degree;

    return pol;
}

dpol_t *
dpol_new(const double *coefs,
         int32_t       degree)
{
    dpol_t *pol;

    if (!coefs)
        return NULL;

    if (degree < 0)
        degree = -1;

    pol = (dpol_t *)malloc(sizeof(dpol_t));
    if (!pol)
        return NULL;

    pol->coefs = coefs;
    pol->degree = degree;

    return pol;
}

void
rpol_delete(rpol_t *pol)
{
    if (!pol)
        return;

    free((void *)pol->coefs);
    free(pol);
}

int32_t
rpol_degree_get(const rpol_t *pol)
{
    if (!pol)
    {
        fprintf(stderr, "[pol] polynomial undefined");
        return -1;
    }

    return pol->degree;
}

rational_t
rpol_coef_get(rpol_t *pol,
              int32_t n)
{
    if (!pol)
    {
        fprintf(stderr, "[pol] polynomial undefined");
        return rat_new(0, 1, 0);
    }

    if ((n < 0) ||
        (n > pol->degree))
    {
        fprintf(stderr, "[pol] wrong index");
        return rat_new(0, 1, 0);
    }

    return pol->coefs[n];
}

rational_t
rpol_monom_get(const rpol_t     *pol,
               int32_t           n,
               int32_t           smooth,
               const rational_t *val)
{
    rational_t res;
    int        j;

    if (!pol)
    {
        fprintf(stderr, "[pol] polynomial undefined");
        return rat_new(0, 1, 0);
    }

    if ((n < 0) ||
        (n > pol->degree))
    {
        fprintf(stderr, "[pol] wrong index");
        return rat_new(0, 1, 0);
    }

    res = pol->coefs[n];

    j = n - smooth - 1;
    if (j < 0)
        return res;
    do
    {
        res  =rat_mul(&res, val);
    } while (--j >= 0);

    return res;
}

rational_t
rpol_value_get(const rpol_t     *pol,
               const rational_t *x)
{
    rational_t        val;
    const rational_t *tmp;
    int32_t           i;

    if (!pol)
    {
        fprintf(stderr, "[pol] polynomial undefined");
        return rat_new(0, 1, 0);
    }

    if (!x)
    {
        fprintf(stderr, "[pol] value undefined");
        return rat_new(0, 1, 0);
    }

    i = pol->degree;
    if (i == 0)
        return rat_new(1, 1, 0);
    tmp = pol->coefs + pol->degree;
    val = *tmp;
    do
    {
        val = rat_mul(&val, x);
        val = rat_add(&val, --tmp);
    } while (--i > 0);

    return val;
}

double
dpol_value_get(const dpol_t *pol,
               double        x)
{
    double        val;
    const double *tmp;
    int32_t       i;

    if (!pol)
    {
        fprintf(stderr, "[pol] polynomial undefined");
        return 0.0;
    }

    i = pol->degree;
    if (i == 0)
        return 1.0;
    tmp = pol->coefs + pol->degree;
    val = *tmp;
    do
    {
        val *= x;
        val += *--tmp;
    } while (--i > 0);

    return val;
}

rational_t
rpol_integrate(const rpol_t *pol, int32_t a, int32_t b)
{
    const rational_t *iter;
    rational_t        res_a;
    rational_t        res_b;
    rational_t        ra;
    rational_t        rb;
    int32_t           i;
    int32_t           j;

    ra = rat_new(a, 1, 0);
    rb = rat_new(b, 1, 0);
    res_a = rat_new(0, 1, 0);
    res_b = rat_new(0, 1, 0);

    iter = pol->coefs;
    for (i = 0; i <= pol->degree; i++, iter++)
    {
        rational_t tmp_a;
        rational_t tmp_b;

        tmp_a = rat_new(i + 1, 1, 0);
        tmp_a = rat_div(iter, &tmp_a);
        tmp_b = tmp_a;

        for (j = 0; j <= i; j++)
        {
            tmp_a = rat_mul(&tmp_a, &ra);
            tmp_b = rat_mul(&tmp_b, &rb);
        }
        res_a = rat_add(&res_a, &tmp_a);
        res_b = rat_add(&res_b, &tmp_b);
    }

    return rat_sub(&res_b, &res_a);
}

static int64_t *
rpol_int_mul(int64_t *vals, int32_t num, int64_t coef)
{
    int64_t *coefs;
    int32_t  i;

    coefs = (int64_t *)malloc(sizeof(int64_t) * (num + 1));
    if (!coefs)
        return NULL;

    coefs[0] = - coef * vals[0];
    for (i = 1; i < num; i++)
    {
        coefs[i] = vals[i - 1] - coef * vals[i];
    }
    coefs[num] = vals[num - 1];

    return coefs;
}

rpol_t *
rpol_lagrange_base_get(int64_t *x, int32_t num, int32_t i)
{
    rpol_t     *pol;
    rational_t *coefs;
    int64_t    *nums;
    int64_t    *n;
    int64_t     den;
    int32_t     size;
    int32_t     j;

    nums = (int64_t *)malloc(sizeof(int64_t));
    if (!nums)
        return NULL;

    coefs = (rational_t *)malloc(sizeof(rational_t) * num);
    if (!coefs)
    {
        free(nums);
        return NULL;
    }

    *nums = 1;
    den = 1;
    size = 1;
    for (j = 0; j < num; j++)
    {
        if (i != j)
        {
            n = rpol_int_mul(nums, size, x[j]);
            free(nums);
            if (!n)
            {
                free(coefs);
                return NULL;
            }
            nums = n;
            den *= x[i] - x[j];
            size++;
        }
    }

    for (j = 0; j < num; j++)
        coefs[j] = rat_new(nums[j], den, 1);
    free(nums);

    pol = rpol_new(coefs, num - 1);
    if (!pol)
    {
        free(coefs);
        return NULL;
    }

    return pol;
}

/* get T_n, n = 0,...,r - 1 */
rpol_t **
sw_rpol_chebychev_pols_get(int32_t r)
{
    rpol_t    **res;
    rational_t *coefs;
    int32_t     i;

    res = (rpol_t **)malloc((r + 1) * sizeof(rpol_t *));
    if (!res)
        return NULL;

    coefs = (rational_t *)malloc(sizeof(rational_t));
    coefs[0] = rat_new(1, 1, 0);
    res[0] = rpol_new(coefs, 0);

    coefs = (rational_t *)malloc(2 * sizeof(rational_t));
    coefs[0] = rat_new(0, 1, 0);
    coefs[1] = rat_new(1, 1, 0);
    res[1] = rpol_new(coefs, 1);

    for (i = 2; i <= r; i++)
    {
        int32_t j;

        coefs = (rational_t *)malloc((i + 1) * sizeof(rational_t));
        coefs[0] = rat_new(0,1,0);
        for (j = 1; j <= i; j++)
        {
            rational_t deux;
            rational_t rat;

            deux = rat_new(2, 1, 0);
            rat = rpol_coef_get(res[i - 1], j - 1);
            coefs[j] = rat_mul(&deux, &rat);
        }
        for (j = 0; j <= i - 2; j++)
        {
            rational_t rat;

            rat = rpol_coef_get(res[i - 2], j);
            coefs[j] = rat_sub(coefs + j, &rat);
        }
        res[i] = rpol_new(coefs, i);
    }

    return res;
}

void
rpol_disp(const rpol_t *pol, uint8_t end)
{
    const rational_t *iter;
    rational_t        tmp;
    int32_t           n;
    int32_t           first = 1;

    if (pol->degree < 0)
    {
        printf("0");
        if (end)
            printf("\n");
        return;
    }

    n = 0;
    iter = pol->coefs;
    do
    {
        if (rat_is_equal(iter, 0))
        {
            iter++;
            continue;
        }

        switch (n)
        {
         case 0:
             rat_disp(iter, 0);
             break;
         case 1:
             tmp = rat_new(0, 1, 0);
             if (rat_greater(iter, &tmp))
             {
                 if (!first)
                     printf(" + ");
             }
             else
                 printf(" - ");
             tmp = rat_abs(iter);
             if (rat_is_equal(&tmp, 1))
                 printf("x");
             else
             {
                 rat_disp(&tmp, 0);
                 printf("x");
             }
             break;
         default:
             tmp = rat_new(0, 1, 0);
             if (rat_greater(iter, &tmp))
             {
                 if (!first)
                     printf(" + ");
             }
             else
                 printf(" - ");
             tmp = rat_abs(iter);
             if (rat_is_equal(&tmp, 1))
                 printf("x^%d", (int)n);
             else
             {
                 rat_disp(&tmp, 0);
                 printf("x^%d", (int)n);
             }
             break;
        }
        iter++;
        first = 0;
    } while (++n <= pol->degree);

    if (end)
        printf("\n");
}

void
dpol_disp(const dpol_t *pol, uint8_t end)
{
    const double *iter;
    double        tmp;
    int32_t       n;

    if (pol->degree < 0)
    {
        printf("0");
        if (end)
            printf("\n");
        return;
    }

    n = 0;
    iter = pol->coefs;
    do
    {
        if (*iter == 0.0)
        {
            iter++;
            continue;
        }

        switch (n)
        {
         case 0:
             printf("%f", *iter);
             break;
         case 1:
             tmp = fabs(*iter);
             if (tmp == 1.0)
                 printf("x");
             else
                 printf("%fx", tmp);
             break;
         default:
             tmp = fabs(*iter);
             if (tmp == 1.0)
                 printf("x^%d", (int)n);
             else
                 printf("%fx^%d", tmp, (int)n);
             break;
        }
        if (n != pol->degree)
        {
            if (*(iter + 1) > 0.0)
                printf(" + ");
            else
                printf(" - ");
        }
        iter++;
    } while (++n <= pol->degree);

    if (end)
        printf("\n");
}
