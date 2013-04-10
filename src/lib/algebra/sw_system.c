#include <stdlib.h>
#include <stdio.h>

#include "sw_system.h"


/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

/**
 * @cond SW_LOCAL
 */

struct system
{
    rational_t *matrix;
    rational_t *vector;
    rational_t *solution;
    int32_t     order;
};

static int32_t
system_pivot_get(const system_t *s,
                 int32_t         k)
{
    rational_t *iter;
    rational_t  max;
    int         l;
    int         i;

    l = k;
    iter = s->matrix + k * (s->order + 1);
    max = rat_abs(iter);
    for (i = k + 1, iter += s->order; i < s->order; i++, iter += s->order)
    {
        if (rat_greater(iter, &max))
        {
            l = i;
            max = rat_abs(iter);
        }
    }

    return l;
}

static void
system_lines_swap(const system_t *s,
                  int32_t         k,
                  int32_t         l)
{
    rational_t  tmp;
    rational_t *iter1;
    rational_t *iter2;
    int32_t     j;

    iter1 = s->matrix + k * (s->order + 1);
    iter2 = s->matrix + l * s->order + k;
    for (j = k; j < s->order; j++, iter1++, iter2++)
    {
        tmp = *iter1;
        *iter1 = *iter2;
        *iter2 = tmp;
    }
    tmp = s->vector[k];
    s->vector[k] = s->vector[l];
    s->vector[l] = tmp;
}

static void
system_pivote(const system_t *s,
              int32_t         k)
{
    rational_t pivot;
    rational_t val;
    int32_t    i;
    int32_t    j;

    for (i = k + 1; i < s->order; i++)
    {
        pivot = rat_div(&s->matrix[i * s->order + k], &s->matrix[k * s->order + k]);
        val = rat_mul(&pivot, &s->vector[k]);

        s->vector[i] = rat_sub(&s->vector[i], &val);
        for (j = k + 1; j < s->order; j++)
        {
            val = rat_mul(&pivot, &s->matrix[k * s->order + j]);
            s->matrix[i * s->order + j] = rat_sub(&s->matrix[i * s->order + j], &val);
        }
    }
}

static void
system_triangularize(const system_t *s)
{
    int32_t k;
    int32_t l;

    for (k = 0; k < s->order; k++)
    {
        l = system_pivot_get(s, k);
        if (l > k)
        {
            system_lines_swap(s, k, l);
        }
        system_pivote(s, k);
    }
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

system_t *
system_new(rational_t *matrix,
           rational_t *vector,
           int32_t     order)
{
    system_t *s;

    if (!matrix || !vector || (order <= 0))
        return NULL;

    s = (system_t *)malloc(sizeof(system_t));
    if (!s)
        return NULL;

    s->matrix = matrix;
    s->vector = vector;
    s->solution = NULL;
    s->order = order;

    return s;
}

void
system_delete(system_t *s)
{
    if (!s)
        return;

    free(s->matrix);
    free(s->vector);
    if (s->solution)
        free(s->solution);
    free(s);
}

void
system_solve(system_t *s)
{
    rational_t *iter_v;
    rational_t  tmp;
    int32_t     i;
    int32_t     j;

    s->solution = (rational_t *)malloc(sizeof(rational_t) * s->order);
    if (!s->solution)
        return;

    system_triangularize(s);

    iter_v = s->vector + s->order - 1;
    for (i = s->order - 1; i >= 0; i--, iter_v--)
    {
        tmp = *iter_v;
        j = i + 1;
        for (j = i + 1; j < s->order; j++)
        {
            rational_t m;

            m = rat_mul(&(s->matrix[i * s->order + j]), &s->solution[j]);
            tmp = rat_sub(&tmp, &m);
        }
        s->solution[i] = rat_div(&tmp, &(s->matrix[i * (s->order + 1)]));
    }
}

const rational_t *
system_solution_get(const system_t *s)
{
    if (!s || !s->solution)
        return NULL;

    return s->solution;
}
