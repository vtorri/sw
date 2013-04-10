#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "sw.h"


/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

/**
 * @cond SW_LOCAL
 */


/* inherit from scale_fct_base_t */

struct scale_fct
{
    scale_fct_base_t scale_fct_base;

    int32_t          order;
    weights_type_t   type;
    spline_t        *spline;

    void (*proj_periodic_backward)(const scale_fct_t *sf,
                                   int32_t            scale,
                                   const double      *ps,
                                   double            *f2);

    void (*proj_dirichlet_backward)(const scale_fct_t *sf,
                                    int32_t            scale,
                                    const double      *ps,
                                    double            *f2);
};


/* inherit from scale_fct_base_t */

struct scale_fct_dual
{
    scale_fct_base_t scale_fct_base;

    int32_t          order;
    int32_t          order_dual;
    weights_type_t   type;

    weights_t       *weights;

    void (*proj_periodic_forward)(const scale_fct_dual_t *sfd,
                                  int32_t                 scale,
                                  const double           *function,
                                  double                 *ps);

    void (*proj_dirichlet_forward)(const scale_fct_dual_t *sfd,
                                   int32_t                 scale,
                                   const double           *function,
                                   double                 *ps);
};

/* static rational_t * */
/* filter_values_get(const filter_t *f, */
/*                    int32_t         scale) */
/*
  { */
/*   rational_t  tmp; */
/*   rational_t *coefs; */
/*   rational_t *val1; */
/*   rational_t *val2; */
/*   int32_t     period; */
/*   int32_t     j; */
/*   int32_t     k; */
/*   int32_t     n; */

/*   period = f->N2 - f->N1; */
/*   coefs = (rational_t *)malloc(sizeof(rational_t) * (f->N2 - f->N1 + 1)); */
/*   if (!coefs) */
/*     return NULL; */

/*   /\* the coefs must be doubled *\/ */
/*   tmp = rat_new(2, 1, 0); */
/*   for (k = 0; k < f->N2 - f->N1 + 1; k++) */
/*     coefs[k] = rat_mul(&f->filter_rat[k], &tmp); */

/*   val1 = (rational_t *)malloc(sizeof(rational_t) * (period << scale)); */
/*   if (!val1)
     { */
/*     free(coefs); */
/*     return NULL; */
/*   } */

/*   for (k = 0; k < (period << scale); k++)
     { */
/*     val1[k] = rat_new(0, 1, 0); */
/*   } */

/*   val2 = (rational_t *)malloc(sizeof(rational_t) * (period << scale)); */
/*   if (!val2)
     { */
/*     free(val1); */
/*     free(coefs); */
/*     return NULL; */
/*   } */

/*   /\* s_0 *\/ */
/*   for (k = 0; k < period; k++)
     { */
/*     int kk; */

/*     kk = k << scale; */
/*     if (k == -f->N1) */
/*       val1[kk] = rat_new(1, 1, 0); */
/*   } */

/*   /\* s_j, j  = 1 ... scale *\/ */
/*   for (j = 1; j <= scale; j++)
     { */
/*     for (k = 0; k < (period << j); k++)
       { */
/*       int kk; */

/*       kk = k << (scale - j); */
/*       val2[kk] = rat_new(0, 1, 0); */
/*       for (n = 0; n < (period << (j - 1)); n++)
         { */
/*         if (((k - 2 * n) >= f->N1) && */
/*             ((k - 2 * n) <= f->N2))
            { */
/*           tmp = rat_mul(coefs + k - 2 * n - f->N1, val1 + (n << (scale - j + 1))); */
/*           val2[kk] = rat_add(val2 + kk, &tmp); */
/*         } */
/*       } */
/*     } */
/*     memcpy(val1, val2, sizeof(rational_t) * (period << scale)); */
/*   } */

/*   free(val1); */
/*   free(coefs); */

/*   return val2; */
/* } */


/*
 * scale function - base class methods
 */

static scale_fct_base_t
scale_fct_base_new(int32_t     N1,
                   int32_t     N2,
                   rational_t *filter)
{
    scale_fct_base_t  sf;
    int32_t           i;

    memset(&sf, 0, sizeof(scale_fct_base_t));

    if (!filter)
        return sf;

    sf.filter = (double *)malloc(sizeof(double) * (N2 - N1 + 1));
    if (!sf.filter)
        return sf;

    for (i = N1; i <= N2; i++)
    {
        sf.filter[i - N1] = rat_double_get(filter + i - N1);
    }

    sf.N1 = N1;
    sf.N2 = N2;
    sf.filter_rat = filter;

    return sf;
}

static void
scale_fct_base_del(scale_fct_base_t *sf)
{
    if (!sf)
        return;

    if (sf->filter_rat)
    {
        free(sf->filter_rat);
        sf->filter_rat = NULL;
    }

    if (sf->filter)
    {
        free(sf->filter);
        sf->filter = NULL;
    }
}


/*
 * scale function methods
 */

static rational_t *
scale_fct_filter_get(int32_t  order)
{
    rational_t *filter;
    int64_t     puiss;
    int32_t     i;

    filter = (rational_t *)malloc(sizeof(rational_t) * (order + 1));
    if (!filter)
        return NULL;

    puiss = 1 << order;

    for (i = 0; i <= order; i++)
        filter[i] = rat_new(sw_binomial(order, i), puiss, 1);

    return filter;
}

static void
scale_fct_inner_product_lagrange_periodic_backward(const scale_fct_t *sf, int32_t scale, const double *ps, double *f2)
{
    double  val;
    int32_t i;
    int32_t k;

/*   for (i = 0; i < (1 << scale); i++)
     { */
/*     val = 0.0; */
/*     for (k = 0; k < (1 << scale); k++)
       { */
/*       val += ps[k] * scale_fct_value_get(sf, i - k); */
/*     } */
/*     f2[i] = val; */
/*   } */

    /*
     * optimisation of the code above : loop on k such that
     * i-k in the support of phi
     * with scale = 5, 7 millions less of call of scale_fct_value_get...
     */

    for (i = 0; i < (1 << scale); i++)
    {
        val = 0.0;
        for (k = i - ((const scale_fct_base_t *)sf)->N2 + 1; k < i - ((const scale_fct_base_t *)sf)->N1; k++)
        {
            if ((k >= 0) && (k < (1 << scale)))
            {
                val += ps[k] * scale_fct_value_get(sf, i - k);
            }
        }
        f2[i] = val;
    }
}

static void
scale_fct_inner_product_lagrange_dirichlet_backward(const scale_fct_t *sf, int32_t scale, const double *ps, double *f2)
{
    double  val;
    int32_t i;
    int32_t k;

/*   for (i = -(1 << scale); i < (1 << scale); i++)
     { */
/*     val = 0.0; */
/*     for (k = -(1 << scale); k < (1 << scale); k++)
       { */
/*       val += ps[k + (1 << scale)] * scale_fct_value_get(sf, i - k); */
/*     } */
/*     f2[i + (1 << scale)] = val; */
/*   } */

    /*
     * optimisation of the code above : loop on k such that
     * i-k in the support of phi
     * see scale_fct_inner_product_lagrange_periodic_backward()
     */

    for (i = -(1 << scale); i < (1 << scale); i++)
    {
        val = 0.0;
        for (k = i - ((const scale_fct_base_t *)sf)->N2 + 1; k < i - ((const scale_fct_base_t *)sf)->N1; k++)
        {
            if ((k >= -(1 << scale)) && (k < (1 << scale)))
            {
                val += ps[k + (1 << scale)] * scale_fct_value_get(sf, i - k);
            }
        }
        f2[i + (1 << scale)] = val;
    }
}


/*
 * dual scale function methods
 */

static rational_t *
scale_fct_dual_filter_get(int32_t  order,
                          int32_t  order_dual)
{
    rational_t *filter;
    int64_t     puiss;
    int32_t     order_total;
    int32_t     i;
    int32_t     j;
    int32_t     l;
    int32_t     n;

    filter = (rational_t *)malloc(sizeof(rational_t) * (order + (order_dual << 1) - 1));
    if (!filter)
        return NULL;

    puiss = 1 << order_dual;
    order_total = order + order_dual;

    for (i = 0; i < (order + (order_dual << 1) - 1); i++)
        filter[i] = rat_new(0, 1, 0);

    for (j = -(order_dual >> 1); j <= (order_dual - (order_dual >> 1)); j++)
    {
        rational_t coef1;

        coef1 = rat_new(sw_binomial(order_dual, j + (order_dual >> 1)), puiss, 1);
        for (n = 0; n < (order_total >> 1); n++)
        {
            rational_t coef2;

            coef2 = rat_new(sw_binomial((order_total >> 1) - 1 + n, n), (1 << (n << 1)), 1);
            coef2 = rat_mul(&coef2, &coef1);

            for (l = -n; l <= n; l++)
            {
                rational_t coef3;

                i = j + l + (order_total >> 1) + (order_dual >> 1) - 1;
                if ((l & 1) == 1)
                    coef3 = rat_new(-sw_binomial(n << 1, n + l), 1, 0);
                else
                    coef3 = rat_new(sw_binomial(n << 1, n + l), 1, 0);

                coef3 = rat_mul(&coef2, &coef3);
                filter[i] = rat_add(&filter[i], &coef3);
            }
        }
    }

    return filter;
}

static int
mod(int64_t v, int64_t b)
{
    if (v >= 0) return (v % b);
    else
    {
        while (v < 0) v +=b;
        return v;
    }
}

static void
scale_fct_dual_inner_product_lagrange_periodic_forward(const scale_fct_dual_t *sfd, int32_t scale, const double *function, double *ps)
{
    const double *w;
    int64_t       k;
    int32_t       degree;

    w = weights_get(sfd->weights, &degree);
    degree--;

    for (k = 0; k < (1 << scale); k++)
    {
        double  val = 0.0;
        int32_t n;

        for (n = -(degree >> 1); n <= degree - (degree >> 1); n++)
        {
            int32_t j;

            /* j = (k + n)  & ((1 << scale) - 1); */
            j = mod(k + n, 1 << scale);
            val += function[j] * w[n + (degree >> 1)];
        }
        ps[k] = val;
    }
}

static void
scale_fct_dual_inner_product_sweldens_periodic_forward(const scale_fct_dual_t *sfd, int32_t scale, const double *function, double *ps)
{
#if 0
    const double *w;
    int64_t       k;
    int32_t       size;

    w = weights_get(sfd->weights, &size);

    for (k = 0; k < (1 << scale); k++)
    {
        double  val = 0.0;
        int32_t n;

        for (n = -(degree >> 1); n <= degree - (degree >> 1); n++)
        {
            int32_t j;

            /* j = (k + n)  & ((1 << scale) - 1); */
            j = mod(k + n, 1 << scale);
            val += function[j] * w[n + (degree >> 1)];
        }
        ps[k] = val;
    }
#endif
}

static void
scale_fct_dual_inner_product_lagrange_dirichlet_forward(const scale_fct_dual_t *sfd, int32_t scale, const double *function, double *ps)
{
    const double *w;
    int64_t       k;
    int32_t       size;
    int32_t       degree;

    w = weights_get(sfd->weights, &degree);
    degree--;

    size = 1 << (scale + 1);
    for (k = 0; k < size; k++)
    {
        double  val = 0.0;
        int32_t n;

        for (n = -(degree >> 1); n <= degree - (degree >> 1); n++)
        {
            int32_t j;

            j = k + n;
            if ((j >= 0) && (j < size))
                val += function[j] * w[n + (degree >> 1)];
        }
        ps[k] = val;
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

/*
 * scale function methods
 */

scale_fct_t *
scale_fct_new(int32_t order)
{
    scale_fct_t *sf;

    if (order <= 0)
        return NULL;

    sf = (scale_fct_t *)malloc(sizeof(scale_fct_t));
    if (!sf)
        return NULL;

    sf->scale_fct_base = scale_fct_base_new(-(order >> 1),
                                            order - (order >> 1),
                                            scale_fct_filter_get (order));
    sf->order = order;
    sf->spline = spline_new(order);
    if (!sf->spline)
    {
        scale_fct_base_del(&sf->scale_fct_base);
        free(sf);
        return NULL;
    }

    return sf;
}

void
scale_fct_del(scale_fct_t *sf)
{
    if (!sf)
        return;

    scale_fct_base_del(&sf->scale_fct_base);
    spline_del(sf->spline);
    free(sf);
}

uint8_t
scale_fct_type_set(scale_fct_t   *sf,
                   weights_type_t type)
{
    if (!sf)
        return 0;

    switch (type)
    {
     case SW_WEIGHTS_TYPE_LAGRANGE:
     case SW_WEIGHTS_TYPE_SWELDENS:
         sf->proj_periodic_backward = scale_fct_inner_product_lagrange_periodic_backward;
         sf->proj_dirichlet_backward = scale_fct_inner_product_lagrange_dirichlet_backward;
         break;
     default:
         return 0;
    }

    sf->type = type;

    return 1;
}

weights_type_t
scale_fct_type_get(const scale_fct_t *sf)
{
    if (!sf)
        return SW_WEIGHTS_TYPE_ERROR;

    return sf->type;
}

const spline_t *
scale_fct_spline_get(const scale_fct_t *sf)
{
    if (!sf)
        return NULL;

    return sf->spline;
}

rational_t
scale_fct_value_rat_get(const scale_fct_t *sf,
                        const rational_t  *val)
{
    if (!sf)
        return rat_new(0, 1, 0);

    return spline_value_rat_get(sf->spline, val);
}

double
scale_fct_value_get(const scale_fct_t *sf,
                    double             val)
{
    return spline_value_get(sf->spline, val);
}

void
scale_fct_proj_periodic_backward(const scale_fct_t *sf,
                                 int32_t            scale,
                                 const double      *ps,
                                 double            *f2)
{
    if (sf)
        sf->proj_periodic_backward(sf, scale, ps, f2);
}

void
scale_fct_proj_dirichlet_backward(const scale_fct_t *sf,
                                  int32_t            scale,
                                  const double      *ps,
                                  double            *f2)
{
    if (sf)
        sf->proj_dirichlet_backward(sf, scale, ps, f2);
}

/* dual scale function methods */

scale_fct_dual_t *
scale_fct_dual_new(int32_t order,
                   int32_t order_dual)
{
    scale_fct_dual_t *sf;

    if ((order <= 0) || (order_dual <= 0))
        return NULL;

    if (((order + order_dual) & 1) == 1)
        return NULL;

    sf = (scale_fct_dual_t *)malloc(sizeof(scale_fct_dual_t));
    if (!sf)
        return NULL;

    sf->scale_fct_base = scale_fct_base_new(-(order_dual >> 1) - ((order + order_dual) >> 1) + 1,
                                            -(order_dual >> 1) + ((order + order_dual) >> 1) + order_dual - 1,
                                            scale_fct_dual_filter_get(order, order_dual));
    sf->order = order;
    sf->order_dual = order_dual;

    return sf;
}

void
scale_fct_dual_del(scale_fct_dual_t *sf)
{
    if (!sf)
        return;

    if (sf->weights) weights_del(sf->weights);
    scale_fct_base_del(&sf->scale_fct_base);
    free(sf);
}

uint8_t
scale_fct_dual_type_set(scale_fct_dual_t *sf,
                        weights_type_t    type,
                        ...)
{
    va_list va;

    if (!sf)
        return 0;

    va_start(va, type);

    switch (type)
    {
     case SW_WEIGHTS_TYPE_LAGRANGE:
     {
         int32_t degree;

         degree = va_arg(va, int32_t);
         sf->weights = sw_weights_lagrange_new(&sf->scale_fct_base, degree);
         if (!sf->weights)
         {
             va_end(va);
             return 0;
         }
         sf->proj_periodic_forward = scale_fct_dual_inner_product_lagrange_periodic_forward;
         sf->proj_dirichlet_forward = scale_fct_dual_inner_product_lagrange_dirichlet_forward;
         break;
     }
     case SW_WEIGHTS_TYPE_SWELDENS:
     {
         int32_t r;
         int32_t s;

         r = va_arg(va, int32_t);
         s = va_arg(va, int32_t);
         sf->weights = sw_weights_sweldens_new(&sf->scale_fct_base, r, s);
         if (!sf->weights)
         {
             va_end(va);
             return 0;
         }
         break;
     }
     default:
         va_end(va);
         return 0;
    }

    va_end(va);

    sf->type = type;

    return 1;
}

weights_type_t
scale_fct_dual_type_get(const scale_fct_dual_t *sf)
{
    if (!sf)
        return SW_WEIGHTS_TYPE_ERROR;

    return sf->type;
}

void
scale_fct_dual_proj_periodic_forward(const scale_fct_dual_t *sfd,
                                     int32_t                 scale,
                                     const double           *function,
                                     double                 *ps)
{
    if (sfd)
        sfd->proj_periodic_forward(sfd, scale, function, ps);
}

void
scale_fct_dual_proj_dirichlet_forward(const scale_fct_dual_t *sfd,
                                      int32_t                 scale,
                                      const double           *function,
                                      double                 *ps)
{
    if (sfd)
        sfd->proj_dirichlet_forward(sfd, scale, function, ps);
}

const double *
scale_fct_dual_weights_get(const scale_fct_dual_t *sfd,
                           int32_t                *size)
{
    return weights_get(sfd->weights, size);
}
