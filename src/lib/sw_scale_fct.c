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


/* inherit from sw_scale_fct_base_t */

struct sw_scale_fct_s
{
    sw_scale_fct_base_t scale_fct_base;

    int32_t             order;
    sw_weights_type_t   type;
    sw_spline_t        *spline;

    void (*proj_periodic_backward)(const sw_scale_fct_t *sf,
                                   int32_t            scale,
                                   const double      *ps,
                                   double            *f2);

    void (*proj_dirichlet_backward)(const sw_scale_fct_t *sf,
                                    int32_t            scale,
                                    const double      *ps,
                                    double            *f2);
};


/* inherit from sw_scale_fct_base_t */

struct sw_scale_fct_dual_s
{
    sw_scale_fct_base_t scale_fct_base;

    int32_t             order;
    int32_t             order_dual;
    sw_weights_type_t   type;

    sw_weights_t          *weights;

    void (*proj_periodic_forward)(const sw_scale_fct_dual_t *sfd,
                                  int32_t                 scale,
                                  const double           *function,
                                  double                 *ps);

    void (*proj_dirichlet_forward)(const sw_scale_fct_dual_t *sfd,
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

static sw_scale_fct_base_t
sw_scale_fct_base_new(int32_t     N1,
		      int32_t     N2,
		      rational_t *filter)
{
    sw_scale_fct_base_t  sf;
    int32_t              i;

    memset(&sf, 0, sizeof(sw_scale_fct_base_t));

    if (!filter)
        return sf;

    sf.filter = (double *)calloc(1, sizeof(double) * (N2 - N1 + 1));
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
sw_scale_fct_base_del(sw_scale_fct_base_t *sf)
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
sw_scale_fct_filter_get(int32_t  order)
{
    rational_t *filter;
    int64_t     puiss;
    int32_t     i;

    filter = (rational_t *)calloc(1, sizeof(rational_t) * (order + 1));
    if (!filter)
        return NULL;

    puiss = 1 << order;

    for (i = 0; i <= order; i++)
        filter[i] = rat_new(sw_binomial(order, i), puiss, 1);

    return filter;
}

static void
sw_scale_fct_inner_product_periodic_backward(const sw_scale_fct_t *sf, int32_t scale, const double *ps, double *f2)
{
    double  val;
    int32_t i;
    int32_t k;

/*   for (i = 0; i < (1 << scale); i++)
     { */
/*     val = 0.0; */
/*     for (k = 0; k < (1 << scale); k++)
       { */
/*       val += ps[k] * sw_scale_fct_value_get(sf, i - k); */
/*     } */
/*     f2[i] = val; */
/*   } */

    /*
     * optimisation of the code above : loop on k such that
     * i-k in the support of phi
     * with scale = 5, 7 millions less of call of sw_scale_fct_value_get...
     */

    for (i = 0; i < (1 << scale); i++)
    {
        val = 0.0;
        for (k = i - ((const sw_scale_fct_base_t *)sf)->N2 + 1; k < i - ((const sw_scale_fct_base_t *)sf)->N1; k++)
        {
            if ((k >= 0) && (k < (1 << scale)))
            {
                val += ps[k] * sw_scale_fct_value_get(sf, i - k);
            }
        }
        f2[i] = val;
    }
}

static void
sw_scale_fct_inner_product_dirichlet_backward(const sw_scale_fct_t *sf, int32_t scale, const double *ps, double *f2)
{
    double  val;
    int32_t i;
    int32_t k;

/*   for (i = -(1 << scale); i < (1 << scale); i++)
     { */
/*     val = 0.0; */
/*     for (k = -(1 << scale); k < (1 << scale); k++)
       { */
/*       val += ps[k + (1 << scale)] * sw_scale_fct_value_get(sf, i - k); */
/*     } */
/*     f2[i + (1 << scale)] = val; */
/*   } */

    /*
     * optimisation of the code above : loop on k such that
     * i-k in the support of phi
     * see sw_scale_fct_inner_product_lagrange_periodic_backward()
     */

    for (i = -(1 << scale); i < (1 << scale); i++)
    {
        val = 0.0;
        for (k = i - ((const sw_scale_fct_base_t *)sf)->N2 + 1; k < i - ((const sw_scale_fct_base_t *)sf)->N1; k++)
        {
            if ((k >= -(1 << scale)) && (k < (1 << scale)))
            {
                val += ps[k + (1 << scale)] * sw_scale_fct_value_get(sf, i - k);
            }
        }
        f2[i + (1 << scale)] = val;
    }
}


/*
 * dual scale function methods
 */

static rational_t *
sw_scale_fct_dual_filter_get(int32_t order,
			     int32_t order_dual)
{
    rational_t *filter;
    int64_t     puiss;
    int32_t     order_total;
    int32_t     i;
    int32_t     j;
    int32_t     l;
    int32_t     n;

    filter = (rational_t *)calloc(1, sizeof(rational_t) * (order + (order_dual << 1) - 1));
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
sw_scale_fct_dual_inner_product_lagrange_periodic_forward(const sw_scale_fct_dual_t *sfd, int32_t scale, const double *function, double *ps)
{
    const double *w;
    int64_t       k;
    int32_t       degree;

    w = sw_weights_get(sfd->weights, &degree);
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
sw_scale_fct_dual_inner_product_sweldens_periodic_forward(const sw_scale_fct_dual_t *sfd, int32_t scale, const double *function, double *ps)
{
#if 0
    const double *w;
    int64_t       k;
    int32_t       size;

    w = sw_weights_get(sfd->weights, &size);

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
sw_scale_fct_dual_inner_product_lagrange_dirichlet_forward(const sw_scale_fct_dual_t *sfd, int32_t scale, const double *function, double *ps)
{
    const double *w;
    int64_t       k;
    int32_t       size;
    int32_t       degree;

    w = sw_weights_get(sfd->weights, &degree);
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

sw_scale_fct_t *
sw_scale_fct_new(int32_t order)
{
    sw_scale_fct_t *sf;

    if (order <= 0)
        return NULL;

    sf = (sw_scale_fct_t *)calloc(1, sizeof(sw_scale_fct_t));
    if (!sf)
        return NULL;

    sf->scale_fct_base = sw_scale_fct_base_new(-(order >> 1),
					       order - (order >> 1),
					       sw_scale_fct_filter_get(order));
    sf->order = order;
    sf->spline = sw_spline_new(order);
    if (!sf->spline)
    {
        sw_scale_fct_base_del(&sf->scale_fct_base);
        free(sf);
        return NULL;
    }

    sf->proj_periodic_backward = sw_scale_fct_inner_product_periodic_backward;
    sf->proj_dirichlet_backward = sw_scale_fct_inner_product_dirichlet_backward;

    return sf;
}

void
sw_scale_fct_del(sw_scale_fct_t *sf)
{
    if (!sf)
        return;

    sw_scale_fct_base_del(&sf->scale_fct_base);
    sw_spline_del(sf->spline);
    free(sf);
}

const sw_spline_t *
sw_scale_fct_spline_get(const sw_scale_fct_t *sf)
{
    if (!sf)
        return NULL;

    return sf->spline;
}

rational_t
sw_scale_fct_value_rat_get(const sw_scale_fct_t *sf,
			   const rational_t     *val)
{
    if (!sf)
        return rat_new(0, 1, 0);

    return sw_spline_value_rat_get(sf->spline, val);
}

double
sw_scale_fct_value_get(const sw_scale_fct_t *sf,
		       double             val)
{
    return sw_spline_value_get(sf->spline, val);
}

void
sw_scale_fct_proj_periodic_backward(const sw_scale_fct_t *sf,
				    int32_t               scale,
				    const double         *ps,
				    double               *f2)
{
    if (sf)
        sf->proj_periodic_backward(sf, scale, ps, f2);
}

void
sw_scale_fct_proj_dirichlet_backward(const sw_scale_fct_t *sf,
				     int32_t               scale,
				     const double      *ps,
				     double            *f2)
{
    if (sf)
        sf->proj_dirichlet_backward(sf, scale, ps, f2);
}

/* dual scale function methods */

sw_scale_fct_dual_t *
sw_scale_fct_dual_new(int32_t order,
		      int32_t order_dual)
{
    sw_scale_fct_dual_t *sf;

    if ((order <= 0) || (order_dual <= 0))
        return NULL;

    if (((order + order_dual) & 1) == 1)
        return NULL;

    sf = (sw_scale_fct_dual_t *)calloc(1, sizeof(sw_scale_fct_dual_t));
    if (!sf)
        return NULL;

    sf->scale_fct_base = sw_scale_fct_base_new(-(order_dual >> 1) - ((order + order_dual) >> 1) + 1,
                                            -(order_dual >> 1) + ((order + order_dual) >> 1) + order_dual - 1,
                                            sw_scale_fct_dual_filter_get(order, order_dual));
    sf->order = order;
    sf->order_dual = order_dual;

    return sf;
}

void
sw_scale_fct_dual_del(sw_scale_fct_dual_t *sf)
{
    if (!sf)
        return;

    if (sf->weights) sw_weights_del(sf->weights);
    sw_scale_fct_base_del(&sf->scale_fct_base);
    free(sf);
}

uint8_t
sw_scale_fct_dual_type_set(sw_scale_fct_dual_t *sf,
			   sw_weights_type_t       type,
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
         sf->proj_periodic_forward = sw_scale_fct_dual_inner_product_lagrange_periodic_forward;
         sf->proj_dirichlet_forward = sw_scale_fct_dual_inner_product_lagrange_dirichlet_forward;
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
    printf(" $$ %p\n", sf->weights);

    return 1;
}

sw_weights_type_t
sw_scale_fct_dual_type_get(const sw_scale_fct_dual_t *sf)
{
    if (!sf)
        return SW_WEIGHTS_TYPE_ERROR;

    return sf->type;
}

void
sw_scale_fct_dual_proj_periodic_forward(const sw_scale_fct_dual_t *sfd,
					int32_t                    scale,
					const double              *function,
					double                     *ps)
{
    if (sfd)
        sfd->proj_periodic_forward(sfd, scale, function, ps);
}

void
sw_scale_fct_dual_proj_dirichlet_forward(const sw_scale_fct_dual_t *sfd,
					 int32_t                    scale,
					 const double              *function,
					 double                    *ps)
{
    if (sfd)
        sfd->proj_dirichlet_forward(sfd, scale, function, ps);
}

const double *
sw_scale_fct_dual_weights_get(const sw_scale_fct_dual_t *sfd,
			      int32_t                   *size)
{
  printf(" ** %p\n", sfd);
  printf(" ** %p\n", sfd->weights);
    return sw_weights_get(sfd->weights, size);
}
