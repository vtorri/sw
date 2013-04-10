#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <check.h>

#include "sw_rational.h"
#include "sw_polynomial.h"
#include "sw_spline.h"
#include "sw_scale_fct.h"
#include "sw_weights.h"
#include "sw_wavelet.h"
#include "sw_mra.h"


static double *
gaussian_x (double lambda, double x0, int32_t scale)
{
    double *f;
    int32_t size;
    int32_t i;

    size = 1 << scale;
    f = sw_new(size);
    if (!f)
        return NULL;

    for (i = 0; i < size; i++)
    {
        double x;

        x = (double)i / (double)size;
        f[i] = exp (-lambda * (x - x0) * (x - x0));
    }

    return f;
}

static double *
gaussian_v (double lambda, int32_t scale)
{
    double *f;
    int32_t size;
    int32_t i;

    size = 1 << (scale + 1);
    f = sw_new(size);
    if (!f)
        return NULL;

    for (i = 0; i < size; i++)
    {
        double x;

        x = -1 + 2.0 * (double)i / (double)size;
        f[i] = exp (-lambda * x * x);
    }

    return f;
}

static double *
step_x (int32_t scale)
{
    double *f;
    int32_t size;
    int32_t i;

    size = 1 << scale;
    f = sw_new(size);
    if (!f)
        return NULL;

    for (i = 0; i < size; i++)
    {
        double x;

        x = (double)i / (double)size;
        if ((x <= 0.25) || (x > 0.75))
            f[i] = 0.0;
        else
            f[i] = 1.0;
    }

    return f;
}

static double *
step_v (int32_t scale)
{
    double *f;
    int32_t size;
    int32_t i;

    size = 1 << (scale + 1);
    f = sw_new(size);
    if (!f)
        return NULL;

    for (i = 0; i < size; i++)
    {
        double x;

        x = -1 + 2.0 * (double)i / (double)size;
        if ((x <= 0.5) || (x > 0.5))
            f[i] = 0.0;
        else
            f[i] = 1.0;
    }

    return f;
}

static double *
gaussian_2d (double lambda, int32_t scale)
{
    double *f;
    int32_t size_x;
    int32_t size_v;
    int32_t i;
    int32_t j;

    size_x = 1 << scale;
    size_v = 1 << (scale + 1);
    f = sw_new(size_x * size_v);
    if (!f)
        return NULL;

    for (i = 0; i < size_x; i++)
    {
        double x;

        x = (double)i / (double)size_x - 0.5;
        for (j = 0; j < size_v; j++)
	{
            double v;

            v = -1.0 + 2.0 * (double)j / (double)size_v;
            f[i * size_v + j] = exp (-lambda * x * x * v * v);
	}
    }

    return f;
}


START_TEST(sw_test_rational_new)
{
    rational_t a;

    a = rat_new (5, 1, 1);
    fail_if((a.num != 5) || (a.den != 1));

    a = rat_new (1, -3, 1);
    fail_if((a.num != -1) || (a.den != 3));

    a = rat_new (0, 4, 1);
    fail_if((a.num != 0) || (a.den != 1));
}
END_TEST

START_TEST(sw_test_rational_mul)
{
    rational_t a;

    a = rat_new (5, 1, 1);
    fail_if((a.num != 5) || (a.den != 1));

    a = rat_new (1, -3, 1);
    fail_if((a.num != -1) || (a.den != 3));

    a = rat_new (0, 4, 1);
    fail_if((a.num != 0) || (a.den != 1));
}
END_TEST

START_TEST(sw_test_rational_div)
{
    rational_t a;

    a = rat_new (5, 1, 1);
    fail_if((a.num != 5) || (a.den != 1));

    a = rat_new (1, -3, 1);
    fail_if((a.num != -1) || (a.den != 3));

    a = rat_new (0, 4, 1);
    fail_if((a.num != 0) || (a.den != 1));
}
END_TEST


START_TEST(sw_test_rpol_int)
{
    rational_t  integral;
    rational_t *coefs;
    rpol_t     *pol;
    int32_t     degree;

    degree = 3;
    coefs = (rational_t *)malloc (sizeof (rational_t) * (degree + 1));
    fail_if(coefs == NULL);

    coefs[0] = rat_new (0, 1, 0);
    coefs[1] = rat_new (2, 1, 0);
    coefs[2] = rat_new (1, 1, 0);
    coefs[3] = rat_new (4, 1, 0);

    pol = rpol_new (coefs, degree);
    fail_if(pol == NULL);

    integral = rpol_integrate (pol, -2, 3);
    fail_if((integral.num != 245) || (integral.den != 3));

    rpol_delete (pol);
}
END_TEST


START_TEST(sw_test_spline_1)
{
    spline_t         *spline;
    const rational_t *coefs;
    int64_t           x1;
    int64_t           x2;
    int32_t           order;

    order = 1;

    spline = spline_new (order);
    fail_if(spline == NULL);

    /* support of the spline */
    x1 = spline_x1_get (spline);
    x2 = spline_x2_get (spline);

    fail_if(x1 != 0);
    fail_if(x2 != 1);

    /* coefficients of the spline */
    coefs = spline_rat_coef_get (spline);
    fail_if(coefs == NULL);
    fail_if(rat_is_equal (&coefs[0], 1) == 0);

    spline_del (spline);
}
END_TEST


START_TEST(sw_test_spline_2)
{
    spline_t         *spline;
    const rational_t *coefs;
    int64_t           x1;
    int64_t           x2;
    int32_t           order;

    order = 2;

    spline = spline_new (order);
    fail_if(spline == NULL);

    /* support of the spline */
    x1 = spline_x1_get (spline);
    x2 = spline_x2_get (spline);

    fail_if(x1 != -1);
    fail_if(x2 != 1);

    /* coefficients of the spline */
    coefs = spline_rat_coef_get (spline);
    fail_if(coefs == NULL);
    fail_if(rat_is_equal (&coefs[0], 1)  == 0);
    fail_if(rat_is_equal (&coefs[1], 1)  == 0);
    fail_if(rat_is_equal (&coefs[2], 1)  == 0);
    fail_if(rat_is_equal (&coefs[3], -1) == 0);

    spline_del (spline);
}
END_TEST


START_TEST(sw_test_spline_3)
{
    rational_t        half;
    spline_t         *spline;
    const rational_t *coefs;
    int64_t           x1;
    int64_t           x2;
    int32_t           order;

    order = 3;

    spline = spline_new (order);
    fail_if(spline == NULL);

    /* support of the spline */
    x1 = spline_x1_get (spline);
    x2 = spline_x2_get (spline);

    fail_if(x1 != -1);
    fail_if(x2 != 2);

    /* coefficients of the spline */
    coefs = spline_rat_coef_get (spline);

    half = rat_new (1, 2, 0);
    fail_if(coefs == NULL);
    fail_if(rat_is_equal_rat (&coefs[0], &half)  == 0);
    fail_if(rat_is_equal     (&coefs[1], 1)  == 0);
    fail_if(rat_is_equal_rat (&coefs[2], &half)  == 0);
    fail_if(rat_is_equal_rat (&coefs[3], &half)  == 0);
    fail_if(rat_is_equal     (&coefs[4], 1)  == 0);
    fail_if(rat_is_equal     (&coefs[5], -1)  == 0);
    fail_if(rat_is_equal     (&coefs[6], 2)  == 0);
    fail_if(rat_is_equal     (&coefs[7], -2)  == 0);
    fail_if(rat_is_equal_rat (&coefs[8], &half)  == 0);

    spline_del (spline);
}
END_TEST


START_TEST(sw_test_integration_1)
{
    spline_t         *spline;
    const rational_t *coefs;
    double            integral;
    int64_t           x1;
    int64_t           x2;
    int32_t           order;

    order = 1;

    spline = spline_new (order);
    fail_if(spline == NULL);

    /* support of the spline */
    x1 = spline_x1_get (spline);
    x2 = spline_x2_get (spline);

    fail_if(x1 != 0);
    fail_if(x2 != 1);

    integral = spline_integral_value_get (spline, (double)x1, (double)x2);
    fail_if (fabs(integral - 1.0) > 1e-15);

    integral = spline_integral_value_get (spline, 0.2, 0.3);
    fail_if (fabs(integral - 0.1) > 1e-15);

    integral = spline_integral_value_get (spline, -1.0, (double)x2);
    fail_if (fabs(integral - 1.0) > 1e-15);

    integral = spline_integral_value_get (spline, -1.0, 2.0);
    fail_if (fabs(integral - 1.0) > 1e-15);

    spline_del (spline);
}
END_TEST


START_TEST(sw_test_integration_2)
{
    spline_t         *spline;
    const rational_t *coefs;
    double            integral;
    int64_t           x1;
    int64_t           x2;
    int32_t           order;

    order = 2;

    spline = spline_new (order);
    fail_if(spline == NULL);

    /* support of the spline */
    x1 = spline_x1_get (spline);
    x2 = spline_x2_get (spline);

    fail_if(x1 != -1);
    fail_if(x2 != 1);

    integral = spline_integral_value_get (spline, (double)x1, (double)x2);
    fail_if (fabs(integral - 1.0) > 1e-15);

    integral = spline_integral_value_get (spline, -0.8, 0.5);
    fail_if (fabs(integral - 0.855) > 1e-15);

    integral = spline_integral_value_get (spline, -2.0, (double)x2);
    fail_if (fabs(integral - 1.0) > 1e-15);

    integral = spline_integral_value_get (spline, -2.0, 2.0);
    fail_if (fabs(integral - 1.0) > 1e-15);

    spline_del (spline);
}
END_TEST


START_TEST(sw_test_integration_3)
{
    spline_t         *spline;
    const rational_t *coefs;
    double            integral;
    int64_t           x1;
    int64_t           x2;
    int32_t           order;

    order = 3;

    spline = spline_new (order);
    fail_if(spline == NULL);

    /* support of the spline */
    x1 = spline_x1_get (spline);
    x2 = spline_x2_get (spline);

    fail_if(x1 != -1);
    fail_if(x2 != 2);

    integral = spline_integral_value_get (spline, (double)x1, (double)x2);
    fail_if (fabs(integral - 1.0) > 1e-15);

    integral = spline_integral_value_get (spline, -0.5, 1.5);
    fail_if (fabs(integral - 23.0/24.0) > 1e-15);

    integral = spline_integral_value_get (spline, -2.0, (double)x2);
    fail_if (fabs(integral - 1.0) > 1e-15);

    integral = spline_integral_value_get (spline, -2.0, 3.0);
    fail_if (fabs(integral - 1.0) > 1e-15);

    spline_del (spline);
}
END_TEST


START_TEST(sw_test_scale_function_1)
{
    rational_t              half;
    scale_fct_t            *sf;
    scale_fct_base_t const *sfb;
    const rational_t       *filter;
    int32_t                 order;
    int32_t                 N1;
    int32_t                 N2;

    order = 1;

    sf = scale_fct_new (order, SW_WEIGHTS_TYPE_LAGRANGE);
    fail_if(sf == NULL);

    sfb = (scale_fct_base_t const *)sf;

    /* support of the filter */
    N1 = scale_fct_base_N1_get (sfb);
    N2 = scale_fct_base_N2_get (sfb);

    fail_if(N1 != 0);
    fail_if(N2 != 1);

    /* coefficients of the filter */
    half = rat_new (1, 2, 0);
    filter = scale_fct_base_filter_rat_get (sfb);
    fail_if(filter == NULL);
    fail_if(rat_is_equal_rat (&filter[0], &half) == 0);
    fail_if(rat_is_equal_rat (&filter[1], &half) == 0);

    scale_fct_del (sf);
}
END_TEST


START_TEST(sw_test_scale_function_2)
{
    rational_t              half;
    rational_t              quarter;
    scale_fct_t            *sf;
    scale_fct_base_t const *sfb;
    const rational_t       *filter;
    int32_t                 order;
    int32_t                 N1;
    int32_t                 N2;

    order = 2;

    sf = scale_fct_new (order, SW_WEIGHTS_TYPE_LAGRANGE);
    fail_if(sf == NULL);

    sfb = (scale_fct_base_t const *)sf;

    /* support of the filter */
    N1 = scale_fct_base_N1_get (sfb);
    N2 = scale_fct_base_N2_get (sfb);

    fail_if(N1 != -1);
    fail_if(N2 != 1);

    /* coefficients of the filter */
    half = rat_new (1, 2, 0);
    quarter = rat_new (1, 4, 0);
    filter = scale_fct_base_filter_rat_get (sfb);
    fail_if(filter == NULL);
    fail_if(rat_is_equal_rat (&filter[0], &quarter) == 0);
    fail_if(rat_is_equal_rat (&filter[1], &half) == 0);
    fail_if(rat_is_equal_rat (&filter[2], &quarter) == 0);

    scale_fct_del (sf);
}
END_TEST


START_TEST(sw_test_scale_function_inner_product)
{
    scale_fct_t            *sf;
    scale_fct_dual_t       *sfd;
    weights_t              *weights;
    double                 *f1;
    double                 *f2;
    double                 *ps;
    int32_t                 order;
    int32_t                 order_dual;
    int32_t                 degree;
    int32_t                 scale;
    int32_t                 size;
    int32_t                 i;

    order = 3;
    order_dual = 3;
    degree = 9;
    scale = 8;
    size = 1 << scale;

    sf = scale_fct_new (order, SW_WEIGHTS_TYPE_LAGRANGE);
    fail_if(sf == NULL);

    sfd = scale_fct_dual_new (order, order_dual, SW_WEIGHTS_TYPE_LAGRANGE, degree);
    fail_if(sf == NULL);

 
    {
        const rational_t *filter;
        int32_t N1;
        int32_t N2;

        N1 = scale_fct_base_N1_get((scale_fct_base_t const *)sf);
        N2 = scale_fct_base_N2_get((scale_fct_base_t const *)sf);
        printf ("bornes : %d %d\n", N1, N2);
        filter = scale_fct_base_filter_rat_get((scale_fct_base_t const *)sf);
        for (i = N1; i <= N2; i++)
            rat_disp(&filter[i - N1], 1);
        printf ("\n");
        printf ("\n");

   
        {
            rational_t x;
            rational_t y1;
            rational_t y2;
            rational_t tmp1;
            rational_t tmp2;
            rational_t tmp3;
            rational_t res;
            int k;

            x = rat_new(1, 4, 0);
            y1 = scale_fct_value_rat_get (sf, &x);

            tmp1 = rat_new(2, 1, 0);
            tmp1 = rat_mul(&x, &tmp1);
            y2 = rat_new(0, 1, 0);
            for (k = N1; k <= N2; k++)
            {
                tmp2 = rat_new(k, 1, 0);
                tmp2 = rat_sub(&tmp1, &tmp2);
                tmp3 = scale_fct_value_rat_get (sf, &tmp2);
                tmp3 = rat_mul(&tmp3, &filter[k - N1]);
                y2 = rat_add(&y2, &tmp3);
            }
            tmp2 = rat_new(2, 1, 0);
            y2 = rat_mul(&y2, &tmp2);
            rat_disp(&y1, 1);
            rat_disp(&y2, 1);
            printf("\n");
        }

        N1 = scale_fct_base_N1_get((scale_fct_base_t const *)sfd);
        N2 = scale_fct_base_N2_get((scale_fct_base_t const *)sfd);
        printf ("bornes : %d %d\n", N1, N2);
        filter = scale_fct_base_filter_rat_get((scale_fct_base_t const *)sfd);
        for (i = N1; i <= N2; i++)
            rat_disp(&filter[i - N1], 1);
        printf ("\n");
    }

    f1 = sw_new(size);
    f2 = sw_new(size);
    ps = sw_new(size);

    for (i = 0; i < size; i++)
    {
        double x;

        x = ((double)i / (double)size - 0.5);
        f1[i] = exp (-300 * x * x);
    }

    scale_fct_dual_proj_periodic_forward(sfd, scale, f1, ps);
    scale_fct_proj_periodic_backward(sf, scale, ps, f2);

    printf ("error : %f\n", sw_error(f1, f2, size));

 
    {
        FILE *f;

        f = fopen("data.data", "wb");

        for (i = 0; i < size; i++)
        {
            fprintf (f, "%f %f %f\n", f1[i], ps[i], f2[i]);
        }

        fclose (f);
    }

    sw_free (ps);
    sw_free (f2);
    sw_free (f1);
    scale_fct_dual_del (sfd);
    scale_fct_del (sf);
}
END_TEST


START_TEST(sw_test_projection)
{
    mra_t  *mra;
    double *gaussian;
    double *f;
    double *ps;
    int32_t order;
    int32_t order_dual;
    int32_t scale;
    int32_t degree;

    order = 6;
    order_dual = 2;
    scale = 8;
    degree = 9;

    mra = mra_new (order, order_dual, 3, scale, SW_WEIGHTS_TYPE_LAGRANGE, degree);
    fail_if(mra == NULL);

    /* X direction */

    gaussian = gaussian_x (80, 0.5, scale);
    fail_if(gaussian == NULL);

    f = sw_new(1 << scale);
    fail_if(f == NULL);

    ps = sw_new(1 << scale);
    fail_if(ps == NULL);

    mra_proj_x_forward (mra, gaussian, ps);
    mra_proj_x_backward (mra, ps, f);

    printf ("error projection x : %e\n", error (gaussian, f, 1 << scale));

    sw_free (ps);
    sw_free (f);
    sw_free (gaussian);

    /* V direction */

    gaussian = gaussian_v (80, scale);
    fail_if(gaussian == NULL);

    f = sw_new(1 << (scale + 1));
    fail_if(f == NULL);

    ps = sw_new(1 << (scale + 1));
    fail_if(ps == NULL);

    mra_proj_v_forward (mra, gaussian, ps);
    mra_proj_v_backward (mra, ps, f);

    printf ("error projection v : %e\n", error (gaussian, f, 1 << scale));

    sw_free (ps);
    sw_free (f);
    sw_free (gaussian);

    mra_del (mra);
}
END_TEST


START_TEST(sw_test_projection_2d)
{
    mra_t  *mra;
    double *gaussian;
    double *ps;
    double *res;
    double *tmp1;
    double *tmp2;
    int32_t order;
    int32_t order_dual;
    int32_t scale;
    int32_t degree;
    int32_t size_x;
    int32_t size_v;

    order = 6;
    order_dual = 2;
    scale = 8;
    degree = 9;

    size_x = 1 << scale;
    size_v = 1 << (scale + 1);

    mra = mra_new (order, order_dual, 3, scale, SW_WEIGHTS_TYPE_LAGRANGE, degree);
    fail_if(mra == NULL);

    gaussian = gaussian_2d (80, scale);
    fail_if(gaussian == NULL);

    ps = sw_new(size_x * size_v);
    fail_if(ps == NULL);

    res = sw_new(size_x * size_v);
    fail_if(res == NULL);

    tmp1 = sw_new(size_x);
    fail_if(tmp1 == NULL);

    tmp2 = sw_new(size_x);
    fail_if(tmp2 == NULL);

    mra_proj_2d_forward (mra, gaussian, ps, tmp1, tmp2);
    mra_proj_2d_backward (mra, ps, res, tmp1, tmp2);

    printf ("error projection 2d : %e\n", error (gaussian, res, size_x * size_v));

    sw_free (tmp2);
    sw_free (tmp1);
    sw_free (res);
    sw_free (ps);
    sw_free (gaussian);

    mra_del (mra);
}
END_TEST


START_TEST(sw_test_wavelet_1_3)
{
    rational_t            rat;
    scale_fct_t          *sf;
    scale_fct_dual_t     *sfd;
    wavelet_t            *w;
    wavelet_base_t const *wb;
    const rational_t     *filter;
    int32_t               order;
    int32_t               order_dual;
    int32_t               N1;
    int32_t               N2;

    order = 1;
    order_dual = 3;

    sf = scale_fct_new (order, SW_WEIGHTS_TYPE_LAGRANGE);
    fail_if(sf == NULL);

    sfd = scale_fct_dual_new (order, order_dual, SW_WEIGHTS_TYPE_LAGRANGE, 9);
    fail_if(sfd == NULL);

    w = wavelet_new (sf, sfd);

    wb = (wavelet_base_t const *)w;

    /* support of the filter */
    N1 = wavelet_base_N1_get (wb);
    N2 = wavelet_base_N2_get (wb);

    fail_if(N1 != -2);
    fail_if(N2 != 3);

    /* coefficients of the filter */
    filter = wavelet_base_filter_rat_get (wb);
    fail_if(filter == NULL);
    rat = rat_new (-1, 16, 0);
    fail_if(rat_is_equal_rat (&filter[0], &rat) == 0);
    rat = rat_new (-1, 16, 0);
    fail_if(rat_is_equal_rat (&filter[1], &rat) == 0);
    rat = rat_new (1, 2, 0);
    fail_if(rat_is_equal_rat (&filter[2], &rat) == 0);
    rat = rat_new (-1, 2, 0);
    fail_if(rat_is_equal_rat (&filter[3], &rat) == 0);
    rat = rat_new (1, 16, 0);
    fail_if(rat_is_equal_rat (&filter[4], &rat) == 0);
    rat = rat_new (1, 16, 0);
    fail_if(rat_is_equal_rat (&filter[5], &rat) == 0);

    printf ("x1: %d\n", wavelet_x1_get(w));
    printf ("x2: %d\n", wavelet_x2_get(w));

 
    {
        FILE *f;
        int i;

        f = fopen("w_1_3.dat", "wb");

        for (i = 0; i < 512; i++)
        {
            double x;

            x = wavelet_x1_get(w) + (wavelet_x2_get(w) - wavelet_x1_get(w)) * (double)i / 512.0;
            fprintf (f, "%f %f\n", x, wavelet_value_get(w, x));
        }

        fclose (f);
    }

    wavelet_del (w);
    scale_fct_dual_del (sfd);
    scale_fct_del (sf);
}
END_TEST


START_TEST(sw_test_wavelet_dual_1_3)
{
    rational_t            rat;
    scale_fct_t          *sf;
    scale_fct_dual_t     *sfd;
    wavelet_dual_t       *wd;
    wavelet_base_t const *wb;
    const rational_t     *filter;
    int32_t               order;
    int32_t               order_dual;
    int32_t               N1;
    int32_t               N2;

    order = 1;
    order_dual = 3;

    sf = scale_fct_new (order, SW_WEIGHTS_TYPE_LAGRANGE);
    fail_if(sf == NULL);

    sfd = scale_fct_dual_new (order, order_dual, SW_WEIGHTS_TYPE_LAGRANGE, 9);
    fail_if(sfd == NULL);

    wd = wavelet_dual_new (sf, sfd);

    wb = (wavelet_base_t const *)wd;

    /* support of the filter */
    N1 = wavelet_base_N1_get (wb);
    N2 = wavelet_base_N2_get (wb);

    fail_if(N1 != 0);
    fail_if(N2 != 1);

    /* coefficients of the filter */
    filter = wavelet_base_filter_rat_get (wb);
    rat_disp(&filter[0], 1);
    rat_disp(&filter[1], 1);
    fail_if(filter == NULL);
    rat = rat_new (1, 2, 0);
    fail_if(rat_is_equal_rat (&filter[0], &rat) == 0);
    rat = rat_new (-1, 2, 0);
    fail_if(rat_is_equal_rat (&filter[1], &rat) == 0);

    printf ("x1: %d\n", wavelet_dual_x1_get(wd));
    printf ("x2: %d\n", wavelet_dual_x2_get(wd));

    /*  
        { */
    /*     FILE *f; */
    /*     int i; */

    /*     f = fopen("wd_1_3.dat", "wb"); */

    /*     for (i = 0; i < 512; i++)
           { */
    /*       double x; */

    /*       x = wavelet_dual_x1_get(wd) + (wavelet_dual_x2_get(wd) - wavelet_dual_x1_get(wd)) * (double)i / 512.0; */
    /*       fprintf (f, "%f %f\n", x, wavelet_dual_value_get(wd, x)); */
    /*     } */

    /*     fclose (f); */
    /*   } */

    wavelet_dual_del (wd);
    scale_fct_dual_del (sfd);
    scale_fct_del (sf);
}
END_TEST


START_TEST(sw_test_fwt_x_1d_gaussian)
{
    mra_t  *mra;
    double *gaussian;
    double *ps_fine;
    double *ps_coarse;
    double *ps_w;
    double *f;
    int32_t order;
    int32_t order_dual;
    int32_t scale_fine;
    int32_t scale_coarse;
    int32_t degree;
    int32_t size_x;

    order = 2;
    order_dual = 2;
    scale_fine = 5;
    scale_coarse = 4;
    degree = 9;

    size_x = 1 << scale_fine;

    mra = mra_new (order, order_dual, scale_coarse, scale_fine, SW_WEIGHTS_TYPE_LAGRANGE, degree);
    fail_if(mra == NULL);

    gaussian = gaussian_x (80, 0.5, scale_fine);
    fail_if(gaussian == NULL);

    ps_fine = sw_new(size_x);
    fail_if(ps_fine == NULL);

    f = sw_new(size_x);
    fail_if(ps_fine == NULL);

    ps_coarse = sw_new(1 << scale_coarse);
    fail_if(ps_coarse == NULL);

    ps_w = sw_new((1 << scale_fine) - (1 << scale_coarse));
    fail_if(ps_w == NULL);

    mra_proj_x_forward (mra, gaussian, ps_fine);
    mra_fwt_x_forward (mra, ps_fine, ps_coarse, ps_w);

    printf ("size ps_w : %d\n", ((1 << scale_fine) - (1 << scale_coarse)));

 
    {
        int i, j, index;
        double val, val1;

        printf (" ## scale fine : \n");
        for (i = 0; i < (1 << scale_fine); i++)
            printf (" * %d %f\n", i, ps_fine[i]);

        printf ("\n ## scale coarse : \n");
        for (i = 0; i < (1 << scale_coarse); i++)
            printf (" * %f\n", ps_coarse[i]);

        printf ("\n ## ondeletttes : \n");
        printf ("\n");
        for (j = scale_fine - 1, index = 0; j >= scale_coarse; index += (1 << j), j--)
     
        {
            printf (" ## %d\n", j);
            for (i = 0; i < (1 << j); i++)
	 
            {
                printf (" * %f\n", ps_w[index + i]);
            }
        }
        printf ("\nresultat\n");

        mra_fwt_x_backward (mra, ps_coarse, ps_w, f);
        for (i = 0, val = 0.0, val1 = 0.0; i < (1 << scale_fine); i++)
     
        {
            val += fabs(f[i] - ps_fine[i]) * fabs(f[i] - ps_fine[i]);
            val1 += fabs(f[i] ) * fabs(f[i]);
            printf (" * %02d %+9f %+9f %f %f\n", i, f[i], ps_fine[i], fabs(f[i] - ps_fine[i]), fabs(100.0 * f[i] / ps_fine[i]));
        }
        printf ("err fwt : %f\n", sqrt(val) / sqrt(val1));

        mra_proj_x_backward (mra, ps_fine, f);
        printf ("error projection x : %e\n", error (gaussian, f, 1 << scale_fine));
    }

    sw_free (ps_w);
    sw_free (ps_coarse);
    sw_free (f);
    sw_free (ps_fine);
    sw_free (gaussian);

    mra_del (mra);
}
END_TEST


START_TEST(sw_test_fwt_x_1d_step)
{
    mra_t  *mra;
    double *step;
    double *ps_fine;
    double *ps_coarse;
    double *ps_w;
    double *f;
    int32_t order;
    int32_t order_dual;
    int32_t scale_fine;
    int32_t scale_coarse;
    int32_t degree;
    int32_t size_x;

    order = 1;
    order_dual = 1;
    scale_fine = 5;
    scale_coarse = 4;
    degree = 9;

    size_x = 1 << scale_fine;

    mra = mra_new (order, order_dual, scale_coarse, scale_fine, SW_WEIGHTS_TYPE_LAGRANGE, degree);
    fail_if(mra == NULL);

    step = step_x (scale_fine);
    fail_if(step == NULL);

    ps_fine = sw_new(size_x);
    fail_if(ps_fine == NULL);

    f = sw_new(size_x);
    fail_if(ps_fine == NULL);

    ps_coarse = sw_new(1 << scale_coarse);
    fail_if(ps_coarse == NULL);

    ps_w = sw_new((1 << scale_fine) - (1 << scale_coarse));
    fail_if(ps_w == NULL);

    mra_proj_x_forward (mra, step, ps_fine);
    mra_fwt_x_forward (mra, ps_fine, ps_coarse, ps_w);

    printf ("size ps_w : %d\n", ((1 << scale_fine) - (1 << scale_coarse)));

 
    {
        int i, j, index;
        double val, val1;

        printf (" ## scale fine : \n");
        for (i = 0; i < (1 << scale_fine); i++)
            printf (" * %d %f\n", i, ps_fine[i]);

        printf ("\n ## scale coarse : \n");
        for (i = 0; i < (1 << scale_coarse); i++)
            printf (" * %f\n", ps_coarse[i]);

        printf ("\n ## ondeletttes : \n");
        printf ("\n");
        for (j = scale_fine - 1, index = 0; j >= scale_coarse; index += (1 << j), j--)
     
        {
            printf (" ## %d\n", j);
            for (i = 0; i < (1 << j); i++)
	 
            {
                printf (" * %f\n", ps_w[index + i]);
            }
        }
        printf ("\nresultat\n");

        mra_fwt_x_backward (mra, ps_coarse, ps_w, f);
        for (i = 0, val = 0.0, val1 = 0.0; i < (1 << scale_fine); i++)
     
        {
            val += fabs(f[i] - ps_fine[i]) * fabs(f[i] - ps_fine[i]);
            val1 += fabs(f[i] ) * fabs(f[i]);
            printf (" * %02d %+9f %+9f %f %f\n", i, f[i], ps_fine[i], fabs(f[i] - ps_fine[i]), fabs(100.0 * f[i] / ps_fine[i]));
        }
        printf ("err fwt : %f\n", sqrt(val) / sqrt(val1));

        mra_proj_x_backward (mra, ps_fine, f);
        printf ("error projection x : %e\n", error (step, f, 1 << scale_fine));
    }

    sw_free (ps_w);
    sw_free (ps_coarse);
    sw_free (f);
    sw_free (ps_fine);
    sw_free (step);

    mra_del (mra);
}
END_TEST


/* START_TEST(sw_test_fwt_v_1d) */
/*
  { */
/*   mra_t  *mra; */
/*   double *step; */
/*   double *ps_fine; */
/*   double *ps_coarse; */
/*   double *ps_w; */
/*   double *f; */
/*   int32_t order; */
/*   int32_t order_dual; */
/*   int32_t scale_fine; */
/*   int32_t scale_coarse; */
/*   int32_t degree; */
/*   int32_t size_v; */

/*   order = 2; */
/*   order_dual = 2; */
/*   scale_fine = 5; */
/*   scale_coarse = 4; */
/*   degree = 9; */

/*   size_v = 1 << (scale_fine + 1); */

/*   mra = mra_new (order, order_dual, scale_coarse, scale_fine, SW_WEIGHTS_TYPE_LAGRANGE, degree); */
/*   fail_if(mra == NULL); */

/*   step = step_v (scale_fine); */
/*   fail_if(gaussian == NULL); */

/*   ps_fine = (double *)malloc (sizeof (double) * size_v); */
/*   fail_if(ps_fine == NULL); */

/*   f = (double *)malloc (sizeof (double) * size_v); */
/*   fail_if(ps_fine == NULL); */

/*   ps_coarse = (double *)malloc (sizeof (double) * (1 << (scale_coarse + 1))); */
/*   fail_if(ps_coarse == NULL); */

/*   ps_w = (double *)malloc (sizeof (double) * ((1 << (scale_fine + 1)) - (1 << (scale_coarse + 1)))); */
/*   fail_if(ps_w == NULL); */

/*   mra_proj_v_forward (mra, gaussian, ps_fine); */
/*   mra_fwt_v_forward (mra, ps_fine, ps_coarse, ps_w); */

/*   printf ("size ps_w : %d\n", ((1 << (scale_fine + 1)) - (1 << (scale_coarse + 1)))); */

/*  
    { */
/*     int i, j, index; */
/*     double val, val1; */

/*     printf (" ## scale fine : \n"); */
/*     for (i = 0; i < size_v; i++) */
/*       printf (" * %d %f\n", i, ps_fine[i]); */

/*     printf ("\n ## scale coarse : \n"); */
/*     for (i = 0; i < (1 << (scale_coarse + 1)); i++) */
/*       printf (" * %f\n", ps_coarse[i]); */

/*     printf ("\n ## ondeletttes : \n"); */
/*     printf ("\n"); */
/*     for (j = scale_fine - 1, index = 0; j >= scale_coarse; index += (1 << j), j--) */
/*      
	{ */
/* 	printf (" ## %d\n", j); */
/* 	for (i = 0; i < (1 << j); i++) */
/* 	 
	 { */
/* 	    printf (" * %f\n", ps_w[index + i]); */
/* 	  } */
/*       } */
/*     printf ("\nresultat\n"); */

/*     mra_fwt_v_backward (mra, ps_coarse, ps_w, f); */
/*     for (i = 0, val = 0.0, val1 = 0.0; i < size_v; i++) */
/*      
	{ */
/* 	val += fabs(f[i] - ps_fine[i]) * fabs(f[i] - ps_fine[i]); */
/* 	val1 += fabs(f[i] ) * fabs(f[i]); */
/* 	printf (" * %02d %+9f %+9f %f %f\n", i, f[i], ps_fine[i], fabs(f[i] - ps_fine[i]), fabs(100.0 * f[i] / ps_fine[i])); */
/*       } */
/*     printf ("err fwt : %f\n", sqrt(val) / sqrt(val1)); */

/*     mra_proj_v_backward (mra, ps_fine, f); */
/*     printf ("error projection x : %e\n", error (gaussian, f, 1 << scale_fine)); */
/*   } */

/*   free (ps_w); */
/*   free (ps_coarse); */
/*   free (f); */
/*   free (ps_fine); */
/*   free (step); */

/*   mra_del (mra); */
/* } */
/* END_TEST */


Suite *
sw_suite(void)
{
    Suite *s;
    TCase *tc;

    s = suite_create("Sw");

    tc = tcase_create("Sw_Rational");
    tcase_add_test(tc, sw_test_rational_new);
    tcase_add_test(tc, sw_test_rational_mul);
    tcase_add_test(tc, sw_test_rational_div);
    suite_add_tcase(s, tc);

    tc = tcase_create("Sw Polynomial");
    tcase_add_test(tc, sw_test_rpol_int);
    suite_add_tcase(s, tc);

    tc = tcase_create("Sw Spline");
    tcase_add_test(tc, sw_test_spline_1);
    tcase_add_test(tc, sw_test_spline_2);
    tcase_add_test(tc, sw_test_spline_3);
    tcase_add_test(tc, sw_test_integration_1);
    tcase_add_test(tc, sw_test_integration_2);
    tcase_add_test(tc, sw_test_integration_3);
    suite_add_tcase(s, tc);

    tc = tcase_create("Sw Scale Function");
    tcase_add_test(tc, sw_test_scale_function_1);
    tcase_add_test(tc, sw_test_scale_function_2);
    tcase_add_test(tc, sw_test_scale_function_inner_product);
    suite_add_tcase(s, tc);

    tc = tcase_create("Sw Projection");
    tcase_add_test(tc, sw_test_projection);
    tcase_add_test(tc, sw_test_projection_2d);
    suite_add_tcase(s, tc);

    tc = tcase_create("Sw Wavelet");
    tcase_add_test(tc, sw_test_wavelet_1_3);
    tcase_add_test(tc, sw_test_wavelet_dual_1_3);
    suite_add_tcase(s, tc);

    tc = tcase_create("Sw FWT Forward");
    tcase_add_test(tc, sw_test_fwt_x_1d_gaussian);
    tcase_add_test(tc, sw_test_fwt_x_1d_step);
    /*    tcase_add_test(tc, sw_test_fwt_v_1d); */
    suite_add_tcase(s, tc);

    return s;
}

int
main(void)
{
    Suite   *s;
    SRunner *sr;
    int      failed_count;

    s = sw_suite();
    sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    failed_count = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (failed_count == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
