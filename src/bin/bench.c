#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "sw.h"

#include "algebra/sw_polynomial.h"

static double gauss(double l, double x)
{
    return exp (-l * x * x);
}

static void
test_ip_x(int32_t order, int32_t order_dual, int32_t scale, int32_t degree)
{
    sw_scale_fct_t      *sf;
    sw_scale_fct_dual_t *sfd;
    double           *f1;
    double           *f2;
    double           *ps;
    int32_t           size;
    int32_t           i;

    size = 1 << scale;

    sf = sw_scale_fct_new(order);
    sfd = sw_scale_fct_dual_new(order, order_dual);
    sw_scale_fct_dual_type_set(sfd, SW_WEIGHTS_TYPE_LAGRANGE, degree);

    f1 = sw_new(size);
    f2 = sw_new(size);
    ps = sw_new(size);

    for (i = 0; i < size; i++)
    {
        double x;

        x = ((double)i / (double)size - 0.5);
        f1[i] = gauss(100, x);
    }

    sw_scale_fct_dual_proj_periodic_forward (sfd, scale, f1, ps);
    sw_scale_fct_proj_periodic_backward(sf, scale, ps, f2);

    printf ("proj x  : %e\n", sw_error(f1, f2, size));

    sw_free(ps);
    sw_free(f2);
    sw_free(f1);
    sw_scale_fct_del (sf);
    sw_scale_fct_dual_del (sfd);
}

static void
test_mra(int32_t order, int32_t order_dual, int32_t scale, int32_t degree)
{
    sw_mra_t                  *mra;
    const sw_scale_fct_t      *sf;
    const sw_scale_fct_dual_t *sfd;
    double                 *f1;
    double                 *f2;
    double                 *ps;
    int32_t                 size;
    int32_t                 i;

    size = 1 << scale;

    mra = sw_mra_new (order, order_dual, 3, scale, SW_WEIGHTS_TYPE_LAGRANGE, degree);
    sf = sw_mra_scale_fct_get(mra);
    sfd = sw_mra_scale_fct_dual_get(mra);

    f1 = sw_new(size);
    f2 = sw_new(size);
    ps = sw_new(size);

    for (i = 0; i < size; i++)
    {
        double x;

        x = ((double)i / (double)size - 0.5);
        f1[i] = gauss(100, x);
    }

    sw_scale_fct_dual_proj_periodic_forward (sfd, scale, f1, ps);
    sw_scale_fct_proj_periodic_backward(sf, scale, ps, f2);

    printf ("mra per : %e\n", sw_error(f1, f2, size));

    sw_free(ps);
    sw_free(f2);
    sw_free(f1);
    sw_mra_del (mra);
}

int main()
{
    test_ip_x(6, 2, 8, 9);
    test_mra(6, 2, 8, 9);

    return 0;
}
