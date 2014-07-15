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
