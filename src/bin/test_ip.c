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

#include "sw_rational.h"
#include "sw_polynomial.h"
#include "sw_spline.h"
#include "sw_scale_fct.h"
#include "sw_weights.h"

void
test_ip_x (int32_t order, int32_t order_dual, int32_t scale, int32_t degree)
{
    scale_fct_t            *sf;
    scale_fct_dual_t       *sfd;
    const scale_fct_base_t *sfb;
    double                 *f1;
    double                 *f2;
    double                 *ps;
    double                  error;
    double                  error1;
    double                  error2;
    int32_t                 size;
    int32_t                 i;

    size = 1 << scale;

    sf = scale_fct_new (order, SW_WEIGHTS_TYPE_LAGRANGE);
    sfd = scale_fct_dual_new (order, order_dual, SW_WEIGHTS_TYPE_LAGRANGE);
    scale_fct_dual_lagrange_data_set (sfd, degree);
    sfb = (scale_fct_base_t const *)sfd;

    f1 = (double *)malloc (sizeof (double) * size);
    f2 = (double *)malloc (sizeof (double) * size);
    ps = (double *)malloc (sizeof (double) * size);

    for (i = 0; i < size; i++)
    {
        double x;

        x = ((double)i / (double)size - 0.5);
        f1[i] = exp (-300 * x * x);
    }

    scale_fct_dual_proj_periodic_forward (sfd, scale, f1, ps);
    scale_fct_proj_periodic_backward(sf, scale, ps, f2);

    error = 0.0;
    error1 = 0.0;
    error2 = 0.0;
    for (i = 0; i < size; i++)
    {
        error1 += f1[i] * f1[i];
        error2 += f2[i] * f2[i];
        error += (f1[i] - f2[i]) * (f1[i] - f2[i]);
    }
    printf ("error en x : %e\n", 2.0 * sqrt(error) / (sqrt(error1) + sqrt(error2)));

    free (f1);
    free (f2);
    free (ps);
    scale_fct_del (sf);
    scale_fct_dual_del (sfd);
}

void
test_ip_v (int32_t order, int32_t order_dual, int32_t scale, int32_t degree)
{
    scale_fct_t            *sf;
    scale_fct_dual_t       *sfd;
    const scale_fct_base_t *sfb;
    double                 *f1;
    double                 *f2;
    double                 *ps;
    double                  error;
    double                  error1;
    double                  error2;
    int32_t                 size;
    int32_t                 i;

    size = 1 << (scale + 1);

    sf = scale_fct_new (order, SW_WEIGHTS_TYPE_LAGRANGE);
    sfd = scale_fct_dual_new (order, order_dual, SW_WEIGHTS_TYPE_LAGRANGE);
    scale_fct_dual_lagrange_data_set (sfd, degree);
    sfb = (scale_fct_base_t const *)sfd;

    f1 = (double *)malloc (sizeof (double) * size);
    f2 = (double *)malloc (sizeof (double) * size);
    ps = (double *)malloc (sizeof (double) * size);

    for (i = -size / 2; i < size / 2; i++)
    {
        double x;

        x = 2.0 * (double)i / (double)size;
        f1[i + size / 2] = exp (-80 * x * x);
    }

    scale_fct_dual_proj_dirichlet_forward (sfd, scale, f1, ps);
    scale_fct_proj_dirichlet_backward(sf, scale, ps, f2);

    error = 0.0;
    error1 = 0.0;
    error2 = 0.0;
    for (i = 0; i < size; i++) {
        error1 += f1[i] * f1[i];
        error2 += f2[i] * f2[i];
        error += (f1[i] - f2[i]) * (f1[i] - f2[i]);
    }
    printf ("error en v : %e\n", 2.0 * sqrt(error) / (sqrt(error1) + sqrt(error2)));

    free (f1);
    free (f2);
    free (ps);
    scale_fct_del (sf);
    scale_fct_dual_del (sfd);
}

int
main ()
{

    test_ip_x (6, 2, 8, 9);
    test_ip_v (6, 2, 8, 9);


    /*   { */
    /*     rational_t *filter; */
    /*     int32_t N1; */
    /*     int32_t N2; */

    /*     N1 = scale_fct_base_N1_get((scale_fct_base_t const *)sf); */
    /*     N2 = scale_fct_base_N2_get((scale_fct_base_t const *)sf); */
    /*     printf ("bornes : %d %d\n", N1, N2); */
    /*     filter = scale_fct_base_filter_rat_get((scale_fct_base_t const *)sf); */
    /*     for (i = N1; i <= N2; i++) */
    /*       rat_disp(&filter[i - N1], 1); */
    /*     printf ("\n"); */
    /*     printf ("\n"); */

    /*     N1 = scale_fct_base_N1_get((scale_fct_base_t const *)sfd); */
    /*     N2 = scale_fct_base_N2_get((scale_fct_base_t const *)sfd); */
    /*     printf ("bornes : %d %d\n", N1, N2); */
    /*     filter = scale_fct_base_filter_rat_get((scale_fct_base_t const *)sfd); */
    /*     for (i = N1; i <= N2; i++) */
    /*       rat_disp(&filter[i - N1], 1); */
    /*     printf ("\n"); */
    /*   } */

    /*   { */
    /*     FILE *f; */

    /*     f = fopen("dataf.dat", "wb"); */

    /*     for (i = 0; i < size; i++) { */
    /*       fprintf (f, "%f %f %f\n", f1[i], ps[i], f2[i]); */
    /*     } */

    /*     fclose (f); */
    /*   } */

    return 0;
}
