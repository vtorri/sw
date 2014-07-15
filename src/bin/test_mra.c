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

#include "sw_mra2.h"

int
main ()
{
    mra_t *mra;
    double                 *f1;
    double                 *f2;
    double                  error;
    double                  error1;
    double                  error2;
    int32_t                 order;
    int32_t                 order_dual;
    int32_t                 degree;
    int32_t                 scale;
    int32_t                 size;
    int32_t                 i;

    printf (" * 0\n");
    order = 6;
    order_dual = 2;
    degree = 9;
    scale = 8;
    size = 1 << scale;

    mra = mra_new (order, order_dual, 3, scale, SW_WEIGHTS_TYPE_LAGRANGE, degree);
    if (!mra)
        return EXIT_FAILURE;

    f1 = (double *)malloc (sizeof (double) * size);
    f2 = (double *)malloc (sizeof (double) * size);

    for (i = 0; i < size; i++)
    {
        double x;

        x = ((double)i / (double)size - 0.5);
        f1[i] = exp (-300 * x * x);
    }

    mra_function_x_set (mra, f1);
    mra_proj_x_forward (mra);
    mra_proj_x_backward (mra);
    mra_function_x_get (mra, f2);

    error = 0.0;
    error1 = 0.0;
    error2 = 0.0;
    for (i = 0; i < size; i++)
    {
        error1 += f1[i] * f1[i];
        error2 += f2[i] * f2[i];
        error += (f1[i] - f2[i]) * (f1[i] - f2[i]);
    }
    printf ("error : %e\n", 2.0 * sqrt(error) / (sqrt(error1) + sqrt(error2)));

    {
        FILE *f;

        f = fopen("data.data", "wb");

        for (i = 0; i < size; i++)
        {
            fprintf (f, "%f %f\n", f1[i], f2[i]);
        }

        fclose (f);
    }

    free (f1);
    free (f2);
    mra_del (mra);

    return 0;
}
