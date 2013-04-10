
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
