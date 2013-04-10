#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "sw.h"

/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/


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

double *
sw_new(size_t nbr)
{
    return (double *)malloc(nbr * sizeof(double));
}

void
sw_free(double *f)
{
    if (f) free(f);
}

double
sw_error(double *f1, double *f2, size_t size)
{
    double e1 = 0;
    double e2 = 0;
    double e3 = 0;
    size_t i;

    for (i = 0; i < size; i++)
    {
        e1 += f1[i] * f1[i];
        e2 += f2[i] * f2[i];
        e3 += (f1[i] - f2[i]) * (f1[i] - f2[i]);
    }

    return 2.0 * sqrt(e3) / (sqrt(e1) + sqrt (e2));
}

int64_t
sw_binomial(uint32_t n, uint32_t p)
{
    int64_t num;
    int64_t den;
    int64_t m;
    int64_t i;

    num = 1;
    den = 1;
    m = (p < (n >> 1)) ? p : n - p;
    for (i = 0; i < m; i++)
    {
        num *= (n - i);
        den *= (i + 1);
    }

    return num / den;
}

double
sw_time_get(void)
{
   struct timeval timev;

   gettimeofday(&timev, NULL);
   return (double)timev.tv_sec + (((double)timev.tv_usec) / 1000000);
}
