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
