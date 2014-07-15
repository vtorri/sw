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
#include <stdio.h>

#include "sw.h"


/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

/**
 * @cond SW_LOCAL
 */

struct sw_wavelet_s
{
    sw_wavelet_base_t     wb;
    const sw_scale_fct_t *sf;
    int32_t               x1;
    int32_t               x2;
};

struct sw_wavelet_dual_s
{
    sw_wavelet_base_t wb;
    int32_t           x1;
    int32_t           x2;
};


static sw_wavelet_base_t
sw_wavelet_base_new (const sw_scale_fct_base_t *sfb)
{
    sw_wavelet_base_t w;
    const rational_t *filter;
    int32_t           i;

    w.N1 = 1 - sw_scale_fct_base_N2_get (sfb);
    w.N2 = 1 - sw_scale_fct_base_N1_get (sfb);

    filter = sw_scale_fct_base_filter_rat_get (sfb);

    w.filter_rat = (rational_t *)malloc (sizeof (rational_t) * (w.N2 - w.N1 + 1));
    if (!w.filter_rat)
        return w;

    w.filter = (double *)malloc (sizeof (double) * (w.N2 - w.N1 + 1));
    if (!w.filter)
    {
        free (w.filter_rat);
        return w;
    }

    for (i = w.N1; i <= w.N2; i++)
    {
        if ((i & 1) == 1)
            w.filter_rat[i - w.N1] = rat_opp (&filter[1 - i - sw_scale_fct_base_N1_get (sfb)]);
        else
            w.filter_rat[i - w.N1] = filter[1 - i - sw_scale_fct_base_N1_get (sfb)];
    }

    for (i = w.N1; i <= w.N2; i++)
    {
        w.filter[i - w.N1] = rat_double_get (&w.filter_rat[i - w.N1]);
    }

    return w;
}

static void
sw_wavelet_base_del (sw_wavelet_base_t *w)
{
    if (!w)
        return;

    if (w->filter)
        free (w->filter);
    if (w->filter_rat)
        free (w->filter_rat);
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
 * wavelet function methods
 */

sw_wavelet_t *
sw_wavelet_new(const sw_scale_fct_t *sf, const sw_scale_fct_dual_t *sfd)
{
    sw_wavelet_t *w;
    int32_t       N1;
    int32_t       N2;
    int32_t       N1d;
    int32_t       N2d;

    if (!sf || !sfd)
        return NULL;

    w = (sw_wavelet_t *)malloc (sizeof (sw_wavelet_t));
    if (!w)
        return NULL;

    N1  = sw_scale_fct_base_N1_get ((const sw_scale_fct_base_t *)sf);
    N2  = sw_scale_fct_base_N2_get ((const sw_scale_fct_base_t *)sf);
    N1d = sw_scale_fct_base_N1_get ((const sw_scale_fct_base_t *)sfd);
    N2d = sw_scale_fct_base_N2_get ((const sw_scale_fct_base_t *)sfd);

    w->wb = sw_wavelet_base_new ((const sw_scale_fct_base_t *)sfd);
    w->sf = sf;
    w->x1 = (N1 - N2d + 1) >> 1;
    w->x2 = (N2 - N1d + 1) >> 1;

    return w;
}

void
sw_wavelet_del(sw_wavelet_t *w)
{
    if (!w)
        return;

    sw_wavelet_base_del (&w->wb);
    free (w);
}

int32_t
sw_wavelet_x1_get(const sw_wavelet_t *w)
{
    return w->x1;
}

int32_t
sw_wavelet_x2_get(const sw_wavelet_t *w)
{
    return w->x2;
}

rational_t
sw_wavelet_value_rat_get(const sw_wavelet_t *w,
                         const rational_t   *val)
{
    rational_t tmp;
    rational_t res;
    int32_t    i;

    tmp = rat_new (2, 1, 0);
    tmp = rat_mul (val, &tmp);

    for (i = ((sw_wavelet_base_t *)w)->N1, res = rat_new (0, 1, 0); i <= ((sw_wavelet_base_t *)w)->N2; i++)
    {
        rational_t x;

        x = rat_new (-i, 1, 0);
        x = rat_add (&tmp, &x);

        x = sw_scale_fct_value_rat_get (w->sf, &x);
        x = rat_mul (&x, &((sw_wavelet_base_t *)w)->filter_rat[i - ((sw_wavelet_base_t *)w)->N1]);
        res = rat_add (&res, &x);
    }

    tmp = rat_new (2, 1, 0);

    return rat_mul (&res, &tmp);
}

double
sw_wavelet_value_get(const sw_wavelet_t *w,
                     double              val)
{
    const double *filter;
    double  res;
    int32_t N1;
    int32_t N2;
    int32_t i;

    N1 = ((sw_wavelet_base_t *)w)->N1;
    N2 = ((sw_wavelet_base_t *)w)->N2;
    filter = ((sw_wavelet_base_t *)w)->filter;

    for (i = N1, res = 0.0; i <= N2; i++)
    {
        double x;

        x = filter[i - N1] * sw_scale_fct_value_get (w->sf, 2.0 * val - (double)i);
        res += x;
    }

    return 2.0 * res;
}


/*
 * dual wavelet function methods
 */

sw_wavelet_dual_t *
sw_wavelet_dual_new(const sw_scale_fct_t *sf, const sw_scale_fct_dual_t *sfd)
{
    sw_wavelet_dual_t *wd;
    int32_t            N1;
    int32_t            N2;
    int32_t            N1d;
    int32_t            N2d;

    if (!sf)
        return NULL;

    wd = (sw_wavelet_dual_t *)malloc (sizeof (sw_wavelet_dual_t));
    if (!wd)
        return NULL;

    N1  = sw_scale_fct_base_N1_get ((const sw_scale_fct_base_t *)sf);
    N2  = sw_scale_fct_base_N2_get ((const sw_scale_fct_base_t *)sf);
    N1d = sw_scale_fct_base_N1_get ((const sw_scale_fct_base_t *)sfd);
    N2d = sw_scale_fct_base_N2_get ((const sw_scale_fct_base_t *)sfd);

    wd->wb = sw_wavelet_base_new ((const sw_scale_fct_base_t *)sf);
    wd->x1 = (N1d - N2 + 1) >> 1;
    wd->x2 = (N2d - N1 + 1) >> 1;

    return wd;
}

void
sw_wavelet_dual_del(sw_wavelet_dual_t *wd)
{
    if (!wd)
        return;

    sw_wavelet_base_del (&wd->wb);
    free (wd);
}

int32_t
sw_wavelet_dual_x1_get(const sw_wavelet_dual_t *wd)
{
    return wd->x1;
}

int32_t
sw_wavelet_dual_x2_get(const sw_wavelet_dual_t *wd)
{
    return wd->x2;
}
