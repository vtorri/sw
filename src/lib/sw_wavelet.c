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

struct wavelet
{
    wavelet_base_t     wb;
    const scale_fct_t *sf;
    int32_t            x1;
    int32_t            x2;
};

struct wavelet_dual
{
    wavelet_base_t wb;
    int32_t        x1;
    int32_t        x2;
};


static wavelet_base_t
wavelet_base_new (const scale_fct_base_t *sfb)
{
    wavelet_base_t    w;
    const rational_t *filter;
    int32_t           i;

    w.N1 = 1 - scale_fct_base_N2_get (sfb);
    w.N2 = 1 - scale_fct_base_N1_get (sfb);

    filter = scale_fct_base_filter_rat_get (sfb);

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
            w.filter_rat[i - w.N1] = rat_opp (&filter[1 - i - scale_fct_base_N1_get (sfb)]);
        else
            w.filter_rat[i - w.N1] = filter[1 - i - scale_fct_base_N1_get (sfb)];
    }

    for (i = w.N1; i <= w.N2; i++)
    {
        w.filter[i - w.N1] = rat_double_get (&w.filter_rat[i - w.N1]);
    }

    return w;
}

static void
wavelet_base_del (wavelet_base_t *w)
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

wavelet_t *
wavelet_new (const scale_fct_t *sf, const scale_fct_dual_t *sfd)
{
    wavelet_t *w;
    int32_t    N1;
    int32_t    N2;
    int32_t    N1d;
    int32_t    N2d;

    if (!sf || !sfd)
        return NULL;

    w = (wavelet_t *)malloc (sizeof (wavelet_t));
    if (!w)
        return NULL;

    N1  = scale_fct_base_N1_get ((const scale_fct_base_t *)sf);
    N2  = scale_fct_base_N2_get ((const scale_fct_base_t *)sf);
    N1d = scale_fct_base_N1_get ((const scale_fct_base_t *)sfd);
    N2d = scale_fct_base_N2_get ((const scale_fct_base_t *)sfd);

    w->wb = wavelet_base_new ((const scale_fct_base_t *)sfd);
    w->sf = sf;
    w->x1 = (N1 - N2d + 1) >> 1;
    w->x2 = (N2 - N1d + 1) >> 1;

    return w;
}

void
wavelet_del (wavelet_t *w)
{
    if (!w)
        return;

    wavelet_base_del (&w->wb);
    free (w);
}

int32_t
wavelet_x1_get (const wavelet_t *w)
{
    return w->x1;
}

int32_t
wavelet_x2_get (const wavelet_t *w)
{
    return w->x2;
}

rational_t
wavelet_value_rat_get (const wavelet_t  *w,
                       const rational_t *val)
{
    rational_t tmp;
    rational_t res;
    int32_t    i;

    tmp = rat_new (2, 1, 0);
    tmp = rat_mul (val, &tmp);

    for (i = ((wavelet_base_t *)w)->N1, res = rat_new (0, 1, 0); i <= ((wavelet_base_t *)w)->N2; i++)
    {
        rational_t x;

        x = rat_new (-i, 1, 0);
        x = rat_add (&tmp, &x);

        x = scale_fct_value_rat_get (w->sf, &x);
        x = rat_mul (&x, &((wavelet_base_t *)w)->filter_rat[i - ((wavelet_base_t *)w)->N1]);
        res = rat_add (&res, &x);
    }

    tmp = rat_new (2, 1, 0);

    return rat_mul (&res, &tmp);
}

double
wavelet_value_get (const wavelet_t *w,
                   double           val)
{
    const double *filter;
    double  res;
    int32_t N1;
    int32_t N2;
    int32_t i;

    N1 = ((wavelet_base_t *)w)->N1;
    N2 = ((wavelet_base_t *)w)->N2;
    filter = ((wavelet_base_t *)w)->filter;

    for (i = N1, res = 0.0; i <= N2; i++)
    {
        double x;

        x = filter[i - N1] * scale_fct_value_get (w->sf, 2.0 * val - (double)i);
        res += x;
    }

    return 2.0 * res;
}


/*
 * dual wavelet function methods
 */

wavelet_dual_t *
wavelet_dual_new (const scale_fct_t *sf, const scale_fct_dual_t *sfd)
{
    wavelet_dual_t *wd;
    int32_t         N1;
    int32_t         N2;
    int32_t         N1d;
    int32_t         N2d;

    if (!sf)
        return NULL;

    wd = (wavelet_dual_t *)malloc (sizeof (wavelet_dual_t));
    if (!wd)
        return NULL;

    N1  = scale_fct_base_N1_get ((const scale_fct_base_t *)sf);
    N2  = scale_fct_base_N2_get ((const scale_fct_base_t *)sf);
    N1d = scale_fct_base_N1_get ((const scale_fct_base_t *)sfd);
    N2d = scale_fct_base_N2_get ((const scale_fct_base_t *)sfd);

    wd->wb = wavelet_base_new ((const scale_fct_base_t *)sf);
    wd->x1 = (N1d - N2 + 1) >> 1;
    wd->x2 = (N2d - N1 + 1) >> 1;

    return wd;
}

void
wavelet_dual_del (wavelet_dual_t *wd)
{
    if (!wd)
        return;

    wavelet_base_del (&wd->wb);
    free (wd);
}

int32_t
wavelet_dual_x1_get (const wavelet_dual_t *wd)
{
    return wd->x1;
}

int32_t
wavelet_dual_x2_get (const wavelet_dual_t *wd)
{
    return wd->x2;
}
