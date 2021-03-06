#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "tr_rk.h"

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

double
tr_rk4(tr_fct_t f, double y, double h)
{
    double k1;
    double k2;
    double k3;
    double k4;

    k1 = h * f(y);
    k2 = h * f(y + 0.5 * k1);
    k3 = h * f(y + 0.5 * k2);
    k4 = h * f(y + k3);

    return y + k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
}
