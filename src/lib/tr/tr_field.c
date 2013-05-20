#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <math.h>

/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

/**
 * @cond SW_LOCAL
 */

static double _tr_Vmax = 1.0;
static double _tr_K = 2.0;
static double _tr_kd = 0.01;
static double _tr_beta = 4.0;
static double _tr_forward = 0;

/**
 * @endcond SW_LOCAL
 */


/******************************************************************************/
/*                                                                            */
/*                                   GLOBAL                                   */
/*                                                                            */
/******************************************************************************/

double
tr_field(double x)
{
  double res;

  if (_tr_beta == 0.0)
    res = _tr_Vmax / 2.0 - _tr_kd * x;
  else
    res = (_tr_Vmax * pow(x, _tr_beta)) / (pow(_tr_K, _tr_beta) + pow(x, _tr_beta)) - _tr_kd * x;

  if (!_tr_forward)
    res = -res;

  return res;
}

/******************************************************************************/
/*                                                                            */
/*                                    API                                     */
/*                                                                            */
/******************************************************************************/

void
tr_field_data_set(double Vmax, double K, double kd, double beta, unsigned int forward)
{
  _tr_Vmax = Vmax;
  _tr_K = K;
  _tr_kd = kd;
  _tr_beta = beta;
  _tr_forward = forward;
}
