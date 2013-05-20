#ifndef TR_RK_H
#define TR_RK_H


#include <stdint.h>

/**
 * @file tr_rk.h
 * @brief Runge-Kutta routines
 *
 * These functions implements Runge-Kutta routines, threaed or not.
 */

/**
 * @typedef tr_fct_t
 * @brief Typedef of a numerical function.
 */
typedef double (*tr_fct_t)(double x);

/**
 * @brief Implement the classical RK4 algorithm.
 *
 * @param f The function which defines the diffeerential equation.
 * @param y The initial point to start from.
 * @param h The step.
 * @return The approximated solution.
 *
 * This function returns the solution of dz(x)/dx = @p f(z(x)) on
 * [0,@p h], with z(0) = @p y, that is, it returns z(@p h), using the 4th
 * order Runge-Kutta algorithm.
 */
double tr_rk4(tr_fct_t f, double y, double h);


#endif /* TR_RK_H */
