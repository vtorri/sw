#ifndef __SW_RK_H__
#define __SW_RK_H__


#include <stdint.h>

/**
 * @file sw_rk.h
 * @brief Runge-Kutta routines
 *
 * These functions implements Runge-Kutta routines, threaed or not.
 */

/**
 * @typedef sw_fct_t
 * @brief Typedef of a numerical function.
 */
typedef double (*sw_fct_t)(double x);

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
SAPI double sw_rk4(sw_fct_t f, double y, double h);


#endif /* __SW_RK_H__ */
