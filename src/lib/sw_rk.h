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
