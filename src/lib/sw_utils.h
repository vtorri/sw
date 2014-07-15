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

#ifndef SW_UTILS_H
#define SW_UTILS_H


#include <stdint.h>

/**
 * @file sw_utils.h
 * @brief Utilities routines
 *
 * These functions implements routines to ease the use of the library
 * (allocation, error, etc...).
 */

/**
 * @brief Allocate an array of double.
 *
 * @param nbr The number of doubles.
 * @return The allocated array.
 *
 * This function just returns the allocation of @p nbr doubles.
 */
SAPI double *sw_new(size_t nbr);

/**
 * @brief Free the given array of double.
 *
 * @param t The array of doubles.
 *
 * This function just frees @p t. If @p t is @c NULL, this function
 * does nothing.
 */
SAPI void sw_free(double *t);

/**
 * @brief Compute the relative error between two functions.
 *
 * @param f1 An array of double.
 * @param f2 An array of double.
 * @param size The size of the arrays.
 * @return The relative error.
 *
 * This function computes the relative error of @p f1 and @p f2. The
 * two arrays must have the same size @p size. No check is done on @p f1
 * or @p f2.
 */
SAPI double sw_error(double *f1, double *f2, size_t size);

/**
 * @brief Return the binomial of two integers.
 *
 * @param n First integer.
 * @param p Second integer.
 * @return The binomial value.
 *
 * This function returns the binomial of (n,p). No test is done on @p n
 * and @p p. They must be greater or equal than 0.
 */
SAPI int64_t sw_binomial(uint32_t n, uint32_t p);

/**
 * @brief Return the current time.
 *
 * @return The current time, in seconds
 *
 * This function returns the current time in seconds, as a double.
 */
SAPI double sw_time_get(void);


#endif /* SW_UTILS_H */
