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

#ifndef SW_SPLINE_H
#define SW_SPLINE_H


/**
 * @file sw_spline.h
 * @brief Spline routines
 *
 * These functions implements routines to manage spline functions
 * (creation, destruction, retrieval of the characteristics, ...).
 */


#include <algebra/sw_rational.h>


/**
 * @typedef sw_spline_t
 * @brief Opaque type for a spline function.
 */
typedef struct spline_s sw_spline_t;

/**
 * @brief Create a new spline function of the given order.
 *
 * @param order The order of the spline.
 * @return The new spline object.
 *
 * This function creates a spline object of order @p order and returns
 * a pointer to a newly created spline. If an error occurs, @c NULL is
 * returned. When this spline is not used anymore, free its memory
 * with sw_spline_del().
 */
SAPI sw_spline_t *sw_spline_new(int32_t order);

/**
 * @brief Free the memory of the given spline.
 *
 * @param spline The spline to free.
 *
 * This function frees the memory of @p spline. @p spline must have
 * been created with sw_spline_new(). If @p spline is @c NULL, this
 * function does nothing.
 */
SAPI void sw_spline_del(sw_spline_t *spline);

/**
 * @brief Return the order of the given spline function.
 *
 * @param spline The spline.
 * @return The order of the spline function
 *
 * This function returns the order of the spline @p spline. If
 * @p spline is @c NULL, this function returns 0.
 */
SAPI int32_t sw_spline_order_get(const sw_spline_t *spline);

/**
 * @brief Return the minimal value of the support of the given
 * spline.
 *
 * @param spline The spline.
 * @return The minimal value of the support of the spline.
 *
 * This function returns the minimal value of the support of the
 * spline @p spline. If @p spline is @c NULL, this function returns
 * 0.
 */
SAPI int64_t sw_spline_x1_get(const sw_spline_t *spline);

/**
 * @brief Return the maximal value of the support of the given
 * spline.
 *
 * @param spline The spline.
 * @return The maximal value of the support of the spline.
 *
 * This function returns the maximal value of the support of the
 * spline @p spline. If @p spline is @c NULL, this function returns
 * 0.
 */
SAPI int64_t sw_spline_x2_get(const sw_spline_t *spline);

/**
 * @brief Return the rational coefficients of the polynoms of the
 * given spline.
 *
 * @param spline The spline to get the coefficients.
 * @return The rational coefficients of the polynoms of the given spline.
 *
 * This function returns the rational coefficients of the polynoms of
 * the spline @p spline. If @p spline is @c NULL, this function
 * returns @c NULL.
 */
SAPI const rational_t *sw_spline_rat_coef_get(const sw_spline_t *spline);

/**
 * @brief Return the real coefficients of the polynoms of the given
 * spline.
 *
 * @param spline The spline to get the coefficients.
 * @return The real coefficients of the polynoms of the given spline.
 *
 * This function returns the real coefficients of the polynoms of the
 * spline @p spline. If @p spline is @c NULL, this function returns
 * @c NULL.
 */
SAPI const double *sw_spline_coef_get(const sw_spline_t *spline);

/**
 * @brief Return the rational value of the given spline at the given
 * rational absciss.
 *
 * @param spline The spline.
 * @param val The rational absciss.
 * @return The rational value of the spline at the given absciss.
 *
 * This function return the rational value of the spline @p spline at
 * the rational absciss @p val.
 */
SAPI rational_t sw_spline_value_rat_get(const sw_spline_t *spline, const rational_t *val);

/**
 * @brief Return the real value of the given spline at the given real
 * absciss.
 *
 * @param spline The spline.
 * @param val The real absciss.
 * @return The real value of the spline at the given absciss.
 *
 * This function return the rational value of the spline @p spline at
 * the rational absciss @p val.
 */
SAPI double sw_spline_value_get(const sw_spline_t *spline, double val);

/**
 * @brief Compute the value of the integral of the given spline
 * function between two points.
 *
 * @param spline The spline function.
 * @param x1 The lower bound of the integration interval.
 * @param x2 The upper bound of the integration interval.
 * @return The value of the integral.
 *
 * This function computes the value of the integral of the spline
 * function @p spline between @p x1 and @p x2.
 */
SAPI double sw_spline_integral_value_get(const sw_spline_t *spline, double x1, double x2);

/**
 * @brief Compute the rational value of the integral of the given spline
 * function between two points.
 *
 * @param spline The spline function.
 * @param x1 The lower bound of the integration interval.
 * @param x2 The upper bound of the integration interval.
 * @return The rational value of the integral.
 *
 * This function computes the rational value of the integral of the spline
 * function @p spline between @p x1 and @p x2.
 */
SAPI rational_t sw_spline_integral_value_rat_get(const sw_spline_t *spline, int x1, int x2);


#endif /* SW_SPLINE_H */
