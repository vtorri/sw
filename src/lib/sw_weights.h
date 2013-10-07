#ifndef SW_WEIGHTS_H
#define SW_WEIGHTS_H


/**
 * @file sw_weights.h
 * @brief Weights routines
 *
 * These functions implements routines to manage weights use in the
 * quadrature formulas. Two kind of quadrature formulas are
 * implemented : one using Lagrange polynomials, the other using the
 * Sweldens one.
 */


#include "sw_forward.h"
#include "sw_scale_fct.h"


/**
 * @typedef sw_weights_t
 * @brief Opaque type for weights.
 */
typedef struct sw_weights_s sw_weights_t;


/**
 * @brief Create a new object that manage weights of the given type
 * and given scale function.
 *
 * @param type The type of quadrature formula.
 * @param sfb The scale function (base class).
 * @return The new weights object.
 *
 * This function creates an object managing the weights of a
 * quadrature formula. The type of the quadrature formula is given by
 * @p type. The scale function used in the quadrature formula is set
 * by @p sfb. It is the base class of a scale function. If @p sfb or
 * an error in memory allocation, @c NULL is returned. Other, a new
 * weights object is returned. When the weights are not used anymore,
 * free its memory with sw_weights_del().
 */
SAPI sw_weights_t *sw_weights_new(sw_weights_type_t          type,
				  const sw_scale_fct_base_t *sfb);

/**
 * @brief Free the memory of the given weights.
 *
 * @param w The weights to free.
 *
 * This function frees the memory of @p w. @p w must have
 * been created with sw_weights_new(). If @p w is @c NULL, this
 * function does nothing.
 */
SAPI void sw_weights_del(sw_weights_t *w);

/**
 * @brief Set the data of the given weights in the Lagrange case.
 *
 * @param w The weights to set.
 * @param degree The degree of the Lagrange polynomial to use.
 *
 * This function set the Lagrange data to the weights @p w. In that
 * case, the data is just the degree @p degree of the Lagrange
 * polynomial used i the quadrature formula. If @p w is @c NULL or if
 * the type of @p w is not @ref SW_WEIGHTS_TYPE_LAGRANGE, then this
 * function does nothing. Otherwise, it set the internal data of the
 * weights object and compute the weights used for the quadrature
 * formula.
 */
SAPI void sw_weights_lagrange_data_set(sw_weights_t *w,
				       int32_t       degree);

/**
 * @brief Set the data of the given weights using the Sweldens algorithm.
 *
 * @param w The weights to set.
 * @param order The order.
 * @param scale The scale.
 * @param tau The offset.
 *
 * @todo Fix the doc.
 *
 * This function set the Lagrange data to the weights @p w. In that
 * case, the data is just the degree @p degree of the Lagrange
 * polynomial used i the quadrature formula. If @p w is @c NULL or if
 * the type of @p w is not @ref SW_WEIGHTS_TYPE_LAGRANGE, then this
 * function does nothing. Otherwise, it set the internal data of the
 * weights object and compute the weights used for the quadrature
 * formula.
 */
SAPI void sw_weights_sweldens_data_set(sw_weights_t     *w,
				       int32_t           order,
				       int32_t           scale,
				       const rational_t *tau);


SAPI sw_weights_t *
sw_weights_lagrange_new(const sw_scale_fct_base_t *sfb, int32_t degree);

SAPI sw_weights_t *
sw_weights_sweldens_new(const sw_scale_fct_base_t *sfb, int32_t r, int32_t s);

/**
 * @brief Return the data of the given weights, according to the type
 * of quadrature.
 *
 * @param w The weights to retrieve the data from.
 * @param weights_size The size of the weights array.
 * @return The array of weights.
 *
 * This function retrieve the data needed by the computation of the
 * quadrature formula: it returns the array of weights stored in @p w
 * and store the size of that array in @p weights_size. If @p w or
 * @p weights_size are @c NULL, @c NULL is returned. The returned
 * array must not be freed.
 *
 * Note that in the Lagrange case,  @p weights_size is the degree of
 * the Lagrange polynomial + 1.
 */
SAPI const double *sw_weights_get(const sw_weights_t *w,
				  int32_t            *weights_size);


#endif /* SW_WEIGHTS_H */
