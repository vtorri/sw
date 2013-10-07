#ifndef SW_SCALE_FCT_H
#define SW_SCALE_FCT_H


/**
 * @file sw_scale_fct.h
 * @brief Scale functions routines
 *
 * These functions implements routines to manage scale functions
 * (creation, destruction, retrieval of the characteristics).
 */

#include <stdlib.h>

#include "algebra/sw_rational.h"
#include "sw_forward.h"
#include "sw_spline.h"
#include "sw_weights.h"


/* scale function - base class methods */

#include "sw_scale_fct_base.x"

/* scale function methods */

/**
 * @brief Create a new scale function of the given order.
 *
 * @param order The order of the spline.
 * @return The new scale function object.
 *
 * This function creates a scale function object. As first part of the
 * bi-orthogonal wavelet, it is mainly a spline function. The order of
 * the spline is given by @p order. If @p order is non positive, then
 * @c NULL is returned. The filter associated to that scale function
 * is retrieved with the base class methods. When the scale function
 * is not used anymore, its memory must be freed with
 * scale_fct_del().
 */
SAPI sw_scale_fct_t *sw_scale_fct_new(int32_t order);

/**
 * @brief Free the memory of the given scale function.
 *
 * @param sf The scale function to free.
 *
 * This function frees the memory of the scale function @p sf. @p sf
 * must have been created with sw_scale_fct_new(). If @p sf is @c NULL,
 * this function does nothing.
 */
SAPI void sw_scale_fct_del(sw_scale_fct_t *sf);

/**
 * @brief Set the type of weights to use for the given scale function.
 *
 * @param sf The scale function.
 * @param type The type of weights.
 * @return 1 on success, 0 otherwise.
 *
 * This function sets the type of weights @p type for the scale function
 * @p sf. If @p sf is @c NULL or @p type is invalid, 0 is returned,
 * otherwise, 1 is returned.
 */
SAPI uint8_t sw_scale_fct_type_set(sw_scale_fct_t   *sf,
				   sw_weights_type_t type);

/**
 * @brief Get the type of weights for the given scale function.
 *
 * @param sf The scale function.
 * @return The type of weights.
 *
 * This function gets the type of weights for the scale function
 * @p sf. If @p sf is @c NULL, #SW_WEIGHTS_TYPE_ERROR is returned,
 * otherwise, the type of weights is returned.
 */
SAPI sw_weights_type_t sw_scale_fct_type_get(const sw_scale_fct_t *sf);

/**
 * @brief Return the spline function associated to the given scale function.
 *
 * @param sf The scale function.
 * @return The spline function.
 *
 * This function returns the spline function associated to @p sf. if
 * @p sf is @c NULL, this function returns @c NULL.
 */
SAPI const sw_spline_t *sw_scale_fct_spline_get(const sw_scale_fct_t *sf);

/**
 * @brief Return the rational value of the given scale function at the given
 * rational absciss.
 *
 * @param sf The scale function.
 * @param val The rational absciss.
 * @return The rational value of the scale function at the given absciss.
 *
 * This function returns the rational value of the scale function @p sf at
 * the rational absciss @p val.
 */
SAPI rational_t sw_scale_fct_value_rat_get(const sw_scale_fct_t *sf,
					   const rational_t     *val);

/**
 * @brief Return the real value of the given scale function at the given real
 * absciss.
 *
 * @param sf The scale function.
 * @param val The real absciss.
 * @return The real value of the scale function at the given absciss.
 *
 * This function return the rational value of the scale function @p sf at
 * the rational absciss @p val.
 */
SAPI double sw_scale_fct_value_get(const sw_scale_fct_t *sf,
				   double                val);

/**
 * @brief Apply the backward projection of a function onto the vector space of
 * fine details using periodic functions.
 *
 * @param sf The scale function.
 * @param scale The scale of the details.
 * @param ps The inner products.
 * @param function The function to get.
 *
 * This function applies the backward periodic projection on the
 * vector space of the fine details whose scale function is @p sf. The
 * scale of the details is given by @p scale. The inner products used
 * to build the function are given by @p ps. @p ps must be an array of
 * size @c 2 power @p scale. The result is stored in @p function, which also
 * must be an array of size @c 2 power @p scale.
 */
SAPI void sw_scale_fct_proj_periodic_backward(const sw_scale_fct_t *sf,
					      int32_t               scale,
					      const double         *ps,
					      double               *function);

/**
 * @brief Apply the backward projection of a function onto the vector space of
 * fine details using dirichlet boundaries.
 *
 * @param sf The scale function.
 * @param scale The scale of the details.
 * @param ps The inner products.
 * @param function The function to get.
 *
 * This function applies the backward dirichlet projection on the
 * vector space of the fine details whose scale function is @p sf. The
 * scale of the details is given by @p scale. The inner products used
 * to build the function are given by @p ps. @p ps must be an array of
 * size @c 2 power(@p scale + 1). The result is stored in @p function,
 * which also must be an array of size @c 2 power(@p scale + 1).
 */
SAPI void sw_scale_fct_proj_dirichlet_backward(const sw_scale_fct_t *sf,
					       int32_t               scale,
					       const double         *ps,
					       double               *function);

/* dual scale function methods */

/**
 * @brief Create a new dual scale function of the given orders.
 *
 * @param order The order of the scale function.
 * @param order_dual The dual order of the scale function.
 * @return The new dual scale function object.
 *
 * This function creates a dual scale function object. The orders of
 * that scale functions are given by @p order and @p order_dual. If
 * @p order or @p order_dual are non positive, then @c NULL is
 * returned. The filter associated to that scale function is retrieved
 * with the base class methods. When the scale function is not used
 * anymore, its memory must be freed with scale_fct_dual_del().
 */
SAPI sw_scale_fct_dual_t *sw_scale_fct_dual_new(int32_t        order,
						int32_t        order_dual);

/**
 * @brief Free the memory of the given dual scale function.
 *
 * @param sfd The dual scale function to free.
 *
 * This function frees the memory of the dual scale function @p sfd. @p sfd
 * must have been created with sw_scale_fct_dual_new(). If @p sfd is @c NULL,
 * this function does nothing.
 */
SAPI void sw_scale_fct_dual_del(sw_scale_fct_dual_t *sfd);

/**
 * @brief Set the type of weights to use for the given dual scale function.
 *
 * @param sf The dual scale function.
 * @param type The type of weights.
 * @return 1 on success, 0 otherwise.
 *
 * This function sets the type of weights @p type for the dual scale function
 * @p sf. If @p sf is @c NULL or @p type is invalid, 0 is returned,
 * otherwise, 1 is returned.
 *
 * @todo comment va_arg
 */
SAPI uint8_t sw_scale_fct_dual_type_set(sw_scale_fct_dual_t   *sf,
					sw_weights_type_t      type,
					...);

/**
 * @brief Get the type of weights for the given dual scale function.
 *
 * @param sf The dual scale function.
 * @return The type of weights.
 *
 * This function gets the type of weights for the dual scale function
 * @p sf. If @p sf is @c NULL, #SW_WEIGHTS_TYPE_ERROR is returned,
 * otherwise, the type of weights is returned.
 */
SAPI sw_weights_type_t sw_scale_fct_dual_type_get(const sw_scale_fct_dual_t *sf);

/**
 * @brief Apply the forward projection of a function onto the vector space of
 * fine details using periodic functions.
 *
 * @param sfd The scale function.
 * @param scale The scale of the details.
 * @param function The function to analyse.
 * @param ps The inner products to get.
 *
 * This function applies the forward periodic projection on the vector
 * space of the fine details whose dual scale function is @p sfd. The
 * scale of the details is given by @p scale. The function to analyse
 * is given by @p ps. @p function must be an array of size @c 2 power
 * @p scale. The computed inner products are stored in @p ps, which
 * also must be an array of size @c 2 power @p scale.
 */
SAPI void sw_scale_fct_dual_proj_periodic_forward(const sw_scale_fct_dual_t *sfd,
						  int32_t                    scale,
						  const double              *function,
						  double                    *ps);

/**
 * @brief Apply the forward projection of a function onto the vector space of
 * fine details using dirichlet boundaries.
 *
 * @param sfd The scale function.
 * @param scale The scale of the details.
 * @param function The function to analyse.
 * @param ps The inner products to get.
 *
 * This function applies the forward dirichlet projection on the vector
 * space of the fine details whose dual scale function is @p sfd. The
 * scale of the details is given by @p scale. The function to analyse
 * is given by @p ps. @p function must be an array of size @c 2 power
 * (@p scale + 1). The computed inner products are stored in @p ps,
 * which also must be an array of size @c 2 power(@p scale + 1).
 */
SAPI void sw_scale_fct_dual_proj_dirichlet_forward(const sw_scale_fct_dual_t *sfd,
						   int32_t                    scale,
						   const double              *function,
						   double                    *ps);

/**
 * @brief Return the data of the given weights, according to the type
 * of quadrature.
 *
 * @param sfd The dual scale function to et the weights from.
 * @param size The size of the weights array.
 * @return The array of weights.
 *
 * This function retrieve the data needed by the computation of the
 * quadrature formula: it returns the array of weights stored in @p sfd
 * and store the size of that array in @p size. If @p sfd or
 * @p weights_size are @c NULL, @c NULL is returned. The returned
 * array must not be freed.
 *
 * Note that in the Lagrange case,  @p weights_size is the degree of
 * the Lagrange polynomial + 1.
 */
SAPI const double *sw_scale_fct_dual_weights_get(const sw_scale_fct_dual_t *sfd,
						 int32_t                   *size);


#endif /* SW_SCALE_FCT_H */
