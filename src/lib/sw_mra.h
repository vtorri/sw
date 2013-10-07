#ifndef SW_MRA_H
#define SW_MRA_H


/**
 * @file sw_mra.h
 * @brief Multi-Resolution Analysis routines
 *
 * These functions implements routines to manage multi-resolution
 * analysis (MRA) functions (creation, destruction, retrieval of the
 * characteristics, projection forward and backward).
 */


#include "sw_scale_fct.h"
#include "sw_wavelet.h"


/**
 * @typedef sw_mra_t
 * @brief Opaque type for a MRA.
 */
typedef struct sw_mra_s sw_mra_t;

/**
 * @struct sw_mra_s
 * @brief Type for a MRA.
 */
struct sw_mra_s
{
    sw_scale_fct_t      *scale_fct;          /**< scale function */
    sw_scale_fct_dual_t *scale_fct_dual;     /**< dual scale function */
    sw_wavelet_t        *wavelet;            /**< wavelet function */
    sw_wavelet_dual_t   *wavelet_dual;       /**< dual wavelet function */

    int32_t              scale_coarse;       /**< coarse scale */
    int32_t              scale_fine;         /**< fine scale */

    int32_t              lambda_fine_x_inf;  /**< inf bound of the support of the function in the x direction */
    int32_t              lambda_fine_x_sup;  /**< sup bound of the support of the function in the x direction */
    int32_t              size_x;             /**< number of points taken in the x direction */

    int32_t              lambda_fine_v_inf;  /**< inf bound of the support of the function in the v direction */
    int32_t              lambda_fine_v_sup;  /**< sup bound of the support of the function in the v direction */
    int32_t              size_v;             /**< number of points taken in the v direction */
};


#include "sw_mra.x"


/**
 * @brief Create a new mra of the given orders and scales.
 *
 * @param order The order of the spline function.
 * @param order_dual The dual order of the dual spline function.
 * @param scale_coarse The coarse scale used in the MRA.
 * @param scale_fine The fine scale used in the MRA.
 * @param type The type of the quadrature formula to une in the projection algorithm.
 * @return The new MRA object.
 *
 * This function creates a MRA object of orders @p order and
 * @p order_dual and returns a pointer to a newly created MRA. The
 * coarse scale and fine scale of the MRA are respectively set by
 * @p scale_coarse and @p scale_fine. @p type sets the type of the
 * algorithm used in sw_mra_proj_x_forward(), sw_mra_proj_v_forward(),
 * sw_mra_proj_x_backward() and sw_mra_proj_v_backward(). If an error
 * occurs, @c NULL is returned. When this MRA is not used anymore,
 * free its memory with sw_mra_del().
 *
 * When @p type is @ref SW_WEIGHTS_TYPE_LAGRANGE, the quadrature
 * formula use a Lagrange interpolation. Hence the next parameter to
 * give is the degree of the Lagrange polynomial to use.
 */
SAPI sw_mra_t *sw_mra_new(int32_t order,
			  int32_t order_dual,
			  int32_t scale_coarse,
			  int32_t scale_fine,
			  sw_weights_type_t type,
			  ...);

/**
 * @brief Free the memory of the given MRA.
 *
 * @param mra The MRA to free.
 *
 * This function frees the memory of @p mra. @p mra must have
 * been created with sw_mra_new(). If @p mra is @c NULL, this
 * function does nothing.
 */
SAPI void sw_mra_del(sw_mra_t *mra);

/**
 * @brief Return the scale function associated to an MRA.
 *
 * @param mra The MRA.
 * @return The scale function.
 *
 * This function returns the scale function associated to @p mra. If
 * @p mra is @c NULL, then @c NULL is returned. The returned value
 * must not be freed.
 */
SAPI const sw_scale_fct_t *sw_mra_scale_fct_get(const sw_mra_t *mra);

/**
 * @brief Return the dual scale function associated to an MRA.
 *
 * @param mra The MRA.
 * @return The dual scale function.
 *
 * This function returns the dual scale function associated to
 * @p mra. If @p mra is @c NULL, then @c NULL is returned. The
 * returned value must not be freed.
 */
SAPI const sw_scale_fct_dual_t *sw_mra_scale_fct_dual_get(const sw_mra_t *mra);

/**
 * @brief Return the size in the X direction at the fine scale.
 *
 * @param mra The MRA.
 * @return The size in the X direction at the fine scale.
 *
 * This function returns the size in the X direction at the fine
 * scale, that is 2 power(the fine scale). If @p mra is NULL, @c 0 is
 * returned.
 */
static __inline__ int32_t sw_mra_size_x_get(const sw_mra_t *mra);

/**
 * @brief Return the size in the V direction at the fine scale.
 *
 * @param mra The MRA.
 * @return The size in the V direction at the fine scale.
 *
 * This function returns the size in the V direction at the fine
 * scale, that is 2 power(the fine scale + 1). If @p mra is NULL,
 * @c 0 is returned.
 */
static __inline__ int32_t sw_mra_size_v_get(const sw_mra_t *mra);

/**
 * @brief Return the value of the minimal bound of arrays in the X direction.
 *
 * @param mra The MRA.
 * @return The value of the minimal bound of arrays in the X direction.
 *
 * This function returns the minimal bound of arrays in the X
 * direction. Actually, it returns 0 !
 */
static __inline__ int32_t sw_mra_size_inf_x_get(const sw_mra_t *mra);

/**
 * @brief Return the value of the minimal bound of arrays in the V direction.
 *
 * @param mra The MRA.
 * @return The value of the minimal bound of arrays in the V direction.
 *
 * This function returns the minimal bound of arrays in the V
 * direction. Actually, it returns - @c 2 power the fine scale !
 */
static __inline__ int32_t sw_mra_size_inf_v_get(const sw_mra_t *mra);

/**
 * @brief Return the value of the maximal bound of arrays in the X direction.
 *
 * @param mra The MRA.
 * @return The value of the maximal bound of arrays in the X direction.
 *
 * This function returns the maximal bound of arrays in the X
 * direction. Actually, it returns @c 2 power the fine scale !
 */
static __inline__ int32_t sw_mra_size_sup_x_get(const sw_mra_t *mra);

/**
 * @brief Return the value of the maximal bound of arrays in the V direction.
 *
 * @param mra The MRA.
 * @return The value of the maximal bound of arrays in the V direction.
 *
 * This function returns the maximal bound of arrays in the V
 * direction. Actually, it returns @c 2 power the fine scale !
 */
static __inline__ int32_t sw_mra_size_sup_v_get(const sw_mra_t *mra);

/**
 * @brief Apply the forward projection on the fine scale in the X direction.
 *
 * @param mra The MRA.
 * @param function The function to analyse.
 * @param ps The inner products to get.
 *
 * This function applies the forward projection on a function in the X
 * direction of @p mra, using the quadrature formula set in
 * sw_mra_new(). The function to analyse is @p function and must be a
 * buffer of size @c 2 power the fine scale set in sw_mra_new(). The
 * inner products are stored in @p ps, which must also be a buffer of
 * the same size.
 */
SAPI void sw_mra_proj_x_forward(const sw_mra_t  *mra,
				const double    *function,
				double          *ps);

/**
 * @brief Apply the backward projection on the fine scale in the X direction.
 *
 * @param mra The MRA.
 * @param ps The inner products.
 * @param function The function to get.
 *
 * This function applies the backward projection on inner products in
 * the X direction of @p mra, using the quadrature formula set in
 * sw_mra_new(). The inner products used to build the function are stored
 * in @p ps, which must be a buffer of size @c 2 power the fine scale
 * set in sw_mra_new(). The function computed from those inner products
 * is then stored in @p function, which must also be a buffer of the
 * same size.
 */
SAPI void sw_mra_proj_x_backward(const sw_mra_t  *mra,
				 const double    *ps,
				 double          *function);

/**
 * @brief Apply the forward projection on the fine scale in the V direction.
 *
 * @param mra The MRA.
 * @param function The function to analyse.
 * @param ps The inner products to get.
 *
 * This function applies the forward projection on a function in the V
 * direction of @p mra, using the quadrature formula set in
 * sw_mra_new(). The function to analyse is @p function and must be a
 * buffer of size @c 2 power the fine scale + 1 set in sw_mra_new(). The
 * inner products are stored in @p ps, which must also be a buffer of
 * the same size.
 */
SAPI void sw_mra_proj_v_forward(const sw_mra_t *mra,
				const double   *function,
				double         *ps);

/**
 * @brief Apply the backward projection on the fine scale in the V direction.
 *
 * @param mra The MRA.
 * @param ps The inner products.
 * @param function The function to get.
 *
 * This function applies the backward projection on inner products in
 * the V direction of @p mra, using the quadrature formula set in
 * sw_mra_new(). The inner products used to build the function are stored
 * in @p ps, which must be a buffer of size @c 2 power the fine scale
 * + 1 set in sw_mra_new(). The function computed from those inner
 * products is then stored in @p function, which must also be a buffer
 * of the same size.
 */
SAPI void sw_mra_proj_v_backward(const sw_mra_t *mra,
				 const double   *ps,
				 double         *function);

SAPI void sw_mra_proj_2d_x_forward(const sw_mra_t *mra,
				   const double   *function,
				   double         *ps,
				   double         *tmp1,
				   double         *tmp2);

SAPI void sw_mra_proj_2d_x_backward(const sw_mra_t *mra,
				    const double   *ps,
				    double         *function,
				    double         *tmp1,
				    double         *tmp2);

SAPI void sw_mra_proj_2d_v_forward(const sw_mra_t *mra,
				   const double   *function,
				   double         *ps);

SAPI void sw_mra_proj_2d_v_backward(const sw_mra_t *mra,
				    const double   *ps,
				    double         *function);

/**
 * @brief Apply the forward projection on the fine scale in the V direction.
 *
 * @param mra The MRA.
 * @param function The function to analyse.
 * @param ps The inner products to get.
 * @param tmp1 A buffer for temporary computations.
 * @param tmp2 A buffer for temporary computations.
 *
 * This function applies the forward projection on a function in 2
 * dimensions, using the quadrature formula set in sw_mra_new(). The
 * function to analyse is @p function and must be a buffer of size_x *
 * size_v, where size_x is equal to @c 2 power the fine scale set in
 * sw_mra_new() and size_v is equal to @c 2 power the fine scale + 1. The
 * inner products are stored in @p ps, which must also be a buffer of
 * the same size than @p function. @p tmp1 and @p tmp2 are both
 * pointers to buffers of size size_x and are used for internal
 * computations.
 */
SAPI void sw_mra_proj_2d_forward(const sw_mra_t *mra,
				 const double   *function,
				 double         *ps,
				 double         *tmp1,
				 double         *tmp2);

/**
 * @brief Apply the backward projection on the fine scale in the V direction.
 *
 * @param mra The MRA.
 * @param ps The inner products.
 * @param function The function to get.
 * @param tmp1 A buffer for temporary computations.
 * @param tmp2 A buffer for temporary computations.
 *
 * This function applies the backward projection on a function in 2
 * dimensions, using the quadrature formula set in sw_mra_new(). The
 * inner products used to build the function are stored in @p ps and
 * must be a buffer of size_x * size_v, where size_x is equal to @c 2
 * power the fine scale set in sw_mra_new() and size_v is equal to @c 2
 * power the fine scale + 1. The function computed from those inner
 * products is then stored in @p function must also be a buffer of the
 * same size than @p ps. @p tmp1 and @p tmp2 are both pointers to
 * buffers of size size_x and are used for internal computations.
 */
SAPI void sw_mra_proj_2d_backward(const sw_mra_t *mra,
				  const double   *ps,
				  double         *function,
				  double         *tmp1,
				  double         *tmp2);

/* double *sw_mra_scale_fct_ps_get(const sw_mra_t *mra, double coef); */

SAPI void sw_mra_advection_v(const sw_mra_t *mra,
			     double          coef,
			     const double   *ps_in,
			     double         *ps_out);


#endif /* SW_MRA_H */
