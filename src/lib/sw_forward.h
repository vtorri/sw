#ifndef SW_FORWARD_H
#define SW_FORWARD_H


/**
 * @file sw_forward.h
 * @brief Forward declarations.
 *
 * This file contains the forward declarations of the weights and the
 * scale functions.
 */


/**
 * @typedef weights_type_t
 * @brief Enumeration of the type of quadrature formulas.
 */
typedef enum
{
    SW_WEIGHTS_TYPE_ERROR,    /**< Error */
    SW_WEIGHTS_TYPE_LAGRANGE, /**< Lagrange approximation */
    SW_WEIGHTS_TYPE_SWELDENS  /**< Sweldens approximation */
} weights_type_t;

/**
 * @typedef scale_fct_base_t
 * @brief Opaque type for a base class of a scale function.
 */
typedef struct scale_fct_base scale_fct_base_t;

/**
 * @typedef scale_fct_t
 * @brief Opaque type for a scale function.
 */
typedef struct scale_fct scale_fct_t;

/**
 * @typedef scale_fct_dual_t
 * @brief Opaque type for a dual scale function.
 */
typedef struct scale_fct_dual scale_fct_dual_t;


#endif /* SW_FORWARD_H */
