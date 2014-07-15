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

/**
 * @struct sw_scale_fct_base
 * @brief Type for a base class of a scale function.
 */
struct sw_scale_fct_base_s
{
    rational_t *filter_rat; /**< filter with rational values */
    double     *filter;     /**< filter with real values */
    int32_t     N1;         /**< inferior bound of the support of the filter */
    int32_t     N2;         /**< superior bound of the support of the filter */
};

/**
 * @brief Return the inferior bound of the support of the scale function filter.
 *
 * @param sf The scale function.
 * @return The value of the inferior bound of the support of the filter.
 *
 * This function returns the value of the inferior bound of the
 * support of the filter of the scale function @p sf. @p sf must be a
 * scale function created by sw_scale_fct_new() or sw_scale_fct_dual_new().
 */
static __inline__ int32_t
sw_scale_fct_base_N1_get(const sw_scale_fct_base_t *sf)
{
    if (sf)
        return sf->N1;

    return 0;
}

/**
 * @brief Return the superior bound of the support of the scale function filter.
 *
 * @param sf The scale function.
 * @return The value of the superior bound of the support of the filter.
 *
 * This function returns the value of the superior bound of the
 * support of the filter of the scale function @p sf. @p sf must be a
 * scale function created by sw_scale_fct_new() or sw_scale_fct_dual_new().
 */
static __inline__ int32_t
sw_scale_fct_base_N2_get(const sw_scale_fct_base_t *sf)
{
    if (sf)
        return sf->N2;

    return 0;
}

/**
 * @brief Return the rational values of the scale function filter.
 *
 * @param sf The scale function.
 * @return The rational values of the filter.
 *
 * This function returns the rational values of the filter of the
 * scale function @p sf. The size of the array is (order + 1) if @p sf
 * is a scale function and is (order + (order_dual * 2) - 1) if @p sf
 * is a dual scale function. In both cases, it can be computed with
 * the values returned by sw_scale_fct_base_N1_get() and
 * sw_scale_fct_base_N2_get(). @p sf must be a scale function created by
 * sw_scale_fct_new() or sw_scale_fct_dual_new().
 */
static __inline__ const rational_t *
sw_scale_fct_base_filter_rat_get(const sw_scale_fct_base_t *sf)
{
    if (sf)
        return sf->filter_rat;

    return NULL;
}

/**
 * @brief Return the real values of the filter.
 *
 * @param sf The scale function.
 * @return The real values of the filter.
 *
 * This function returns the real values of the filter of the
 * scale function @p sf. The size of the array is (order + 1) if @p sf
 * is a scale function and is (order + (order_dual * 2) - 1) if @p sf
 * is a dual scale function. @p sf must be a scale function created by
 * sw_scale_fct_new() or sw_scale_fct_dual_new().
 */
static __inline__ const double *
sw_scale_fct_base_filter_get(const sw_scale_fct_base_t *sf)
{
    if (sf)
        return sf->filter;

    return NULL;
}
