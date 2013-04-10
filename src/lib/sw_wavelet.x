
/**
 * @struct wavelet_base
 * @brief Type for a base class of a wavelet function.
 */
struct wavelet_base
{
    rational_t *filter_rat; /**< filter with rational values */
    double     *filter;     /**< filter with real values */
    int32_t     N1;         /**< inferior bound of the support of the filter */
    int32_t     N2;         /**< superior bound of the support of the filter */
};

/**
 * @brief Return the inferior bound of the support of the wavelet filter.
 *
 * @param wb The wavelet function.
 * @return The value of the inferior bound of the support of the wavelet filter.
 *
 * This function returns the value of the inferior bound of the
 * support of the filter of the wavelet function @p wb. @p wb must be a
 * wavelet function created by wavelet_new() or wavelet_dual_new().
 */
static __inline__ int32_t
wavelet_base_N1_get(const wavelet_base_t *wb)
{
    return wb->N1;
}

/**
 * @brief Return the superior bound of the support of the wavelet filter.
 *
 * @param wb The wavelet function.
 * @return The value of the superior bound of the support of the wavelet filter.
 *
 * This function returns the value of the superior bound of the
 * support of the filter of the wavelet function @p wb. @p wb must be a
 * wavelet function created by wavelet_new() or wavelet_dual_new().
 */
static __inline__ int32_t
wavelet_base_N2_get(const wavelet_base_t *wb)
{
    return wb->N2;
}

/**
 * @brief Return the rational values of the wavelet filter.
 *
 * @param wb The scale function.
 * @return The rational values of the filter.
 *
 * This function returns the rational values of the filter of the
 * wavelet function @p wb. The size of the array is (order + (order_dual * 2) - 1) if @p wb
 * is a wavelet function and is (order + 1) if @p wb
 * is a dual wavelet function. In both cases, it can be computed with
 * the values returned by wavelet_base_N1_get() and
 * wavelet_base_N2_get(). @p wb must be a wavelet function created by
 * wavelet_new() or wavelet_dual_new().
 */
static __inline__ const rational_t *
wavelet_base_filter_rat_get(const wavelet_base_t *wb)
{
    return wb->filter_rat;
}

/**
 * @brief Return the real values of the wavelet filter.
 *
 * @param wb The scale function.
 * @return The real values of the filter.
 *
 * This function returns the real values of the filter of the
 * wavelet function @p wb. The size of the array is (order + (order_dual * 2) - 1) if @p wb
 * is a wavelet function and is (order + 1) if @p wb
 * is a dual wavelet function. In both cases, it can be computed with
 * the values returned by wavelet_base_N1_get() and
 * wavelet_base_N2_get(). @p wb must be a wavelet function created by
 * wavelet_new() or wavelet_dual_new().
 */
static __inline__ const double *
wavelet_base_filter_get(const wavelet_base_t *wb)
{
    return wb->filter;
}
