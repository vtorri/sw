
static __inline__ int32_t
mra_size_x_get(const mra_t *mra)
{
    if (mra)
        return mra->size_x;

    return 0;
}

static __inline__ int32_t
mra_size_v_get(const mra_t *mra)
{
    if (mra)
        return mra->size_v;

    return 0;
}

static __inline__ int32_t
mra_size_inf_x_get(const mra_t *mra)
{
    if (mra)
        return mra->lambda_fine_x_inf;

    return 0;
}

static __inline__ int32_t
mra_size_inf_v_get(const mra_t *mra)
{
    if (mra)
        return mra->lambda_fine_v_inf;

    return 0;
}

static __inline__ int32_t
mra_size_sup_x_get(const mra_t *mra)
{
    if (mra)
        return mra->lambda_fine_x_sup;

    return 0;
}

static __inline__ int32_t
mra_size_sup_v_get(const mra_t *mra)
{
    if (mra)
        return mra->lambda_fine_v_sup;

    return 0;
}
