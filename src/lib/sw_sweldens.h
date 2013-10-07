#ifndef SW_SWELDENS_H
#define SW_SWELDENS_H


#include "sw_scale_fct.h"


typedef struct sw_sweldens_s sw_sweldens_t;


SAPI sw_sweldens_t *sw_sweldens_new(const sw_scale_fct_base_t *sfb,
				    int32_t                    order,
				    int32_t                    scale,
				    const rational_t          *tau);

SAPI void        sw_sweldens_delete(sw_sweldens_t *w);

SAPI int32_t     sw_sweldens_order_get(const sw_sweldens_t *s);

SAPI rational_t *sw_sweldens_absciss_get(const sw_sweldens_t *s);

SAPI rational_t *sw_sweldens_weights_get(const sw_sweldens_t *s);


#endif /* SW_SWELDENS_H */
