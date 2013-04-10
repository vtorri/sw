#ifndef SW_SWELDENS_H
#define SW_SWELDENS_H


#include "sw_scale_fct.h"


typedef struct sweldens sweldens_t;


SAPI sweldens_t *sweldens_new(const scale_fct_base_t *sfb,
                              int32_t                 order,
                              int32_t                 scale,
                              const rational_t       *tau);

SAPI void        sweldens_delete(sweldens_t *w);

SAPI int32_t     sweldens_order_get(const sweldens_t *s);

SAPI rational_t *sweldens_absciss_get(const sweldens_t *s);

SAPI rational_t *sweldens_weights_get(const sweldens_t *s);


#endif /* SW_SWELDENS_H */
