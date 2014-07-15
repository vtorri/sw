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
