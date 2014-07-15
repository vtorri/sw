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

#ifndef SW_WAVELET_H
#define SW_WAVELET_H


/**
 * @file sw_wavelet.h
 * @brief Wavelet functions routines
 *
 * These functions implements routines to manage wavelet functions
 * (creation, destruction, retrieval of the characteristics).
 */


#include "sw_scale_fct.h"


/**
 * @typedef sw_wavelet_base_t
 * @brief Opaque type for a base class of a wavelet function.
 */
typedef struct sw_wavelet_base_s sw_wavelet_base_t;

/**
 * @typedef sw_wavelet_t
 * @brief Opaque type for a wavelet function.
 */
typedef struct sw_wavelet_s sw_wavelet_t;

/**
 * @typedef sw_wavelet_dual_t
 * @brief Opaque type for a dual wavelet function.
 */
typedef struct sw_wavelet_dual_s sw_wavelet_dual_t;


/* wavelet function - base class methods */

#include "sw_wavelet.x"


/*
 * wavelet function methods
 */

/**
 * @brief Create a new wavelet function from the given scale functions.
 *
 * @param sf The scale function.
 * @param sfd The dual scale function.
 * @return The new wavelet function object.
 *
 * This function creates a wavelet function built as a linear
 * combination of translated-dilated of @p sf, using the filter
 * coefficient of @p sdf. If @p sf or @p sfd are @c NULL, @c NULL is
 * returned. The filter and bounds of the wavelet function are
 * retrieved using the wavelet base class. When the returned wavelet is
 * not needed anymore, it must be freed with sw_wavelet_del().
 */
SAPI sw_wavelet_t *sw_wavelet_new(const sw_scale_fct_t      *sf,
                                  const sw_scale_fct_dual_t *sfd);

/**
 * @brief Free the memory of the given wavelet function.
 *
 * @param w The wavelet function to free.
 *
 * This function frees the memory of the wavelet function @p w. @p w
 * must have been created with sw_wavelet_new(). If @p w is @c NULL,
 * this function does nothing.
 */
SAPI void sw_wavelet_del(sw_wavelet_t *w);

/**
 * @brief Return the inferior bound of the support of the given wavelet function.
 *
 * @param w The wavelet function.
 * @return The inferior bound of the support of the wavelet function.
 *
 * This function returns the inferior bound of the support of the
 * wavelet function @p w. No test are done on @p w. If @p w is not
 * valid, a seg fault or an undefined behavior might occur.
 */
SAPI int32_t sw_wavelet_x1_get(const sw_wavelet_t *w);

/**
 * @brief Return the superior bound of the support of the given wavelet function.
 *
 * @param w The wavelet function.
 * @return The superior bound of the support of the wavelet function.
 *
 * This function returns the superior bound of the support of the
 * wavelet function @p w. No test are done on @p w. If @p w is not
 * valid, a seg fault or an undefined behavior might occur.
 */
SAPI int32_t sw_wavelet_x2_get(const sw_wavelet_t *w);

/**
 * @brief Return the rational value of the given wavelet function at the given
 * rational absciss.
 *
 * @param wd The wavelet function.
 * @param val The rational absciss.
 * @return The rational value of the wavelet function at the given absciss.
 *
 * This function returns the rational value of the wavelet function @p wd at
 * the rational absciss @p val.
 */
SAPI rational_t sw_wavelet_value_rat_get(const sw_wavelet_t  *wd,
                                         const rational_t    *val);

/**
 * @brief Return the real value of the given scale function at the given real
 * absciss.
 *
 * @param wd The scale function.
 * @param val The real absciss.
 * @return The real value of the scale function at the given absciss.
 *
 * This function return the rational value of the scale function @p wd at
 * the rational absciss @p val.
 */
SAPI double sw_wavelet_value_get(const sw_wavelet_t *wd,
                                 double              val);

/*
 * dual wavelet function methods
 */

/**
 * @brief Create a new dual wavelet function from the given scale functions.
 *
 * @param sf The scale function.
 * @param sfd The dual scale function.
 * @return The new wavelet dual function object.
 *
 * This function creates a dual wavelet function built as a linear
 * combination of translated-dilated of @p sfd, using the filter
 * coefficient of @p sd. If @p sf or @p sfd are @c NULL, @c NULL is
 * returned. The filter and bounds of the wavelet function are
 * retrieved using the wavelet base class. When the returned wavelet is
 * not needed anymore, it must be freed with sw_wavelet_del().
 */
SAPI sw_wavelet_dual_t *sw_wavelet_dual_new(const sw_scale_fct_t      *sf,
                                            const sw_scale_fct_dual_t *sfd);

/**
 * @brief Free the memory of the given dual wavelet function.
 *
 * @param wd The dual wavelet function to free.
 *
 * This function frees the memory of the dual wavelet function @p wd. @p wd
 * must have been created with sw_wavelet_dual_new(). If @p wd is @c NULL,
 * this function does nothing.
 */
SAPI void sw_wavelet_dual_del(sw_wavelet_dual_t *wd);

/**
 * @brief Return the inferior bound of the support of the given dual wavelet function.
 *
 * @param wd The dual wavelet function.
 * @return The inferior bound of the support of the wavelet function.
 *
 * This function returns the inferior bound of the support of the
 * dual wavelet function @p wd. No test are done on @p wd. If @p wd is not
 * valid, a seg fault or an undefined behavior might occur.
 */
SAPI int32_t sw_wavelet_dual_x1_get(const sw_wavelet_dual_t *wd);

/**
 * @brief Return the superior bound of the support of the given dual wavelet function.
 *
 * @param wd The dual wavelet function.
 * @return The superior bound of the support of the wavelet function.
 *
 * This function returns the superior bound of the support of the
 * dual wavelet function @p wd. No test are done on @p wd. If @p wd is not
 * valid, a seg fault or an undefined behavior might occur.
 */
SAPI int32_t sw_wavelet_dual_x2_get(const sw_wavelet_dual_t *wd);


#endif /* SW_WAVELET_H */
