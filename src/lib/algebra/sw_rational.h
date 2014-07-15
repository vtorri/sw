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

#ifndef SW_RATIONAL_H
#define SW_RATIONAL_H


#include <stdint.h>

/**
 * @typedef rational_t
 * @brief The rational type
 */
typedef struct _rational_t rational_t;

/**
 * @struct _rational_t
 * @brief The rational type
 */
struct _rational_t
{
  int64_t num; /**< The numerator */
  int64_t den; /**< The denominator */
};

/**
 * @brief Return a new rational.
 *
 * @param num The numerator.
 * @param den The denominator.
 * @param simpl Whether the fraction is simplified or not.
 *
 * This function returns a new rational with @p num and @p den for
 * respectively the numerator and the denominator. To simplify the
 * fraction, pass 1 to @p simpl, otherwise pass 0.
 */
rational_t rat_new             (int64_t num, int64_t den, uint8_t simpl);
double     rat_double_get      (const rational_t *r);

rational_t rat_add             (const rational_t *r1, const rational_t *r2);
rational_t rat_sub             (const rational_t *r1, const rational_t *r2);
rational_t rat_opp             (const rational_t *r);
rational_t rat_mul             (const rational_t *r1, const rational_t *r2);
rational_t rat_div             (const rational_t *r1, const rational_t *r2);

rational_t rat_abs             (const rational_t *r);

uint8_t    rat_greater         (const rational_t *r1, const rational_t *r2);
uint8_t    rat_greater_or_equal(const rational_t *r1, const rational_t *r2);

uint8_t    rat_lesser          (const rational_t *r1, const rational_t *r2);
uint8_t    rat_lesser_or_equal (const rational_t *r1, const rational_t *r2);
uint8_t    rat_is_equal        (const rational_t *r,  int64_t val);

uint8_t    rat_is_equal_rat    (const rational_t *r1,  const rational_t *r2);

void       rat_disp            (const rational_t *r, uint8_t end);


#endif /* SW_RATIONAL_H */
