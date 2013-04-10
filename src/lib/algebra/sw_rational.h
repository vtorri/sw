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
