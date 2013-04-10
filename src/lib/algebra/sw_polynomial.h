#ifndef SW_POLYNOMIAL_H
#define SW_POLYNOMIAL_H


#include <stdint.h>

#include "sw_rational.h"


/**
 * @typedef rpol_t
 * @brief The polynomial type with rational coefficients
 */
typedef struct rpolynomial rpol_t;

/**
 * @typedef dpol_t
 * @brief The polynomial type with double coefficients
 */
typedef struct dpolynomial dpol_t;


/**
 * @struct rpolynomial
 * @brief The polynomial type with rational coefficients
 */
struct rpolynomial
{
  const rational_t *coefs;  /**< the coefficients */
  int32_t           degree; /**< the degree */
};

/**
 * @struct dpolynomial
 * @brief The polynomial type with double coefficients
 */
struct dpolynomial
{
  const double *coefs;  /**< the coefficients */
  int32_t       degree; /**< the degree */
};


rpol_t    *rpol_new(const rational_t *coefs,
                    int32_t           degree);

dpol_t    *dpol_new(const double *coefs,
                    int32_t       degree);

void       rpol_delete(rpol_t *pol);

int32_t    rpol_degree_get(const rpol_t *pol);

rational_t rpol_coef_get(rpol_t *pol,
                         int32_t n);

rational_t rpol_monom_get(const rpol_t     *pol,
                          int32_t           n,
                          int32_t           smooth,
                          const rational_t *val);

rational_t rpol_value_get(const rpol_t     *pol,
                          const rational_t *x);

double     dpol_value_get(const dpol_t *pol,
                          double        x);

rational_t rpol_integrate(const rpol_t *pol, int32_t a, int32_t b);

rpol_t    *rpol_lagrange_base_get(int64_t *x, int32_t num, int32_t i);

rpol_t   **sw_rpol_chebychev_pols_get(int32_t r);

void       rpol_disp(const rpol_t *pol, uint8_t end);

void       dpol_disp(const dpol_t *pol, uint8_t end);


#endif /* SW_POLYNOMIAL_H */
