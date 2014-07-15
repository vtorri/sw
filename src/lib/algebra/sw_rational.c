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

#include <stdlib.h>
#include <stdio.h>

#include "sw_rational.h"


/******************************************************************************/
/*                                                                            */
/*                                   LOCAL                                    */
/*                                                                            */
/******************************************************************************/

/**
 * @cond SW_LOCAL
 */

static int64_t
pgcd(int64_t p, int64_t q)
{
    int64_t a = p;
    int64_t b = q;
    int64_t r;

    while (b != 0)
    {
        r = a % b;
        a = b;
        b = r;
    }
    return a;
}

static int64_t
ppcm(int64_t p, int64_t q)
{
    int64_t d = pgcd(p, q);
    d = p / d;

    return d * q;
}

static void
rat_simplify(rational_t *r)
{
    int64_t d;

    if (r->num == 0)
    {
        r->den = 1;
        return;
    }

    d = pgcd(r->num, r->den);
    if (d != 1)
    {
        r->num /= d;
        r->den /= d;
    }
    if (r->den < 0)
    {
        r->num = -r->num;
        r->den = -r->den;
    }
}

/**
 * @endcond SW_LOCAL
 */


/******************************************************************************/
/*                                                                            */
/*                                   GLOBAL                                   */
/*                                                                            */
/******************************************************************************/


/******************************************************************************/
/*                                                                            */
/*                                    API                                     */
/*                                                                            */
/******************************************************************************/


rational_t
rat_new(int64_t num, int64_t den, uint8_t simpl)
{
    rational_t r = { 0, 1 };

    if (den == 0)
    {
        fprintf(stderr, "[rat] ERROR: rat_new : denominator NULL\n");
        abort();
        return r;
    }

    r.num = num;
    r.den = den;
    if (simpl)
        rat_simplify(&r);

    return r;
}

double
rat_double_get(const rational_t *r)
{
    return (double)r->num / (double)r->den;
}

rational_t
rat_add(const rational_t *r1, const rational_t *r2)
{
    rational_t  res;
    int64_t   den;
    int64_t   d;
    int64_t   num;

    if (r1->num == 0)
        return *r2;

    if (r2->num == 0)
        return *r1;

    den = ppcm(r1->den, r2->den);
    d = pgcd(r1->den, r2->den);
    num = r1->num * (r2->den / d) + r2->num * (r1->den / d);
    res = rat_new(num, den, 1);

    return res;
}

rational_t
rat_sub(const rational_t *r1, const rational_t *r2)
{
    rational_t res;
    rational_t r;

    r = rat_new(-r2->num, r2->den, 0);
    res = rat_add(r1, &r);

    return res;
}

rational_t
rat_opp(const rational_t *r)
{
    return rat_new(-r->num, r->den, 0);
}

rational_t
rat_mul(const rational_t *r1, const rational_t *r2)
{
    rational_t res = { 0, 1 };
    rational_t tmp1;
    rational_t tmp2;
    int64_t    num;
    int64_t    den;

    if ((r1->num == 0) || (r2->num == 0))
        return res;

    tmp1 = rat_new(r1->num, r2->den, 1);
    tmp2 = rat_new(r2->num, r1->den, 1);

    num = tmp1.num * tmp2.num;
    den = tmp1.den * tmp2.den;

/*   if (den == 0) { */
/*     printf("\n\n * ERROR : %lld, %lld, %lld\n\n", tmp1.den, tmp2.den, tmp1.den * tmp2.den); */
/*   } */

    res = rat_new(num, den, 0);

    return res;
}

rational_t
rat_div(const rational_t *r1, const rational_t *r2)
{
    rational_t res = { 0, 1 };
    rational_t r;

    if (r2->den == 0)
        return res;

    if (r2->num < 0)
        r = rat_new(-r2->den, -r2->num, 0);
    else
        r = rat_new(r2->den, r2->num, 0);

    res = rat_mul(r1, &r);
    return res;
}

rational_t
rat_abs(const rational_t *r)
{
    return rat_new((r->num > 0) ? r->num : -r->num, r->den, 0);
}

uint8_t
rat_greater(const rational_t *r1, const rational_t *r2)
{

    return ((r1->num * r2->den) > (r1->den * r2->num));
}

uint8_t
rat_greater_or_equal(const rational_t *r1, const rational_t *r2)
{

    return ((r1->num * r2->den) >= (r1->den * r2->num));
}

uint8_t
rat_lesser(const rational_t *r1, const rational_t *r2)
{

    return ((r1->num * r2->den) < (r1->den * r2->num));
}

uint8_t
rat_lesser_or_equal(const rational_t *r1, const rational_t *r2)
{

    return ((r1->num * r2->den) <= (r1->den * r2->num));
}

uint8_t
rat_is_equal(const rational_t *r,  int64_t val)
{
    return ((r->den == 1) && (r->num == val));
}

uint8_t
rat_is_equal_rat(const rational_t *r1,  const rational_t *r2)
{
    return ((r1->num * r2->den) == (r1->den * r2->num));
}

void
rat_disp(const rational_t *r, uint8_t end)
{
    printf("%2d", (int)r->num);
    if (r->den != 1)
        printf("/%d", (int)r->den);
    if (end)
        printf("\n");
}
