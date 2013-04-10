#ifndef SW_SYSTEM_H
#define SW_SYSTEM_H


#include "sw_rational.h"


typedef struct system system_t;


system_t         *system_new(rational_t *matrix,
                             rational_t *vector,
                             int32_t     order);

void              system_delete(system_t *s);

void              system_solve(system_t *s);

const rational_t *system_solution_get(const system_t *s);


#endif /* SW_SYSTEM_H */
