#include <arm_neon.h>

#ifndef BASEFIELD_H
#define BASEFIELD_H

typedef struct {
	poly64x2_t p0, p1;
} pmullres;

void bf_print_elem_expr(poly64x2_t p);

poly64x2_t bf_rand_elem();

pmullres bf_pmull32(poly64x2_t a, poly64x2_t b);

#endif
