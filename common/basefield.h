#include <arm_neon.h>

#ifndef BASEFIELD_H
#define BASEFIELD_H

//Formal interface

typedef struct {
	poly64x2_t p0, p1;
} pmullres;

void bf_print_expr(poly64x2_t p);

void bf_print_hex(poly64x2_t p);

poly64x2_t bf_rand_elem();

poly64x2_t bf_add(poly64x2_t a, poly64x2_t b);

pmullres bf_pmull(poly64x2_t a, poly64x2_t b);

pmullres bf_psquare(poly64x2_t a);

poly64x2_t bf_red(pmullres c);

poly64x2_t bf_inv(poly64x2_t a);

//Implementation alternatives:

pmullres bf_pmull32(poly64x2_t a, poly64x2_t b);

poly64x2_t bf_red_generic(pmullres c);

poly64x2_t fermat_inv(poly64x2_t a);

#endif
