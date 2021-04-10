#include <arm_neon.h>

#ifndef BASEFIELD_H
#define BASEFIELD_H

//Formal interface

typedef struct {
	poly64x2_t p0, p1;
} bf_polyx2;

void bf_print_expr(poly64x2_t p);

void bf_print_hex(poly64x2_t p);

poly64x2_t bf_rand_elem();

poly64x2_t bf_add(poly64x2_t a, poly64x2_t b);

bf_polyx2 bf_pmull(poly64x2_t a, poly64x2_t b);

bf_polyx2 bf_psquare(poly64x2_t a);

poly64x2_t bf_red(bf_polyx2 c);

poly64x2_t bf_red_psquare(bf_polyx2 c);

poly64x2_t bf_inv(poly64x2_t a);

//Implementation alternatives:

bf_polyx2 bf_pmull32(poly64x2_t a, poly64x2_t b);

bf_polyx2 bf_pmull64(poly64x2_t a, poly64x2_t b);

poly64x2_t bf_red_generic(bf_polyx2 c);

poly64x2_t bf_red_formula(bf_polyx2 c);

poly64x2_t bf_red_neon(bf_polyx2 c);

poly64x2_t bf_red_psquare_formula(bf_polyx2 c);

poly64x2_t bf_red_psquare_neon(bf_polyx2 c);

poly64x2_t bf_fermat_inv(poly64x2_t a);

poly64x2_t bf_addition_chain_inv(poly64x2_t a);

#endif
