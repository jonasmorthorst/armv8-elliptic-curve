#include <arm_neon.h>

#ifndef BASEFIELD_H
#define BASEFIELD_H

//Formal interface

void bf_print_expr(poly64x2_t p);

void bf_print_expr_nl(poly64x2_t p);

void bf_print_hex(poly64x2_t p);

void bf_print_hex_nl(poly64x2_t p);

poly64x2_t bf_create_elem(uint64_t l, uint64_t h);

poly64x2_t bf_rand_elem();

poly64x2_t bf_add(poly64x2_t a, poly64x2_t b);

poly64x2x2_t bf_pmull(poly64x2_t a, poly64x2_t b);

poly64x2x2_t bf_psquare(poly64x2_t a);

poly64x2_t bf_red(poly64x2x2_t c);

poly64x2_t bf_red_psquare(poly64x2x2_t c);

poly64x2_t bf_inv(poly64x2_t a);

//Implementation alternatives:

poly64x2x2_t bf_pmull32(poly64x2_t a, poly64x2_t b);

poly64x2x2_t bf_pmull64(poly64x2_t a, poly64x2_t b);

poly64x2_t bf_red_generic(poly64x2x2_t c);

poly64x2_t bf_red_formula(poly64x2x2_t c);

poly64x2_t bf_red_neon(poly64x2x2_t c);

poly64x2_t bf_red_psquare_formula(poly64x2x2_t c);

poly64x2_t bf_red_psquare_neon(poly64x2x2_t c);

poly64x2_t bf_fermat_inv(poly64x2_t a);

poly64x2_t bf_addchain_inv(poly64x2_t a);

poly64x2_t bf_addchain_lookup_inv(poly64x2_t a);

#endif
