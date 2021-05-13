#include <arm_neon.h>
#include "ec.h"

#ifndef EC_SCALARMULL_H
#define EC_SCALARMULL_H

typedef uint64_t elt[2];

ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_single_lproj(ec_point_lproj P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_double(ec_point_lproj P, uint64x2x2_t k, ec_point_lproj Q, uint64x2x2_t l);

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k);

void linear_pass(ec_point_laffine *P1, ec_point_laffine *P2, ec_point_laffine* table, uint64_t index1, uint64_t index2, uint64_t l);

ec_point_laffine ec_scalarmull_single_endo_w3_randaccess(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w4_randaccess(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w6_randaccess(ec_point_laffine P, uint64x2x2_t k);

void scmul_wreg(signed char *naf, int *len, elt k, int n, int w);

void precompute(ec_point_laffine P, ec_point_laffine* table);

void precompute_w3(ec_point_laffine P, ec_point_laffine* table);

void precompute_w4(ec_point_laffine P, ec_point_laffine* table);

void precompute_w6(ec_point_laffine P, ec_point_laffine* table);

void ec_to_naf(uint64x2_t k, uint64_t w, signed char* naf);

void ec_print_naf(signed char *naf, uint64_t l);

void print_table(ec_point_laffine* table);

#endif
