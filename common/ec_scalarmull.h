#include <arm_neon.h>
#include "ec.h"
#include <stdio.h>

#ifndef EC_SCALARMULL_H
#define EC_SCALARMULL_H

#include "utils.h"

extern void lin_pass_w3();

extern void lin_pass_w4();

extern void lin_pass_w5();

extern void lin_pass_w6();

extern void csel_asm();

ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_single_lproj(ec_point_lproj P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_double(ec_point_lproj P, uint64x2x2_t k, ec_point_lproj Q, uint64x2x2_t l);

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w3_randaccess(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w4_randaccess(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k);

ec_point_laffine ec_scalarmull_single_endo_w6_randaccess(ec_point_laffine P, uint64x2x2_t k);

void precompute_w3(ec_point_laffine P, ec_point_laffine* table);

void precompute_w4(ec_point_laffine P, ec_point_laffine* table);

void precompute_w5(ec_point_laffine P, ec_point_laffine* table);

void precompute_w6(ec_point_laffine P, ec_point_laffine* table);

void reg_rec(uint64x2_t k, uint64_t w, signed char* rec, uint64_t l);

void ec_print_rec(signed char *rec, uint64_t l);

#endif
