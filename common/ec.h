#include "extensionfield.h"

#ifndef EC_H
#define EC_H

#define A {{{0, 0}, {1, 0}}}
#define B {{{0x2E6D944FA54DE7E5, 0x59C8202CB9E6E0AE}, {0,0}}}

#define GEN {{{{0x344038B63FBA32DE, 0x203B6A93395E0432}, {0x5396E0681AA10E0D, 0x78E51FD0C310696D}}}, {{{0xDEB59C6137074B50, 0x5BD7653482085F55}, {0x24568FA5A1033946, 0x7F90D98B1589A17F}}} , {{{1,0}, {0,0}}}}
#define INFTY {{{{1, 0}, {0, 0}}}, {{{1, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}}}

#define SUBGROUP_ORDER {{{0x877DABA2A44750A5, 0xDAC40D1195270779},{0xFFFFFFFFFFFFFFFF, 0x1FFFFFFFFFFFFFFF}}}
#define TRACE 0xD792EA76691524E3

typedef struct {
	ef_elem x, l, z;
} ec_point_lproj;

typedef struct {
	ef_elem x, l;
} ec_point_laffine;

typedef struct {
	uint64x2_t k1, k2;
	uint64_t k1_sign, k2_sign;
} ec_split_scalar;

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem l, ef_elem z);

ec_point_laffine ec_create_point_laffine(ef_elem x, ef_elem l);

void ec_print_expr(ec_point_lproj P);

void ec_print_expr_laffine(ec_point_laffine P);

void ec_print_hex(ec_point_lproj P);

void ec_print_hex_laffine(ec_point_laffine P);

ec_point_lproj ec_laffine_to_lproj(ec_point_laffine P);

ec_point_laffine ec_lproj_to_laffine(ec_point_lproj P);

uint64_t ec_is_on_curve(ec_point_lproj P);

uint64_t ec_equal_point_lproj(ec_point_lproj P, ec_point_lproj Q);

uint64_t ec_equal_point_mixed(ec_point_laffine P, ec_point_lproj Q);

poly64x2x2_t ec_rand_scalar();

ec_point_lproj ec_rand_point_lproj();

ec_point_laffine ec_rand_point_laffine();

ec_point_lproj ec_neg(ec_point_lproj P);

ec_point_laffine ec_neg_laffine(ec_point_laffine P);

ec_point_lproj ec_add(ec_point_lproj P1, ec_point_lproj P2);

ec_point_lproj ec_add_mixed(ec_point_laffine P, ec_point_lproj Q);

ec_point_lproj ec_double(ec_point_lproj P);

ec_point_lproj ec_double_then_add(ec_point_laffine P, ec_point_lproj Q);

ec_point_lproj ec_double_then_addtwo(ec_point_laffine P1, ec_point_laffine P2, ec_point_lproj Q);

ec_point_laffine ec_endo_laffine(ec_point_laffine P);

ec_split_scalar ec_scalar_decomp(uint64x2x2_t k);

#endif
