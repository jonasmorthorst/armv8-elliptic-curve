#include "extensionfield.h"

#ifndef EC_H
#define EC_H

#define A {{{0, 0}, {1, 0}}}
#define B {{{0x2E6D944FA54DE7E5, 0x59C8202CB9E6E0AE}, {0,0}}}

#define GEN {{{{0x344038B63FBA32DE, 0x203B6A93395E0432}, {0x5396E0681AA10E0D, 0x78E51FD0C310696D}}}, {{{0xDEB59C6137074B50, 0x5BD7653482085F55}, {0x24568FA5A1033946, 0x7F90D98B1589A17F}}} , {{{1,0}, {0,0}}}}
#define INFTY {{{{1, 0}, {0, 0}}}, {{{1, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}}}

#define SUBGROUP_ORDER {{{0x877DABA2A44750A5, 0xDAC40D1195270779},{0xFFFFFFFFFFFFFFFF, 0x1FFFFFFFFFFFFFFF}}}

typedef struct {
	ef_elem x, l, z;
} ec_point_lproj;

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem l, ef_elem z);

void ec_print_expr(ec_point_lproj P);

void ec_print_hex(ec_point_lproj P);

uint64_t ec_is_on_curve(ec_point_lproj P);

uint64_t ec_equal_point_lproj(ec_point_lproj P, ec_point_lproj Q);

poly64x2x2_t ec_rand_scalar();

ec_point_lproj ec_rand_point_lproj();

ec_point_lproj ec_neg(ec_point_lproj P);

ec_point_lproj ec_add(ec_point_lproj P1, ec_point_lproj P2);

ec_point_lproj ec_double(ec_point_lproj P);

#endif
