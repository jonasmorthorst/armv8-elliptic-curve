#include "extensionfield.h"

#ifndef ECPOINT_H
#define ECPOINT_H

#define A {{0, 0}, {1, 0}}
#define B {{0x2E6D944FA54DE7E5, 0x59C8202CB9E6E0AE}, {0,0}}

#define GEN {{{0x344038B63FBA32DE, 0x203B6A93395E0432}, {0x5396E0681AA10E0D, 0x78E51FD0C310696D}}, {{0xDEB59C6137074B50, 0x5BD7653482085F55}, {0x24568FA5A1033946, 0x7F90D98B1589A17F}} , {{1,0}, {0,0}}}

#define ORDER {{0x877DABA2A44750A5, 0xDAC40D1195270779},{0xFFFFFFFFFFFFFFFF, 0x1FFFFFFFFFFFFFFF}}

typedef struct {
	ef_elem x, l, z;
} ec_point_lproj;

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem l, ef_elem z);

void ec_print_point_lproj_expr(ec_point_lproj p);

void ec_print_point_lproj_hex(ec_point_lproj p);

ec_point_lproj ec_rand_point_lproj();

ec_point_lproj ec_point_lproj_scalar_mull(uint64_t k, ec_point_lproj p);

ec_point_lproj ec_point_lproj_double_scalar_mull(uint64_t k, ec_point_lproj p1, uint64_t j, ec_point_lproj p2);

ec_point_lproj ec_point_lproj_neg(ec_point_lproj p);

ec_point_lproj ec_point_lproj_add(ec_point_lproj p1, ec_point_lproj p2);

ec_point_lproj ec_point_lproj_double(ec_point_lproj p);

#endif
