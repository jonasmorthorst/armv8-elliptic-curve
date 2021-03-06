#include "extensionfield_interleaved.h"

#ifndef EC_H
#define EC_H

#define A {{{0, 1}, {0, 0}}}
#define B {{{0x2E6D944FA54DE7E5, 0}, {0x59C8202CB9E6E0AE,0}}}

#define GEN {{{{0x344038B63FBA32DE, 0x5396E0681AA10E0D}, {0x203B6A93395E0432, 0x78E51FD0C310696D}}}, {{{0xDEB59C6137074B50, 0x24568FA5A1033946}, {0x5BD7653482085F55, 0x7F90D98B1589A17F}}} , {{{1,0}, {0,0}}}}
#define INFTY {{{{1, 0}, {0, 0}}}, {{{1, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}}}

#define SUBGROUP_ORDER {{{0x877DABA2A44750A5, 0xDAC40D1195270779},{0xFFFFFFFFFFFFFFFF, 0x1FFFFFFFFFFFFFFF}}}

typedef struct {
	ef_intrl_elem x, l, z;
} ec_point_lproj;

typedef struct {
	ef_intrl_elem x, l;
} ec_point_laffine;

ec_point_lproj ec_create_point_lproj(ef_intrl_elem x, ef_intrl_elem l, ef_intrl_elem z);

ec_point_laffine ec_create_point_laffine(ef_intrl_elem x, ef_intrl_elem l);

void ec_print_expr(ec_point_lproj P);

void ec_print_expr_laffine(ec_point_laffine P);

void ec_print_hex(ec_point_lproj P);

void ec_print_hex_laffine(ec_point_laffine P);

ec_point_lproj ec_laffine_to_lproj(ec_point_laffine P);

ec_point_laffine ec_lproj_to_laffine(ec_point_lproj P);

uint64_t ec_is_on_curve(ec_point_lproj P);

uint64_t ec_is_on_curve_laffine(ec_point_laffine P);

uint64_t ec_equal_point_lproj(ec_point_lproj P, ec_point_lproj Q);

uint64_t ec_equal_point_mixed(ec_point_laffine P, ec_point_lproj Q);

uint64_t ec_equal_point_laffine(ec_point_laffine P, ec_point_laffine Q);

uint64x2x2_t ec_rand_scalar();

ec_point_lproj ec_rand_point_lproj();

ec_point_laffine ec_rand_point_laffine();

static inline ec_point_lproj ec_neg(ec_point_lproj P) {
	P.l = ef_intrl_add(P.l, P.z);
	return P;
}

static inline ec_point_laffine ec_neg_laffine(ec_point_laffine P) {
	P.l.val[0][0] ^= 1;
	return P;
}

ec_point_lproj ec_add(ec_point_lproj P1, ec_point_lproj P2);

ec_point_lproj ec_add_unchecked(ec_point_lproj P1, ec_point_lproj P2);

ec_point_lproj ec_add_mixed(ec_point_laffine P, ec_point_lproj Q);

ec_point_lproj ec_add_mixed_unchecked(ec_point_laffine P, ec_point_lproj Q);

ec_point_lproj ec_add_laffine_unchecked(ec_point_laffine P, ec_point_laffine Q);

ec_point_lproj ec_double(ec_point_lproj P);

ec_point_lproj ec_double_mixed(ec_point_laffine P);

ec_point_lproj ec_double_alt(ec_point_lproj P);

ec_point_lproj ec_double_then_add(ec_point_laffine P, ec_point_lproj Q);

ec_point_lproj ec_double_then_addtwo(ec_point_laffine P1, ec_point_laffine P2, ec_point_lproj Q);

static inline ec_point_laffine ec_endo_laffine(ec_point_laffine P) {
	/*P.x.val[0] = bf_add(P.x.val[0], P.x.val[1]);
	P.l.val[0] = bf_add(P.l.val[0], P.l.val[1]);
	P.l.val[1] = bf_add(P.l.val[1], (poly64x2_t) {1,0});*/

	//x_1u + (x0+x1)
	poly64x2_t t = vextq_p64(P.x.val[0], P.x.val[0], 1);
	P.x.val[0][0] ^= t[0];
	t = vextq_p64(P.x.val[1], P.x.val[1], 1);
	P.x.val[1][0] ^= t[0];

	//(l_1+1)u + (l0+l1)
	t[0] = 1;
	t = vextq_p64(P.l.val[0], t, 1);
	P.l.val[0] = (poly64x2_t) veorq_u64((uint64x2_t) P.l.val[0], (uint64x2_t) t);
	t = vextq_p64(P.l.val[1], P.l.val[1], 1);
	P.l.val[1][0] ^= t[0];

	return P;
}

#endif
