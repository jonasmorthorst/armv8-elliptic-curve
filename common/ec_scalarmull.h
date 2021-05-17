#include <arm_neon.h>
#include "ec.h"
#include <stdio.h>

#ifndef EC_SCALARMULL_H
#define EC_SCALARMULL_H

#include "utils.h"

typedef uint64_t elt[2];

extern void lin_pass();

extern void c_sel();

ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_single_lproj(ec_point_lproj P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_double(ec_point_lproj P, uint64x2x2_t k, ec_point_lproj Q, uint64x2x2_t l);

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k);

static inline void linear_pass(ec_point_laffine *P1, ec_point_laffine *P2, ec_point_laffine* table, uint64_t index1, uint64_t index2, uint64_t l) {
	uint64_t val, new_ptr;
	ec_point_laffine item, P1_tmp, P2_tmp;

	for (uint64_t i = 0; i < l; i++) {
		// c_sel(index1, i, &P1_tmp, &table[i]);
		// c_sel(index2, i, &P2_tmp, &table[i]);
		CSEL(val, index1, i, P1_tmp, table[i], new_ptr, typeof(ec_point_laffine));
		CSEL(val, index2, i, P2_tmp, table[i], new_ptr, typeof(ec_point_laffine));
	}

	*P1 = P1_tmp;
	*P2 = ec_endo_laffine(P2_tmp);
}

static inline void linear_pass_inline_asm(ec_point_laffine *P, ec_point_laffine* table, uint64_t index, uint64_t l) {
	uint64_t val, new_ptr;
	asm volatile ("SUBS %0, %2, #0;" \
								"CSEL %1, %3, %1, EQ;" \
								"ADD %3, %3, #64;" \
								"SUBS %0, %2, #1;" \
								"CSEL %1, %3, %1, EQ;" \
								"ADD %3, %3, #64;" \
								"SUBS %0, %2, #2;" \
								"CSEL %1, %3, %1, EQ;" \
								"ADD %3, %3, #64;" \
								"SUBS %0, %2, #3;" \
								"CSEL %1, %3, %1, EQ;" \
								"ADD %3, %3, #64;" \
								"SUBS %0, %2, #4;" \
								"CSEL %1, %3, %1, EQ;" \
								"ADD %3, %3, #64;" \
								"SUBS %0, %2, #5;" \
								"CSEL %1, %3, %1, EQ;" \
								"ADD %3, %3, #64;" \
								"SUBS %0, %2, #6;" \
								"CSEL %1, %3, %1, EQ;" \
								"ADD %3, %3, #64;" \
								"SUBS %0, %2, #7;" \
								"CSEL %1, %3, %1, EQ;" \
								: "=r" ( val ), "=r" ( new_ptr ), "+r" (index), "+r" (table)
								);

	*P = *((ec_point_laffine*)new_ptr);
}

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
