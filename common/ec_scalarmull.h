#include <arm_neon.h>
#include "ec.h"
#include <stdio.h>

#ifndef EC_SCALARMULL_H
#define EC_SCALARMULL_H

#include "utils.h"

typedef uint64_t elt[2];

ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_single_lproj(ec_point_lproj P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_double(ec_point_lproj P, uint64x2x2_t k, ec_point_lproj Q, uint64x2x2_t l);

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k);

static inline void linear_pass(ec_point_laffine *P1, ec_point_laffine *P2, ec_point_laffine* table, uint64_t index1, uint64_t index2, uint64_t l) {
	uint64_t val, new_ptr;
	ec_point_laffine item, P1_tmp, P2_tmp;

	for (uint64_t i = 0; i < l; i++) {
		CSEL(val, index1, i, P1_tmp, table[i], new_ptr, typeof(ec_point_laffine));
		CSEL(val, index2, i, P2_tmp, table[i], new_ptr, typeof(ec_point_laffine));
	}

	*P1 = P1_tmp;
	*P2 = ec_endo_laffine(P2_tmp);
}

void linear_pass_new1(ec_point_laffine *P, ec_point_laffine* table, uint64_t index, uint64_t l);

// static inline void linear_pass_new2(ec_point_laffine *P, ec_point_laffine* table, uint64_t index, uint64_t l) {
// 	uint64_t val, new_ptr;
// 	ec_point_laffine P_tmp;
//
// 	asm volatile ("SUBS %0, %2, #0;" \
// 								"CSEL %1, %3, %1, EQ;" \
// 								"SUBS %0, %2, #1;" \
// 								"CSEL %1, %4, %1, EQ;" \
// 								"SUBS %0, %2, #2;" \
// 								"CSEL %1, %5, %1, EQ;" \
// 								"SUBS %0, %2, #3;" \
// 								"CSEL %1, %6, %1, EQ;" \
// 								"SUBS %0, %2, #4;" \
// 								"CSEL %1, %7, %1, EQ;" \
// 								"SUBS %0, %2, #5;" \
// 								"CSEL %1, %8, %1, EQ;" \
// 								"SUBS %0, %2, #6;" \
// 								"CSEL %1, %9, %1, EQ;" \
// 								"SUBS %0, %2, #7;" \
// 								"CSEL %1, %10, %1, EQ;" \
// 								: "+r" ( val ), "+r" ( new_ptr ) \
// 								: "r" (index), "r" (&table[0]), "r" (&table[1]), "r" (&table[2]), "r" (&table[3]), "r" (&table[4]), "r" (&table[5]), "r" (&table[6]), "r" (&table[7]) \
// 								);
// 	*P = *((ec_point_laffine*)new_ptr);
// 	// asm volatile ("SUBS %0, %2, #0;" \
// 	// 							"CSEL %1, %4, %1, EQ;" \
// 	// 							"SUBS %0, %2, #1;" \
// 	// 							"CSEL %1, %5, %1, EQ;" \
// 	// 							"SUBS %0, %2, #2;" \
// 	// 							"CSEL %1, %6, %1, EQ;" \
// 	// 							"SUBS %0, %2, #3;" \
// 	// 							"CSEL %1, %7, %1, EQ;" \
// 	// 							"SUBS %0, %2, #4;" \
// 	// 							"CSEL %1, %8, %1, EQ;" \
// 	// 							"SUBS %0, %2, #5;" \
// 	// 							"CSEL %1, %9, %1, EQ;" \
// 	// 							"SUBS %0, %2, #6;" \
// 	// 							"CSEL %1, %10, %1, EQ;" \
// 	// 							"SUBS %0, %2, #7;" \
// 	// 							"CSEL %1, %11, %1, EQ;" \
// 	// 							: "+r" ( val ), "+r" ( new_ptr2 ) \
// 	// 							: "r" (index2), "r" (&P2_tmp), "r" (&table[0]), "r" (&table[1]), "r" (&table[2]), "r" (&table[3]), "r" (&table[4]), "r" (&table[5]), "r" (&table[6]), "r" (&table[7]) \
// 	// 							);
// 	// P2_tmp = *((ec_point_laffine*)new_ptr2);
//
// 	// *P2 = ec_endo_laffine(P2_tmp);
// }

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
