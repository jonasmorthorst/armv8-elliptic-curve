#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/extensionfield.h"
#include "common/ec.h"
#include "common/ec_scalarmull.h"
#include "common/setup.h"
#include "common/utils.h"

void experiment() {
	uint64x2_t k1 = { 3, 0 };
	uint64x2_t k2 = { 0, 0 };
	k2[0] = k2[0]+1;

	ec_naf naf_k1 = ec_to_naf(k1);
	ec_naf naf_k2 = ec_to_naf(k2);

	// ec_print_naf(naf_k1);
	// ec_print_naf(naf_k2);


	// printf("%hhd\n", naf_k1.val[64]);

	ec_point_laffine P = ec_lproj_to_laffine((ec_point_lproj) GEN);

	//K1 naf 1  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -13
	//K2 naf 1  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15  -15

	uint64x2x2_t k1_digit = (uint64x2x2_t) {{{1, 0}, {0, 0}}};
	uint64x2x2_t k2_digit = (uint64x2x2_t) {{{1, 0}, {0, 0}}};

	ec_point_laffine P1 = ec_lproj_to_laffine(ec_scalarmull_single(P, k1_digit));
	ec_point_laffine P2 = ec_lproj_to_laffine(ec_scalarmull_single(P, k2_digit));
	P2 = ec_endo_laffine(P2);

	ec_point_lproj Q = ec_add_mixed(P1, ec_laffine_to_lproj(P2));


	for (int i = 63; i > 0; i--) {
		//Iteration i k1_digit=-15  k2_digit=-15
		Q = ec_double(ec_double(ec_double(Q)));
		k1_digit = (uint64x2x2_t) {{{15, 0}, {0, 0}}};
		k2_digit = (uint64x2x2_t) {{{15, 0}, {0, 0}}};

		P1 = ec_lproj_to_laffine(ec_scalarmull_single(P, k1_digit));
		P2 = ec_lproj_to_laffine(ec_scalarmull_single(P, k2_digit));
		P1 = ec_neg_laffine(P1);
		P2 = ec_neg_laffine(ec_endo_laffine(P2));

		Q = ec_double_then_addtwo(P1, P2, Q);

		printf("Q after iteration i=%d\n", i);
		ec_print_hex(Q);
	}

	//Iteration i=0 k1_digit=-13  k2_digit=-15
	Q = ec_double(ec_double(ec_double(Q)));
	k1_digit = (uint64x2x2_t) {{{13, 0}, {0, 0}}};
	k2_digit = (uint64x2x2_t) {{{15, 0}, {0, 0}}};

	P1 = ec_lproj_to_laffine(ec_scalarmull_single(P, k1_digit));
	P2 = ec_lproj_to_laffine(ec_scalarmull_single(P, k2_digit));
	P1 = ec_neg_laffine(P1);
	P2 = ec_neg_laffine(ec_endo_laffine(P2));

	Q = ec_double_then_addtwo(P1, P2, Q);

	printf("Q after iteration i=%d\n", 0);
	ec_print_hex(Q);

	// After loop
	uint64x2x2_t one = (uint64x2x2_t) {{{1, 0}, {0, 0}}};
	ec_point_lproj cP = ec_scalarmull_single(P, one);
	Q = ec_add(Q, ec_neg(cP));

	printf("Returning Q \n");
	ec_print_hex(Q);
}


int main() {
	init_components();

	experiment();

	printf("\n -------------------- \n");

	// uint64x2_t k = (uint64x2_t) {123, 124};
	//
	// ec_naf result = ec_to_naf(k);
	// printf("%s\n", "Our naf");
	// ec_print_naf(result);
	//
	// elt k1_naf = {k[0], k[1]};
	// signed char naf[65];
	// int len = 0;
	// int order_len = 253;
	// int w = 5;
	//
	// scmul_wreg(naf, &len, k1_naf, order_len, w);
	// printf("\n\n\n%s\n", "Diego's naf");
	// printf("Len %d\n", len);
	//
	// ec_print_naf_arr(naf);

	//uint64x2x2_t k = (uint64x2x2_t) {{{1, 1}, {1, 1}}};

	uint64x2x2_t k = (uint64x2x2_t) {{{3, 0}, {0, 0}}};
	ec_split_scalar decomp = ec_scalar_decomp(k);



	ec_naf result = ec_to_naf(bf_create_elem(decomp.k1[0], decomp.k1[1]));
	// printf("%s\n", "Our naf");
	// ec_print_naf(result);

	result = ec_to_naf(bf_create_elem(decomp.k2[0], decomp.k2[1]));
	// printf("\n");
	// ec_print_naf(result);

	elt k1_naf = {decomp.k1[0], decomp.k1[1]};
	signed char naf[100];
	int len = 0;
	int order_len = 253;
	int w = 5;

	scmul_wreg(naf, &len, k1_naf, order_len, w);
	// printf("\n\n\n%s\n", "Diego's naf");
	// printf("Len %d\n", len);

	// ec_print_naf_arr(naf);

	elt k2_naf = {decomp.k2[0], decomp.k2[1]};
	scmul_wreg(naf, &len, k2_naf, order_len, w);
	// printf("\n");
	//
	// ec_print_naf_arr(naf);



	//Arrange
	//uint64x2x2_t k = (uint64x2x2_t) {{{2243156791409331652485, 0}, {0, 0}}};
	// ef_elem PX = ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3));
	// ef_elem PL = ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC));
	// ef_elem PZ = ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000));
	// ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 78632917462800214 * GEN

	ec_point_lproj P = (ec_point_lproj) GEN;

	ec_point_lproj expected = ec_scalarmull_single(ec_lproj_to_laffine(P), k);

	//Act
	ec_point_laffine actual = ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), k);

	printf("Expected: \n");
	ec_print_expr(expected);

	printf("\n\nActual: \n");
	ec_print_expr_laffine(actual);

	printf("\nEqual: %lu", ec_equal_point_laffine(ec_lproj_to_laffine(expected), actual));
	printf("\nOn curve: %lu \n\n", ec_is_on_curve(ec_laffine_to_lproj(actual)));
	printf("\nOn curve: %lu \n\n", ec_is_on_curve(expected));


	// ef_elem PX = ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266));
	// ef_elem PL = ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490));
	// ef_elem PZ = ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001));
	// ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //99921481365893197563 * GEN
	//
	// ec_point_lproj table[4];
	// precompute(ec_lproj_to_laffine(P), table);


	// ec_naf result = ec_to_naf(bf_create_elem(153881, 0));
	// ec_naf result = ec_to_naf(bf_create_elem(35236236236, 2478623785632));
	// ec_print_naf(result);

	//printf("%p\n", &result);

	//153881

	// k = 3
	// m = 8
	//
	// i = 5
	// while (37 > 8)
	//
	// k0 = 1
	// k1 = -5
	// k2 = -3
	// k3 = 5
	// k4 = -3
	//
	// n = 5

	//sub_res_0 = 153880


	// i  |  ki  |  k
	// 0  |  -6  |  4540453076705
	// 1  |  -7  |

	// 1001 1010


	// 24411958 | 14699749183737301352

	// 1011101000111111100110110

// 01000011 11101001 10000000 00000000 00000000 00000000 00000000 00000001
//    1011101000111111100110 | 0001100110000000000000000000000000000000000000000000000000000000
//    1011101000111111100110	 1101100110000000000000000000000000000000000000000000000000000000
	// after division

	// 3051494 | 15672526703249326381

	dispose_components();
	return 0;
}
