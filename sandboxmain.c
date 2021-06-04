#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/extensionfield.h"
#include "common/ec.h"
#include "common/ec_scalarmull.h"
#include "common/setup.h"
#include "common/utils.h"

int main() {
	init_components();

	ec_point_laffine P = ec_rand_point_laffine();
	ec_point_laffine table[8];
	precompute_w5(P, table);

	ec_point_laffine P1;
	ec_point_laffine P2;
	lin_pass_w5(&P1, &P2, &table, 4, 7);

	P2 = ec_endo_laffine(P2);
	ec_point_laffine expected1 = table[4];
	ec_point_laffine expected2 = ec_endo_laffine(table[7]);

	uint64_t equal1 = ec_equal_point_laffine(expected1, P1);
	printf("equal1: %lu\n", equal1);
	uint64_t equal2 = ec_equal_point_laffine(expected2, P2);
	printf("equal2: %lu\n", equal2);




	// uint64x2x2_t k = (uint64x2x2_t) {{{253256326376, 457436346236}, {124525, 11352535}}};
	//
	// ec_point_laffine P = ec_rand_point_laffine();
	//
	// ec_point_lproj expected = ec_scalarmull_single(P, k);
	// //Act
	// ec_point_laffine actual = ec_scalarmull_single_endo_w5_randaccess(P, k);
	//
	// //Assert
	// uint64_t equal = ec_equal_point_lproj(expected, ec_laffine_to_lproj(actual));
	// uint64_t on_curve = ec_is_on_curve_laffine(actual);
	//
	// printf("Equal: %lu\n", equal);
	// printf("On curve: %lu\n", on_curve);



	// ec_point_laffine P = ec_rand_point_laffine();
	// ec_point_laffine table[8];
	// precompute(P, table);
	//
	// ec_point_laffine P1;
	//
	// // uint64_t r1 = (uint64_t) &table;
	// // uint64_t r2 = (uint64_t) &table[1];
	// //
	// printf("%p\n", &P1);
	// printf("%p\n", table);

	// ec_print_hex_laffine(P1);


	// ec_print_hex_laffine(P1);
	// ec_print_hex_laffine(table[7]);


	// linear_pass_new1(&P1, table, 7, 8);
	// lin_pass(&P, &table, 0);
	// ec_point_laffine expected1 = table[0];
	//
	// uint64_t equal1 = ec_equal_point_laffine(expected1, P);
	// printf("equal1: %lu\n", equal1);
	//
	// uint64_t on_curve = ec_is_on_curve_laffine(table[7]);
	// printf("On curve: %lu\n", on_curve);

	//
	// ec_point_laffine P1;
	// ec_point_laffine P2;
	//
	// ec_print_hex_laffine(P1);
	//
	//
	// linear_pass_new(&P1, table, 5, 8);
	// linear_pass_new(&P2, table, 1, 8);
	// ec_point_laffine expected1 = table[5];
	// ec_point_laffine expected2 = table[1];
	//
	// uint64_t equal1 = ec_equal_point_laffine(expected1, P1);
	// ec_print_hex_laffine(expected1);
	// ec_print_hex_laffine(P1);
	// printf("equal1: %lu\n", equal1);
	// uint64_t equal2 = ec_equal_point_laffine(expected2, P2);
	// printf("Equal2: %lu\n", equal2);







	// poly64x2x2_t a = (poly64x2x2_t) {{{0, 0}, {9223372036854775808U, 0}}};
	// poly64x2_t ared = bf_red(a);
	// poly64x2_t aredpsquare = bf_red_psquare_formula(a);
	// poly64x2_t aredlazyformula = bf_red_lazy_formula(a);
	// poly64x2_t aredlazy = bf_red_lazy(a);
	// bf_print_expr_nl(ared);
	// bf_print_expr_nl(aredpsquare);
	// bf_print_expr_nl(aredlazyformula);
	// bf_print_expr_nl(aredlazy);

	// signed char k1_digit = -15;
	// uint64_t k1_sign = ((unsigned char)k1_digit >> 7);
	//
	// uint64_t zero = 0;
	// signed char k1_val = k1_digit;
	// k1_val ^= (zero - k1_sign);
	// k1_val += k1_sign;
	//
	// printf("k1 sign: %lu\n", k1_sign);
	// printf("k1 digit: %hhd\n", k1_digit);
	// printf("k1 abs value: %hhd\n", k1_val);


	// printf("\n -------------------- \n");

// 	k: 17543080488742735143, 14748808945444813383, 8374720340618217181, 1767656777874844946
// x: p0: 1fc34c09780b0839||0dbe04483b1838ea p1: 4be7046f6821989b||dbfe7553e9fd2369
//  l: p0: 284c3c7af24c119f||02352b17f9da5d0a p1: 1cb1c6318570bd4a||7e49abaf1ffd8b9f

	// c1=1 k1sign=1
	//uint64x2x2_t k = (uint64x2x2_t) {{{431572255129370311, 5384670855758030008}, {4353603234575721204, 773725909533223690}}};

	// c1,2=1 k1sign=1
	// uint64x2x2_t k = (uint64x2x2_t) {{{23523616, 436346236}, {239852385867236, 4262357238457}}};
	//k1 sign 0 k2 sign 1 c1: 0 c2: 0
	// uint64x2x2_t k = (uint64x2x2_t) {{{253256326376, 457436346236}, {22385867236, 427238457}}};
	//
	// uint64x2x2_t k = (uint64x2x2_t) {{{1, 0}, {0, 0}}};
	//
	// ef_elem PX = ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3));
	// ef_elem PL = ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC));
	// ef_elem PZ = ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000));
	// ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 78632917462800214 * GEN
	//
	// ec_point_lproj expected = ec_scalarmull_single(ec_lproj_to_laffine(P), k);
	// //Act
	// ec_point_laffine actual = ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), k);
	//
	// //Assert
	// uint64_t equal = ec_equal_point_lproj(expected, ec_laffine_to_lproj(actual));
	// uint64_t on_curve = ec_is_on_curve_laffine(actual);
	//
	// printf("Equal: %lu\n", equal);
	// printf("On curve: %lu\n", on_curve);


	// uint64x2x2_t k = ec_rand_scalar();
	// ec_point_laffine P = ec_rand_point_laffine();
	//
	// ec_point_lproj expected = ec_scalarmull_single(P, k);
	// //Act
	// ec_point_lproj actual = ec_scalarmull_single_endo_w3_randaccess(P, k);
	//
	// //Assert
	// uint64_t equal = ec_equal_point_lproj(expected, actual);
	// uint64_t on_curve = ec_is_on_curve(actual);



	// ec_split_scalar decomp = ec_scalar_decomp(k);
	//
	//
	// //
	// ec_naf result = ec_to_naf(decomp.k1, 6);
	// printf("%s\n", "Our naf");
	// ec_print_naf(result, 6);
	//
	// elt k1_naf = {decomp.k1[0], decomp.k1[1]};
	// signed char naf[200];
	// int len = 0;
	// int order_len = 253;
	// int w = 2;
	//
	// scmul_wreg(naf, &len, k1_naf, order_len, w);
	// printf("\n\n\n%s\n", "Diego's naf");
	// printf("Len %d\n", len);
	//
	// ec_print_naf_arr(naf, 52);

	// uint64x2x2_t k = (uint64x2x2_t) {{{1, 1}, {1, 1}}};
	//
	// // uint64x2x2_t k = (uint64x2x2_t) {{{3, 0}, {0, 0}}};
	// ec_split_scalar decomp = ec_scalar_decomp(k);
	//
	//
	//
	// ec_naf result = ec_to_naf(bf_create_elem(decomp.k1[0], decomp.k1[1]));
	// // printf("%s\n", "Our naf");
	// // ec_print_naf(result);
	//
	// result = ec_to_naf(bf_create_elem(decomp.k2[0], decomp.k2[1]));
	// // printf("\n");
	// // ec_print_naf(result);
	//
	// elt k1_naf = {decomp.k1[0], decomp.k1[1]};
	// signed char naf[100];
	// int len = 0;
	// int order_len = 253;
	// int w = 5;
	//
	// scmul_wreg(naf, &len, k1_naf, order_len, w);
	// // printf("\n\n\n%s\n", "Diego's naf");
	// // printf("Len %d\n", len);
	//
	// // ec_print_naf_arr(naf);
	//
	// elt k2_naf = {decomp.k2[0], decomp.k2[1]};
	// scmul_wreg(naf, &len, k2_naf, order_len, w);
	// // printf("\n");
	// //
	// // ec_print_naf_arr(naf);
	//
	// //Arrange
	// //uint64x2x2_t k = (uint64x2x2_t) {{{2243156791409331652485, 0}, {0, 0}}};
	// // ef_elem PX = ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3));
	// // ef_elem PL = ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC));
	// // ef_elem PZ = ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000));
	// // ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 78632917462800214 * GEN
	//
	// ec_point_lproj P = (ec_point_lproj) GEN;
	//
	// ec_point_lproj expected = ec_scalarmull_single(ec_lproj_to_laffine(P), k);
	//
	// //Act
	// ec_point_laffine actual = ec_scalarmull_single_endo_w5_randaccess(ec_lproj_to_laffine(P), k);
	//
	// printf("Expected: \n");
	// ec_print_expr(expected);
	//
	// printf("\n\nActual: \n");
	// ec_print_expr_laffine(actual);
	//
	// printf("\nEqual: %lu", ec_equal_point_laffine(ec_lproj_to_laffine(expected), actual));
	// printf("\nOn curve: %lu \n\n", ec_is_on_curve(ec_laffine_to_lproj(actual)));
	// printf("\nOn curve: %lu \n\n", ec_is_on_curve(expected));


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
