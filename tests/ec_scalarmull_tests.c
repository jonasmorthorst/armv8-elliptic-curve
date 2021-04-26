#include <stdio.h>

#include "ec_scalarmull_tests.h"
#include "../common/ec_scalarmull.h"
#include "../common/ec.h"
#include "../common/extensionfield.h"
#include "../common/basefield.h"

void ec_scalarmull_single_test_example(test_ctr *ctr) {
	//Arrange
	ef_elem k = ef_create_elem(bf_create_elem(1984, 0), bf_create_elem(0, 0));

	ef_elem EX = ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25));
	ef_elem EL = ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB));
	ef_elem EZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0));
	ec_point_lproj expected = ec_create_point_lproj(EX, EL, EZ); // expected = 1984 * GEN

	//Act
	ec_point_lproj actual = ec_scalarmull_single((ec_point_lproj)GEN, k);

	//Assert
	uint64_t equal = ec_equal_point_lproj(expected, actual);
	uint64_t on_curve = ec_is_on_curve(actual);
	assert_true(equal && on_curve, ctr, "ec: ec_scalar_mull_test_example FAILED");
}

void ec_scalarmull_single_test_linearity(test_ctr *ctr) {
	//Arrange
	ef_elem k = ef_create_elem(bf_create_elem(7876556743120, 65714569742132121), bf_create_elem(11, 0));

	ef_elem PX = ef_create_elem(bf_create_elem(0X51441C4EE272FE55, 0X2DB9775DAEDDE550), bf_create_elem(0X12DD1A65F1D5B480, 0X6CB9034E20AD0EEB));
	ef_elem PL = ef_create_elem(bf_create_elem(0X7AA01DC08C73455A, 0X51F2DF8B2F5FA18C), bf_create_elem(0X6730BC49B9A98F41, 0X57DEBE6DBCE321DE));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X80, 0), bf_create_elem(0X0000000200000000, 0X2000000000));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //(order-1) * GEN

	ef_elem QX = ef_create_elem(bf_create_elem(0XB17380115E6C4F69, 0XE1DA273142A52B3), bf_create_elem(0X529D0294EB09720D, 0X26690AF5D1CC9E78));
	ef_elem QL = ef_create_elem(bf_create_elem(0XA75AABA91400402C, 0X2EB5E7A7EE3E6DEA), bf_create_elem(0X82C9609FFC5E6794, 0X19B8A04EA04D0DF));
	ef_elem QZ = ef_create_elem(bf_create_elem(0X8150BD7B681CEB67, 0X1773D8F08AD460A5), bf_create_elem(0X64E2D12FB4409DDD, 0X391936225DE0A796));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //14329 * GEN

	//Act
	ec_point_lproj add_first = ec_scalarmull_single(ec_add(P, Q), k);
	ec_point_lproj add_after = ec_add(ec_scalarmull_single(P, k), ec_scalarmull_single(Q, k));

	//Assert
	uint64_t correct = ec_equal_point_lproj(add_first, add_after) && ec_is_on_curve(add_first);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_linearity FAILED");
}

void ec_scalarmull_single_test_negation_order_indifference(test_ctr *ctr) {
	//Arrange
	ef_elem k = ef_create_elem(bf_create_elem(76591954, 7695159674259757), bf_create_elem(95124743611111, 56214));
	ef_elem PX = ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4));
	ef_elem PL = ef_create_elem(bf_create_elem(0XB55AF6853BAA916A, 0X368B1A5434DF7331), bf_create_elem(0X4B3628231E3A83C2, 0XCA5DED000C90027));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //2243156791409331652485 * P

	//Act
	ec_point_lproj neg_first = ec_scalarmull_single(ec_neg(P), k);
	ec_point_lproj neg_last = ec_neg(ec_scalarmull_single(P, k));

	//Assert
	uint64_t correct = ec_equal_point_lproj(neg_first, neg_last);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_negation_order_indifference FAILED");
}

void ec_scalarmull_single_test_order_of_gen_is_order_of_subgroup(test_ctr *ctr) {
	//Arrange
	poly64x2x2_t order_minus_1 = (poly64x2x2_t) SUBGROUP_ORDER;
	order_minus_1.val[0][0]--;

	//Act
	ec_point_lproj gen_inv = ec_scalarmull_single((ec_point_lproj) GEN, order_minus_1);
	ec_point_lproj gen_plus_gen_inv = ec_add((ec_point_lproj) GEN, gen_inv);

	//Assert
	uint64_t correct = ec_equal_point_lproj(gen_plus_gen_inv, (ec_point_lproj) INFTY);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_order_of_gen_is_order_of_subgroup FAILED");
}

void ec_scalarmull_single_test_final_k_at_once_same_as_factor_one_at_a_time(test_ctr *ctr) {
	//Arrange
	ef_elem k1 = ef_create_elem(bf_create_elem(0XB56483101AD77613, 0X123A67F2366799C), bf_create_elem(0, 0));
	ef_elem k2 = ef_create_elem(bf_create_elem(0XEEFFEE6752A38174, 0X7544), bf_create_elem(0, 0));
	ef_elem k = ef_create_elem(bf_create_elem(0XD6D3C77A003A139C, 0XACF7CBD9DC139043), bf_create_elem(0X99A09D4FE2DACA55, 0X85)); //k=k1*k2

	ef_elem PX = ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266));
	ef_elem PL = ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //99921481365893197563 * GEN

	//Act
	ec_point_lproj k1P = ec_scalarmull_single(P, k1);
	ec_point_lproj k1k2P = ec_scalarmull_single(k1P, k2);
	ec_point_lproj kP = ec_scalarmull_single(P, k);

	//Assert
	uint64_t correct = ec_equal_point_lproj(kP, k1k2P) && ec_is_on_curve(k1P) && ec_is_on_curve(k1k2P);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_final_k_at_once_same_as_factor_one_at_a_time FAILED");
}

void ec_scalarmull_single_test_k_one_is_identity(test_ctr *ctr) {
	//Arrange
	ef_elem k = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));

	ef_elem PX = ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3));
	ef_elem PL = ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 78632917462800214 * GEN

	//Act
	ec_point_lproj result = ec_scalarmull_single(P, k);

	//Assert
	uint64_t correct = ec_equal_point_lproj(result, P) && ec_is_on_curve(result);
	assert_true(correct, ctr, "ec: ec_scalarmull_single_test_k_one_is_identity FAILED");
}

void ec_scalarmull_double_test_example(test_ctr *ctr) {
	//Arrange
	ef_elem k1 = ef_create_elem(bf_create_elem(0X4168B72, 0), bf_create_elem(0, 0));
	ef_elem k2 = ef_create_elem(bf_create_elem(0X3F28CB71575AF, 0), bf_create_elem(0, 0));

	ef_elem PX = ef_create_elem(bf_create_elem(0X277FDD6F541DF784, 0X3EEE553504EC41D7), bf_create_elem(0X911D456F43BC03E9, 0X73784D7E0FF24D9D));
	ef_elem PL = ef_create_elem(bf_create_elem(0X7002B4AE8DC97103, 0X758EF7198CCF5839), bf_create_elem(0X9830E5CA1DC22270, 0X35183A0FF5A7A2E1));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X2000002000, 0), bf_create_elem(0X2000002000, 0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 133713371337 * GEN

	ef_elem QX = ef_create_elem(bf_create_elem(0XE20478404AAA1B7F, 0X17A4B025EE434AAD), bf_create_elem(0X3AF8DD6E5CE149C3, 0X4DEC004752AFECF9));
	ef_elem QL = ef_create_elem(bf_create_elem(0XA7CFD4415C7C3FE0, 0X142DB23E47B34EF1), bf_create_elem(0X586F58E463F3EEE0, 0X4E8DA1EAC6215EE2));
	ef_elem QZ = ef_create_elem(bf_create_elem(0X210, 0), bf_create_elem(0, 0X40000000));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //Q = 949494 * GEN

	ef_elem EX = ef_create_elem(bf_create_elem(0X54932C0CDE84DD0C, 0X2920B6D886EA380B), bf_create_elem(0XE16091FE1767CA8D, 0X386448A745715798));
	ef_elem EL = ef_create_elem(bf_create_elem(0X974AAEB6BCA972D6, 0X23C05EFC539C393B), bf_create_elem(0X80434D5FE9F454A8, 0X1BB6F877D66F45F5));
	ef_elem EZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
	ec_point_lproj expected = ec_create_point_lproj(EX, EL, EZ);

	//Act
	ec_point_lproj actual = ec_scalarmull_double(P, k1, Q, k2);

	//Assert
	uint64_t correct = ec_equal_point_lproj(expected, actual) && ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(actual);
	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_example FAILED");
}

void ec_scalarmull_double_test_linearity(test_ctr *ctr) {
	//Arrange
	ef_elem k1 = ef_create_elem(bf_create_elem(0XABCDEF123, 0XABCDEF123), bf_create_elem(0XABCDEF123, 0XABCDEF123));
	ef_elem k2 = ef_create_elem(bf_create_elem(0XBCADBCADBCAD, 0XFFFFFFFF1), bf_create_elem(0X77CC77, 0XAAA));

	ef_elem P1X = ef_create_elem(bf_create_elem(0X277FDD6F541DF784, 0X3EEE553504EC41D7), bf_create_elem(0X911D456F43BC03E9, 0X73784D7E0FF24D9D));
	ef_elem P1L = ef_create_elem(bf_create_elem(0X7002B4AE8DC97103, 0X758EF7198CCF5839), bf_create_elem(0X9830E5CA1DC22270, 0X35183A0FF5A7A2E1));
	ef_elem P1Z = ef_create_elem(bf_create_elem(0X2000002000, 0), bf_create_elem(0X2000002000, 0));
	ec_point_lproj P1 = ec_create_point_lproj(P1X, P1L, P1Z); //P1 = 133713371337 * GEN

	ef_elem P2X = ef_create_elem(bf_create_elem(0XE232234BCDA1F7A9, 0X72E20184752C18C3), bf_create_elem(0XA249666BD031EF41,0X78BFF5D2CAFC6FC2));
	ef_elem P2L = ef_create_elem(bf_create_elem(0X1064B65732340CD9, 0X63C8BCD9FBA1773F), bf_create_elem(0X128743886A3EC579,0X6E8E5EC1D3091E7F));
	ef_elem P2Z = ef_create_elem(bf_create_elem(0X212FB30BAC886248, 0X797A835808F9BF57), bf_create_elem(0X8CA3223B872913CE,0X4B0C008D15FE6EB1));
	ec_point_lproj P2 = ec_create_point_lproj(P2X, P2L, P2Z); //P2 = 78632917462800214 * 2 * GEN

	ef_elem P3X = ef_create_elem(bf_create_elem(0XE20478404AAA1B7F, 0X17A4B025EE434AAD), bf_create_elem(0X3AF8DD6E5CE149C3, 0X4DEC004752AFECF9));
	ef_elem P3L = ef_create_elem(bf_create_elem(0XA7CFD4415C7C3FE0, 0X142DB23E47B34EF1), bf_create_elem(0X586F58E463F3EEE0, 0X4E8DA1EAC6215EE2));
	ef_elem P3Z = ef_create_elem(bf_create_elem(0X210, 0), bf_create_elem(0, 0X40000000));
	ec_point_lproj P3 = ec_create_point_lproj(P3X, P3L, P3Z); //P3 = 949494 * GEN

	ef_elem P4X = ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4));
	ef_elem P4L = ef_create_elem(bf_create_elem(0XB55AF6853BAA91FA, 0X368B1A5434DF7331), bf_create_elem(0X4B36A8231E3A83C2, 0XCA5DED000C90427));
	ef_elem P4Z = ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400));
	ec_point_lproj P4 = ec_create_point_lproj(P4X, P4L, P4Z); //P4 = -2243156791409331652485 * GEN

	//Act
	ec_point_lproj result_double = ec_scalarmull_double(ec_add(P1, P2), k1, ec_add(P3, P4), k2);
	ec_point_lproj result_single = ec_add(ec_add(ec_add(ec_scalarmull_single(P1, k1), ec_scalarmull_single(P2, k1)), ec_scalarmull_single(P3, k2)), ec_scalarmull_single(P4, k2));

	//Assert
	uint64_t correct = ec_equal_point_lproj(result_double, result_single) && ec_is_on_curve(result_double);
	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_linearity FAILED");
}

void ec_scalarmull_double_test_commutative(test_ctr *ctr) {
	//Arrange
	ef_elem k1 = ef_create_elem(bf_create_elem(0XAAA, 0XBBB), bf_create_elem(0XCCC, 0XDDD));
	ef_elem k2 = ef_create_elem(bf_create_elem(0XEEE, 0XFFF), bf_create_elem(0X111, 0X222));

	ef_elem PX = ef_create_elem(bf_create_elem(0X277FDD6F541DF784, 0X3EEE553504EC41D7), bf_create_elem(0X911D456F43BC03E9, 0X73784D7E0FF24D9D));
	ef_elem PL = ef_create_elem(bf_create_elem(0X7002B4AE8DC97103, 0X758EF7198CCF5839), bf_create_elem(0X9830E5CA1DC22270, 0X35183A0FF5A7A2E1));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X2000002000, 0), bf_create_elem(0X2000002000, 0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 133713371337 * GEN

	ef_elem QX = ef_create_elem(bf_create_elem(0XE232234BCDA1F7A9, 0X72E20184752C18C3), bf_create_elem(0XA249666BD031EF41,0X78BFF5D2CAFC6FC2));
	ef_elem QL = ef_create_elem(bf_create_elem(0X1064B65732340CD9, 0X63C8BCD9FBA1773F), bf_create_elem(0X128743886A3EC579,0X6E8E5EC1D3091E7F));
	ef_elem QZ = ef_create_elem(bf_create_elem(0X212FB30BAC886248, 0X797A835808F9BF57), bf_create_elem(0X8CA3223B872913CE,0X4B0C008D15FE6EB1));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //Q = 78632917462800214 * 2 * GEN

	//Act
	ec_point_lproj k1P_plus_k2Q = ec_scalarmull_double(P, k1, Q, k2);
	ec_point_lproj k2Q_plus_k1P = ec_scalarmull_double(Q, k2, P, k1);

	//Assert
	uint64_t correct = ec_equal_point_lproj(k1P_plus_k2Q, k2Q_plus_k1P) && ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(k1P_plus_k2Q);
	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_commutative FAILED");
}

void ec_scalarmull_double_test_negation_order_indifference(test_ctr *ctr) {
	//Arrange
	ef_elem k1 = ef_create_elem(bf_create_elem(0XA7A7A, 0XB8B8B), bf_create_elem(0XC9C9C, 0XD0D0D));
	ef_elem k2 = ef_create_elem(bf_create_elem(0XE7E7E, 0XF8F8F), bf_create_elem(0X19191, 0X20202));

	ef_elem PX = ef_create_elem(bf_create_elem(0X7674C426F68A7C0D, 0X26C3E68569307393), bf_create_elem(0X9BFA0D5F1CB2BB3F, 0X53889FE5B08254D3));
	ef_elem PL = ef_create_elem(bf_create_elem(0X4F88EF9F49D18A5E, 0X5C7C38B577B3EAF4), bf_create_elem(0XCDD4DCBE486CC880, 0X18FEF6543ECA3ABC));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X20000000000004, 0), bf_create_elem(0X8000000000000000, 0X80000));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //78632917462800214 * GEN

	ef_elem QX = ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25));
	ef_elem QL = ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB));
	ef_elem QZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , z = u + z

	//Act
	ec_point_lproj neg_first = ec_scalarmull_double(ec_neg(P), k1, ec_neg(Q), k2);
	ec_point_lproj neg_after = ec_neg(ec_scalarmull_double(P, k1, Q, k2));

	//Assert
	uint64_t correct = ec_equal_point_lproj(neg_first, neg_after) && ec_is_on_curve(neg_first);
	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_negation_order_indifference FAILED");
}

void ec_scalarmull_double_test_k_ones_is_add(test_ctr *ctr) {
	//Arrange
	ef_elem k = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));

	ef_elem PX = ef_create_elem(bf_create_elem(0XD2C27333EFC0AE61, 0X4306673487679D76), bf_create_elem(0X909BEC5477E860BB, 0X480D39C8A1B98266));
	ef_elem PL = ef_create_elem(bf_create_elem(0XF84FB0B45D95FC31, 0X24C3FF4B68C78BE3), bf_create_elem(0X963FE2DA0544E1A4, 0X17B6B0A1380A490));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X100, 0), bf_create_elem(0X8000000000000000, 0X4000000000000001));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 99921481365893197563 * GEN

	ef_elem QX = ef_create_elem(bf_create_elem(0XB17380115E6C4F69, 0XE1DA273142A52B3), bf_create_elem(0X529D0294EB09720D, 0X26690AF5D1CC9E78));
	ef_elem QL = ef_create_elem(bf_create_elem(0XA75AABA91400402C, 0X2EB5E7A7EE3E6DEA), bf_create_elem(0X82C9609FFC5E6794, 0X19B8A04EA04D0DF));
	ef_elem QZ = ef_create_elem(bf_create_elem(0X8150BD7B681CEB67, 0X1773D8F08AD460A5), bf_create_elem(0X64E2D12FB4409DDD, 0X391936225DE0A796));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //Q = 14329 * GEN

	//Act
	ec_point_lproj scalarmull_result = ec_scalarmull_double(P, k, Q, k);
	ec_point_lproj add_result = ec_add(P, Q);

	//Assert
	uint64_t correct = ec_equal_point_lproj(scalarmull_result, add_result);
	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_k_ones_is_add FAILED");
}

void ec_scalarmull_double_test_crosscheck_scalarmull_single(test_ctr *ctr) {
	//Arrange
	ef_elem k1 = ef_create_elem(bf_create_elem(12345, 6789), bf_create_elem(10111213, 141516));
	ef_elem k2 = ef_create_elem(bf_create_elem(24690, 13578), bf_create_elem(20222426, 283032)); //k2 = 2 * k1

	ef_elem PX = ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4));
	ef_elem PL = ef_create_elem(bf_create_elem(0XB55AF6853BAA91FA, 0X368B1A5434DF7331), bf_create_elem(0X4B36A8231E3A83C2, 0XCA5DED000C90427));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = -2243156791409331652485 * GEN

	//Act
	ec_point_lproj double_result = ec_scalarmull_double(P, k1, P, k1);
	ec_point_lproj single_result = ec_scalarmull_single(P, k2);

	//Assert
	uint64_t correct = ec_equal_point_lproj(double_result, single_result);
	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_crosscheck_scalarmull_single FAILED");
}

void ec_scalarmull_double_test_point_and_neg_cancel(test_ctr *ctr) {
	//Arrange
	ef_elem k = ef_create_elem(bf_create_elem(999213654, 2084685670), bf_create_elem(32123142, 6462462412477));

	ef_elem PX = ef_create_elem(bf_create_elem(0X277FDD6F541DF784, 0X3EEE553504EC41D7), bf_create_elem(0X911D456F43BC03E9, 0X73784D7E0FF24D9D));
	ef_elem PL = ef_create_elem(bf_create_elem(0X7002B4AE8DC97103, 0X758EF7198CCF5839), bf_create_elem(0X9830E5CA1DC22270, 0X35183A0FF5A7A2E1));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X2000002000, 0), bf_create_elem(0X2000002000, 0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 133713371337 * GEN

	//Act
	ec_point_lproj result = ec_scalarmull_double(P, k, ec_neg(P), k);

	//Assert
	uint64_t correct = ec_equal_point_lproj((ec_point_lproj) INFTY, result);
	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_point_and_neg_cancel FAILED");
}

void ec_scalarmull_double_test_point_and_neg_interfere(test_ctr *ctr) {
	//Arrange
	ef_elem k1 = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
	ef_elem k2 = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(0, 0));

	ef_elem PX = ef_create_elem(bf_create_elem(0X2A955F2AB87B63BD, 0X7B5F4418CA97D542), bf_create_elem(0XDD55878DFFE87C62, 0X1CE397B5452EAAD4));
	ef_elem PL = ef_create_elem(bf_create_elem(0XB55AF6853BAA916A, 0X368B1A5434DF7331), bf_create_elem(0X4B3628231E3A83C2, 0XCA5DED000C90027));
	ef_elem PZ = ef_create_elem(bf_create_elem(0X90, 0), bf_create_elem(0X0000800000000000, 0X400));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //2243156791409331652485 * P

	//Act
	ec_point_lproj P_neg = ec_neg(P);
	ec_point_lproj result = ec_scalarmull_double(P, k1, P_neg, k2);

	//Assert
	uint64_t correct = ec_equal_point_lproj(P_neg, result);
	assert_true(correct, ctr, "ec: ec_scalarmull_double_test_point_and_neg_interfere FAILED");
}

void ec_scalarmull_tests(test_ctr *ctr) {
	ec_scalarmull_single_test_example(ctr);
	ec_scalarmull_single_test_linearity(ctr);
	ec_scalarmull_single_test_negation_order_indifference(ctr);
	ec_scalarmull_single_test_order_of_gen_is_order_of_subgroup(ctr);
	ec_scalarmull_single_test_final_k_at_once_same_as_factor_one_at_a_time(ctr);
	ec_scalarmull_single_test_k_one_is_identity(ctr);

	ec_scalarmull_double_test_example(ctr);
	ec_scalarmull_double_test_linearity(ctr);
	ec_scalarmull_double_test_commutative(ctr);
	ec_scalarmull_double_test_negation_order_indifference(ctr);
	ec_scalarmull_double_test_k_ones_is_add(ctr);
	ec_scalarmull_double_test_crosscheck_scalarmull_single(ctr);
	ec_scalarmull_double_test_point_and_neg_cancel(ctr);
	ec_scalarmull_double_test_point_and_neg_interfere(ctr);
}
