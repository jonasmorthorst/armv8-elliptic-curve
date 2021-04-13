#include <stdio.h>

#include "extensionfield_tests.h"
#include "../common/extensionfield.h"
#include "../common/utils.h"

void ef_create_elem_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {16785442, 721476};
	poly64x2_t a1 = {78099554548664, 6547959942615};
	
	//Act
	ef_elem a = ef_create_elem(a0, a1);
	
	//Assert
	uint64_t correct = equal_poly64x2(a.p0, a0) & equal_poly64x2(a.p1, a1);
	assert_true(correct, ctr, "extensionfield: ef_create_elem_test_example FAILED");
}

void ef_add_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {31, 7};
	poly64x2_t a1 = {0, 64};
	ef_elem a = ef_create_elem(a0, a1);
	
	poly64x2_t b0 = {9, 13};
	poly64x2_t b1 = {1, 64};
	ef_elem b = ef_create_elem(b0, b1);
	
	poly64x2_t e0 = {22, 10};
	poly64x2_t e1 = {1, 0};
	ef_elem expected = ef_create_elem(e0, e1);
	
	//Act
	ef_elem actual = ef_add(a, b);
	
	//Assert
	uint64_t correct = equal_ef_elem(expected, actual);
	assert_true(correct, ctr, "extensionfield: ef_add_test_example FAILED");
}

void ef_add_test_doubling_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {6751410153, 112308523};
	poly64x2_t a1 = {998742121496553186, 718656530};
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t zero0 = {0,0};
	ef_elem zero = ef_create_elem(zero0, zero0);
	
	//Act
	ef_elem actual = ef_add(a, a);
	
	//Assert
	uint64_t is_zero = equal_ef_elem(zero, actual);
	assert_true(is_zero, ctr, "extensionfield: ef_add_test_doubling_is_zero FAILED");
}

void ef_add_test_zero_is_identity(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {233207851, 6715398730024};
	poly64x2_t a1 = {67651597459761, 12435676};
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t zero0 = {0,0};
	ef_elem zero = ef_create_elem(zero0, zero0);
	
	//Act
	ef_elem result = ef_add(a, zero);
	
	//Assert
	uint64_t correct = equal_ef_elem(a, result);
	assert_true(correct, ctr, "extensionfield: ef_add_test_zero_is_identity FAILED");
}

void ef_add_test_associative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {7651232241, 2233445678};
	poly64x2_t a1 = {101, 59874};
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t b0 = {149073,8569319135};
	poly64x2_t b1 = {2287569446512,433318594};
	ef_elem b = ef_create_elem(b0, b1);
	poly64x2_t c0 = {1786054612, 786156412};
	poly64x2_t c1 = {54364258769123, 3521913758};
	ef_elem c = ef_create_elem(c0, c1);
	
	//Act
	ef_elem a_plus_b_first = ef_add(ef_add(a, b), c);
	ef_elem b_plus_c_first = ef_add(a, ef_add(b, c));
	
	//Assert
	uint64_t equal = equal_ef_elem(a_plus_b_first, b_plus_c_first);
	assert_true(equal, ctr, "extensionfield: ef_add_test_associative FAILED");
}

void ef_add_test_commutative(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {9785313565, 7656548341100};
	poly64x2_t a1 = {94271531642, 2176487654};
	ef_elem a = ef_create_elem(a0, a1);
	poly64x2_t b0 = {78695976431414, 5677568658888};
	poly64x2_t b1 = {77012835631, 2176456423658};
	ef_elem b = ef_create_elem(b0, b1);
	
	//Act
	ef_elem a_plus_b = ef_add(a, b);
	ef_elem b_plus_a = ef_add(b, a);
	
	//Assert
	uint64_t equal = equal_ef_elem(a_plus_b, b_plus_a);
	assert_true(equal, ctr, "extensionfield: ef_add_test_commutative FAILED");
}
	
void ef_mull_test_example(test_ctr *ctr) {
	
}

void ef_mull_test_associative(test_ctr *ctr) {
	
}

void ef_mull_test_associative_rnd(test_ctr *ctr) {
	
}

void ef_mull_test_commutative(test_ctr *ctr) {
	
}

void ef_mull_test_commutative_rnd(test_ctr *ctr) {
	
}

void ef_mull_test_one_is_identity(test_ctr *ctr) {
	
}

void ef_mull_test_zero_is_zero(test_ctr *ctr) {
	
}
	
void ef_square_test_example(test_ctr *ctr) {
	
}

void ef_square_test_every_possible_term(test_ctr *ctr) {
	
}

void ef_square_test_one_is_one(test_ctr *ctr) {
	
}

void ef_square_test_zero_is_zero(test_ctr *ctr) {
	
}

void ef_square_test_crosscheck_pmull(test_ctr *ctr) {
	
}

void ef_square_test_freshmans_dream(test_ctr *ctr) {
	
}

void ef_square_test_freshmans_dream_rnd(test_ctr *ctr) {
	
}

void ef_inv_test_example(test_ctr *ctr) {
	
}

void ef_inv_test_inverse_of_one_is_one(test_ctr *ctr) {
	
}

void ef_inv_zero_outputs_zero(test_ctr *ctr) {
	
}

void ef_inv_test_inverse_of_inverse_is_original(test_ctr *ctr) {
	
}

void ef_inv_test_inverse_of_inverse_is_original_rnd(test_ctr *ctr) {
	
}

void ef_inv_test_prod_of_inverses_is_inverse_of_prod(test_ctr *ctr) {
	
}
	
void ef_inv_test_prod_of_inverses_is_inverse_of_prod_rnd(test_ctr *ctr) {
	
}

void ef_inv_test_prod_with_inv_is_one(test_ctr *ctr) {
	
}

void ef_inv_test_prod_with_inv_is_one_rnd(test_ctr *ctr) {
	
}

void extensionfield_tests(test_ctr *ctr) {
	ef_create_elem_test_example(ctr);
	
	ef_add_test_example(ctr);
	ef_add_test_doubling_is_zero(ctr);
	ef_add_test_zero_is_identity(ctr);
	ef_add_test_associative(ctr);
	ef_add_test_commutative(ctr);
	
	ef_mull_test_example(ctr);
	ef_mull_test_associative(ctr);
	ef_mull_test_associative_rnd(ctr);
	ef_mull_test_commutative(ctr);
	ef_mull_test_commutative_rnd(ctr);
	ef_mull_test_one_is_identity(ctr);
	ef_mull_test_zero_is_zero(ctr);
	
	ef_square_test_example(ctr);
	ef_square_test_every_possible_term(ctr);
	ef_square_test_one_is_one(ctr);
	ef_square_test_zero_is_zero(ctr);
	ef_square_test_crosscheck_pmull(ctr);
	ef_square_test_freshmans_dream(ctr);
	ef_square_test_freshmans_dream_rnd(ctr);
	
	ef_inv_test_example(ctr);
	ef_inv_test_inverse_of_one_is_one(ctr);
	ef_inv_zero_outputs_zero(ctr);
	ef_inv_test_inverse_of_inverse_is_original(ctr);
	ef_inv_test_inverse_of_inverse_is_original_rnd(ctr);
	ef_inv_test_prod_of_inverses_is_inverse_of_prod(ctr);
	ef_inv_test_prod_of_inverses_is_inverse_of_prod_rnd(ctr);
	ef_inv_test_prod_with_inv_is_one(ctr);
	ef_inv_test_prod_with_inv_is_one_rnd(ctr);
}
