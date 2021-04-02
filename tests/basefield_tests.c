#include "basefield_tests.h"
#include "../common/basefield.h"

void bf_add_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {4, 756};
	poly64x2_t b = {15, 2};
	poly64x2_t expected = {11, 758};
	
	//Act
	poly64x2_t actual = bf_add(a, b);
	
	//Assert
	uint64_t equal = (actual[0] == expected[0]) && (actual[1] == expected[1]); 
	assert_true(equal, ctr, "basefield: bf_add_test_example FAILED");
}

void bf_add_doubling_is_zero(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {638540023, 896234759148761U};
	
	//Act
	poly64x2_t result = bf_add(a, a);
	
	//Assert
	uint64_t is_zero = (result[0] == 0) && (result[1] == 0);
	assert_true(is_zero, ctr, "basefield: bf_add_doubling_is_zero FAILED");
}

void bf_add_zero_is_identity(test_ctr *ctr) {
	//Arrange
	poly64x2_t a = {129548752, 8236754276};
	poly64x2_t zero = {0,0};
	
	//Act
	poly64x2_t result = bf_add(a, zero);
	
	//Assert
	uint64_t is_identity = (result[0] == a[0]) && (result[1] == a[1]);
	assert_true(is_identity, ctr, "basefield: bf_add_zero_is_identity FAILED");
}

void basefield_tests(test_ctr *ctr) {
	bf_add_test_example(ctr);
	bf_add_doubling_is_zero(ctr);
	bf_add_zero_is_identity(ctr);
}
