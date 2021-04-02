#include "util_tests.h"
#include "../common/utils.h"

void compare_doubles_test_equal(test_ctr *ctr) {
	//Arrange
	double a = 9.3;
	double b = 9.3;
	double errmargin = 0.0000001;
	
	//Act & Assert
	uint64_t are_equal = compare_doubles(a, b, errmargin) == 0;
	assert_true(are_equal, ctr, "util: compare_doubles_test_equal FAILED");
}

void compare_doubles_test_approxequal(test_ctr *ctr) {
	//Arrange
	double a = -9.3;
	double b = -9.30000001;
	double errmargin = 0.0000001;
	
	//Act & Assert
	uint64_t are_equal = compare_doubles(a, b, errmargin) == 0;
	assert_true(are_equal, ctr, "util: compare_doubles_test_approxequal FAILED");
}

void compare_doubles_test_notequal(test_ctr *ctr) {
	//Arrange
	double a = 9.3;
	double b = 9.3000002;
	double errmargin = 0.0000001;
	
	//Act & Assert
	uint64_t are_equal = compare_doubles(a, b, errmargin) == 0;
	assert_false(are_equal, ctr, "util: compare_doubles_test_notequal FAILED");
}

void compare_doubles_test_lessthan(test_ctr *ctr) {
	//Arrange
	double a = 9.3;
	double b = 9.3000002;
	double errmargin = 0.0000001;
	
	//Act & Assert
	uint64_t b_greater = compare_doubles(a, b, errmargin) == -1;
	assert_true(b_greater, ctr, "util: compare_doubles_test_lessthan FAILED");
}

void compare_doubles_test_greaterthan(test_ctr *ctr) {
	//Arrange
	double a = -13.111111;
	double b = -13.111112;
	double errmargin = 0.0000001;
	
	//Act & Assert
	uint64_t a_greater = compare_doubles(a, b, errmargin) == 1;
	assert_true(a_greater, ctr, "util: compare_doubles_test_greaterthan FAILED");
}

void average_test(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = { 15, 2, 9};
	double expected = (15.0+2.0+9.0)/3.0;
	double errmargin = 0.0000001;
	
	//Act
	double actual = average(times, 3);
	
	//Assert
	uint64_t in_err_interval = compare_doubles(expected, actual, errmargin) == 0;
	assert_true(in_err_interval, ctr, "util: average_test FAILED");
}

void median_test_even(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = {2, 15, 17, 18, 23, 36};
	double expected = 17.5;
	double errmargin = 0.0000001;
	
	//Act
	double actual = median(times, 6);
	
	//Assert
	uint64_t in_err_interval = compare_doubles(expected, actual, errmargin) == 0;
	assert_true(in_err_interval, ctr, "util: median_test_even FAILED");
}

void median_test_odd(test_ctr *ctr) {
	//Arrange
	uint64_t times[] = {2, 15, 17, 18, 20, 23, 36};
	double expected = 18;
	double errmargin = 0.0000001;
	
	//Act
	double actual = median(times, 7);
	
	//Assert
	uint64_t in_err_interval = compare_doubles(expected, actual, errmargin) == 0;
	assert_true(in_err_interval, ctr, "util: median_test_odd FAILED");
}

void util_tests(test_ctr *ctr) {
	compare_doubles_test_equal(ctr);
	compare_doubles_test_approxequal(ctr);
	compare_doubles_test_notequal(ctr);
	compare_doubles_test_lessthan(ctr);
	compare_doubles_test_greaterthan(ctr);
	average_test(ctr);
	median_test_even(ctr);
	median_test_odd(ctr);
}
