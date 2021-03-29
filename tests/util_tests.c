#include "util_tests.h"
#include "../common/utils.h"

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

void util_tests(test_ctr *ctr) {
	average_test(ctr);
}
