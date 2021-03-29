#include <stdio.h>

#include "tests/test_tool.h"

#include "tests/benchmark_tool_tests.h"
#include "tests/util_tests.h"

int main() {
	test_ctr *ctr = init_ctr();
	
	benchmark_tool_tests(ctr);
	util_tests(ctr);
	
	print_results(ctr);
	return 0;
}
