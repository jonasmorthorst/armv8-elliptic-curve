#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/extensionfield.h"
#include "common/setup.h"
#include "common/utils.h"

int main() {
	init_components();
	
	uint64x2_t result = mult_u64(18446744073709551615U, 18446744073709551615U);
	printf("%lu, %lu \n", result[0], result[1]);
	
	dispose_components();
	return 0;
}
