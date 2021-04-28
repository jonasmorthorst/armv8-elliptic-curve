#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/extensionfield.h"
#include "common/ec.h"
#include "common/setup.h"
#include "common/utils.h"

int main() {
	init_components();
	
	uint64x2x2_t k = (uint64x2x2_t) {{{18446744073709551615U, 0}, {0,0}}};
	ec_split_scalar res = ec_scalar_decomp(k);
	printf("k1: %lu, %lu, sign: %lu\n", res.k1[0], res.k1[1], res.k1_sign);
	printf("k2: %lu, %lu, sign: %lu\n", res.k2[0], res.k2[1], res.k2_sign);
	
	dispose_components();
	return 0;
}
