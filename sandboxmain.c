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

	ec_naf result = ec_to_naf(bf_create_elem(153881, 0));
	//ec_naf result = ec_to_naf(bf_create_elem(153881, 1562365363));

	ec_print_naf(result);

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
