#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"

int main() {
	pmullres a;
	poly64x2_t a0 = {1, 4611686018427387905U};
	bf_print_expr(a0);
	poly64x2_t a1 = {2305843009213693952U,0};
	bf_print_expr(a1);
	a.p0 = a0;
	a.p1 = a1;
	poly64x2_t ared = bf_red(a);
	bf_print_expr(ared);
	
	return 0;
}
