#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/utils.h"

int main() {
	poly64x2_t a = {7, 4611686018427387904}; //z^126 + z^2+z+1
	bf_polyx2 a_squared = bf_psquare(a);
	
	poly64x2_t reduced_neon = bf_red_psquare_neon(a_squared);
	printf("neon: ");
	bf_print_expr(reduced_neon);
	poly64x2_t reduced_formula = bf_red_psquare_formula(a_squared);
	printf("formula: ");
	bf_print_expr(reduced_formula);

	return 0;
}
