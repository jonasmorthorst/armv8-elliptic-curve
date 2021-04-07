#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/utils.h"

int main() {
	poly64x2_t a0 = {0,0};//{9223372036854775808U, 9223372036854775808U}; //z^127 + z^63
	poly64x2_t a1 = {9223372036854775808U, 0}; //z^128 + z^192
	bf_polyx2 a = concat_bf_poly(a0, a1);
	printf("a0: ");
	bf_print_expr(a0);
	printf("a1: ");
	bf_print_expr(a1);
	poly64x2_t a_redneon = bf_red_neon(a);
	printf("neon: ");
	bf_print_expr(a_redneon);
	poly64x2_t a_redformula = bf_red(a);
	printf("formula: ");
	bf_print_expr(a_redformula);
	
	return 0;
}
