#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/utils.h"

int main() {
	/*poly64x2_t a0 = {0, 9223372036854775808U};
	poly64x2_t a1 = {0,0};
	bf_polyx2 a = concat_bf_poly(a0, a1); 
	
	bf_print_expr(a0);
	bf_print_expr(bf_red_formula(a));
	printf("\n");
	
	poly64x2_t b0 = {0, 0};
	poly64x2_t b1 = {1,0};
	bf_polyx2 b = concat_bf_poly(b0, b1); 
	
	printf("z^128\n");
	bf_print_expr(bf_red_formula(b));
	printf("\n");
	
	poly64x2_t c0 = {1,0};
	poly64x2_t c1 = {0, 1729382256910270464U};
	bf_polyx2 c = concat_bf_poly(c0, c1);
	
	printf("z^252 + z^251 + 1\n");
	bf_print_expr(bf_red_formula(c));
	printf("\n");*/
	
	//Arrange
	poly64x2_t a = {98314422310010283, 87159124330183471};
	poly64x2_t b = {3310980040412344311, 1120087461344321001};
	poly64x2_t c = {2276490184612864, 334652613944236613};
	
	//Act
	bf_polyx2 ab;
	ab.p0[0] = 0;
	ab.p0[1] = 0;
	ab.p1[1] = 0;
	ab.p1[0] = 9223372036854775808U; //z^191 so problem is in bit 191
	poly64x2_t ab_generic = bf_red_generic(ab);
	poly64x2_t ab_formula = bf_red_formula(ab);
	printf("ab_generic: ");
	bf_print_expr(ab_generic);
	printf("\n");
	printf("ab_formula: ");
	bf_print_expr(ab_formula);
	printf("\n");
	printf("diff: ");
	bf_print_expr(bf_add(ab_generic, ab_formula));
	printf("\n");
	/*poly64x2_t bc = bf_red(bf_pmull(b,c));
	poly64x2_t ab_times_c = bf_red(bf_pmull(ab, c));
	poly64x2_t a_times_bc = bf_red(bf_pmull(a, bc));*/

	
	return 0;
}
