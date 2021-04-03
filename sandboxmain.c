#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"

int main() {
	poly64x2_t a = {2, 0};
	poly64x2_t inva = bf_inv(a);
	
	bf_print_expr(a);
	bf_print_expr(inva);
	
	return 0;
}
