#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/extensionfield.h"
#include "common/setup.h"
#include "common/utils.h"

int main() {
	init_components();
	
	ef_print_expr_nl(ef_rand_elem());
	poly64x2_t one = {1, 0};
	ef_print_expr_nl(ef_create_elem(one, one));
	
	dispose_components();
	return 0;
}
