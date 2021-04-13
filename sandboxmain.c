#include <arm_neon.h>
#include <stdio.h>

#include "common/basefield.h"
#include "common/extensionfield.h"
#include "common/utils.h"

int main() {
	ef_print_expr_nl(ef_rand_elem());
	poly64x2_t one = {1, 0};
	ef_print_expr_nl(create_elem(one, one));
	return 0;
}
