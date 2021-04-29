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
	
	ef_elem EX = ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U));
	ef_elem EL = ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U));
	ef_elem EZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
	ec_point_lproj E = ec_create_point_lproj(EX, EL, EZ); //E = 12345 * GEN

	uint64x2x2_t k = (uint64x2x2_t) {{{12345, 0}, {0,0}}};
	ec_split_scalar decomp = ec_scalar_decomp(k);
	printf("k1: %lu, %lu, sign: %lu\n", decomp.k1[0], decomp.k1[1], decomp.k1_sign);
	printf("k2: %lu, %lu, sign: %lu\n", decomp.k2[0], decomp.k2[1], decomp.k2_sign);
	ec_point_lproj P_proj = (ec_point_lproj) GEN;
	ec_point_laffine P_affine = ec_lproj_to_laffine(P_proj);
	ec_point_laffine Q_affine = ec_endo_laffine(P_affine);
	ec_point_lproj Q_proj = ec_laffine_to_lproj(Q_affine);
	if(decomp.k1_sign) {
		P_proj = ec_neg(P_proj);
		P_affine = ec_lproj_to_laffine(P_proj);
	}
	if(decomp.k2_sign) {
		Q_proj = ec_neg(Q_proj);
		Q_affine = ec_lproj_to_laffine(Q_proj);
	}
	uint64x2x2_t k1 = {{ decomp.k1, {0, 0}}};
	uint64x2x2_t k2 = {{ decomp.k2, {0, 0}}};

	ec_point_lproj R_proj = ec_scalarmull_double(P_proj, k1, Q_proj, k2);
	ec_point_laffine R_affine = ec_lproj_to_laffine(R_proj);
	printf("R.x: %lu, %lu, %lu, %lu\n", R_affine.x.val[0][0], R_affine.x.val[0][1], R_affine.x.val[1][0], R_affine.x.val[1][1]);
	printf("R.l: %lu, %lu, %lu, %lu\n", R_affine.l.val[0][0], R_affine.l.val[0][1], R_affine.l.val[1][0], R_affine.l.val[1][1]);

	printf("Are they equal? %lu\n", ec_equal_point_lproj(R_proj, E));

	dispose_components();
	return 0;
}
