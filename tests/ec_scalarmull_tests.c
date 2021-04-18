#include <stdio.h>

#include "ec_scalarmull_tests.h"
#include "../common/ec_scalarmull.h"

void ec_scalar_mull_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t a0 = {16785442, 721476};
	poly64x2_t a1 = {78099554548664, 6547959942615};

	ef_elem x = ef_create_elem(a0, a1);

	poly64x2_t b0 = {972573, 353523623};
	poly64x2_t b1 = {23523988509283, 2435335};

	ef_elem y = ef_create_elem(b0, b1);

	poly64x2_t c0 = {54689374, 4584};
	poly64x2_t c1 = {34548734, 5648574685};

	ef_elem z = ef_create_elem(c0, c1);
	ec_point_lproj P = ec_create_point_lproj(x, y, z);


	poly64x2_t a02 = {31245135, 1353};
	poly64x2_t a12 = {153535, 35135};

	ef_elem x2 = ef_create_elem(a02, a12);

	poly64x2_t b02 = {135135, 135};
	poly64x2_t b12 = {315135135, 135135};

	ef_elem y2 = ef_create_elem(b02, b12);

	poly64x2_t c02 = {54689374, 4584};
	poly64x2_t c12 = {34, 314134134};

	ef_elem z2 = ef_create_elem(c02, c12);
	ec_point_lproj P2 = ec_create_point_lproj(x2, y2, z2);

	//Act
	//ec_point_lproj added = ec_point_lproj_add(p, p);
	ec_point_lproj Q = ec_scalarmull_single(P, 135235681);

	ec_print_expr(Q);

	//Assert
	// uint64_t correct = equal_poly64x2(a.p0, a0) & equal_poly64x2(a.p1, a1);
	// assert_true(correct, ctr, "extensionfield: ef_create_elem_test_example FAILED");
}

void ec_scalarmull_tests(test_ctr *ctr) {
	
}
