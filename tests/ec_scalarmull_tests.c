#include <stdio.h>

#include "ec_scalarmull_tests.h"
#include "../common/ec_scalarmull.h"
#include "../common/ec.h"
#include "../common/extensionfield.h"
#include "../common/basefield.h"

void ec_scalar_mull_test_example(test_ctr *ctr) {
	//Arrange
	poly64x2_t k0 = {1984, 0};
	poly64x2_t k1 = {0, 0};
	ef_elem k = ef_create_elem(k0, k1);

	ef_elem QX = ef_create_elem(bf_create_elem(0X2CFAFD7ACC5AFCFF, 0X2D234D04135BB6AC), bf_create_elem(0XB061B3D61FBCA71D, 0X6EE833ECBA9D25));
	ef_elem QL = ef_create_elem(bf_create_elem(0XEB821D89C63B9871, 0X6E5C83A975E6B141), bf_create_elem(0XB2CC95280AEB5B47, 0X73EE26ACBE0918AB));
	ef_elem QZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0));
	ec_point_lproj expected = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN

	//Act
	ec_point_lproj actual = ec_scalarmull_single((ec_point_lproj)GEN, k);

	//Assert
	uint64_t equal = ec_equal_point_lproj(expected, actual);
	uint64_t on_curve = ec_is_on_curve(actual);
	assert_true(equal && on_curve, ctr, "ec: ec_scalar_mull_test_example FAILED");
}

void ec_scalarmull_tests(test_ctr *ctr) {
	ec_scalar_mull_test_example(ctr);
}
