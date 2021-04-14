#include <stdio.h>

#include "ecpoint_tests.h"
#include "../common/ecpoint.h"
#include "../common/utils.h"

void ec_create_point_lproj_test_example(test_ctr *ctr) {
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

	//Act
	ec_point_lproj p = ec_create_point_lproj(x, y, z);

  //ec_print_point_lproj_expr(p);

	//Assert
	// uint64_t correct = equal_poly64x2(a.p0, a0) & equal_poly64x2(a.p1, a1);
	// assert_true(correct, ctr, "extensionfield: ef_create_elem_test_example FAILED");
}

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
	ec_point_lproj p = ec_create_point_lproj(x, y, z);


	poly64x2_t a02 = {31245135, 1353};
	poly64x2_t a12 = {153535, 35135};

  ef_elem x2 = ef_create_elem(a02, a12);

  poly64x2_t b02 = {135135, 135};
	poly64x2_t b12 = {315135135, 135135};

  ef_elem y2 = ef_create_elem(b02, b12);

	poly64x2_t c02 = {54689374, 4584};
	poly64x2_t c12 = {34, 314134134};

  ef_elem z2 = ef_create_elem(c02, c12);
	ec_point_lproj p2 = ec_create_point_lproj(x2, y2, z2);

	//Act
	//ec_point_lproj added = ec_point_lproj_add(p, p);
	//ec_point_lproj q = ec_point_lproj_scalar_mull(135235681, p);

  //ec_print_point_lproj_expr(q);

	//Assert
	// uint64_t correct = equal_poly64x2(a.p0, a0) & equal_poly64x2(a.p1, a1);
	// assert_true(correct, ctr, "extensionfield: ef_create_elem_test_example FAILED");
}

void ec_add_infinity_test(test_ctr *ctr) {
	poly64x2_t p_one = {1, 0};
  poly64x2_t p_zero = {0, 0};

  ef_elem ef_one = ef_create_elem(p_one, p_zero);
  ef_elem ef_zero = ef_create_elem(p_zero, p_zero);

  ec_point_lproj infinity = ec_create_point_lproj(ef_one, ef_one, ef_zero);

	ec_point_lproj rand_point = ec_rand_point_lproj();

	ec_point_lproj res = ec_point_lproj_add(infinity, rand_point);

	ec_print_point_lproj_expr(res);
}

void test_generator_on_curve(test_ctr *ctr) {
	//Arrange
	ec_point_lproj gen_elem = (ec_point_lproj) GEN;
	
	//Act
	ef_elem lhs = ef_mull(ef_add(ef_square(gen_elem.l), ef_add(ef_mull(gen_elem.l, gen_elem.z), ef_mull((ef_elem) A, ef_square(gen_elem.z)))), ef_square(gen_elem.x)); //(L^2 + LZ + AZ^2)X^2
	ef_elem rhs = ef_add(ef_square(ef_square(gen_elem.x)), ef_mull((ef_elem) B, ef_square(ef_square(gen_elem.z)))); //X^4 + BZ^4
	
	//Assert
	uint64_t correct = equal_ef_elem(lhs, rhs);
	assert_true(correct, ctr, "ecpoint: test_generator_on_curve FAILED");
}

void ec_rand_point_lproj_test_on_curve(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange & Act
		ec_point_lproj p = ec_rand_point_lproj();
		
		ef_elem lhs = ef_mull(ef_add(ef_square(p.l), ef_add(ef_mull(p.l, p.z), ef_mull((ef_elem) A, ef_square(p.z)))), ef_square(p.x)); //(L^2 + LZ + AZ^2)X^2
		ef_elem rhs = ef_add(ef_square(ef_square(p.x)), ef_mull((ef_elem) B, ef_square(ef_square(p.z)))); //X^4 + BZ^4
		
		//Assert
		correct = equal_ef_elem(lhs, rhs);
		if(!correct) {
			printf("p: ");
			ec_print_point_lproj_hex(p);
			break;
		}
	}
	assert_true(correct, ctr, "ecpoint: ec_rand_point_lproj_test_on_curve FAILED");
}

void ecpoint_tests(test_ctr *ctr) {
	ec_create_point_lproj_test_example(ctr);
	ec_add_infinity_test(ctr);
	ec_scalar_mull_test_example(ctr);
	test_generator_on_curve(ctr);
	ec_rand_point_lproj_test_on_curve(ctr);
}
