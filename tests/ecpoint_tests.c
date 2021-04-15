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
	ec_point_lproj q = ec_point_lproj_scalar_mull(135235681, p);

  ec_print_point_lproj_expr(q);

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

void ec_point_lproj_equal_test_equivalent_example(test_ctr *ctr) {
	//Arrange
	ef_elem PX = ef_create_elem(bf_create_elem(8580737671568810207U, 2871255672283442416U), bf_create_elem(14038621212797386049U, 7102795227656941388U));
	ef_elem PL = ef_create_elem(bf_create_elem(8047436918354421319U, 1946646612861660792U), bf_create_elem(3809596628914525439U, 6366755232823307697U));
	ef_elem PZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0,0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 1984 * GEN
	
	ef_elem QX = ef_create_elem(bf_create_elem(3241181585702845695U, 3252528035791615660U), bf_create_elem(12709637355653080861U, 31217557150801189U));
	ef_elem QL = ef_create_elem(bf_create_elem(16970158823458969713U, 7952375805880217921U), bf_create_elem(12883836633214573383U, 8353656882183346347U));
	ef_elem QZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , Q.z = u + z
	
	//Act & Assert
	uint64_t equal = ec_point_lproj_equal(P, Q);
	uint64_t on_curve = ec_is_on_curve(P) && ec_is_on_curve(Q);
	assert_true(equal && on_curve, ctr, "ecpoint: ec_point_lproj_equal_test_equivalent_example FAILED");
}

void ec_point_lproj_equal_test_notequivalent_example(test_ctr *ctr) {
	//Arrange
	ef_elem PX = ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U));
	ef_elem PL = ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U));
	ef_elem PZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 12345 * GEN
	
	ef_elem QX = ef_create_elem(bf_create_elem(3241181585702845695U, 3252528035791615660U), bf_create_elem(12709637355653080861U, 31217557150801189U));
	ef_elem QL = ef_create_elem(bf_create_elem(16970158823458969713U, 7952375805880217921U), bf_create_elem(12883836633214573383U, 8353656882183346347U));
	ef_elem QZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , Q.z = u + z
	
	//Act & Assert
	uint64_t equal = ec_point_lproj_equal(P, Q);
	assert_false(equal, ctr, "ecpoint: ec_point_lproj_equal_test_notequivalent_example FAILED");
}

void ec_is_on_curve_test_generator_on_curve(test_ctr *ctr) {
	//Arrange, Act & Assert
	uint64_t on_curve = ec_is_on_curve((ec_point_lproj) GEN);
	assert_true(on_curve, ctr, "ecpoint: test_generator_on_curve FAILED");
}

void ec_is_on_curve_test_point_not_on_curve(test_ctr *ctr) {
	//Arrange
	ef_elem PX = ef_create_elem(bf_create_elem(6574758758697437212U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U));
	ef_elem PL = ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U));
	ef_elem PZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); 
	
	//Act & Assert
	uint64_t on_curve = ec_is_on_curve((ec_point_lproj) P);
	assert_false(on_curve, ctr, "ecpoint: test_generator_on_curve FAILED");
}


void ec_rand_point_lproj_test_on_curve(test_ctr *ctr) {
	uint64_t correct = 1;
	for(int i = 0; i < 10; i++) {
		//Arrange & Act
		ec_point_lproj P = ec_rand_point_lproj();
		
		//Assert
		correct = ec_is_on_curve(P);
		if(!correct) {
			printf("p: ");
			ec_print_point_lproj_hex(P);
			break;
		}
	}
	assert_true(correct, ctr, "ecpoint: ec_rand_point_lproj_test_on_curve FAILED");
}

void ec_point_lproj_add_test_example(test_ctr *ctr) {
	//Arrange
	ef_elem PX = ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U));
	ef_elem PL = ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U));
	ef_elem PZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 12345 * GEN
	
	ef_elem QX = ef_create_elem(bf_create_elem(3241181585702845695U, 3252528035791615660U), bf_create_elem(12709637355653080861U, 31217557150801189U));
	ef_elem QL = ef_create_elem(bf_create_elem(16970158823458969713U, 7952375805880217921U), bf_create_elem(12883836633214573383U, 8353656882183346347U));
	ef_elem QZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0));
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , z = u + z
	
	ef_elem eX = ef_create_elem(bf_create_elem(11714707759024576024U, 1425171228285341038U), bf_create_elem(8911788759205811087U, 1880974208631851155U));
	ef_elem eL = ef_create_elem(bf_create_elem(11786282500562796401U, 494288005942721389U), bf_create_elem(13537247985355034756U, 7807710330745742853U));
	ef_elem eZ = ef_create_elem(bf_create_elem(10401784470963653389U, 3885991217892589471U), bf_create_elem(7010465213269421323U, 7477155044425948780U));
	ec_point_lproj expected = ec_create_point_lproj(eX, eL, eZ);
	ec_print_point_lproj_expr(expected);
	
	//Act
	ec_point_lproj actual = ec_point_lproj_add(P, Q);
	
	//Assert
	uint64_t equal = ec_point_lproj_equal(expected, actual);
	uint64_t on_curve = ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(actual);
	assert_true(equal && on_curve, ctr, "ecpoint: ec_point_lproj_add_test_example FAILED");
}

void ec_point_lproj_add_test_associative(test_ctr *ctr) {
	//Arrange
	ef_elem PX = ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U));
	ef_elem PL = ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U));
	ef_elem PZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 12345 * GEN
	
	ef_elem QX = ef_create_elem(bf_create_elem(5130258657669722291U, 2683433950625433362U), bf_create_elem(15861652668403718055U, 2280238350963310U));
	ef_elem QL = ef_create_elem(bf_create_elem(13076644468273807311U, 8504358646598259325U), bf_create_elem(16381575749271831681U, 1938714279827322046U));
	ef_elem QZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0)); 
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); //Q = 674559848943297 * GEN
	
	ef_elem RX = ef_create_elem(bf_create_elem(3241181585702845695U, 3252528035791615660U), bf_create_elem(12709637355653080861U, 31217557150801189U));
	ef_elem RL = ef_create_elem(bf_create_elem(16970158823458969713U, 7952375805880217921U), bf_create_elem(12883836633214573383U, 8353656882183346347U));
	ef_elem RZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0));
	ec_point_lproj R = ec_create_point_lproj(RX, RL, RZ); // R = 1984 * GEN , z = u + z
	
	//Act
	ec_point_lproj P_plus_Q = ec_point_lproj_add(P, Q);
	ec_point_lproj P_plus_Q_first = ec_point_lproj_add(P_plus_Q, R);
	ec_point_lproj Q_plus_R = ec_point_lproj_add(Q, R);
	ec_point_lproj Q_plus_R_first = ec_point_lproj_add(P, Q_plus_R);
	
	//Assert
	uint64_t equal = ec_point_lproj_equal(P_plus_Q_first, Q_plus_R_first);
	uint64_t on_curve = ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(R) && ec_is_on_curve(P_plus_Q) && ec_is_on_curve(P_plus_Q_first) && ec_is_on_curve(Q_plus_R);
	assert_true(equal && on_curve, ctr, "ecpoint ec_point_lproj_add_test_associative FAILED");
}

void ec_point_lproj_add_test_commutative(test_ctr *ctr) {
	//Arrange
	ef_elem PX = ef_create_elem(bf_create_elem(6574758758697437213U, 9029938885770167062U), bf_create_elem(2238988517843169761U, 2100850367268144132U));
	ef_elem PL = ef_create_elem(bf_create_elem(6503207926092968936U, 3219845272070329794U), bf_create_elem(10930336994397273494U, 8947327516344479714U));
	ef_elem PZ = ef_create_elem(bf_create_elem(1, 0), bf_create_elem(0, 0));
	ec_point_lproj P = ec_create_point_lproj(PX, PL, PZ); //P = 12345 * GEN
	
	ef_elem QX = ef_create_elem(bf_create_elem(3241181585702845695U, 3252528035791615660U), bf_create_elem(12709637355653080861U, 31217557150801189U));
	ef_elem QL = ef_create_elem(bf_create_elem(16970158823458969713U, 7952375805880217921U), bf_create_elem(12883836633214573383U, 8353656882183346347U));
	ef_elem QZ = ef_create_elem(bf_create_elem(2, 0), bf_create_elem(1, 0)); 
	ec_point_lproj Q = ec_create_point_lproj(QX, QL, QZ); // Q = 1984 * GEN , z = u + z
	
	//Act
	ec_point_lproj P_plus_Q = ec_point_lproj_add(P, Q);
	ec_point_lproj Q_plus_P = ec_point_lproj_add(Q, P);
	
	//Assert
	uint64_t equal = ec_point_lproj_equal(P_plus_Q, Q_plus_P);
	uint64_t on_curve = ec_is_on_curve(P) && ec_is_on_curve(Q) && ec_is_on_curve(P_plus_Q);
	assert_true(equal && on_curve, ctr, "ecpoint: ec_point_lproj_add_test_commutative FAILED");
}

void ecpoint_tests(test_ctr *ctr) {
	//ec_create_point_lproj_test_example(ctr);
	//ec_add_infinity_test(ctr);
	//ec_scalar_mull_test_example(ctr);
	ec_point_lproj_equal_test_equivalent_example(ctr);
	ec_point_lproj_equal_test_notequivalent_example(ctr);
	ec_is_on_curve_test_generator_on_curve(ctr);
	ec_is_on_curve_test_point_not_on_curve(ctr);
	ec_rand_point_lproj_test_on_curve(ctr);
	ec_point_lproj_add_test_example(ctr);
	ec_point_lproj_add_test_associative(ctr);
	ec_point_lproj_add_test_commutative(ctr);
}
