#include <stdio.h>

#include "ec.h"
#include "utils.h"
#include "ec_scalarmull.h"

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem l, ef_elem z) {
	ec_point_lproj P;
	P.x = x;
	P.l = l;
	P.z = z;
	return P;
}

ec_point_laffine ec_create_point_laffine(ef_elem x, ef_elem l) {
	ec_point_laffine P;
	P.x = x;
	P.l = l;
	return P;
}

ec_point_lproj ec_point_laffine_to_point_lproj(ec_point_laffine P) {
	ec_point_lproj R;
	R.x = P.x;
	R.l = P.l;
	R.z = (ef_elem) {{{1, 0}, {0, 0}}};
	return R;
}

ec_point_laffine ec_point_lproj_to_point_laffine(ec_point_lproj P) {
	ef_elem Z_inv = ef_inv(P.z);
	ec_point_laffine R;
	R.x = ef_mull(P.x, Z_inv);
	R.l = ef_mull(P.l, Z_inv);
	return R;
}

void ec_print_expr(ec_point_lproj P) {
	printf("x: ");
	ef_print_expr_nl(P.x);
	printf(" l: ");
	ef_print_expr_nl(P.l);
	printf(" z: ");
	ef_print_expr_nl(P.z);
}

void ec_print_hex(ec_point_lproj P) {
	printf("x: ");
	ef_print_hex_nl(P.x);
	printf(" l: ");
	ef_print_hex_nl(P.l);
	printf(" z: ");
	ef_print_hex_nl(P.z);
}

void ec_print_hex_laffine(ec_point_laffine P) {
	printf("x: ");
	ef_print_hex_nl(P.x);
	printf(" l: ");
	ef_print_hex_nl(P.l);
}

uint64_t ec_is_on_curve(ec_point_lproj P) {
	ef_elem lhs = ef_mull(ef_add(ef_square(P.l), ef_add(ef_mull(P.l, P.z), ef_mull((ef_elem) A, ef_square(P.z)))), ef_square(P.x)); //(L^2 + LZ + AZ^2)X^2
	ef_elem rhs = ef_add(ef_square(ef_square(P.x)), ef_mull((ef_elem) B, ef_square(ef_square(P.z)))); //X^4 + BZ^4
	return equal_ef_elem(lhs, rhs);
}

uint64_t ec_equal_point_lproj(ec_point_lproj P, ec_point_lproj Q) {
	ef_elem zP_inv = ef_inv(P.z);
	ef_elem xP_normalized = ef_mull(P.x, zP_inv);
	ef_elem lP_normalized = ef_mull(P.l, zP_inv);
	ef_elem zQ_inv = ef_inv(Q.z);
	ef_elem xQ_normalized = ef_mull(Q.x, zQ_inv);
	ef_elem lQ_normalized = ef_mull(Q.l, zQ_inv);

	return equal_ef_elem(xP_normalized, xQ_normalized) && equal_ef_elem(lP_normalized, lQ_normalized);
}

uint64_t ec_equal_point_mixed(ec_point_laffine P, ec_point_lproj Q) {
	ec_point_laffine Q_affine = ec_point_lproj_to_point_laffine(Q);
	return equal_ef_elem(P.x, Q_affine.x) && equal_ef_elem(P.l, Q_affine.l);
}

// Generate random number in range [1, ORDER-1]
poly64x2x2_t ec_rand_scalar() {
	poly64x2x2_t order = (poly64x2x2_t) SUBGROUP_ORDER;
	poly64x2x2_t k;

	int in_range = 0;
	while (!in_range) {
		poly64_t a0 = rand_uint64();
		poly64_t a1 = rand_uint64();
		poly64_t a2 = rand_uint64();
		poly64_t a3 = rand_uint64();

		if (a0 > order.val[1][1]) continue;
		if (a0 == order.val[1][1] && a1 > order.val[1][0]) continue;
		if (a0 == order.val[1][1] && a1 == order.val[1][0] && a2 > order.val[0][1]) continue;
		if (a0 == order.val[1][1] && a1 == order.val[1][0] && a2 == order.val[0][1] && a3 >= order.val[0][0]) continue;

		in_range = 1;

		poly64x2_t p1 = { a0, a1 };
		poly64x2_t p2 = { a2, a3 };

		k.val[0] = p1;
		k.val[1] = p2;
	}

	return k;
}

ec_point_lproj ec_rand_point_lproj() {
	poly64x2x2_t k = ec_rand_scalar();

	return ec_scalarmull_single_lproj((ec_point_lproj) GEN, k);
}

ec_point_laffine ec_rand_point_laffine() {
	ec_point_laffine L;
	ec_point_lproj P = ec_rand_point_lproj();

	L.x = ef_mull(P.x, ef_inv(P.z));
	L.l = ef_mull(P.l, ef_inv(P.z));

	return L;
}

ec_point_lproj ec_neg(ec_point_lproj P) {
	P.l = ef_add(P.l, P.z);

	return P;
}

// Non constant implementation atm.
ec_point_lproj ec_add(ec_point_lproj P, ec_point_lproj Q) {
	if(ec_equal_point_lproj(P, (ec_point_lproj) INFTY)) {
		return Q;
	}
	if(ec_equal_point_lproj(Q, (ec_point_lproj) INFTY)) {
		return P;
	}
	if(ec_equal_point_lproj(P, ec_neg(Q))) {
		return (ec_point_lproj) INFTY;
	}
	if(ec_equal_point_lproj(P, Q)) {
		return ec_double(P);
	}

	ef_elem u = ef_add(ef_mull(P.l, Q.z), ef_mull(Q.l, P.z)); // U = L_P * Z_Q + L_Q * Z_P
	ef_elem w1 = ef_mull(P.x, Q.z); //W1 = X_P * Z_Q
	ef_elem w2 = ef_mull(Q.x, P.z); //W2 = X_Q * Z_P
	ef_elem v = ef_square(ef_add(w1, w2)); //V = (X_P * Z_Q + X_Q * Z_P)^2
	ef_elem w3 = ef_mull(u, w2); //W3 = U * X_Q * Z_P
	ef_elem w4 = ef_mull(u, ef_mull(v, Q.z)); //W4 = U * V * Z_Q
	ec_point_lproj R;
	R.x = ef_mull(u, ef_mull(w1, w3));
	R.l = ef_add(ef_square(ef_add(w3, v)), ef_mull(w4, ef_add(P.l, P.z)));
	R.z = ef_mull(w4, P.z);
	return R;
}

ec_point_lproj ec_add_mixed(ec_point_laffine P, ec_point_lproj Q) {
	if(ec_equal_point_lproj(Q, (ec_point_lproj) INFTY)) {
		return ec_point_laffine_to_point_lproj(P);
	}
	if(ec_equal_point_mixed(P, ec_neg(Q))) {
		return (ec_point_lproj) INFTY;
	}
	if(ec_equal_point_mixed(P, Q)) {
		return ec_double(Q);
	}

	ef_elem E = ef_add(ef_mull(P.l, Q.z), Q.l); //A = L_P * Z_Q + L_Q
	ef_elem F = ef_mull(P.x, Q.z); //X_P * Z_Q
	ef_elem G = ef_square(ef_add(F, Q.x)); //B = (X_P * Z_Q + X_Q)^2
	ef_elem H = ef_mull(E, Q.x); //A * Q.x
	//P.l.val[0] = bf_add(P.l.val[0], (poly64x2_t) {1, 0});
	ec_point_lproj R;
	R.x = ef_mull(ef_mull(E, F), H); //A * (X_P * Z_Q) * A * Q.x
	R.z = ef_mull(ef_mull(E, G), Q.z); //A * B * Z_Q
	R.l = ef_add(ef_square(ef_add(G, H)), ef_mull(R.z, ef_add(P.l, (ef_elem) {{{1, 0}, {0, 0}}}))); //(G+H)^2 + R.z + R.z * P.l
	return R;
}

/*
 * I think I have found a proof for why we don't need to check that P = -P here.
 * It seems that our projective lambda coords do not allow finite points of order 2:
 * 2*P = infty <-> P = -P <-> (x, l, z) = (x, l+z, z) => z=0
 *
 * For the case of the point at infty,
 * a quick calculation shows that this alg outputs infty again!
 */
ec_point_lproj ec_double(ec_point_lproj P) {
	ef_elem u = ef_mull(P.l, P.z); //U = L_P * Z_P
	ef_elem v = ef_square(P.z); //V = Z_P^2
	ef_elem w = ef_add(ef_square(P.l), ef_add(u, ef_mull_A(v))); //W = L_P^2 + (L_P * Z_P) + A * Z_P^2
	ec_point_lproj R;
	R.x = ef_square(w);
	R.z = ef_mull(w, v);
	R.l = ef_add(ef_add(ef_add(ef_square(ef_mull(P.x, P.z)), R.x), ef_mull(w, u)), R.z);
	return R;
}

ec_point_lproj ec_double_then_add(ec_point_laffine P, ec_point_lproj Q) {
	ef_elem LQ_sqr = ef_square(Q.l);
	ef_elem ZQ_sqr = ef_square(Q.z);
	ef_elem T = ef_add(ef_add(LQ_sqr, ef_mull(Q.l, Q.z)), ef_mull_A(ZQ_sqr));

	ef_elem one = (ef_elem) {{{1, 0}, {0, 0}}};

	ef_elem E = ef_add(ef_mull(ef_square(Q.x), ZQ_sqr), ef_mull(T,ef_add(LQ_sqr, ef_mull(ef_add(ef_add((ef_elem)A, one), P.l), ZQ_sqr))));
	ef_elem F = ef_mull(P.x, ZQ_sqr);
	ef_elem G = ef_square(ef_add(F, T));

	ec_point_lproj R;
	R.x = ef_mull(F, ef_square(E));
	R.z = ef_mull(ef_mull(E, G), ZQ_sqr);
	R.l = ef_add(ef_mull(T, ef_square(ef_add(E, G))), ef_mull(ef_add(P.l, one), R.z));
	return R;
}

ec_point_lproj ec_double_then_addtwo(ec_point_laffine P1, ec_point_laffine P2, ec_point_lproj Q) {
	ef_elem one = (ef_elem) {{{1, 0}, {0, 0}}};
	ef_elem LP1_plus_1 = ef_add(P1.l, one);
	ef_elem LQ_sqr = ef_square(Q.l);
	ef_elem ZQ_sqr = ef_square(Q.z);
	ef_elem T = ef_add(ef_add(LQ_sqr, ef_mull(Q.l, Q.z)), ef_mull_A(ZQ_sqr)); //L_Q^2 + L_Q*Z_Q + A*Z_Q^2
	ef_elem E = ef_add(ef_mull(ef_square(Q.x), ZQ_sqr), ef_mull(T, ef_add(LQ_sqr, ef_mull(ef_add((ef_elem) A, LP1_plus_1), ZQ_sqr)))); //X_Q^2*Z_Q^2 + T * (L_Q^2 + (a + L_P1 + 1)*Z_Q^2)
	ef_elem F = ef_mull(P1.x, ZQ_sqr); //X_P1 * Z_Q^2
	ef_elem G = ef_square(ef_add(F, T)); //(F+T)^2
	ef_elem H = ef_mull(ef_square(E), F); //X_2Q+P1 = F * E^2
	ef_elem I = ef_mull(ef_mull(E, G), ZQ_sqr); //Z_2Q+P1 = E * G * Z_Q^2
	ef_elem J = ef_add(ef_mull(ef_add(LP1_plus_1, P2.l), I), ef_mull(T, ef_square(ef_add(E, G)))); //(L_P1 + L_P2 + 1)*I+T*(E+G)^2
	ef_elem K = ef_mull(P2.x, I); // X_P2 * I
	ef_elem L = ef_square(ef_add(H, K)); //(H+K)^2
	ef_elem M = ef_mull(H, J); //H * J
	ec_point_lproj R;
	R.x = ef_mull(ef_mull(J, K), M);
	R.z = ef_mull(ef_mull(I, J), L);
	R.l = ef_add(ef_square(ef_add(L, M)), ef_mull(R.z, ef_add(P2.l, one)));
	return R;
}
