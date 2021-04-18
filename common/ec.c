#include <stdio.h>

#include "ec.h"
#include "utils.h"

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem l, ef_elem z) {
	ec_point_lproj P;
	P.x = x;
	P.l = l;
	P.z = z;
	return P;
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

ec_point_lproj ec_rand_point_lproj() {
	return ec_create_point_lproj(ef_rand_elem(), ef_rand_elem(), ef_rand_elem());
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
	ef_elem w = ef_add(ef_square(P.l), ef_add(u, ef_mull((ef_elem) A, v))); //W = L_P^2 + (L_P * Z_P) + A * Z_P^2
	ec_point_lproj R;
	R.x = ef_square(w);
	R.z = ef_mull(w, v);
	R.l = ef_add(ef_add(ef_add(ef_square(ef_mull(P.x, P.z)), R.x), ef_mull(w, u)), R.z);
	return R;
}
