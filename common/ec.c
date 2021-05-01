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

ec_point_lproj ec_laffine_to_lproj(ec_point_laffine P) {
	ec_point_lproj R;
	R.x = P.x;
	R.l = P.l;
	R.z = (ef_elem) {{{1, 0}, {0, 0}}};
	return R;
}

ec_point_laffine ec_lproj_to_laffine(ec_point_lproj P) {
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

void ec_print_expr_laffine(ec_point_laffine P) {
	printf("x: ");
	ef_print_expr_nl(P.x);
	printf(" l: ");
	ef_print_expr_nl(P.l);
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
	ec_point_laffine Q_affine = ec_lproj_to_laffine(Q);
	return equal_ef_elem(P.x, Q_affine.x) && equal_ef_elem(P.l, Q_affine.l);
}

// Generate random number in range [1, ORDER-1]
uint64x2x2_t ec_rand_scalar() {
	uint64x2x2_t order = (uint64x2x2_t) SUBGROUP_ORDER;
	uint64x2x2_t k;

	int in_range = 0;
	while (!in_range) {
		uint64_t a0 = rand_uint64();
		uint64_t a1 = rand_uint64();
		uint64_t a2 = rand_uint64();
		uint64_t a3 = rand_uint64();

		if (a0 > order.val[1][1]) continue;
		if (a0 == order.val[1][1] && a1 > order.val[1][0]) continue;
		if (a0 == order.val[1][1] && a1 == order.val[1][0] && a2 > order.val[0][1]) continue;
		if (a0 == order.val[1][1] && a1 == order.val[1][0] && a2 == order.val[0][1] && a3 >= order.val[0][0]) continue;

		in_range = 1;

		uint64x2_t p1 = { a0, a1 };
		uint64x2_t p2 = { a2, a3 };

		k.val[0] = p1;
		k.val[1] = p2;
	}

	return k;
}

ec_point_lproj ec_rand_point_lproj() {
	uint64x2x2_t k = ec_rand_scalar();

	return ec_scalarmull_single_lproj((ec_point_lproj) GEN, k);
}

ec_point_laffine ec_rand_point_laffine() {
	ec_point_lproj P = ec_rand_point_lproj();
	return ec_lproj_to_laffine(P);
}

ec_point_lproj ec_neg(ec_point_lproj P) {
	P.l = ef_add(P.l, P.z);
	return P;
}

ec_point_laffine ec_neg_laffine(ec_point_laffine P) {
	P.l.val[0] = bf_add(P.l.val[0], (poly64x2_t) {1,0});
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
		return ec_laffine_to_lproj(P);
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

//We will only need this endomorphism for lambda affine coords for scalar mul GLV trick
ec_point_laffine ec_endo_laffine(ec_point_laffine P) {
	P.x.val[0] = bf_add(P.x.val[0], P.x.val[1]);
	P.l.val[0] = bf_add(P.l.val[0], P.l.val[1]);
	P.l.val[1] = bf_add(P.l.val[1], (poly64x2_t) {1,0});
	return P;
}

//c += a	  
#define ADDACC_128(a0, a1, c0, c1)\
	asm ("ADDS %0, %0, %2;"\
				  "ADC %1, %1, %3;"\
				  : "+r" (c0), "+r" (c1)\
				  : "r" (a0), "r" (a1)\
				  );

//c -= a	
#define SUBACC_128(a0, a1, c0, c1)\
	asm ("SUBS %0, %0, %2;" \
		 "SBC %1, %1, %3;"\
		: "+r" (c0), "+r" (c1)\
		: "r" (a0), "r" (a1)\
		);

ec_split_scalar ec_scalar_decomp(uint64x2x2_t k) {
	ec_split_scalar result;
	uint64_t tmp0, tmp1, tmp2, tmp3, sign;
	uint64_t zero = 0;
	
	// Step 1: b1 = k / 2^127 (where / is integer division) (Original comments say it is -k / 2^127, mystery why that is)
	uint64x2_t b1;
	b1[0] = (k.val[1][0] << 1) | (k.val[0][1] >> 63);
	b1[1] = (k.val[1][1] << 1) | (k.val[1][0] >> 63);
	
	//Step 2: b2 = k*trace / 2^254
	// Mult trace by each word in k.
	//All of this can probably be optimized because we can't reach the top word without ext atm.
	//Also there's no reliable 64 bit add instructions with carry :/
	uint64x2_t tk0 = mult_u64((uint64_t) TRACE, k.val[0][0]);
	uint64x2_t tk1 = mult_u64((uint64_t) TRACE, k.val[0][1]);
	uint64x2_t tk2 = mult_u64((uint64_t) TRACE, k.val[1][0]);
	uint64x2_t tk3 = mult_u64((uint64_t) TRACE, k.val[1][1]);
	
	//Add previous [1] with next [0]
	//We are only interested in the last two words of the result, rest will be consumed by division by 2^254.
	asm ("ADDS %[tmp0], %[tk01], %[tk10];"
		 "ADCS %[tmp0], %[tk11], %[tk20];"
		 "ADCS %[tmp0], %[tk21], %[tk30];"
		 "ADC %[tmp1], %[tk31], %[zero];"
		: [tmp0] "+r" (tmp0), [tmp1] "+r" (tmp1)
		: [tk01] "r" (tk0[1]), [tk10] "r" (tk1[0]), [tk11] "r" (tk1[1]), [tk20] "r" (tk2[0]), [tk21] "r" (tk2[1]), [tk30] "r" (tk3[0]), [tk31] "r" (tk3[1]), [zero] "r" (zero)
		);
	
	//Divide k*trace by 2^254
	uint64_t b2 = (tmp0 >> 62) | (tmp1 << 2);
	b2 |= 1; //Mystery or!
	
	//Step 3: b1*t
	uint64x2_t b1_times_t = mult_u64((uint64_t) TRACE, b1[0]);
	uint64x2_t b11_times_t = mult_u64((uint64_t) TRACE, b1[1]);
	b1_times_t[1] += b11_times_t[0]; //He ignores last word of result???
	
	//Step 4: b2*t
	uint64x2_t b2_times_t = mult_u64((uint64_t) TRACE, b2);
	
	//k1 computation
	
	//Step 5: {tmp0, tmp1} = b1*q (q = 2^127)
	tmp0 = 0;
	tmp1 = b1[0] << 63;
	
	//Step 6: {tmp0, tmp1} = b1*q + b2*t
	ADDACC_128(b2_times_t[0], b2_times_t[1], tmp0, tmp1);
		
	//Step 7: {tmp0, tmp1} = b1*q + b2*t - b1
	//Hope operands are in the correct order
	SUBACC_128(b1[0], b1[1], tmp0, tmp1);
	
	//Step 8: Determine sign of k-tmp (0 for positive), just by checking second word?? Maybe if res1=k.val[0][1] you can show that it will not be negative.
	sign = tmp1 < k.val[0][1];
	
	//Step 9: {tmp0, tmp1} = (b1*q + b2*t - b1) - k
	//Opposite order of the paper, why tho?? And we just subtract with the first two words of k like it is nothing.
	SUBACC_128(k.val[0][0], k.val[0][1], tmp0, tmp1);
	
	//Step 10: Now take two's complement if needed and then flip sign, no matter order of ops may need to do two's complement transform.
	tmp0 ^= (zero - sign);
	tmp1 ^= (zero - sign);
	ADDACC_128(sign, zero, tmp0, tmp1);
	result.k1[0] = tmp0;
	result.k1[1] = tmp1;
	result.k1_sign = sign ^ 0x1;
	
	//k2 computation
	
	//Step 11: {tmp2, tmp3} = b1*t + b2
	tmp2 = b1_times_t[0];
	tmp3 = b1_times_t[1];
	ADDACC_128(b2, zero, tmp2, tmp3);
	
	//Step 12: {tmp0, tmp1} = b2*q (q = 2^127)
	tmp0 = 0;
	tmp1 = b2 << 63;
	
	//Step 13: k2 sign (0 for positive)
	sign = tmp1 < tmp3;
	
	//Step 14: {tmp0, tmp1} = b2*q - (b1*t + b2)
	SUBACC_128(tmp2, tmp3, tmp0, tmp1);
	
	//Step 15: Now take two's complement if needed and then flip sign again.
	tmp0 ^= (zero - sign);
	tmp1 ^= (zero - sign);
	ADDACC_128(sign, zero, tmp0, tmp1);
	result.k2[0] = tmp0;
	result.k2[1] = tmp1;
	result.k2_sign = sign ^ 0x1;
	
	return result;
}
