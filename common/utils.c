#include <stdlib.h>
#include <time.h>

#include "utils.h"

void utils_init() {
	// Intializes random number generator
	time_t t;
	srand((unsigned) time(&t));
}

double average(uint64_t nums[], uint64_t len) {
	double sum = 0.0;
	for (int i = 0; i < len; i++) {
		sum += nums[i];
	}
	return sum / len;
}

/* Returns:
 *  0 if diff within [-errmargin, errmargin],
 *  1 if a is greater
 * -1 if b is greater
 */
uint64_t compare_doubles(double a, double b, double errmargin) {
	double diff = a - b;
	if (diff > errmargin) {
		return 1;
	}
	if (diff < -errmargin) {
		return -1;
	}
	return 0;
}

uint64_t equal_poly64x2(poly64x2_t a, poly64x2_t b) {
	return (a[0] == b[0]) && (a[1] == b[1]);
}

uint64_t equal_poly64x2x2(poly64x2x2_t a, poly64x2x2_t b) {
	return (a.val[0][0] == b.val[0][0]) &&
		   (a.val[0][1] == b.val[0][1]) &&
		   (a.val[1][0] == b.val[1][0]) &&
		   (a.val[1][1] == b.val[1][1]);
}

uint64_t equal_ef_elem(ef_elem a, ef_elem b) {
	return (a.val[0][0] == b.val[0][0]) &&
		   (a.val[0][1] == b.val[0][1]) &&
		   (a.val[1][0] == b.val[1][0]) &&
		   (a.val[1][1] == b.val[1][1]);
}

poly64x2x2_t concat_bf_poly(poly64x2_t p0, poly64x2_t p1) {
	poly64x2x2_t p;
	p.val[0] = p0;
	p.val[1] = p1;
	return p;
}

double median(uint64_t sorted_nums[], uint64_t len) {
	if (len % 2 == 0) {
		return (sorted_nums[len / 2 - 1] + sorted_nums[len / 2]) / 2.0;
	}
	return sorted_nums[len / 2];
}

uint64_t rand_uint64() {
	uint64_t r = 0;
	for (int i=0; i<64; i++) {
		r = r*2 + rand()%2;
	}
	return r;
}

ec_point_laffine lproj_to_laffine(ec_point_lproj P) {
	ec_point_laffine L;
	L.x = ef_mull(P.x, ef_inv(P.z));
	L.l = ef_mull(P.l, ef_inv(P.z));

	return L;
}
