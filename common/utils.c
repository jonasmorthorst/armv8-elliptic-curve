#include <stdlib.h>
#include <time.h>

#include "utils.h"

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

uint64_t equal_bf_polyx2(bf_polyx2 a, bf_polyx2 b) {
	return (a.p0[0] == b.p0[0]) &&
		   (a.p0[1] == b.p0[1]) &&
		   (a.p1[0] == b.p1[0]) &&
		   (a.p1[1] == b.p1[1]);
}

uint64_t equal_ef_elem(ef_elem a, ef_elem b) {
	return (a.p0[0] == b.p0[0]) &&
		   (a.p0[1] == b.p0[1]) &&
		   (a.p1[0] == b.p1[0]) &&
		   (a.p1[1] == b.p1[1]);
}

bf_polyx2 concat_bf_poly(poly64x2_t p0, poly64x2_t p1) {
	bf_polyx2 p;
	p.p0 = p0;
	p.p1 = p1;
	return p;
}

double median(uint64_t sorted_nums[], uint64_t len) {
	if (len % 2 == 0) {
		return (sorted_nums[len / 2 - 1] + sorted_nums[len / 2]) / 2.0;
	}
	return sorted_nums[len / 2];
}

uint64_t rand_uint64() {
	time_t t;
   
	// Intializes random number generator
	srand((unsigned) time(&t));
	
	uint64_t r = 0;
	for (int i=0; i<64; i++) {
		r = r*2 + rand()%2;
	}
	return r;
}
