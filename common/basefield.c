#include <stdio.h>

#include "basefield.h"
#include "utils.h"

void bf_print_elem_expr(poly64x2_t p) {
	poly64_t c;
	int wasFirst = 1;
	for (int i=1; i>=0; i--) {
		// 2^63 = the value of the leftmost bit in a word
		c = 9223372036854775808U;
		for (int j = 63; j>=0; j--) {
			poly64_t polybitcopy = c & p[i];
			if (polybitcopy == c) {
				if(!wasFirst) {
					printf("+");
				}
				wasFirst = 0;
				
				if (i == 0 && j == 0) {
					printf("1");
				} else {
					printf("z^%d", (i)*64 +j);
				}
			}
			c /= 2;
		}
	}
	
	if (wasFirst) {
		printf("0");
	}
	printf("\n");
}

poly64x2_t bf_rand_elem() { 
	// 2^63-1 = 01111111...
	long c = 9223372036854775807;
	
	poly64_t p1 = rand_uint64();
	poly64_t p2 = rand_uint64() & c;
	
	poly64x2_t p = {p1, p2};
	
	return p;
}

pmullres bf_pmull32(poly64x2_t a, poly64x2_t b) {
	poly64x2_t t = {0,0};
	pmullres r;
	r.p0 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[0], b[0]));
	r.p1 = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(a[1], b[1]));
	t[1] = (poly64_t) veor_u64((uint64x1_t) b[0], (uint64x1_t) b[1]);
	t[0] = (poly64_t) veor_u64((uint64x1_t) a[0], (uint64x1_t) a[1]);
	t = (poly64x2_t) vreinterpretq_u64_p128(vmull_p64(t[1], t[0]));
	t = (poly64x2_t) veorq_u64((uint64x2_t) t, (uint64x2_t) r.p0);
	t = (poly64x2_t) veorq_u64((uint64x2_t) t, (uint64x2_t) r.p1);
	r.p0[1] = (poly64_t) veor_u64((uint64x1_t) r.p0[1], (uint64x1_t) t[0]);
	r.p1[0] = (poly64_t) veor_u64((uint64x1_t) r.p1[0], (uint64x1_t) t[1]);
	
	return r;
}
