#include <stdio.h>
#include <arm_neon.h>
#include <stdlib.h>
#include <time.h>

void print_poly(poly64x2_t p) {
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

uint64_t rand_uint64() {
  uint64_t r = 0;
  for (int i=0; i<64; i++) {
    r = r*2 + rand()%2;
  }
  return r;
}

poly64x2_t rand_poly() {
	time_t t;
   
	// Intializes random number generator
	srand((unsigned) time(&t));
   
	// 2^63-1 = 01111111...
	long c = 9223372036854775807;
	
	poly64_t p1 = rand_uint64();
	poly64_t p2 = rand_uint64() & c;
	
	poly64x2_t p = {p1, p2};
	
	return p;
}

typedef struct {
	poly64x2_t p0, p1;
} pmullres;

pmullres pmull32(poly64x2_t a, poly64x2_t b) {
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

static inline uint64_t read_pmccntr() {
		uint64_t val;
		asm volatile("mrs %0, pmccntr_el0" : "=r"(val));
		return val;
}

int main() {
	//printf("Cycle counter: %u ", read_pmccntr());
	poly64x2_t a = {11, 0};
	poly64x2_t b = {6, 0};
	pmullres res = pmull32(a, b);
	print_poly(a);
	print_poly(b);
	print_poly(res.p0);
	/*
	print_poly(rand_poly());
	poly64x2_t r = {5, 0};
	print_poly(r);
	printf("\n");
	poly64x2_t castedR = (poly64x2_t) vreinterpretq_u64_p128(5);
	print_poly(castedR);*/
	//poly128_t res = vmull_p64(castedR, castedR);
	return 0;
}
