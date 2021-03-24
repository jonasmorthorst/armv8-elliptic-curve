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
				
				if (i == 1 && j == 0) {
					printf("1");
				} else {
					printf("z^%d", (1-i)*64 +j);
				}
			}
			c /= 2;
		}
	}
	
	if (wasFirst) {
		printf("0");
	}
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
	
	poly64_t p1 = rand_uint64() & c;
	poly64_t p2 = rand_uint64();
	
	poly64x2_t p = {p1, p2};
	
	return p;
}

//poly128_t pmullTheirs64(poly128_t a, poly128_t b) {
//	poly128_t z = 0;
//	vmull_p64
//}



int main() {
	//print_poly(rand_poly());
	poly128_t r = 5;
	//poly64x2_t castedR = (poly64x2_t) vreinterpretq_u64_p128(r);
	//poly128_t res = vmull_p64(castedR, castedR);
	return 0;
}
