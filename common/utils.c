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
