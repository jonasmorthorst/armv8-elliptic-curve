#include <stdio.h>
#include "benchmark_tool.h"
#include "../common/utils.h"

inline uint64_t read_pmccntr() {
		uint64_t val;
		asm volatile("mrs %0, pmccntr_el0" : "=r"(val));
		return val;
}

void insert_sorted(uint64_t time, uint64_t *sorted_times, uint64_t curr_len) {
	for (int i = 0; i < curr_len; i++) {
		if (sorted_times[i] > time) {
			uint64_t tmp = time;
			time = sorted_times[i];
			sorted_times[i] = tmp;
		}
	}
	sorted_times[curr_len] = time;
}
