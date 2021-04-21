#include <stdio.h>

#include "../common/ec.h"
#include "../common/utils.h"
#include "benchmark_ec.h"
#include "benchmark_tool.h"


void benchmark_ec_add() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];

	for(int i = 0; i < num_runs; i++) {
    ec_point_lproj a = ec_rand_point_lproj();
    ec_point_lproj b = ec_rand_point_lproj();

		uint64_t start = read_pmccntr();
		ec_add(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK ec_add\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_double() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];

	for(int i = 0; i < num_runs; i++) {
    ec_point_lproj a = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_double(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK ec_double\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_all() {
	benchmark_ec_add();
	benchmark_ec_double();
}
