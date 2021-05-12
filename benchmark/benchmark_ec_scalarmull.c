#include <stdio.h>

#include "../common/ec_scalarmull.h"
#include "../common/ec.h"
#include "../common/utils.h"
#include "benchmark_ec_scalarmull.h"
#include "benchmark_tool.h"


void benchmark_ec_scalarmull_single() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		ec_point_laffine P = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_lproj R = ec_scalarmull_single(P, k);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(R, sum);
	}
	ec_print_hex(sum);
	printf("BENCHMARK benchmark_ec_scalarmull_single\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_scalarmull_single_endo() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		ec_point_laffine P = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_laffine R = ec_scalarmull_single_endo(P, k);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add_mixed(R, sum);
	}
	ec_print_hex(sum);
	printf("BENCHMARK benchmark_ec_scalarmull_single_endo\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_scalarmull_single_endo_w5_randaccess() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		ec_point_laffine P = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_laffine R = ec_scalarmull_single_endo_w5_randaccess(P, k);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add_mixed(R, sum);
	}
	ec_print_hex(sum);
	printf("BENCHMARK benchmark_ec_scalarmull_single_endo_w5_randaccess\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_scalarmull_single_endo_w6_randaccess() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		ec_point_laffine P = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_laffine R = ec_scalarmull_single_endo_w6_randaccess(P, k);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add_mixed(R, sum);
	}
	ec_print_hex(sum);
	printf("BENCHMARK benchmark_ec_scalarmull_single_endo_w6_randaccess\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_scalarmull_double() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		ec_point_lproj P = ec_rand_point_lproj();
		uint64x2x2_t l = ec_rand_scalar();
		ec_point_lproj Q = ec_rand_point_lproj();
		uint64_t start = read_pmccntr();
		ec_point_lproj R = ec_scalarmull_double(P, k, Q, l);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add(sum, R);
	}
	ec_print_hex(sum);
	printf("BENCHMARK benchmark_ec_scalarmull_double\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_scalarmull_all() {
	benchmark_ec_scalarmull_single();
	benchmark_ec_scalarmull_single_endo();
	benchmark_ec_scalarmull_single_endo_w5_randaccess();
	benchmark_ec_scalarmull_single_endo_w6_randaccess();
	benchmark_ec_scalarmull_double();
}
