#include <stdio.h>
#include <stdlib.h>

#include "../common/ec_scalarmull.h"
#include "../common/ec.h"
#include "../common/utils.h"
#include "benchmark_ec_scalarmull.h"
#include "benchmark_tool.h"


void benchmark_ec_linaer_pass_w5() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_laffine sum = ec_rand_point_laffine();

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();

		ec_point_laffine table[8];
		precompute_w5(P, table);

		uint64_t index1 = rand() % (7 + 1);
		uint64_t index2 = rand() % (7 + 1);
		ec_point_laffine P1 = ec_rand_point_laffine();
		ec_point_laffine P2 = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		lin_pass_w5(&P1, &P2, &table, index1, index2);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_lproj_to_laffine(ec_add_laffine_unchecked(P1, sum));
		sum = ec_lproj_to_laffine(ec_add_laffine_unchecked(P2, sum));
	}
	ec_print_hex_laffine(sum);
	printf("BENCHMARK benchmark_ec_linaer_pass_w5\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_csel_asm() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_laffine sum = ec_rand_point_laffine();

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();
		ec_point_laffine P_neg = ec_neg_laffine(P);
		uint64_t new_ptr, con = 1;
		uint64_t rand_num = rand() % 2;

		uint64_t start = read_pmccntr();
		csel_asm(rand_num, con, &P, &P_neg);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_lproj_to_laffine(ec_add_laffine_unchecked(P, sum));
	}
	ec_print_hex_laffine(sum);
	printf("BENCHMARK benchmark_csel_asm\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_csel_inline_asm() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_laffine sum = ec_rand_point_laffine();

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();
		ec_point_laffine P_neg = ec_neg_laffine(P);
		uint64_t new_ptr, con = 1;
		uint64_t rand_num = rand() % 2;

		uint64_t start = read_pmccntr();
		CSEL(rand_num, con, P, P_neg, new_ptr, typeof(ec_point_laffine));
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_lproj_to_laffine(ec_add_laffine_unchecked(P, sum));
	}
	ec_print_hex_laffine(sum);
	printf("BENCHMARK benchmark_csel_inline_asm\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_cmov() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_laffine sum = ec_rand_point_laffine();

	for(int i = 0; i < num_runs; i++) {
		ec_point_laffine P = ec_rand_point_laffine();
		ec_point_laffine P_neg = ec_neg_laffine(P);
		uint64x1_t old_ptr, new_ptr, tmp, cond;
		cond[0]=1;

		uint64_t rand_num = rand() % 2;

		uint64_t start = read_pmccntr();
		CMOV(tmp, rand_num, cond, P, P_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_lproj_to_laffine(ec_add_laffine_unchecked(P, sum));
	}
	ec_print_hex_laffine(sum);
	printf("BENCHMARK benchmark_cmov\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

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

void benchmark_ec_scalarmull_single_endo_w3_randaccess() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		ec_point_laffine P = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_laffine R = ec_scalarmull_single_endo_w3_randaccess(P, k);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add_mixed(R, sum);
	}
	ec_print_hex(sum);
	printf("BENCHMARK benchmark_ec_scalarmull_single_endo_w3_randaccess\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_ec_scalarmull_single_endo_w4_randaccess() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs];
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		ec_point_laffine P = ec_rand_point_laffine();
		uint64_t start = read_pmccntr();
		ec_point_laffine R = ec_scalarmull_single_endo_w4_randaccess(P, k);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
		sum = ec_add_mixed(R, sum);
	}
	ec_print_hex(sum);
	printf("BENCHMARK benchmark_ec_scalarmull_single_endo_w4_randaccess\n");
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

void benchmark_ec_scalarmull_single_endo_w5_randaccess_time() {
	uint64_t num_runs = 2000;
	unsigned long long start, end, nsec = 0;
	ec_point_lproj sum = ec_rand_point_lproj();

	for(int i = 0; i < num_runs; i++) {
		uint64x2x2_t k = ec_rand_scalar();
		ec_point_laffine P = ec_rand_point_laffine();
		int64_t start = cpu_nseconds();
		ec_point_laffine R = ec_scalarmull_single_endo_w5_randaccess(P, k);
		int64_t end = cpu_nseconds();
		nsec = nsec+(end-start);
		sum = ec_add_mixed(R, sum);
	}
	ec_print_hex(sum);
	printf("BENCHMARK benchmark_ec_scalarmull_single_endo_w5_randaccess_time\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average time: %8lld\n\n", nsec/num_runs);
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
	benchmark_ec_linaer_pass_w5();
	benchmark_csel_asm();
	benchmark_csel_inline_asm();
	benchmark_cmov();
	benchmark_ec_scalarmull_single();
	benchmark_ec_scalarmull_single_endo();
	benchmark_ec_scalarmull_single_endo_w3_randaccess();
	benchmark_ec_scalarmull_single_endo_w4_randaccess();
	benchmark_ec_scalarmull_single_endo_w5_randaccess();
	benchmark_ec_scalarmull_single_endo_w5_randaccess_time();
	benchmark_ec_scalarmull_single_endo_w6_randaccess();
	benchmark_ec_scalarmull_double();
}
