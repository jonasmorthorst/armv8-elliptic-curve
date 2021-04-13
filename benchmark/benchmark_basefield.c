#include <stdio.h>

#include "../common/basefield.h"
#include "../common/utils.h"
#include "benchmark_tool.h"

void benchmark_bf_add() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		uint64_t start = read_pmccntr();
		bf_add(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_add\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_pmull32() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		uint64_t start = read_pmccntr();
		bf_pmull32(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_pmull32\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_pmull64() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		uint64_t start = read_pmccntr();
		bf_pmull64(a, b);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_pmull64\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_psquare() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		bf_psquare(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_psquare\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_red() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		poly64x2_t b = bf_rand_elem();
		bf_polyx2 c = bf_pmull(a,b); //To get more avrg input
		uint64_t start = read_pmccntr();
		bf_red(c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_red\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_red_psquare() {
	uint64_t num_runs = 20000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		bf_polyx2 c = bf_psquare(a); //To get more avrg input
		uint64_t start = read_pmccntr();
		bf_red_psquare(c);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_red_psquare\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_fermat_inv() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		bf_fermat_inv(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_fermat_inv\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_addchain_inv() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		bf_addchain_inv(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_addchain_inv\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_addchain_lookup_inv() {
	uint64_t num_runs = 2000;
	uint64_t times[num_runs]; 
	
	for(int i = 0; i < num_runs; i++) {
		poly64x2_t a = bf_rand_elem();
		uint64_t start = read_pmccntr();
		bf_addchain_lookup_inv(a);
		uint64_t end = read_pmccntr();
		insert_sorted(end-start, times, i);
	}
	printf("BENCHMARK bf_addchain_lookup_inv\n");
	printf("Number of iterations: %lu\n", num_runs);
	printf("Average: %lf\n", average(times, num_runs));
	printf("Median: %lf\n\n", median(times, num_runs));
}

void benchmark_bf_all() {
	bf_init();
	
	benchmark_bf_add();
	benchmark_bf_pmull32();
	benchmark_bf_pmull64();
	benchmark_bf_psquare();
	benchmark_bf_red();
	benchmark_bf_red_psquare();
	benchmark_bf_fermat_inv();
	benchmark_bf_addchain_inv();
	benchmark_bf_addchain_lookup_inv();
	bf_post();
}

