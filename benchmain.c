#include <unistd.h>

#include "benchmark/benchmark_tool.h"
#include "benchmark/benchmark_basefield.h"

int main() {
	nice(-30);
	benchmark_bf_all();
	return 0;
}
