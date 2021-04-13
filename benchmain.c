#include <unistd.h>

#include "benchmark/benchmark_tool.h"
#include "benchmark/benchmark_basefield.h"
#include "benchmark/benchmark_extensionfield.h"

int main() {
	nice(-30);
	benchmark_bf_all();
	benchmark_ef_all();
	return 0;
}
