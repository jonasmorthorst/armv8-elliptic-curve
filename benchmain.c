#include <unistd.h>

#include "benchmark/benchmark_tool.h"
#include "benchmark/benchmark_basefield.h"
#include "benchmark/benchmark_extensionfield.h"
#include "common/setup.h"

int main() {
	init_components();
	
	nice(-30);
	benchmark_bf_all();
	benchmark_ef_all();
	
	dispose_components();
	return 0;
}
