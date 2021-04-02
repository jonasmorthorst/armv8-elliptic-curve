GCCFLAGS = -march=armv8.4-a+crypto

runsandbox:
	gcc -o run_sandbox $(GCCFLAGS) sandboxmain.c common/*.c
	./run_sandbox

runtests:
	gcc -o run_tests $(GCCFLAGS) testmain.c common/*.c tests/*.c benchmark/benchmark_tool.c
	./run_tests
runbench:
	gcc -o run_bench $(GCCFLAGS) benchmain.c common/*.c benchmark/*.c
	./run_bench
