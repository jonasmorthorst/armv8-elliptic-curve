#To compile with gcc: "make runtests gcc=1"
GCCFLAGS = -march=armv8.4-a+crypto

runsandbox:
ifdef gcc
	gcc -o run_sandbox_gcc $(GCCFLAGS) sandboxmain.c common/*.c
	./run_sandbox_gcc
else
	clang -o run_sandbox_clang sandboxmain.c common/*.c
	./run_sandbox_clang
endif
runtests:
ifdef gcc
	gcc -o run_tests_gcc $(GCCFLAGS) testmain.c common/*.c tests/*.c benchmark/benchmark_tool.c
	./run_tests_gcc
else
	clang -o run_tests_clang testmain.c common/*.c tests/*.c benchmark/benchmark_tool.c
	./run_tests_clang
endif
runbench:
ifdef gcc
	gcc -o run_bench_gcc $(GCCFLAGS) benchmain.c common/*.c benchmark/*.c
	./run_bench_gcc
else
	clang -o run_bench_clang benchmain.c common/*.c benchmark/*.c
	./run_bench_clang
endif
