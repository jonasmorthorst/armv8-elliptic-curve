runtests:
	gcc -o runtests -march=armv8-a+crypto testmain.c common/*.c tests/*.c
	./runtests
runbench:
	gcc -o runbench -march=armv8-a+crypto benchmain.c common/*.c benchmark/*.c
	./runbench
