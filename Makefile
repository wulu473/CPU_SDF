main: main.cpp
	@echo "Compiling..."
	g++-6 -O0 -g -Wextra  -Wall -pedantic -std=c++11 \
		$^ -o $@
	@echo "done"

main.cpp: SignedDistance.hpp 2DSDF.hpp

.PHONY: all
all: main test

test: test.cpp
	@echo "Compiling..."
	g++-6 -O0 -g -Wextra  -Wall -pedantic -std=c++11 \
	       	-I/lsc/opt/modules/gcc-6.4.0/boost-1.61.0/include/ \
		$^ -o $@ \
		-lboost_unit_test_framework -L/lsc/opt/modules/gcc-6.4.0/boost-1.61.0/lib  
	@echo "done"

test.cpp: SignedDistance.hpp 2DSDF.hpp

.PHONY: clean
clean:
	rm -rf test main

.PHONY: memtest
memtest:
	valgrind --leak-check=full --track-origins=yes ./main

