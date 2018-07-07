test: test.cpp
	g++-6 -O0 -g -Wextra -o test -Wall -pedantic -std=c++11  $^ -o $@

.PHONY: clean
clean:
	rm test
