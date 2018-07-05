
test: test.C SignedDistance.C 2DSDF.C
	g++-6 -O0 -g -Wextra  -Wall -pedantic -std=c++11  $^ -o $@

.PHONY: clean
clean:
	rm -r test
