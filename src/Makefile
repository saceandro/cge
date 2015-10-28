CXX = g++

%: %.cc
	$(CXX) -I/usr/local/Cellar/boost/1.58.0/include/ -Wall -g -o $@ $< -lm

.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only em.cc

