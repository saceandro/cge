CXX = g++

.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only setting.hh calc_read_prob.hh calc_read_prob.cc em.hh em.cc
%: %.cc
	$(CXX) -Wall -g -o $@ $<
