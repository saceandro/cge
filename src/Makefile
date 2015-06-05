.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)
%: %.cc
	$(CXX) -Wall -g -o $@ $<
