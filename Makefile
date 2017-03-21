# Ludwig Schmidt (ludwigschmidt2@gmail.com) 2017
#
# This makefile is based on http://make.paulandlesley.org/autodep.html .

CXX = clang++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -fPIC
#CXXFLAGS = -std=c++11 -Wall -Wextra -g -fPIC
GTESTDIR = external/googletest/googletest

SRCDIR = src
OBJDIR = obj

.PHONY: clean

clean:
	rm -rf $(OBJDIR)
	rm -f pcst_fast_wrap.cxx
	rm -f _pcst_fast.so
	rm -f pcst_fast.py
	rm -f pcst_fast.pyc
	rm -f pcst_fast_test

run_tests: run_pcst_fast_test

mexfiles: cluster_grid_pcst_mexfile cluster_grid_pcst_binsearch_mexfile

# gtest
$(OBJDIR)/gtest-all.o: $(GTESTDIR)/src/gtest-all.cc
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -I $(GTESTDIR) -c -o $@ $<

$(OBJDIR)/gtest_main.o: $(GTESTDIR)/src/gtest_main.cc
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -c -o $@ $<


PCST_FAST_TEST_OBJS = gtest-all.o gtest_main.o
pcst_fast_test: $(PCST_FAST_TEST_OBJS:%=$(OBJDIR)/%) $(SRCDIR)/pcst_fast_test.cc $(SRCDIR)/pcst_fast.cc
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -c -o $(OBJDIR)/pcst_fast.o $(SRCDIR)/pcst_fast.cc
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -c -o $(OBJDIR)/pcst_fast_test.o $(SRCDIR)/pcst_fast_test.cc
	$(CXX) $(CXXFLAGS) -o $@ $(PCST_FAST_TEST_OBJS:%=$(OBJDIR)/%) $(OBJDIR)/pcst_fast_test.o $(OBJDIR)/pcst_fast.o -pthread

run_pcst_fast_test: pcst_fast_test
	./pcst_fast_test


PCST_FAST_PY_SRC = pcst_fast_pybind.cc
PCST_FAST_PY_SRC_DEPS = $(PCST_FAST_PY_SRC) pcst_fast.h pcst_fast.cc
pcst_fast_py: $(PCST_FAST_PY_SRC_DEPS:%=$(SRCDIR)/%)
	$(CXX) $(CXXFLAGS) -shared -I $(SRCDIR) -I external/pybind11/include `python-config --cflags --ldflags` $(SRCDIR)/pcst_fast_pybind.cc $(SRCDIR)/pcst_fast.cc -o pcst_fast.so

run_pcst_fast_py_test: pcst_fast_py
	python -m pytest src/test_pcst_fast.py
