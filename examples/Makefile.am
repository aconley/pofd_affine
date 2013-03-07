if DOEXAMPLES

OTHERDIR=../src/

bin_PROGRAMS = multi_gauss_example rosenbrock_example \
	single_gauss_example polyfit_example

single_gauss_example_SOURCES = single_gauss_example.cc
single_gauss_example_LDADD = $(OTHERDIR)/libaffine.la

multi_gauss_example_SOURCES = multi_gauss_example.cc
multi_gauss_example_LDADD = $(OTHERDIR)/libaffine.la

polyfit_example_SOURCES = polyfit_example.cc
polyfit_example_LDADD = $(OTHERDIR)/libaffine.la

rosenbrock_example_SOURCES = rosenbrock_example.cc
rosenbrock_example_LDADD = $(OTHERDIR)/libaffine.la

endif

EXTRA_DIST= exampledata/polyexample.txt
