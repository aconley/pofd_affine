if DOTEST

OTHERDIR=../src/

check_PROGRAMS = chain_test ensemble_test queue_test

chain_test_SOURCES = chain_test.cc
chain_test_LDADD = $(OTHERDIR)/libaffine.la

ensemble_test_SOURCES = ensemble_test.cc
ensemble_test_LDADD = $(OTHERDIR)/libaffine.la

queue_test_SOURCES = queue_test.cc
queue_test_LDADD = $(OTHERDIR)/libaffine.la

endif

EXTRA_DIST=testdata/testvec.txt testdata/testinvcov.txt
