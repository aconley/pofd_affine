if DOTEST

OTHERDIR=../src/

check_PROGRAMS = chain_test ensemble_test mpi_copy_test initFile_test \
	specFile_test beam_test model_test queue_test data_test pofd_test

chain_test_SOURCES = chain_test.cc
chain_test_LDADD = $(OTHERDIR)/libaffine.la

ensemble_test_SOURCES = ensemble_test.cc
ensemble_test_LDADD = $(OTHERDIR)/libaffine.la

queue_test_SOURCES = queue_test.cc
queue_test_LDADD = $(OTHERDIR)/libaffine.la

mpi_copy_test_SOURCES = mpi_copy_test.cc
mpi_copy_test_LDADD = $(OTHERDIR)/libpd.la

initFile_test_SOURCES = initFile_test.cc
initFile_test_LDADD = $(OTHERDIR)/libpd.la

specFile_test_SOURCES = specFile_test.cc
specFile_test_LDADD = $(OTHERDIR)/libpd.la

beam_test_SOURCES = beam_test.cc
beam_test_LDADD = $(OTHERDIR)/libpd.la

model_test_SOURCES = model_test.cc
model_test_LDADD = $(OTHERDIR)/libpd.la

pofd_test_SOURCES = pofd_test.cc
pofd_test_LDADD = $(OTHERDIR)/libpd.la

data_test_SOURCES = data_test.cc
data_test_LDADD = $(OTHERDIR)/libpd.la

endif

EXTRA_DIST=testdata/testvec.txt testdata/testinvcov.txt \
	testdata/fiducial_1Dmodel.txt testdata/initfile_test1D.txt \
	testdata/specfile_test1D.txt testdata/specfile_test2D.txt \
	testdata/band1_beam.fits testdata/band2_beam.fits \
	testdata/testmodel2D_band1.fits testdata/testmodel2D_band2.fits \
	testdata/band12_beam_f100.fits
