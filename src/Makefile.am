noinst_LTLIBRARIES = libaffine.la libpd.la

libaffine_la_SOURCES = affineChainSet.cc affineEnsemble.cc affineExcept.cc \
	paramSet.cc proposedStep.cc ran.cc hdf5utils.cc hashbar.cc

libpd_la_SOURCES = utility.cc numberCounts.cc numberCountsKnots.cc \
	numberCountsKnotsSpline.cc beam.cc numberCountsDouble.cc \
	numberCountsDoubleLogNormal.cc doublebeam.cc PD.cc PDDouble.cc \
	PDFactory.cc PDFactoryDouble.cc fitsData.cc fitsDataDouble.cc \
	calcLike.cc calcLikeDouble.cc specFile.cc specFileDouble.cc
libpd_la_LIBADD = libaffine.la

bin_PROGRAMS = pofd_affine_getdNdS pofd_affine_getR pofd_affine_getPD \
	pofd_affine_getLike pofd_affine_mcmc

pofd_affine_getdNdS_SOURCES = pofd_affine_getdNdS.cc
pofd_affine_getdNdS_LDADD = libpd.la

pofd_affine_getR_SOURCES = pofd_affine_getR.cc
pofd_affine_getR_LDADD = libpd.la

pofd_affine_getPD_SOURCES = pofd_affine_getPD.cc 
pofd_affine_getPD_LDADD = libpd.la

pofd_affine_getLike_SOURCES = pofd_affine_getLike.cc 
pofd_affine_getLike_LDADD = libpd.la

pofd_affine_mcmc_SOURCES = pofdMCMC.cc pofdMCMCDouble.cc \
	pofd_affine_mcmc.cc 
pofd_affine_mcmc_LDADD = libpd.la

