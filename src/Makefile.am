noinst_LTLIBRARIES = libaffine.la

libaffine_la_SOURCES = affineChainSet.cc affineEnsemble.cc affineExcept.cc \
	paramSet.cc proposedStep.cc ran.cc

bin_PROGRAMS = pofd_affine_getdNdS pofd_affine_getR pofd_affine_getPD \
	pofd_affine_getLike pofd_affine_mcmc

pofd_affine_getdNdS_SOURCES = numberCounts.cc numberCountsKnots.cc \
	beam.cc utility.cc numberCountsKnotsSpline.cc \
	numberCountsDouble.cc numberCountsDoubleLogNormal.cc doublebeam.cc \
	pofd_affine_getdNdS.cc
pofd_affine_getdNdS_LDADD = libaffine.la

pofd_affine_getR_SOURCES =numberCounts.cc numberCountsKnots.cc beam.cc \
	utility.cc numberCountsKnotsSpline.cc numberCountsDouble.cc \
	numberCountsDoubleLogNormal.cc doublebeam.cc \
	pofd_affine_getR.cc
pofd_affine_getR_LDADD = libaffine.la

pofd_affine_getPD_SOURCES = numberCounts.cc numberCountsKnots.cc 
	beam.cc utility.cc numberCountsKnotsSpline.cc PD.cc PDFactory.cc \
	fitsData.cc numberCountsDouble.cc numberCountsDoubleLogNormal.cc \
	doublebeam.cc ran.cc PDDouble.cc PDFactoryDouble.cc fitsDataDouble.cc \
	pofd_affine_getPD.cc 
pofd_affine_getPD_LDADD = libaffine.la

pofd_affine_getLike_SOURCES = numberCounts.cc numberCountsKnots.cc beam.cc \
	utility.cc numberCountsKnotsSpline.cc PD.cc PDFactory.cc fitsData.cc \
	calcLike.cc numberCountsDouble.cc numberCountsDoubleLogNormal.cc \
	doublebeam.cc PDDouble.cc PDFactoryDouble.cc fitsDataDouble.cc \
	calcLikeDouble.cc specFile.cc specFileDouble.cc \
	pofd_affine_getLike.cc 
pofd_affine_getLike_LDADD = libaffine.la

pofd_affine_mcmc_SOURCES = numberCounts.cc numberCountsKnots.cc beam.cc \
	utility.cc numberCountsKnotsSpline.cc PD.cc PDFactory.cc fitsData.cc \
	calcLike.cc numberCountsDouble.cc numberCountsDoubleLogNormal.cc \
	doublebeam.cc PDDouble.cc PDFactoryDouble.cc fitsDataDouble.cc \
	calcLikeDouble.cc specFile.cc specFileDouble.cc pofdMCMC.cc \
	pofdMCMCDouble.cc \
	pofd_affine_getLike.cc 
pofd_affine_mcmc_LDADD = libaffine.la
