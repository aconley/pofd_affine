if BUILDTELIB
noinst_LTLIBRARIES = libaffine.la

libaffine_la_SOURCES = affineChainSet.cc affineEnsemble.cc affineExcept.cc \
	paramSet.cc proposedStep.cc ran.cc
endif

bin_PROGRAMS = pofd_affine_getdNdS pofd_affine_getR pofd_affine_getPD \
	pofd_affine_getLike pofd_affine_mcmc

pofd_affine_getdNdS_SOURCES =global_settings.h numberCounts.h \
	numberCountsKnots.h numberCountsKnotsSpline.h paramSet.h \
	affineExcept.h beam.h utility.h numberCountsDouble.h \
	numberCountsDoubleLogNormal.h doublebeam.h ran.h \
	numberCounts.cc numberCountsKnots.cc beam.cc utility.cc \
	numberCountsKnotsSpline.cc paramSet.cc affineExcept.cc ran.cc \
	numberCountsDouble.cc numberCountsDoubleLogNormal.cc doublebeam.cc \
	pofd_affine_getdNdS.cc

pofd_affine_getR_SOURCES =global_settings.h numberCounts.h numberCountsKnots.h \
	numberCountsKnotsSpline.h paramSet.h affineExcept.h beam.h utility.h \
	numberCountsDouble.h numberCountsDoubleLogNormal.h doublebeam.h ran.h \
	numberCounts.cc numberCountsKnots.cc beam.cc utility.cc \
	numberCountsKnotsSpline.cc paramSet.cc affineExcept.cc ran.cc \
	numberCountsDouble.cc numberCountsDoubleLogNormal.cc doublebeam.cc \
	pofd_affine_getR.cc

pofd_affine_getPD_SOURCES = global_settings.h numberCounts.h \
	numberCountsKnots.h numberCountsKnotsSpline.h paramSet.h \
	affineExcept.h beam.h utility.h PD.h PDFactory.h fitsData.h \
	numberCountsDouble.h numberCountsDoubleLogNormal.h doublebeam.h \
	PDDouble.h PDFactoryDouble.h fitsDataDouble.h ran.h \
	numberCounts.cc numberCountsKnots.cc beam.cc utility.cc \
	numberCountsKnotsSpline.cc paramSet.cc affineExcept.cc \
	PD.cc PDFactory.cc fitsData.cc numberCountsDouble.cc \
	 numberCountsDoubleLogNormal.cc doublebeam.cc ran.cc \
	PDDouble.cc PDFactoryDouble.cc fitsDataDouble.cc \
	pofd_affine_getPD.cc 

pofd_affine_getLike_SOURCES = global_settings.h numberCounts.h \
	numberCountsKnots.h numberCountsKnotsSpline.h paramSet.h \
	affineExcept.h beam.h utility.h PD.h PDFactory.h fitsData.h \
	calcLike.h numberCountsDouble.h numberCountsDoubleLogNormal.h \
	doublebeam.h PDDouble.h PDFactoryDouble.h fitsDataDouble.h \
	calcLikeDouble.h ran.h specFile.h specFileDouble.h \
	numberCounts.cc numberCountsKnots.cc beam.cc utility.cc \
	numberCountsKnotsSpline.cc paramSet.cc affineExcept.cc ran.cc \
	PD.cc PDFactory.cc fitsData.cc calcLike.cc numberCountsDouble.cc \
	numberCountsDoubleLogNormal.cc doublebeam.cc PDDouble.cc \
	PDFactoryDouble.cc fitsDataDouble.cc calcLikeDouble.cc \
	specFile.cc specFileDouble.cc \
	pofd_affine_getLike.cc 

pofd_affine_mcmc_SOURCES = global_settings.h numberCounts.h \
	numberCountsKnots.h numberCountsKnotsSpline.h paramSet.h \
	affineExcept.h beam.h utility.h PD.h PDFactory.h fitsData.h \
	calcLike.h numberCountsDouble.h numberCountsDoubleLogNormal.h \
	doublebeam.h PDDouble.h PDFactoryDouble.h fitsDataDouble.h \
	calcLikeDouble.h ran.h specFile.h specFileDouble.h \
	affineEnsemble.h affineChainSet.h affineQueue.h proposedStep.h \
	pofdMCMC.h pofdMCMCDouble.h \
	numberCounts.cc numberCountsKnots.cc beam.cc utility.cc \
	numberCountsKnotsSpline.cc paramSet.cc affineExcept.cc ran.cc \
	PD.cc PDFactory.cc fitsData.cc calcLike.cc numberCountsDouble.cc \
	numberCountsDoubleLogNormal.cc doublebeam.cc PDDouble.cc \
	PDFactoryDouble.cc fitsDataDouble.cc calcLikeDouble.cc \
	specFile.cc specFileDouble.cc affineEnsemble.cc affineChainSet.cc \
	proposedStep.cc pofdMCMC.cc pofdMCMCDouble.cc \
	pofd_affine_getLike.cc 
