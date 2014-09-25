//Program to test MPI copying of various classes
//Not all are tested yet
//To some extent this just tests to make sure that all the
// messages that get sent are taken.  It would also be good to 
// test the internal values, but I have only done that in a few
// limited cases so far.

#include<iostream>
#include<string>
#include<cmath>
#include<sstream>
#include<utility>

#include<mpi.h>
#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/affineExcept.h"
#include "../include/beam.h"
#include "../include/doublebeam.h"
#include "../include/paramSet.h"
#include "../include/proposedStep.h"
#include "../include/numberCountsKnotsSpline.h"
#include "../include/numberCountsDoubleLogNormal.h"
#include "../include/fitsData.h"
#include "../include/fitsDataDouble.h"
#include "../include/PDFactory.h"
#include "../include/PDFactoryDouble.h"
#include "../include/calcLike.h"
#include "../include/calcLikeDouble.h"

enum test_enums { STOP=1000000, STARTTESTS=1000001, 
		  TESTSUCCEEDED=1000002, NEXTTEST=1000003,
		  DONE=1000004, SENDUIVAL=1000005, SENDDBLVAL=1000006 };

void master( int argc, char **argv ) {
  //Process arguments
  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help",no_argument,0,'h'},
    {"verbose",no_argument,0,'v'},
    {"version",no_argument,0,'V'},
    {0,0,0,0}
  };
  
  bool earlyexit = false;
  bool verbose   = false;

  char optstring[] = "hvV";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\ttest_copy -- test MPI copying operations"
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\ttest_copy [options]" << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-h, --help" << std::endl;
      std::cerr << "\t\tPrint this message and exit" << std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tOutput informational messages as tests run"
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tPrint the version number and exit" << std::endl;
      earlyexit = true;
      break;
    case 'v' :
      verbose = true;
      break;
    case 'V' :
      std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
		<< std::endl;
      earlyexit = true;
      break;
    }

  int rank, nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int jnk;

  if (rank != 0) {
    std::cerr << "Master should only be run on node 0" << std::endl;
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&jnk, 1, MPI_INT, i, STOP, MPI_COMM_WORLD);
    return;
  }

  if (earlyexit) {
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&jnk, 1, MPI_INT, i, STOP, MPI_COMM_WORLD);
    return;
  }

  if (nproc < 2) {
    std::cerr << "Need at least two processes" << std::endl;
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&jnk, 1, MPI_INT, i, STOP, MPI_COMM_WORLD);
    return;
  }

  //Hardwired test files
  std::string fitsfile1 = "testdata/testmodel2D_band1.fits";
  std::string fitsfile2 = "testdata/testmodel2D_band2.fits";
  std::string psffile1 = "testdata/band1_beam.fits";
  std::string psffile2 = "testdata/band2_beam.fits";
  std::string psffile1_f100 = "testdata/band12_beam_f100.fits[0]";
  std::string psffile2_f100 = "testdata/band12_beam_f100.fits[1]";
    

  MPI_Status Info;

  //Start the testing
  for (int i = 1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, STARTTESTS, MPI_COMM_WORLD);

  //First, test beam
  try {
    if (verbose) std::cout << "Beam test" << std::endl;
    beam bm(psffile1, true); //Test with histogramming

    unsigned int n;
    n = bm.getNPos();
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&n, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);

    n = bm.getNNeg();
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&n, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);

    n = bm.getNHistPos();
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&n, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);

    double val;
    val = bm.getEffectiveArea();
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&val, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);

    val = bm.getPixSize();
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&val, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);

    for (int i = 1; i < nproc; ++i)
      bm.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "Beam test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //Next doublebeam.  Test on filtered beam
  try {
    if (verbose) std::cout << "Doublebeam test" << std::endl;
    doublebeam bm(psffile1_f100, psffile2_f100, true); //Test with histogramming
    for (int i=1; i < nproc; ++i)
      bm.sendSelf(MPI_COMM_WORLD, i);

    unsigned int hassgn;
    for (unsigned int i = 0; i < 4; ++i) {
      hassgn = bm.hasSign(i) ? 1 : 0;
      for (int i = 1; i < nproc; ++i)
	MPI_Send(&hassgn, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);
    }

    unsigned int npix, ishist, nhist;
    for (unsigned int i = 0; i < 4; ++i) {
      hassgn = bm.hasSign(i) ? 1 : 0;
      if (hassgn) {
	npix = bm.getNPix(i);
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&npix, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);
	dblpair pr = bm.getMinMax1(i);
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&pr.first, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&pr.second, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
	pr = bm.getMinMax2(i);
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&pr.first, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&pr.second, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
	ishist = bm.isHistogrammed(i) ? 1 : 0;
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&ishist, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);
	if (ishist) {
	  nhist = bm.getNHist(i);
	  for (int i = 1; i < nproc; ++i)
	    MPI_Send(&nhist, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);
	}
      }
      
    }
    for (int i = 1; i < nproc; ++i)
	MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "Doublebeam test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }

  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //fitsData
  try {
    if (verbose) std::cout << "fitsData test" << std::endl;
    fitsData data(fitsfile1);
    for (int i=1; i < nproc; ++i)
      data.sendSelf(MPI_COMM_WORLD, i);
    
    unsigned int n = data.getN();
    for (int i=1; i < nproc; ++i)
      MPI_Send(&n, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);

    double val = data.getMean();
    for (int i=1; i < nproc; ++i)
      MPI_Send(&val, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
    val = data.getMin();
    for (int i=1; i < nproc; ++i)
      MPI_Send(&val, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
    val = data.getMax();
    for (int i=1; i < nproc; ++i)
      MPI_Send(&val, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);

    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);

    if (verbose) std::cout << "fitsData test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }

  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //fitsDataDouble
  try {
    if (verbose) std::cout << "fitsDataDouble test" << std::endl;
    fitsDataDouble data(fitsfile1, fitsfile2);
    for (int i=1; i < nproc; ++i)
      data.sendSelf(MPI_COMM_WORLD, i);

    unsigned int n = data.getN();
    for (int i=1; i < nproc; ++i)
      MPI_Send(&n, 1, MPI_UNSIGNED, i, SENDUIVAL, MPI_COMM_WORLD);

    std::pair<double,double> pr;
    pr = data.getMean();
    for (int i=1; i < nproc; ++i)
      MPI_Send(&pr.first, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
    for (int i=1; i < nproc; ++i)
      MPI_Send(&pr.second, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
    pr = data.getMin();
    for (int i=1; i < nproc; ++i)
      MPI_Send(&pr.first, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
    for (int i=1; i < nproc; ++i)
      MPI_Send(&pr.second, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
    pr = data.getMax();
    for (int i=1; i < nproc; ++i)
      MPI_Send(&pr.first, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
    for (int i=1; i < nproc; ++i)
      MPI_Send(&pr.second, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);

    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "fitsDataDouble test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }

  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //paramSet
  try {
    if (verbose) std::cout << "paramSet test" << std::endl;
    paramSet pars(9);
    for (int i=1; i < nproc; ++i)
      pars.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "paramSet test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //proposedStep
  try {
    if (verbose) std::cout << "proposedStep test" << std::endl;
    proposedStep pr(5);
    pr.update_idx = 3;
    pr.oldLogLike = 5.0;
    pr.newLogLike = 2.0;
    for (int i=1; i < nproc; ++i)
      pr.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "proposedStep test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //numberCountsKnotsSpline
  try {
    if (verbose) std::cout << "numberCountsKnotsSpline test" << std::endl;
    const float kpos[] = {0.1, 0.3, 1.0};
    const float kval[] = {10, 5, 3};
    paramSet p(3, kval);
    numberCountsKnotsSpline model(3, kpos);
    double nkval = model.getNumberCounts(0.2);
    model.setParams(p);
    for (int i=1; i < nproc; ++i)
      model.sendSelf(MPI_COMM_WORLD, i);
    for (int i=1; i < nproc; ++i)
      MPI_Send(&nkval, 1, MPI_DOUBLE, i, SENDDBLVAL, MPI_COMM_WORLD);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "numberCountsDoubleLogNormal test succeeded" 
			   << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //initFileKnots
  try {
    if (verbose) std::cout << "initFileKnots test" << std::endl;
    initFileKnots ifile;
    for (int i=1; i < nproc; ++i)
      ifile.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "initFileKnots test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);


  //numberCountsDoubleLogNormal
  try {
    if (verbose) std::cout << "numberCountsDoubleLogNormal test" << std::endl;
    numberCountsDoubleLogNormal model(8, 4, 3);
    for (int i=1; i < nproc; ++i)
      model.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "numberCountsDoubleLogNormal test succeeded" 
			   << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //initFileDoubleLogNormal
  try {
    if (verbose) std::cout << "initFileDoubleLogNormal test" << std::endl;
    initFileDoubleLogNormal ifile;
    for (int i=1; i < nproc; ++i)
      ifile.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "initFileDoubleLogNormal test succeeded" 
			   << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //PDFactory
  try {
    if (verbose) std::cout << "PDFactory test" << std::endl;
    PDFactory pfactory(100);
    for (int i=1; i < nproc; ++i)
      pfactory.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "PDFactory test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);


  //PDFactoryDouble
  try {
    if (verbose) std::cout << "PDFactoryDouble test" << std::endl;
    PDFactoryDouble pfactory(110);
    for (int i=1; i < nproc; ++i)
      pfactory.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "PDFactoryDouble test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //calcLikeSingle
  try {
    if (verbose) std::cout << "calcLikeSingle test" << std::endl;
    calcLikeSingle like(50);
    for (int i=1; i < nproc; ++i)
      like.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "calcLikeSingle test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //calcLike
  try {
    if (verbose) std::cout << "calcLike test" << std::endl;
    calcLike like(2048, 130, false);
    for (int i=1; i < nproc; ++i)
      like.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "calcLike test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);


  //calcLikeDoubleSingle
  try {
    if (verbose) std::cout << "calcLikeDoubleSingle test" << std::endl;
    calcLikeDoubleSingle like(12);
    for (int i=1; i < nproc; ++i)
      like.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "calcLikeDoubleSingle test succeeded" 
			   << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  //calcLikeDouble
  try {
    if (verbose) std::cout << "calcLike test" << std::endl;
    calcLikeDouble like(4096, 104, false, false);
    for (int i=1; i < nproc; ++i)
      like.sendSelf(MPI_COMM_WORLD, i);
    for (int i = 1; i < nproc; ++i)
      MPI_Recv(&jnk, 1, MPI_INT, i, TESTSUCCEEDED, MPI_COMM_WORLD, &Info);
    if (verbose) std::cout << "calcLike test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  for (int i=1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, NEXTTEST, MPI_COMM_WORLD);

  for (int i=1; i < nproc; ++i)
    MPI_Recv(&jnk, 1, MPI_INT, i, DONE, MPI_COMM_WORLD, &Info);

  std::cout << "Tests passed" << std::endl;

}

void slave() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int jnk;

  //Wait for message, either starttests or stop
  MPI_Status Info;
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  int this_tag = Info.MPI_TAG;

  if (this_tag == STOP) return;

  if (this_tag != STARTTESTS) {
    std::cerr << "Unexpected tag: " << this_tag << " in slave: "
	      << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //First, beam
  try {
    beam bm;
    bm.receiveCopy(MPI_COMM_WORLD, 0);

    unsigned int n;
    MPI_Recv(&n, 1, MPI_UNSIGNED, 0, SENDUIVAL, MPI_COMM_WORLD, &Info);
    if (n != bm.getNPos())
      throw affineExcept("test_copy", "slave", "beam has wrong Npos");
    MPI_Recv(&n, 1, MPI_UNSIGNED, 0, SENDUIVAL, MPI_COMM_WORLD, &Info);
    if (n != bm.getNNeg())
      throw affineExcept("test_copy", "slave", "beam has wrong Nneg");
    MPI_Recv(&n, 1, MPI_UNSIGNED, 0, SENDUIVAL, MPI_COMM_WORLD, &Info);
    if (n != bm.getNHistPos())
      throw affineExcept("test_copy", "slave", "beam has wrong NHistPos");
    double val;
    MPI_Recv(&val, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    double diff = (bm.getEffectiveArea() - val)/val;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "beam has wrong Eff area");

    MPI_Recv(&val, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    diff = (bm.getPixSize() - val)/val;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "beam has wrong pixel size");

    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }

  //Make sure we get nexttest
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after beam test in slave: "
	      << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //Doublebeam, same story
  try {
    doublebeam bm;
    bm.receiveCopy(MPI_COMM_WORLD, 0);

    unsigned int hassgn;
    for (unsigned int i = 0; i < 4; ++i) {
      MPI_Recv(&hassgn, 1, MPI_UNSIGNED, 0, SENDUIVAL, MPI_COMM_WORLD, &Info);
      if (hassgn != (bm.hasSign(i) ? 1 : 0))
	throw affineExcept("test_copy", "slave", 
			   "beam doesn't match presence of component");
    }
    
    double val, diff;
    dblpair pr;
    unsigned int npix, ishist, nhist;
    for (unsigned int i = 0; i < 4; ++i)
      if (bm.hasSign(i)) {
	MPI_Recv(&npix, 1, MPI_UNSIGNED, 0, SENDUIVAL, MPI_COMM_WORLD, &Info);
	if (npix != bm.getNPix(i))
	  throw affineExcept("test_copy", "slave", 
			     "beam doesn't match number of elements in component");
	
	pr = bm.getMinMax1(i);
	MPI_Recv(&val, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
	diff = fabs(pr.first - val);
	if (diff > 1e-5)
	  throw affineExcept("test_copy", "slave", 
			     "Wrong min in band 1 doublebeam component");
	MPI_Recv(&val, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
	diff = fabs(pr.second - val);
	if (diff > 1e-5)
	  throw affineExcept("test_copy", "slave", 
			     "Wrong max in band 1 doublebeam component");
	pr = bm.getMinMax2(i);
	MPI_Recv(&val, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
	diff = fabs(pr.first - val);
	if (diff > 1e-5)
	  throw affineExcept("test_copy", "slave", 
			     "Wrong min in band 2 doublebeam component");
	MPI_Recv(&val, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
	diff = fabs(pr.second - val);
	if (diff > 1e-5)
	  throw affineExcept("test_copy", "slave", 
			     "Wrong max in band 2 doublebeam component");
	
	MPI_Recv(&ishist, 1, MPI_UNSIGNED, 0, SENDUIVAL, MPI_COMM_WORLD, &Info);
	if (ishist != (bm.isHistogrammed(i) ? 1 : 0))
	  throw affineExcept("test_copy", "slave", 
			     "doublebeam doesn't match histogram presence component");	
	if (bm.isHistogrammed(i)) {
	  MPI_Recv(&nhist, 1, MPI_UNSIGNED, 0, SENDUIVAL,
		   MPI_COMM_WORLD, &Info);
	  if (nhist != bm.getNHist(i))
	    throw affineExcept("test_copy", "slave", 
			       "doublebeam doesn't match number of hist bins");	
	}
      }
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }

  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after doublebeam test in slave: "
	      << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //fitsData
  try {
    fitsData data;
    data.receiveCopy(MPI_COMM_WORLD, 0);

    //Image statistics
    unsigned int n;
    MPI_Recv(&n, 1, MPI_UNSIGNED, 0, SENDUIVAL, MPI_COMM_WORLD, &Info);
    if (n != data.getN())
      throw affineExcept("test_copy", "slave", "fitsData had wrong data size");
    double val;
    MPI_Recv(&val, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    double diff = (data.getMean() - val); //Mean close to zero
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "fitsData had wrong mean");
    MPI_Recv(&val, 1, MPI_DOUBLE, 0 ,SENDDBLVAL, MPI_COMM_WORLD, &Info);
    diff = (data.getMin() - val)/val;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "fitsData had wrong min");
    MPI_Recv(&val, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    diff = (data.getMax() - val)/val;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "fitsData had wrong max");

    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after fitsData test in slave: "
	      << rank << std::endl;
    std::cerr << "Got message code: " << this_tag << std::endl;
    MPI_Finalize();
    return;
  }


  //fitsDataDouble, same story
  try {
    fitsDataDouble data;
    data.receiveCopy(MPI_COMM_WORLD, 0);

    //Image statistics
    unsigned int n;
    MPI_Recv(&n, 1, MPI_UNSIGNED, 0, SENDUIVAL, MPI_COMM_WORLD, &Info);
    if (n != data.getN())
      throw affineExcept("test_copy", "slave", 
			 "fitsDataDouble had wrong data size");
    std::pair<double,double> pr;
    MPI_Recv(&pr.first, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    MPI_Recv(&pr.second, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    double diff = (data.getMean().first - pr.first); //Mean close to zero
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", 
			 "fitsDataDouble had wrong mean1");
    diff = (data.getMean().second - pr.second);
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", 
			 "fitsDataDouble had wrong mean2");

    MPI_Recv(&pr.first, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    MPI_Recv(&pr.second, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    diff = (data.getMin().first - pr.first)/pr.first;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "fitsData had wrong min1");
    diff = (data.getMin().second - pr.second)/pr.second;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "fitsData had wrong min2");

    MPI_Recv(&pr.first, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    MPI_Recv(&pr.second, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    diff = (data.getMax().first - pr.first)/pr.first;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "fitsData had wrong max1");
    diff = (data.getMax().second - pr.second)/pr.second;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", "fitsData had wrong max2");

    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after fitsDataDouble test in slave: "
	      << rank << std::endl;
    std::cerr << "Got message code: " << this_tag << std::endl;
    MPI_Finalize();
    return;
  }

  //paramSet
  try {
    paramSet pars;
    pars.receiveCopy(MPI_COMM_WORLD, 0);
    if (pars.getNParams() != 9)
      throw affineExcept("test_copy", "slave", "paramSet should have 9 params");
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after paramSet test in slave: "
	      << rank << std::endl;
    std::cerr << "Expected: " << NEXTTEST << " got " << this_tag
	      << std::endl;
    MPI_Finalize();
    return;
  }

  //proposedStep
  try {
    proposedStep pr(0);
    pr.receiveCopy(MPI_COMM_WORLD, 0);
    if (pr.oldStep.getNParams() != 5)
      throw affineExcept("test_copy", "slave",
			 "proposedStep.oldStep should have 5 params");
    if (pr.newStep.getNParams() != 5)
      throw affineExcept("test_copy", "slave", 
			 "proposedStep.newStep should have 5 params");
    if (pr.update_idx != 3)
      throw affineExcept("test_copy", "slave", 
			 "proposedStep.update_idx should have 3 params");
    double diff = fabs((pr.oldLogLike - 5.0)/5.0);
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", 
			 "proposedStep.oldLogLike should be 5.0");
    diff = fabs((pr.newLogLike - 2.0)/2.0);
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy", "slave", 
			 "proposedStep.newLogLike should be 5.0");

    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after proposedStep test in slave: "
	      << rank << std::endl;
    std::cerr << "Expected: " << NEXTTEST << " got " << this_tag
	      << std::endl;
    MPI_Finalize();
    return;
  }
  
  //numberCountsKnotsSpline
  try {
    const double kpos[] = {0.1, 0.3, 1.0};
    const double kval[] = {10, 5, 3};
    numberCountsKnotsSpline model;
    double nkval, nkval_this;
    model.receiveCopy(MPI_COMM_WORLD, 0);
    MPI_Recv(&nkval, 1, MPI_DOUBLE, 0, SENDDBLVAL, MPI_COMM_WORLD, &Info);
    if (model.getNKnots() != 3)
      throw affineExcept("test_copy", "slave", 
			 "numberCountsKnotsSpline.getNKnots() should be 3");
    std::vector<double> kposvc;
    model.getKnotPositions(kposvc);
    double dist = (kposvc[0] - kpos[0]) * (kposvc[0] - kpos[0]) +
      (kposvc[1] - kpos[1]) * (kposvc[1] - kpos[1]) +
      (kposvc[2] - kpos[2]) * (kposvc[2] - kpos[2]);
    if (dist > 1e-5)
      throw affineExcept("test_copy", "slave", 
			 "numberCountsKnotsSpline knot positions are off");
    paramSet p(3);
    model.getParams(p);
    dist = (p[0] - kval[0]) * (p[0] - kval[0]) +
      (p[1] - kval[1]) * (p[1] - kval[1]) +
      (p[2] - kval[2]) * (p[2] - kval[2]);
    if (dist > 1e-5)
      throw affineExcept("test_copy", "slave", 
			 "numberCountsKnotsSpline knot values are off");
      
    nkval_this = model.getNumberCounts(0.2);
    if (fabs((nkval_this - nkval)/nkval) > 1e-4)
      throw affineExcept("test_copy", "slave", 
			 "numberCountsKnotsSpline number counts are off");
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after numberCountsKnotsSpline test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //initFileKnots
  try {
    initFileKnots ifile;
    ifile.receiveCopy(MPI_COMM_WORLD, 0);
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after initFileKnots test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //numberCountsDoubleLogNormal
  try {
    numberCountsDoubleLogNormal model;
    model.receiveCopy(MPI_COMM_WORLD, 0);
    if (model.getNKnots() != 8)
      throw affineExcept("test_copy", "slave", 
			 "numberCountsDoubleLogNormal.getNKnots() should be 8");
    if (model.getNSigmas() != 4)
      throw affineExcept("test_copy", "slave", 
			 "numberCountsDoubleLogNormal.getNSigmas() should be 4");
    if (model.getNOffsets() != 3)
      throw affineExcept("test_copy", "slave", 
			 "numberCountsDoubleLogNormal.getNOffsets() should be 3");
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after numberCountsDoubleLogNormal test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //initFileDoubleLogNormal
  try {
    initFileDoubleLogNormal ifile;
    ifile.receiveCopy(MPI_COMM_WORLD, 0);
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after initFileDoubleLogNormal test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //PDFactory
  try {
    PDFactory pfactory;
    pfactory.receiveCopy(MPI_COMM_WORLD, 0);
    if (pfactory.getNInterp() != 100)
      throw affineExcept("test_copy", "slave", 
			 "PDFactory.getNInterp() should be 100");
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after PDFactory test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }


  //PDFactoryDouble
  try {
    PDFactoryDouble pfactory;
    pfactory.receiveCopy(MPI_COMM_WORLD, 0);
    if (pfactory.getNEdge() != 110)
      throw affineExcept("test_copy", "slave", 
			 "PDFactoryDouble.getNEdge() should be 110");
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after PDFactoryDouble test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //calcLikeSingle
  try {
    calcLikeSingle like;
    like.receiveCopy(MPI_COMM_WORLD, 0);
    if (like.getNInterp() != 50)
      throw affineExcept("test_copy", "slave", 
			 "calcLikeSingle.getNInterp() should be 50");
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after calcLikeSingle test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //calcLike
  try {
    calcLike like;
    like.receiveCopy(MPI_COMM_WORLD, 0);
    if (like.getFFTSize() != 2048)
      throw affineExcept("test_copy", "slave", 
			 "calcLike.getFFTSize() should be 2048");
    if (like.getNInterp() != 130) {
      std::stringstream errstr;
      errstr << "calcLike.getNInterp() should be 130 but is "
	     << like.getNInterp();
      throw affineExcept("test_copy", "slave", errstr.str());
    }
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after calcLike test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //calcLikeDoubleSingle
  try {
    calcLikeDoubleSingle like;
    like.receiveCopy(MPI_COMM_WORLD, 0);
    if (like.getNEdge() != 12)
      throw affineExcept("test_copy", "slave", 
			 "calcLikeDoubleSingle.getNEdge() should be 12");
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after calcLikeDoubleSingle test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  //calcLikeDouble
  try {
    calcLikeDouble like;
    like.receiveCopy(MPI_COMM_WORLD, 0);
    if (like.getFFTSize() != 4096)
      throw affineExcept("test_copy", "slave", 
			 "calcLikeDouble.getFFTSize() should be 4096");
    if (like.getNEdge() != 104)
      throw affineExcept("test_copy", "slave", 
			 "calcLike.getNEdge() should be 104");
    if (like.getEdgeInteg())
      throw affineExcept("test_copy", "slave", 
			 "calcLikeDouble.getEdgeInteg() should be false");
    MPI_Send(&jnk, 1, MPI_INT, 0, TESTSUCCEEDED, MPI_COMM_WORLD);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI_Finalize();
    return;
  }
  MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);
  this_tag = Info.MPI_TAG;
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after calcLike test "
	      << "in slave: " << rank << std::endl;
    MPI_Finalize();
    return;
  }

  MPI_Send(&jnk, 1, MPI_INT, 0, DONE, MPI_COMM_WORLD);

  std::cout << "Slave " << rank << " passed all tests" << std::endl;

}

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int rank, nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (nproc == 0) {
    if (rank == 0)
      std::cerr << "Must run on more than one process";
    MPI_Finalize();
    return 1;
  }

  if (rank == 0) master(argc,argv); else slave();

  MPI_Finalize();

  return 0;
  
}
