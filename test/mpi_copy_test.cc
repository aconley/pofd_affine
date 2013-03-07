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

  unsigned int rank, nproc;
  rank = MPI::COMM_WORLD.Get_rank();
  nproc = MPI::COMM_WORLD.Get_size();
  int jnk;

  if (rank != 0) {
    std::cerr << "Master should only be run on node 0" << std::endl;
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,STOP);
    return;
  }

  if (earlyexit) {
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,STOP);
    return;
  }

  if (nproc < 2) {
    std::cerr << "Need at least two processes" << std::endl;
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,STOP);
    return;
  }

  //Hardwired test files
  std::string fitsfile1 = "test/fiducial_2Dsim_PSW.fits";
  std::string fitsfile2 = "test/fiducial_2Dsim_PMW.fits";
  std::string psffile1 = "test/band1_beam.fits";
  std::string psffile2 = "test/band2_beam.fits";
    
  //Start the testing
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,STARTTESTS);

  //First, test beam
  try {
    if (verbose) std::cout << "Beam test" << std::endl;
    beam bm(psffile1, true); //Test with histogramming

    unsigned int n;
    n = bm.getNPos();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&n,1,MPI::UNSIGNED,i,SENDUIVAL);

    n = bm.getNNeg();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&n,1,MPI::UNSIGNED,i,SENDUIVAL);

    double val;
    val = bm.getEffectiveArea();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&val,1,MPI::DOUBLE,i,SENDDBLVAL);

    val = bm.getPixSize();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&val,1,MPI::DOUBLE,i,SENDDBLVAL);

    for (int i=1; i < static_cast<int>(nproc); ++i)
      bm.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "Beam test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //Next doublebeam
  try {
    if (verbose) std::cout << "Doublebeam test" << std::endl;
    doublebeam bm(psffile1,psffile2,true); //Test with histogramming
    for (int i=1; i < static_cast<int>(nproc); ++i)
      bm.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "Doublebeam test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }

  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //fitsData
  try {
    if (verbose) std::cout << "fitsData test" << std::endl;
    fitsData data(fitsfile1);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      data.sendSelf(MPI::COMM_WORLD,i);
    
    unsigned int n = data.getN();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&n,1,MPI::UNSIGNED,i,SENDUIVAL);

    double val = data.getMean();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&val,1,MPI::DOUBLE,i,SENDDBLVAL);
    val = data.getMin();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&val,1,MPI::DOUBLE,i,SENDDBLVAL);
    val = data.getMax();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&val,1,MPI::DOUBLE,i,SENDDBLVAL);

    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);

    if (verbose) std::cout << "fitsData test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }

  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //fitsDataDouble
  try {
    if (verbose) std::cout << "fitsDataDouble test" << std::endl;
    fitsDataDouble data(fitsfile1, fitsfile2);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      data.sendSelf(MPI::COMM_WORLD,i);

    unsigned int n = data.getN();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&n,1,MPI::UNSIGNED,i,SENDUIVAL);

    std::pair<double,double> pr;
    pr = data.getMean();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&pr.first,1,MPI::DOUBLE,i,SENDDBLVAL);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&pr.second,1,MPI::DOUBLE,i,SENDDBLVAL);
    pr = data.getMin();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&pr.first,1,MPI::DOUBLE,i,SENDDBLVAL);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&pr.second,1,MPI::DOUBLE,i,SENDDBLVAL);
    pr = data.getMax();
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&pr.first,1,MPI::DOUBLE,i,SENDDBLVAL);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Send(&pr.second,1,MPI::DOUBLE,i,SENDDBLVAL);

    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "fitsDataDouble test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }

  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //paramSet
  try {
    if (verbose) std::cout << "paramSet test" << std::endl;
    paramSet pars(9);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      pars.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "paramSet test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //proposedStep
  try {
    if (verbose) std::cout << "proposedStep test" << std::endl;
    proposedStep pr(5);
    pr.update_idx = 3;
    pr.oldLogLike = 5.0;
    pr.newLogLike = 2.0;
    for (int i=1; i < static_cast<int>(nproc); ++i)
      pr.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "proposedStep test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //numberCountsKnotsSpline
  try {
    if (verbose) std::cout << "numberCountsKnotsSpline test" << std::endl;
    const double kpos[] = {0.1, 0.3, 1.0};
    const double kval[] = {10, 5, 3};
    paramSet p(3, kval);
    numberCountsKnotsSpline model(3, kpos);
    model.setParams(p);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      model.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "numberCountsDoubleLogNormal test succeeded" 
			   << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //initFileKnots
  try {
    if (verbose) std::cout << "initFileKnots test" << std::endl;
    initFileKnots ifile;
    for (int i=1; i < static_cast<int>(nproc); ++i)
      ifile.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "initFileKnots test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);


  //numberCountsDoubleLogNormal
  try {
    if (verbose) std::cout << "numberCountsDoubleLogNormal test" << std::endl;
    numberCountsDoubleLogNormal model(8,4,3);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      model.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "numberCountsDoubleLogNormal test succeeded" 
			   << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //initFileDoubleLogNormal
  try {
    if (verbose) std::cout << "initFileDoubleLogNormal test" << std::endl;
    initFileDoubleLogNormal ifile;
    for (int i=1; i < static_cast<int>(nproc); ++i)
      ifile.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "initFileDoubleLogNormal test succeeded" 
			   << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //PDFactory
  try {
    if (verbose) std::cout << "PDFactory test" << std::endl;
    PDFactory pfactory(100);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      pfactory.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "PDFactory test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);


  //PDFactoryDouble
  try {
    if (verbose) std::cout << "PDFactoryDouble test" << std::endl;
    PDFactoryDouble pfactory(110);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      pfactory.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "PDFactoryDouble test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //calcLikeSingle
  try {
    if (verbose) std::cout << "calcLikeSingle test" << std::endl;
    calcLikeSingle like(50);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      like.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "calcLikeSingle test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //calcLike
  try {
    if (verbose) std::cout << "calcLike test" << std::endl;
    calcLike like(2048, 130, false);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      like.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "calcLike test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);


  //calcLikeDoubleSingle
  try {
    if (verbose) std::cout << "calcLikeDoubleSingle test" << std::endl;
    calcLikeDoubleSingle like(12);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      like.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "calcLikeDoubleSingle test succeeded" 
			   << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  //calcLikeDouble
  try {
    if (verbose) std::cout << "calcLike test" << std::endl;
    calcLikeDouble like(4096, 104, false, false);
    for (int i=1; i < static_cast<int>(nproc); ++i)
      like.sendSelf(MPI::COMM_WORLD,i);
    for (int i = 1; i < static_cast<int>(nproc); ++i)
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,TESTSUCCEEDED);
    if (verbose) std::cout << "calcLike test succeeded" << std::endl;
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,NEXTTEST);

  for (int i=1; i < static_cast<int>(nproc); ++i)
    MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,i,DONE);

  std::cout << "Tests passed" << std::endl;

}

void slave() {
  unsigned int rank;
  rank = MPI::COMM_WORLD.Get_rank();
  int jnk;

  //Wait for message, either starttests or stop
  MPI::Status Info;
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  int this_tag = Info.Get_tag();

  if (this_tag == STOP) return;

  if (this_tag != STARTTESTS) {
    std::cerr << "Unexpected tag: " << this_tag << " in slave: "
	      << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //First, beam
  try {
    beam bm;
    bm.recieveCopy(MPI::COMM_WORLD,0);

    unsigned int n;
    MPI::COMM_WORLD.Recv(&n,1,MPI::UNSIGNED,0,SENDUIVAL);
    if (n != bm.getNPos())
      throw affineExcept("test_copy","slave", "beam has wrong Npos",
			 34);
    MPI::COMM_WORLD.Recv(&n,1,MPI::UNSIGNED,0,SENDUIVAL);
    if (n != bm.getNNeg())
      throw affineExcept("test_copy","slave", "beam has wrong Nneg",
			 35);
    double val;
    MPI::COMM_WORLD.Recv(&val,1,MPI::DOUBLE,0,SENDDBLVAL);
    double diff = (bm.getEffectiveArea() - val)/val;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "beam has wrong Eff area",
			 36);

    MPI::COMM_WORLD.Recv(&val,1,MPI::DOUBLE,0,SENDDBLVAL);
    diff = (bm.getPixSize() - val)/val;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "beam has wrong pixel size",
			 37);      

    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }

  //Make sure we get nexttest
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after beam test in slave: "
	      << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //Doublebeam, same story
  try {
    doublebeam bm;
    bm.recieveCopy(MPI::COMM_WORLD,0);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after doublebeam test in slave: "
	      << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //fitsData
  try {
    fitsData data;
    data.recieveCopy(MPI::COMM_WORLD,0);

    //Image statistics
    unsigned int n;
    MPI::COMM_WORLD.Recv(&n,1,MPI::UNSIGNED,0,SENDUIVAL);
    if (n != data.getN())
      throw affineExcept("test_copy","slave", "fitsData had wrong data size",
			 23);
    double val;
    MPI::COMM_WORLD.Recv(&val,1,MPI::DOUBLE,0,SENDDBLVAL);
    double diff = (data.getMean() - val); //Mean close to zero
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "fitsData had wrong mean",
			 24);
    MPI::COMM_WORLD.Recv(&val,1,MPI::DOUBLE,0,SENDDBLVAL);
    diff = (data.getMin() - val)/val;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "fitsData had wrong min",
			 25);
    MPI::COMM_WORLD.Recv(&val,1,MPI::DOUBLE,0,SENDDBLVAL);
    diff = (data.getMax() - val)/val;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "fitsData had wrong max",
			 26);

    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after fitsData test in slave: "
	      << rank << std::endl;
    std::cerr << "Got message code: " << this_tag << std::endl;
    MPI::Finalize();
    return;
  }


  //fitsDataDouble, same story
  try {
    fitsDataDouble data;
    data.recieveCopy(MPI::COMM_WORLD,0);

    //Image statistics
    unsigned int n;
    MPI::COMM_WORLD.Recv(&n,1,MPI::UNSIGNED,0,SENDUIVAL);
    if (n != data.getN())
      throw affineExcept("test_copy","slave", 
			 "fitsDataDouble had wrong data size", 27);
    std::pair<double,double> pr;
    MPI::COMM_WORLD.Recv(&pr.first,1,MPI::DOUBLE,0,SENDDBLVAL);
    MPI::COMM_WORLD.Recv(&pr.second,1,MPI::DOUBLE,0,SENDDBLVAL);
    double diff = (data.getMean().first - pr.first); //Mean close to zero
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", 
			 "fitsDataDouble had wrong mean1", 28);
    diff = (data.getMean().second - pr.second);
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", 
			 "fitsDataDouble had wrong mean2", 29);

    MPI::COMM_WORLD.Recv(&pr.first,1,MPI::DOUBLE,0,SENDDBLVAL);
    MPI::COMM_WORLD.Recv(&pr.second,1,MPI::DOUBLE,0,SENDDBLVAL);
    diff = (data.getMin().first - pr.first)/pr.first;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "fitsData had wrong min1",
			 30);
    diff = (data.getMin().second - pr.second)/pr.second;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "fitsData had wrong min2",
			 31);

    MPI::COMM_WORLD.Recv(&pr.first,1,MPI::DOUBLE,0,SENDDBLVAL);
    MPI::COMM_WORLD.Recv(&pr.second,1,MPI::DOUBLE,0,SENDDBLVAL);
    diff = (data.getMax().first - pr.first)/pr.first;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "fitsData had wrong max1",
			 32);
    diff = (data.getMax().second - pr.second)/pr.second;
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", "fitsData had wrong max2",
			 33);

    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after fitsDataDouble test in slave: "
	      << rank << std::endl;
    std::cerr << "Got message code: " << this_tag << std::endl;
    MPI::Finalize();
    return;
  }

  //paramSet
  try {
    paramSet pars;
    pars.recieveCopy(MPI::COMM_WORLD,0);
    if (pars.getNParams() != 9)
      throw affineExcept("test_copy","slave", "paramSet should have 9 params",
			 1);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after paramSet test in slave: "
	      << rank << std::endl;
    std::cerr << "Expected: " << NEXTTEST << " got " << this_tag
	      << std::endl;
    MPI::Finalize();
    return;
  }

  //proposedStep
  try {
    proposedStep pr(0);
    pr.recieveCopy(MPI::COMM_WORLD,0);
    if (pr.oldStep.getNParams() != 5)
      throw affineExcept("test_copy","slave", 
			 "proposedStep.oldStep should have 5 params",1);
    if (pr.newStep.getNParams() != 5)
      throw affineExcept("test_copy","slave", 
			 "proposedStep.newStep should have 5 params",2);
    if (pr.update_idx != 3)
      throw affineExcept("test_copy","slave", 
			 "proposedStep.update_idx should have 3 params", 3);
    double diff = fabs((pr.oldLogLike - 5.0)/5.0);
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", 
			 "proposedStep.oldLogLike should be 5.0", 4);
    diff = fabs((pr.newLogLike - 2.0)/2.0);
    if (fabs(diff) > 1e-5)
      throw affineExcept("test_copy","slave", 
			 "proposedStep.newLogLike should be 5.0", 5);

    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after proposedStep test in slave: "
	      << rank << std::endl;
    std::cerr << "Expected: " << NEXTTEST << " got " << this_tag
	      << std::endl;
    MPI::Finalize();
    return;
  }
  
  //numberCountsKnotsSpline
  try {
    const double kpos[] = {0.1, 0.3, 1.0};
    const double kval[] = {10, 5, 3};
    numberCountsKnotsSpline model;
    model.recieveCopy(MPI::COMM_WORLD,0);
    if (model.getNKnots() != 3)
      throw affineExcept("test_copy","slave", 
			 "numberCountsKnotsSpline.getNKnots() should be 3", 6);
    std::vector<double> kposvc;
    model.getKnotPositions(kposvc);
    double dist = (kposvc[0] - kpos[0]) * (kposvc[0] - kpos[0]) +
      (kposvc[1] - kpos[1]) * (kposvc[1] - kpos[1]) +
      (kposvc[2] - kpos[2]) * (kposvc[2] - kpos[2]);
    if (dist > 1e-5)
      throw affineExcept("test_copy","slave", 
			 "numberCountsKnotsSpline knot positions are off", 21);
    paramSet p(3);
    model.getParams(p);
    dist = (p[0] - kval[0]) * (p[0] - kval[0]) +
      (p[1] - kval[1]) * (p[1] - kval[1]) +
      (p[2] - kval[2]) * (p[2] - kval[2]);
    if (dist > 1e-5)
      throw affineExcept("test_copy","slave", 
			 "numberCountsKnotsSpline knot values are off", 22);
      
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after numberCountsKnotsSpline test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //initFileKnots
  try {
    initFileKnots ifile;
    ifile.recieveCopy(MPI::COMM_WORLD,0);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after initFileKnots test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //numberCountsDoubleLogNormal
  try {
    numberCountsDoubleLogNormal model;
    model.recieveCopy(MPI::COMM_WORLD,0);
    if (model.getNKnots() != 8)
      throw affineExcept("test_copy","slave", 
			 "numberCountsDoubleLogNormal.getNKnots() should be 8",
			 7);
    if (model.getNSigmas() != 4)
      throw affineExcept("test_copy","slave", 
			 "numberCountsDoubleLogNormal.getNSigmas() should be 4",
			 8);
    if (model.getNOffsets() != 3)
      throw affineExcept("test_copy","slave", 
			 "numberCountsDoubleLogNormal.getNOffsets() should be 3",
			 9);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after numberCountsDoubleLogNormal test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //initFileDoubleLogNormal
  try {
    initFileDoubleLogNormal ifile;
    ifile.recieveCopy(MPI::COMM_WORLD,0);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after initFileDoubleLogNormal test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //PDFactory
  try {
    PDFactory pfactory;
    pfactory.recieveCopy(MPI::COMM_WORLD,0);
    if (pfactory.getNInterp() != 100)
      throw affineExcept("test_copy","slave", 
			 "PDFactory.getNInterp() should be 100", 10);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after PDFactory test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }


  //PDFactoryDouble
  try {
    PDFactoryDouble pfactory;
    pfactory.recieveCopy(MPI::COMM_WORLD,0);
    if (pfactory.getNEdge() != 110)
      throw affineExcept("test_copy","slave", 
			 "PDFactoryDouble.getNEdge() should be 110", 11);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after PDFactoryDouble test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //calcLikeSingle
  try {
    calcLikeSingle like;
    like.recieveCopy(MPI::COMM_WORLD,0);
    if (like.getNInterp() != 50)
      throw affineExcept("test_copy","slave", 
			 "calcLikeSingle.getNInterp() should be 50", 12);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after calcLikeSingle test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //calcLike
  try {
    calcLike like;
    like.recieveCopy(MPI::COMM_WORLD,0);
    if (like.getFFTSize() != 2048)
      throw affineExcept("test_copy","slave", 
			 "calcLike.getFFTSize() should be 2048", 13);
    if (like.getNInterp() != 130) {
      std::stringstream errstr;
      errstr << "calcLike.getNInterp() should be 130 but is "
	     << like.getNInterp();
      throw affineExcept("test_copy","slave", errstr.str(), 14);
    }
    if (like.getEdgeFix())
      throw affineExcept("test_copy","slave", 
			 "calcLike.getEdgeFix() should be false", 15);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after calcLike test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //calcLikeDoubleSingle
  try {
    calcLikeDoubleSingle like;
    like.recieveCopy(MPI::COMM_WORLD,0);
    if (like.getNEdge() != 12)
      throw affineExcept("test_copy","slave", 
			 "calcLikeDoubleSingle.getNEdge() should be 12", 16);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after calcLikeDoubleSingle test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }

  //calcLikeDouble
  try {
    calcLikeDouble like;
    like.recieveCopy(MPI::COMM_WORLD,0);
    if (like.getFFTSize() != 4096)
      throw affineExcept("test_copy","slave", 
			 "calcLikeDouble.getFFTSize() should be 4096", 17);
    if (like.getNEdge() != 104)
      throw affineExcept("test_copy","slave", 
			 "calcLike.getNEdge() should be 104", 18);
    if (like.getEdgeFix())
      throw affineExcept("test_copy","slave", 
			 "calcLikeDouble.getEdgeFix() should be false", 19);
    if (like.getEdgeInteg())
      throw affineExcept("test_copy","slave", 
			 "calcLikeDouble.getEdgeInteg() should be false", 20);
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,TESTSUCCEEDED);
  } catch ( const affineExcept& ex ) {
    std::cerr << ex << std::endl;
    MPI::Finalize();
    return;
  }
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);
  this_tag = Info.Get_tag();
  if (this_tag != NEXTTEST) {
    std::cerr << "MPI error encountered after calcLike test "
	      << "in slave: " << rank << std::endl;
    MPI::Finalize();
    return;
  }



  MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,DONE);

  std::cout << "Slave " << rank << " passed all tests" << std::endl;

}

int main(int argc, char **argv) {
  unsigned int rank;
  MPI::Init(argc,argv);
  rank = MPI::COMM_WORLD.Get_rank();

  if (rank == 0) master(argc,argv); else slave();
  MPI::Finalize();
  return 0;
  
}
