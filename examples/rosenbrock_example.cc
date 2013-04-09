#include<string>
#include<iostream>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/affineExcept.h"
#include "../include/affineEnsemble.h"
#include "../include/paramSet.h"

/*!
  \brief Rosenbrock density
*/
class rosenbrockDensity : public affineEnsemble {
private :
  double a1, a2;
public:
  rosenbrockDensity(double, double, unsigned int, 
		    unsigned int, unsigned int);
  void initChains();
  double getLogLike(const paramSet&);
};

rosenbrockDensity::rosenbrockDensity( double A1, double A2, 
				      unsigned int NWALKERS, 
				      unsigned int NPARAMS,
				      unsigned int NSAMPLES ) :
  affineEnsemble(NWALKERS,NPARAMS,NSAMPLES), a1(A1), a2(A2) { }
  

void rosenbrockDensity::initChains() {
  //Just generate random positions from zero to two for
  // initial vectors

  if (rank != 0) return;

  chains.clear();
  chains.addChunk(1);

  paramSet p( getNParams() );
  double logLike;
  unsigned int npar = getNParams();
  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    for (unsigned int j = 0; j < npar; ++j)
      p[j] = 2*rangen.doub();
    logLike = getLogLike(p);
    chains.addNewStep( i, p, logLike );
    naccept[i] = 1;
  }
  chains.setSkipFirst();
}

double rosenbrockDensity::getLogLike(const paramSet& p) {
  double val1, val2;
  val1 = (p[1]-p[0]*p[0]);
  val2 = 1 - p[0];
  return - (a1*val1*val1 + val2*val2)/a2;
}

int main(int argc, char** argv) {

  unsigned int nwalkers, nsamples;
  std::string outfile;
  double a1, a2;
  bool verbose;

  verbose = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help",no_argument,0,'h'},
    {"verbose",no_argument,0,'v'},
    {"Version",no_argument,0,'V'},
    {0,0,0,0}
  };
  char optstring[] = "hvV";

  int rank, nproc;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &nproc);

  if (nproc < 2) {
    if (rank == 0) {
      std::cerr << "Must run on multiple processes" << std::endl;
    }
    MPI::Finalize();
    return 1;
  }

  while ((c = getopt_long(argc,argv,optstring,long_options, 
			  &option_index)) != -1 ) 
    switch(c) {
    case 'h' :
      if (rank == 0) {
	std::cerr << "NAME" << std::endl;
	std::cerr << "\trosenbrock_example -- draws samples from a "
		  << " Rosenbrock " 
		  << std::endl;
	std::cerr << "\t\tdensity." << std::endl;
	std::cerr << std::endl;
	std::cerr << "SYNOPSIS" << std::endl;
	std::cerr << "\trosenbrock_example a1 a2 nwalkers nsamples outfile"
		  << std::endl;
	std::cerr << "DESCRIPTION:" << std::endl;
	std::cerr << "\tDraws samples from a Rosenbrock density using"
		  << std::endl;
	std::cerr << "\tan affine-invariant MCMC code, using nwalkers walkers"
		  << std::endl;
	std::cerr << "\tand generating (approximately) nsamples samples.  The"
		  << std::endl;
	std::cerr << "\tresults are written to outfile.  a1 and a2 are the" 
		  << std::endl;
	std::cerr << "\tparameters of the Rosenbrock density" << std::endl;
	std::cerr << "OPTIONS" << std::endl;
	std::cerr << "\t-h, --help" << std::endl;
	std::cerr << "\t\tOutput this help." << std::endl;
	std::cerr << "\t-v, --verbose" << std::endl;
	std::cerr << "\t\tOuput informational messages as it runs"
		  << std::endl;
	std::cerr << "\t-V, --version" << std::endl;
	std::cerr << "\t\tOuput the version number of mcmc_affine in use"
		  << std::endl;
      }
      MPI_Finalize();
      return 0;
      break;
    case 'v' :
      verbose = true;
      break;
    case 'V' :
      if (rank == 0) {
	std::cerr << "mcmc_affine version number: " << mcmc_affine::version 
		  << std::endl;
      }
      MPI_Finalize();
      return 0;
      break;
    }
  if ( optind >= argc-4 ) {
    if (rank == 0) {
      std::cerr << "Required arguments missing" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  a1       = atof( argv[optind] );
  a2       = atof( argv[optind+1] );
  nwalkers = atoi( argv[optind+2] );
  nsamples = atoi( argv[optind+3] );
  outfile  = std::string( argv[optind+4] );

  if (nwalkers == 0) return 0;
  if (nsamples == 0) return 0;

  //Hardwired cov matrix
  try {
    rosenbrockDensity rd(a1,a2,nwalkers,2,nsamples);
    if (verbose) rd.setVerbose();
    
    if (verbose && rank == 0)
      std::cout << "Entering main loop" << std::endl;
    if (rank == 0) rd.initChains();
    rd.doSteps(rd.getNSteps());
    
    if (rank == 0) {
      std::vector<double> accept;
      rd.getAcceptanceFrac(accept);
      double mnacc = accept[0];
      for (unsigned int i = 1; i < nwalkers; ++i)
	mnacc += accept[i];
      std::cout << "Mean acceptance: " << mnacc / static_cast<double>(nwalkers)
		<< std::endl;

      std::vector<double> acor;
      bool succ = rd.getAcor(acor);
      if (succ) std::cout << "Autocorrelation length: " << acor[0] 
			  << " " << acor[1] << std::endl; 
      else std::cout << "Failed to compute autocorrelation" << std::endl;
  
      //And write
      rd.writeToFile(outfile);
    }
  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    int jnk;
    if (rank == 0)
      for (int i = 1; i < nproc; ++i)
	MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
    MPI_Finalize();
    return 1;
  }
  MPI_Finalize();

  return 0;
}
