//Multi-dimensional Gaussian example

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
  \brief Multi-dimensional Gaussian, centered at 0.5 in all dimensions
*/
class multiGauss : public affineEnsemble {
private:
  double* invCovMatrix; //!< Inverse covariance matrix, nparams by nparams
  double* work1; //!< Working vector
  double* work2; //!< Working vector
public:
  multiGauss(unsigned int, unsigned int, unsigned int, 
	     const std::string&);
  virtual ~multiGauss();

  void initChains();
  double getLogLike(const paramSet&);
 
};

multiGauss::multiGauss( unsigned int NWALKERS, unsigned int NPARAMS,
			unsigned int NSAMPLES, const std::string& filename ) :
  affineEnsemble(NWALKERS,NPARAMS,NSAMPLES) {
  if (NPARAMS > 0) {
    invCovMatrix = new double[NPARAMS*NPARAMS];
    work1 = new double[NPARAMS];
    work2 = new double[NPARAMS];
  } else {
    invCovMatrix=NULL;
    work1 = NULL;
    work2 = NULL;
  }
  
  //Read in file
  std::ifstream ifs(filename.c_str());
  if (!ifs)
    throw affineExcept("multiGauss","multiGauss",
		       "Couldn't open input file",1);
  for (unsigned int i = 0; i < NPARAMS*NPARAMS; ++i)
    ifs >> invCovMatrix[i];
  ifs.close();

}

multiGauss::~multiGauss() {
  if (invCovMatrix != NULL) delete[] invCovMatrix;
  if (work1 != NULL) delete[] work1;
  if (work2 != NULL) delete[] work2;
}

void multiGauss::initChains() {
  //Just generate random positions from zero to one for
  // initial vectors

  if (rank != 0) return;

  //For master
  chains.clear();
  chains.addChunk(1);
    
  paramSet p( getNParams() );
  double logLike;
  unsigned int npar = getNParams();
  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    for (unsigned int j = 0; j < npar; ++j)
      p[j] = rangen.doub();
    logLike = getLogLike(p);
    chains.addNewStep( i, p, logLike );
    naccept[i] = 1;
  }
  chains.setSkipFirst();
}

double multiGauss::getLogLike(const paramSet& p) {
  unsigned int npar = getNParams();

  //Mean subtracted vector; mean is 0.5 in all dim
  for (unsigned int i = 0; i < npar; ++i)
    work1[i] = p[i] - 0.5;

  //First product
  double sum, *ptr;
  for (unsigned int i = 0; i < npar; ++i) {
    ptr = invCovMatrix + npar*i;
    sum = ptr[0]*work1[0];
    for (unsigned int j = 1; j < npar; ++j)
      sum += ptr[j]*work1[j];
    work2[i] = sum;
  }
  //Second
  sum = work1[0]*work2[0];
  for (unsigned int i = 1; i < npar; ++i)
    sum += work1[i]*work2[i];

  return -0.5*sum;
}

int main(int argc, char** argv) {

  unsigned int nwalkers, nsamples;
  std::string covfile, outfile;
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
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (nproc < 2) {
    if (rank == 0) {
      std::cerr << "Must run on multiple processes" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  while ( ( c = getopt_long(argc,argv,optstring,long_options,
                            &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      if (rank == 0) {
	std::cerr << "NAME" << std::endl;
	std::cerr << "\tmulti_gauss_example -- draws samples from a "
		  << " 10-dimensional " 
		  << std::endl;
	std::cerr << "\t\tmulti-variate Gaussian" << std::endl;
	std::cerr << std::endl;
	std::cerr << "SYNOPSIS" << std::endl;
	std::cerr << "\tmulti_gauss_example covfile nwalkers nsamples outfile"
		  << std::endl;
	std::cerr << "DESCRIPTION:" << std::endl;
	std::cerr << "\tDraws samples from a multi-variate Gaussian using"
		  << std::endl;
	std::cerr << "\tan affine-invariant MCMC code, using nwalkers walkers"
		  << std::endl;
	std::cerr << "\tand generating (approximately) nsamples samples.  The"
		  << std::endl;
	std::cerr << "\tresults are written to outfile." << std::endl;
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
  if ( optind >= argc-3 ) {
    if (rank == 0) {
      std::cerr << "Required arguments missing" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  covfile  = std::string(argv[optind]);
  nwalkers = atoi(argv[optind+1]);
  nsamples = atoi(argv[optind+2]);
  outfile  = std::string(argv[optind+3]);

  if (nwalkers == 0 || nsamples == 0) {
    MPI_Finalize();
    return 0;
  }

  //Hardwired cov matrix
  try {
    multiGauss *mg = new multiGauss(nwalkers,10,nsamples,
				    covfile);
    if (verbose) mg->setVerbose();
    
    if (verbose && rank == 0)
      std::cout << "Entering main loop" << std::endl;
    mg->initChains();
    mg->doSteps(mg->getNSteps());
    
    if (rank == 0) {
      std::vector<double> accept;
      mg->getAcceptanceFrac(accept);
      double mnacc = accept[0];
      for (unsigned int i = 1; i < nwalkers; ++i)
	mnacc += accept[i];
      std::cout << "Mean acceptance: " << mnacc / static_cast<double>(nwalkers)
		<< std::endl;
      
      std::vector<double> acor;
      bool succ = mg->getAcor(acor);
      if (succ) std::cout << "Autocorrelation length: " << acor[0] 
			  << " " << acor[1] << std::endl; 
      else std::cout << "Failed to compute autocorrelation" << std::endl;
  
      //Write
      mg->writeToFile(outfile);
    }
    delete mg;

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
