//Multi-dimensional Gaussian example

#include<string>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>

#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/affineExcept.h"
#include "../include/affineEnsemble.h"
#include "../include/paramSet.h"

//The covariance matrix is set up as
// 1) A correlation matrix with the band structure [-0.1, 0.2, 1, 0.2, 0.1]
// 2) Sigmas of [1.0, 0.9, 0.8, 0.7, 0.6]

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
  void getStats(std::vector<float>&, std::vector<float>&) const;
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
    chains.addNewStep(i, p, logLike);
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

void multiGauss::getStats(std::vector<float>& mn, 
			  std::vector<float>& var) const {
  unsigned int npar = getNParams();
  mn.resize(npar);
  var.resize(npar);
  float cmn, cvar, lowlim, uplim;
  for (unsigned int i = 0; i < npar; ++i) {
    chains.getParamStats(i, cmn, cvar, lowlim, uplim);
    mn[i] = cmn;
    var[i] = cvar;
  }
}

///////////////////////////////////

int main(int argc, char** argv) {

  const unsigned int ndim = 5;
  unsigned int nwalkers, nsamples, nburn;
  std::string invcovfile, outfile;
  bool verbose;

  verbose = false;
  nburn = 20;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help",no_argument,0,'h'},
    {"nburn", required_argument, 0, 'n'},
    {"verbose",no_argument,0,'v'},
    {"Version",no_argument,0,'V'},
    {0,0,0,0}
  };
  char optstring[] = "hn:vV";

  int rank, nproc;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  while ( ( c = getopt_long(argc,argv,optstring,long_options,
                            &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      if (rank == 0) {
	std::cerr << "NAME" << std::endl;
	std::cerr << "\tmulti_gauss_example -- draws samples from a "
		  << " 5-dimensional " 
		  << std::endl;
	std::cerr << "\t\tmulti-variate Gaussian" << std::endl;
	std::cerr << std::endl;
	std::cerr << "SYNOPSIS" << std::endl;
	std::cerr << "\tmulti_gauss_example nwalkers nsamples outfile"
		  << std::endl;
	std::cerr << "DESCRIPTION:" << std::endl;
	std::cerr << "\tDraws samples from a multi-variate Gaussian using"
		  << std::endl;
	std::cerr << "\tan affine-invariant MCMC code, using nwalkers walkers"
		  << std::endl;
	std::cerr << "\tand generating (approximately) nsamples samples.  The"
		  << std::endl;
	std::cerr << "\tresults are written to outfile.  The input covariance" 
		  << std::endl;
	std::cerr << "\tmatrix is fixed, and the mean is 0.5 in all " << 
	  "dimensions." << std::endl;
	std::cerr << std::endl;
	std::cerr << "OPTIONS" << std::endl;
	std::cerr << "\t-h, --help" << std::endl;
	std::cerr << "\t\tOutput this help." << std::endl;
	std::cerr << "\t-n, --nburn NBURN" << std::endl;
	std::cerr << "\t\tNumber of burn-in steps to do per walker (def: 20)"
		  << std::endl;
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
    case 'n':
      nburn = atoi(optarg);
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

  if (nproc < 2) {
    if (rank == 0) {
      std::cerr << "Must run on multiple processes" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  if ( optind >= argc-2 ) {
    if (rank == 0) {
      std::cerr << "Required arguments missing" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  nwalkers = atoi(argv[optind]);
  nsamples = atoi(argv[optind+1]);
  outfile  = std::string(argv[optind+2]);

  if (nwalkers == 0 || nsamples == 0) {
    MPI_Finalize();
    return 0;
  }

  invcovfile = "exampledata/test_invcovmat.txt";

  //Hardwired cov matrix
  try {
    multiGauss *mg = new multiGauss(nwalkers, ndim, nsamples,
				    invcovfile);
    if (verbose) mg->setVerbose();
    
    if (verbose && rank == 0)
      std::cout << "Entering main loop" << std::endl;
    mg->initChains();
    mg->doSteps(mg->getNSteps(), nburn);
    
    if (rank == 0) {
      std::vector<float> mn, var;
      mg->getStats(mn, var);
      std::cout << "Results" << std::endl;
      for (unsigned int i = 0; i < mn.size(); ++i)
	std::cout << "  Param " << i << " " << mn[i] << " +- "
		  << sqrt(var[i]) << std::endl;

      std::vector<float> accept;
      mg->getAcceptanceFrac(accept);
      double mnacc = accept[0];
      for (unsigned int i = 1; i < nwalkers; ++i)
	mnacc += accept[i];
      std::cout << "Mean acceptance: " << mnacc / static_cast<double>(nwalkers)
		<< std::endl;
      
      std::vector<float> acor;
      bool succ = mg->getAcor(acor);
      if (succ) {
	std::cout << "Autocorrelation length: " << acor[0];
	for (unsigned int i = 1; i < acor.size(); ++i)
	  std::cout << " " << acor[i];
	std::cout << std::endl;
      } else std::cout << "Failed to compute autocorrelation" << std::endl;
  
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
