//Single-dimensional Gaussian example (very basic test)

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

/*!
  \brief Single-dimensional Gaussian,
*/
class singleGauss : public affineEnsemble {
private:
  double mean; //!< Mean value
  double var; //!< Variance
  double gfac; //!< -1/(2*var)
public:
  singleGauss(double, double, unsigned int, unsigned int,
	      unsigned int, bool);
  virtual ~singleGauss();

  void initChains();
  double getLogLike(const paramSet&);
  void getStats(float&, float&) const;
 
};

singleGauss::singleGauss(double MN, double VAR,
			 unsigned int NWALKERS, unsigned int NSAMPLES, 
			 unsigned int MIN_BURN, bool FIXED_BURN=false) :
  affineEnsemble(NWALKERS, 1, NSAMPLES, 0, MIN_BURN, FIXED_BURN) {
    mean = MN;
    var  = VAR;
    gfac = - 1.0/(2.0*var);
}

singleGauss::~singleGauss() {}

void singleGauss::initChains() {
  //Just generate random positions from mean-sigma*5 to mean+sigma*5
  // initial vectors

  if (rank != 0) return;

  //For master
  chains.clear();
  chains.addChunk(1);
  
  paramSet p(1);
  double logLike;
  double rng = sqrt(var) * 10.0;
  double genmn = mean - 0.5*rng;

  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    p[0] = rng * rangen.doub() + genmn;
    logLike = getLogLike(p);
    chains.addNewStep(i, p, logLike);
    naccept[i] = 1;
  }
  chains.setSkipFirst();
}

double singleGauss::getLogLike(const paramSet& p) {
  double val = p[0] - mean;
  return gfac * val * val;
}

void singleGauss::getStats(float& mn, float& var) const {
  float lowlim, uplim;
  chains.getParamStats(0, mn, var, lowlim, uplim);
}

//////////////////////////////////

int main(int argc, char** argv) {

  unsigned int nwalkers, nsamples, nburn;
  double mean, sigma;
  std::string outfile;
  bool verbose, fixed_burn;

  nburn = 20;
  verbose = false;
  fixed_burn = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help",no_argument,0,'h'},
    {"fixedburn", no_argument, 0, 'f'},
    {"nburn", required_argument, 0, 'n'},
    {"verbose",no_argument,0,'v'},
    {"Version",no_argument,0,'V'},
    {0,0,0,0}
  };
  char optstring[] = "hfn:vV";

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
	std::cerr << "\tsingle_gauss_example -- draws samples from a "
		  << " 1-dimensional Gaussian." 
		  << std::endl;
	std::cerr << std::endl;
	std::cerr << "SYNOPSIS" << std::endl;
	std::cerr << "\tsingle_gauss_example nwalkers nsamples mean sigma "
		  << "outfile" << std::endl;
	std::cerr << "DESCRIPTION:" << std::endl;
	std::cerr << "\tDraws samples from a Gaussian using"
		  << std::endl;
	std::cerr << "\tan affine-invariant MCMC code, using nwalkers walkers"
		  << std::endl;
	std::cerr << "\tand generating (approximately) nsamples samples.  The"
		  << std::endl;
	std::cerr << "\tresults are written to outfile." << std::endl;
	std::cerr << "OPTIONS" << std::endl;
	std::cerr << "\t-h, --help" << std::endl;
	std::cerr << "\t\tOutput this help." << std::endl;
        std::cerr << "\t-f, --fixedburn" << std::endl;
        std::cerr << "\t\tUsed a fixed burn in length before starting main"
                  << " sample," << std::endl;
        std::cerr << "\t\trather than using the autocorrelation."
		  << std::endl;
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
    case 'f':
      fixed_burn = true;
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

  if (optind >= argc - 4) {
    if (rank == 0) {
      std::cerr << "Required arguments missing" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }
  nwalkers = atoi(argv[optind]);
  nsamples = atoi(argv[optind+1]);
  mean     = atof(argv[optind+2]);
  sigma    = atof(argv[optind+3]);
  outfile  = std::string(argv[optind+4]);

  if (nwalkers == 0 || nsamples == 0) {
    MPI_Finalize();
    return 0;
  }

  //Hardwired cov matrix
  try {
    singleGauss *sg = new singleGauss(mean, sigma*sigma,
				      nwalkers, nsamples, nburn,
				      fixed_burn);
    
    if (verbose) sg->setVerbose();
    
    if (verbose && rank == 0)
      std::cout << "Entering main loop" << std::endl;
    sg->sample();

    if (rank == 0) {
      float mn, var;
      sg->getStats(mn, var);
      std::cout << "Mean: " << mn << " sigma: " << sqrt(var) << std::endl;

      std::vector<float> accept;
      sg->getAcceptanceFrac(accept);
      std::cout << "Acceptance fractions:";
      for (unsigned int i = 0; i < nwalkers; ++i)
	std::cout << " " << accept[i];
      std::cout << std::endl;
      
      sg->writeToFile(outfile);
    }

    delete sg;

  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    int jnk;
    if (rank == 0)
      for (unsigned int i = 1; i < nproc; ++i)
	MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
    MPI_Finalize();
    return 1;
  }
  MPI_Finalize();
  return 0;
}
