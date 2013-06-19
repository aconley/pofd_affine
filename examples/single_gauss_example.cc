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
  singleGauss(double, double, unsigned int, unsigned int, unsigned int,
	      double, unsigned int, bool);
  virtual ~singleGauss();

  void initChains();
  void generateInitialPosition(const paramSet&);
  double getLogLike(const paramSet&);
  void getStats(float&, float&) const;
 
};

singleGauss::singleGauss(double MN, double SIGMA,
			 unsigned int NWALKERS, unsigned int NSAMPLES, 
			 unsigned int INIT_STEPS=0, double INIT_TEMP=2.0,
			 unsigned int MIN_BURN=50, bool FIXED_BURN=false) :
  affineEnsemble(NWALKERS, 1, NSAMPLES, INIT_STEPS, INIT_TEMP, MIN_BURN, 
		 FIXED_BURN) {
    mean = MN;
    var  = SIGMA * SIGMA;
    gfac = - 0.5 / var;
}

singleGauss::~singleGauss() {}

void singleGauss::initChains() {
  // For this example, this doesn't do much of anything except
  // call generateInitialPosition around mean
  is_init = true;

  if (rank != 0) return;

  paramSet p(1);
  //Don't quite start in the right place
  p[0] = mean + 2 * sqrt(var) * rangen.gauss(); 
  generateInitialPosition(p);
}

void singleGauss::generateInitialPosition(const paramSet& p) {
  if (rank != 0) return;

  chains.clear();
  chains.addChunk(1);

  //Just generate random positions from mean-sigma*5 to mean+sigma*5
  double rng = sqrt(var) * 10.0;
  double genmn = p[0] - 0.5*rng;
  
  unsigned int nwalk = getNWalkers();
  paramSet p2(1);
  for (unsigned int i = 0; i < nwalk; ++i) {
    p2[0] = rng * rangen.doub() + genmn;
    chains.addNewStep(i, p2, -std::numeric_limits<double>::infinity());
    naccept[i] = 0;
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

  unsigned int nwalkers, nsamples, min_burn, init_steps;
  double mean, sigma, init_temp;
  std::string outfile;
  bool verbose, fixed_burn, write_as_hdf;

  min_burn = 20;
  init_steps = 20;
  init_temp = 2.0;
  verbose = false;
  fixed_burn = false;
  write_as_hdf = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help",no_argument,0,'h'},
    {"fixedburn", no_argument, 0, 'f'},
    {"hdf", no_argument, 0, 'H'},
    {"initsteps", required_argument, 0, 'i'},
    {"inittemp", required_argument, 0, 'I'},
    {"minburn", required_argument, 0, 'm'},
    {"verbose",no_argument,0,'v'},
    {"Version",no_argument,0,'V'},
    {0,0,0,0}
  };
  char optstring[] = "hfHi:I:m:vV";

  int rank, nproc;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  while ((c = getopt_long(argc, argv, optstring, long_options,
			  &option_index)) != -1) 
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
	std::cerr << "\t-H, --hdf" << std::endl;
	std::cerr << "\t\tWrite as HDF5 file rather than text file."
		  << std::endl;
	std::cerr << "\t-i, --initsteps STEPS" << std::endl;
	std::cerr << "\t\tTake this many initial steps per walker, then "
		  << "recenter" << std::endl;
	std::cerr << "\t\taround the best one before starting burn in"
		  << " (def: 20)" << std::endl;
	std::cerr << "\t-I, --inittemp VALUE" << std::endl;
	std::cerr << "\t\tTemperature used during initial steps (def: 2.0)"
		  << std::endl;
	std::cerr << "\t-m, --minburn NSTEPS" << std::endl;
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
    case 'H':
      write_as_hdf = true;
      break;
    case 'i':
      init_steps = atoi(optarg);
      break;
    case 'I':
      init_temp = atof(optarg);
      break;
    case 'm':
      min_burn = atoi(optarg);
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

  if (optind >= argc - 4) {
    if (rank == 0)
      std::cerr << "Required arguments missing" << std::endl;
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

  if (nproc < 2) {
    if (rank == 0) 
      std::cerr << "Must run on multiple processes" << std::endl;
    MPI_Finalize();
    return 1;
  }

  //Hardwired cov matrix
  try {
    singleGauss sg(mean, sigma, nwalkers, nsamples, 
		   init_steps, init_temp, min_burn, fixed_burn);
    if (verbose) {
      sg.setVerbose();
      if (rank == 0)
	std::cout << sg << std::endl;
    }
    
    // Do main loop
    sg.sample(); //Also initializes

    if (rank == 0) {
      float mn, var;
      sg.getStats(mn, var);
      std::cout << "Mean: " << mn << " sigma: " << sqrt(var) << std::endl;

      std::vector<float> accept;
      sg.getAcceptanceFrac(accept);
      std::cout << "Acceptance fractions:";
      for (unsigned int i = 0; i < nwalkers; ++i)
	std::cout << " " << accept[i];
      std::cout << std::endl;
      
      if (write_as_hdf)
	sg.writeToHDF5(outfile);
      else
	sg.writeToFile(outfile);
    }

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
