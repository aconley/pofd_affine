#include<string>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<limits>

#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/affineExcept.h"
#include "../include/affineEnsemble.h"
#include "../include/paramSet.h"
#include "../include/hdf5utils.h"

/*!
  \brief Rosenbrock density
*/
class rosenbrockDensity : public affineEnsemble {
private :
  double a1, a2;
public:
  rosenbrockDensity(double, double, unsigned int, unsigned int,
		    unsigned int, double, unsigned int, bool);
  void initChains();
  void generateInitialPosition(const paramSet&);
  double getLogLike(const paramSet&, bool&);
};

rosenbrockDensity::rosenbrockDensity(double A1, double A2, 
				     unsigned int NWALKERS, 
				     unsigned int NSAMPLES,
				     unsigned int INIT_STEPS=0,
				     double INIT_TEMP=2.0,
				     unsigned int MIN_BURN=50,
				     bool FIXED_BURN=false) :
  affineEnsemble(NWALKERS, 2, NSAMPLES, INIT_STEPS, INIT_TEMP,
		 MIN_BURN, FIXED_BURN), a1(A1), a2(A2) { }
  
void rosenbrockDensity::initChains() {
  is_init = true;
  paramSet p(2);
  p[0] = a1 + 0.25 * rangen.gauss();
  p[1] = a2 + 0.15 * rangen.gauss();
  has_initStep = true;
  initStep = p;
  generateInitialPosition(p);

}

void rosenbrockDensity::generateInitialPosition(const paramSet& p) {
  // Uniform initial values
  if (rank != 0) return;

  chains.clear();
  chains.addChunk(1);

  paramSet p2(2);
  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    p2[0] = 2.0*rangen.doub() - 1.0 + p[0];
    p2[1] = rangen.doub() - 0.5 + p[1];
    chains.addNewStep(i, p2, -std::numeric_limits<double>::infinity());
    naccept[i] = 0;
  }
  chains.setSkipFirst();
}

double rosenbrockDensity::getLogLike(const paramSet& p, bool& rej) {
  double val1, val2;
  val1 = (p[1] -p[0] * p[0]);
  val2 = 1 - p[0];
  if (std::isnan(val1) || std::isinf(val1) ||
      std::isnan(val2) || std::isinf(val2)) rej = true; 
  else rej=false;
  return - (a1 * val1 * val1 + val2 * val2) / a2;
}

int main(int argc, char** argv) {

  unsigned int nwalkers, nsamples, min_burn, init_steps;
  std::string outfile;
  double a1, a2, init_temp;
  bool verbose, fixed_burn, write_as_hdf;

  verbose = false;
  min_burn = 20;
  init_steps = 20;
  init_temp = 2.0;
  fixed_burn = false;
  write_as_hdf = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help",no_argument,0,'h'},
    {"fixedburn", no_argument, 0, 'f'},
    {"initsteps", required_argument, 0, 'i'},
    {"inittemp", required_argument, 0, 'I'},
    {"minburn", required_argument, 0, 'm'},
    {"verbose",no_argument,0,'v'},
    {"Version",no_argument,0,'V'},
    {0,0,0,0}
  };
  char optstring[] = "hfi:I:m:vV";

  int rank, nproc;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
	std::cerr << "\trosenbrock_example nwalkers nsamples a1 a2 outfile"
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
        std::cerr << "\t-f, --fixedburn" << std::endl;
        std::cerr << "\t\tUsed a fixed burn in length before starting main"
                  << " sample," << std::endl;
        std::cerr << "\t\trather than using the autocorrelation."
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
    if (rank == 0) {
      std::cerr << "Required arguments missing" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  nwalkers = atoi(argv[optind]);
  nsamples = atoi(argv[optind + 1]);
  a1       = atof(argv[optind + 2]);
  a2       = atof(argv[optind + 3]);
  outfile  = std::string(argv[optind + 4]);

  if (nwalkers == 0 || nsamples == 0) {
    MPI_Finalize();
    return 0;
  }

  if (nproc < 2) {
    if (rank == 0)
      std::cerr << "Must run on multiple processes" << std::endl;
    MPI::Finalize();
    return 1;
  }

  //Hardwired cov matrix
  try {
    rosenbrockDensity rd(a1, a2, nwalkers, nsamples, init_steps, init_temp,
			 min_burn, fixed_burn);
    rd.setParamName(0, "a1");
    rd.setParamName(1, "a2");
    if (verbose) {
      rd.setVerbose();
      if (rank == 0)
	std::cout << rd << std::endl;
    }

    rd.sample(); //Also does initialization
    
    if (rank == 0) {
      rd.printStatistics();

      std::vector<float> accept;
      rd.getAcceptanceFrac(accept);
      double mnacc = accept[0];
      for (unsigned int i = 1; i < nwalkers; ++i)
	mnacc += accept[i];
      std::cout << "Mean acceptance: " << mnacc / static_cast<double>(nwalkers)
		<< std::endl;

      std::vector<float> acor;
      bool succ = rd.getAcor(acor);
      if (succ) std::cout << "Autocorrelation length: " << acor[0] 
			  << " " << acor[1] << std::endl; 
      else std::cout << "Failed to compute autocorrelation" << std::endl;
  
      hdf5utils::outfiletype ftype;
      ftype = hdf5utils::getOutputFileType(outfile);
      switch(ftype) {
      case hdf5utils::TXT:
	rd.writeToFile(outfile);
	break;
      case hdf5utils::FITS:
	throw affineExcept("rosenbrock_example", "main",
			   "No support for FITS output");
	break;
      case hdf5utils::UNKNOWN:
      case hdf5utils::HDF5:
	rd.writeToHDF5(outfile);
	break;
      }
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
