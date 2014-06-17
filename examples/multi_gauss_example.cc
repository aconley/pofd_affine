//Multi-dimensional Gaussian example

#include<string>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<limits>

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
  multiGauss(const std::string&, unsigned int, unsigned int, unsigned int, 
	     unsigned int, double, unsigned int, bool);
  virtual ~multiGauss();

  void initChains();
  void generateInitialPosition(const paramSet&);
  double getLogLike(const paramSet&, bool& params_rejected);
  bool areParamsValid(const paramSet& p) const;
  void getStats(std::vector<float>&, std::vector<float>&) const;
};

multiGauss::multiGauss(const std::string& filename, unsigned int NWALKERS, 
		       unsigned int NPARAMS, unsigned int NSAMPLES, 
		       unsigned int INIT_STEPS=0, double INIT_TEMP=2.0,
		       unsigned int MIN_BURN=50, bool FIXED_BURN=false) :
  affineEnsemble(NWALKERS, NPARAMS, NSAMPLES, INIT_STEPS, INIT_TEMP, MIN_BURN, 
		 FIXED_BURN) {
  if (NPARAMS > 0) {
    invCovMatrix = new double[NPARAMS * NPARAMS];
    work1 = new double[NPARAMS];
    work2 = new double[NPARAMS];
  } else {
    invCovMatrix = NULL;
    work1 = NULL;
    work2 = NULL;
  }
  
  //Read in file
  std::ifstream ifs(filename.c_str());
  if (!ifs)
    throw affineExcept("multiGauss", "multiGauss", "Couldn't open input file");
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
  // This doesn't do much of anything except call generateInitialPosition 
  // the around the mean
  is_init = true;
  if (rank != 0) return;

  unsigned int npar = getNParams();
  paramSet p(npar);
  for (unsigned int i = 0; i < npar; ++i)
    p[i] = 0.55; //Slightly off
  has_initStep = true;
  initStep = p;
  generateInitialPosition(p);
}

// This isn't really necessary for this simple distribution, but
// as an example of how to use it, reject all points outside [-10, 10]
bool multiGauss::areParamsValid(const paramSet& p) const {
  unsigned int npar = getNParams();
  for (unsigned int i = 0; i < npar; ++i)
    if (fabs(p[i]) > 10) return false;
  return true;
}

// Note we use areParamsValid to reject invalid starting positions
void multiGauss::generateInitialPosition(const paramSet& p) {

  const unsigned int maxiters = 20;

  if (rank != 0) return;

  chains.clear();
  chains.addChunk(1);
    
  unsigned int npar = p.getNParams();
  if (npar != getNParams())
    throw affineExcept("multiGauss", "generateInitialPosition",
		       "Wrong number of params in p");

  if (!areParamsValid(p))
    throw affineExcept("multiGauss", "generateInitialPosition",
		       "Initial position seed is invalid");

  paramSet p2(npar);
  unsigned int nwalk = getNWalkers();
  bool is_valid;
  unsigned int iter;
  for (unsigned int i = 0; i < nwalk; ++i) {
    iter = 0;
    is_valid = false;
    while (!is_valid) {
      if (iter > maxiters)
	throw affineExcept("multiGauss", "generateInitialPosition",
			   "Unable to generate initial position");
      for (unsigned int j = 0; j < npar; ++j)
	p2[j] = rangen.doub() - 0.5 + p[j];
      is_valid = areParamsValid(p2);
      ++iter;
    }
    chains.addNewStep(i, p2, -std::numeric_limits<double>::infinity());
    naccept[i] = 0;
  }
  chains.setSkipFirst();
}


double multiGauss::getLogLike(const paramSet& p, bool& params_rejected) {
  unsigned int npar = getNParams();
  if (npar == 0) {
    params_rejected = true;
    return 0;
  } else params_rejected = false;

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
  unsigned int nwalkers, nsamples, min_burn, init_steps;
  std::string invcovfile, outfile;
  bool verbose, fixed_burn, write_as_hdf;
  double init_temp;

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
    {"hdf", no_argument, 0, 'H'},
    {"inittemp", required_argument, 0, 'I'},
    {"initsteps", required_argument, 0, 'i'},
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

  if (optind >= argc - 2) {
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

  if (nproc < 2) {
    if (rank == 0) {
      std::cerr << "Must run on multiple processes" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  invcovfile = "exampledata/test_invcovmat.txt";

  //Hardwired cov matrix
  try {
    multiGauss mg(invcovfile, nwalkers, ndim, nsamples, init_steps,
		  init_temp, min_burn, fixed_burn);
    if (verbose) {
      mg.setVerbose();
      if (rank == 0)
	std::cout << mg << std::endl;
    }

    mg.setParamName(0, "mn[0]");
    mg.setParamName(1, "mn[1]");
    mg.setParamName(2, "mn[2]");
    mg.setParamName(3, "mn[3]");
    mg.setParamName(4, "mn[4]");

    mg.sample(); //Also initializes
    
    if (rank == 0) {
      mg.printStatistics();

      std::vector<float> accept;
      mg.getAcceptanceFrac(accept);
      double mnacc = accept[0];
      for (unsigned int i = 1; i < nwalkers; ++i)
	mnacc += accept[i];
      std::cout << "Mean acceptance: " << mnacc / static_cast<double>(nwalkers)
		<< std::endl;
      
      std::vector<float> acor;
      bool succ = mg.getAcor(acor);
      if (succ) {
	std::cout << "Autocorrelation length: " << acor[0];
	for (unsigned int i = 1; i < acor.size(); ++i)
	  std::cout << " " << acor[i];
	std::cout << std::endl;
      } else std::cout << "Failed to compute autocorrelation" << std::endl;
  
      //Write
      if (write_as_hdf)
	mg.writeToHDF5(outfile);
      else
	mg.writeToFile(outfile);

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
