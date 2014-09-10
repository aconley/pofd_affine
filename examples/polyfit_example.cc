#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<ios>

#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/affineExcept.h"
#include "../include/affineEnsemble.h"
#include "../include/paramSet.h"
#include "../include/hdf5utils.h"

/*!
  \brief Fit a polynomial to a fixed set of input data

  Has one bonus parameter, the rms of the data points around the
  model.
*/
class polyFit : public affineEnsemble {
private :
  unsigned int ndata;
  double *x;
  double *y;
  double *iyvar;

  double rms; //!< RMS of data around model, computed by getLogLike
public:
  polyFit(const std::string&, unsigned int, unsigned int, 
	  unsigned int, unsigned int, double, unsigned int, bool);
  virtual ~polyFit();
  void initChains();
  void generateInitialPosition(const paramSet&);
  void readFile(const std::string&);
  double getLogLike(const paramSet&, bool&);

  void fillBonusParams(paramSet&, bool);
};

polyFit::polyFit(const std::string& datafile, unsigned int NWALKERS, 
		 unsigned int NPARAMS, unsigned int NSAMPLES,
		 unsigned INIT_STEPS=0, double INIT_TEMP=2.0,
		 unsigned int MIN_BURN=50, bool FIXED_BURN=false) :
  affineEnsemble(NWALKERS, NPARAMS, NSAMPLES, INIT_STEPS, INIT_TEMP,
		 MIN_BURN, FIXED_BURN) {
  ndata = 0;
  x = NULL;
  y = NULL;
  iyvar = NULL;
  readFile(datafile);
}

polyFit::~polyFit() {
  if (x != NULL) delete[] x;
  if (y != NULL) delete[] y;
  if (iyvar != NULL) delete[] iyvar;
}
  
void polyFit::readFile(const std::string& datafile) {
  std::ifstream ifs;
  ifs.open(datafile.c_str());
  if (!ifs) {
    std::stringstream errstr;
    errstr << "Couldn't open input file: " << datafile;
    throw affineExcept("polyFit", "readFile", errstr.str());
  }

  unsigned int nentries = 0;
  std::string line;
  while (std::getline(ifs,line)) ++nentries;
  ifs.clear(); //Clear error flag resulting from last failed read
  ifs.seekg(0,std::ios_base::beg); //Rewind

  //Reallocate space
  if (nentries != ndata) {
    if (x != NULL) delete[] x;
    if (y != NULL) delete[] y;
    if (iyvar != NULL) delete[] iyvar;
    if (nentries > 0) {
      ndata = nentries;
      x = new double[ndata];
      y = new double[ndata];
      iyvar = new double[ndata];
    } else {
      ndata = 0;
      x = y = iyvar = NULL;
    }
  }

  //Read
  for (unsigned int i = 0; i < ndata; ++i)
    ifs >> x[i] >> y[i] >> iyvar[i];
  for (unsigned int i = 0; i < ndata; ++i)
    iyvar[i] = 1.0 / (iyvar[i] * iyvar[i]);
  ifs.close();
}


void polyFit::initChains() {
  is_init = true;
  if (rank != 0) return;

  unsigned int npar = getNParams();
  paramSet p(npar);
  for (unsigned int i = 0; i < npar-1; ++i)
    p[i] = 0.0;
  p[npar - 1] = std::numeric_limits<double>::quiet_NaN();
  has_initStep = true;
  initStep = p;
  generateInitialPosition(p);
}

void polyFit::generateInitialPosition(const paramSet& p) {
  if (rank != 0) return;

  chains.clear();
  chains.addChunk(1);
    
  unsigned int npar = p.getNParams();
  if (npar != getNParams())
    throw affineExcept("polyFit", "generateInitialPosition",
		       "Wrong number of params in p");

  paramSet p2(npar);
  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    for (unsigned int j = 0; j < npar-1; ++j)
      p2[j] = rangen.doub() - 0.5 + p[j];
    p2[npar-1] = std::numeric_limits<double>::quiet_NaN(); // RMS
    chains.addNewStep(i, p2, -std::numeric_limits<double>::infinity());
    naccept[i] = 0;
  }
  chains.setSkipFirst();
}

double polyFit::getLogLike(const paramSet& p, bool& rej) {

  //Evaluate the polynomial at each data point
  //Order starts with smallest order first

  unsigned int npar = getNParams();
  if (npar > 0) rej=false;
  else {
    rej = true;
    return 0;
  }

  double val, xval, deltasq, logLike;
  logLike = 0.0;
  rms = 0.0;
  for (unsigned int i = 0; i < ndata; ++i) {
    xval = x[i];
    val = p[npar-2]; // highest poly coeff
    for (int j=static_cast<int>(npar)-3; j >= 0; --j)
      val = val * xval + p[j];
    deltasq = y[i] - val;
    deltasq *= deltasq;
    logLike -= 0.5 * deltasq * iyvar[i];
    rms += deltasq;
  }

  rms = sqrt(rms / static_cast<double>(ndata));
  return logLike;
}

void polyFit::fillBonusParams(paramSet& p, bool rej) {
  if (rej)
    p[getNParams()-1] = std::numeric_limits<double>::quiet_NaN();
  else
    p[getNParams()-1] = rms;
}

/////////////////////////////////////

int main(int argc, char** argv) {

  const unsigned int nterms = 4;

  bool verbose, fixed_burn;
  unsigned int nwalkers, nsamples, min_burn, init_steps;
  double scalefac, init_temp;
  std::string outfile;

  //Defaults
  verbose     = false;
  fixed_burn  = false;
  init_steps  = 30;
  init_temp   = 2.0;
  min_burn    = 50;
  scalefac    = 2;
  
  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"fixedburn", no_argument, 0, 'f'},
    {"initsteps", required_argument, 0, 'i'},
    {"inittemp", required_argument, 0, 'I'},
    {"minburn", required_argument, 0, 'm'},
    {"scalefac", required_argument,0, 's'},
    {"verbose", no_argument, 0, 'v'},
    {"Version", no_argument, 0, 'V'},
    {0, 0, 0, 0}
  };
  char optstring[] = "hfi:I:m:s:vV";

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
	std::cerr << "\tpolyfit_example -- fits a polynomial to example" 
		  << std::endl;
	std::cerr << "\t\tdata set." << std::endl;
	std::cerr << std::endl;
	std::cerr << "SYNOPSIS" << std::endl;
	std::cerr << "\tpolyfit_example nwalkers nsamples outfile"
		  << std::endl;
	std::cerr << "DESCRIPTION:" << std::endl;
	std::cerr << "\tFits the polynomail using an affine-invariant " 
		  << "MCMC code," << std::endl;
	std::cerr << "\twith nwalkers walkers and generating (approximately)"
		  << " nsamples samples." << std::endl;
	std::cerr << "\tThe RMS of the data around the model is also computed"
		  << std::endl;
	std::cerr << "\tfor each step in the chain.  The chain is output to" 
		  << std::endl;
	std::cerr << "\toutfile in HDF5 format." << std::endl;
	std::cerr << std::endl;
	std::cerr << "\tThe input polynomial is: " << std::endl
		  << "\t\t0.5 + 0.04 x - 0.15 x^2 + 0.3 x^3"
		  << std::endl;
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
		  << " (def: 30)" << std::endl;
	std::cerr << "\t-I, --inittemp VALUE" << std::endl;
	std::cerr << "\t\tTemperature used during initial steps (def: 2.0)"
		  << std::endl;
	std::cerr << "\t-m, --minburn MINBURN" << std::endl;
	std::cerr << "\t\tMinimum number of burn-in steps to do per "
		  << "walker (def: 50)" << std::endl;
	std::cerr << "\t-s, --scalefac SCALEFAC" << std::endl;
	std::cerr << "\t\tScaling factor in MCMC proposal density (def: 2)"
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
    case 'm' :
      min_burn = atoi(optarg);
      break;
    case 's' :
      scalefac = atof(optarg);
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
    if (rank == 0)
      std::cerr << "Required arguments missing" << std::endl;
    MPI_Finalize();
    return 1;
  }

  nwalkers = atoi(argv[optind]);
  nsamples = atoi(argv[optind + 1]);
  outfile = std::string(argv[optind + 2]);
  if (nwalkers == 0 || nsamples == 0) {
    MPI_Finalize();
    return 0;
  }
  if (scalefac <= 0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) scaling factor " << scalefac
		<< std::endl;
    MPI_Finalize();
    return 1;
  }

  if (nproc < 2) {
    if (rank == 0) {
      std::cerr << "Must run on multiple processes" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  try {
    std::string infile="exampledata/polyexample.txt";
    polyFit ply(infile, nwalkers, nterms + 1, nsamples, init_steps,
		init_temp, min_burn, fixed_burn);
    ply.setScalefac(scalefac);

    ply.setParamName(0, "c[0]");
    ply.setParamName(1, "c[1]*x");
    ply.setParamName(2, "c[2]*x^2");
    ply.setParamName(3, "c[3]*x^3");
    ply.setParamName(4, "RMS");
    ply.setParamBonus(4);

    if (verbose) {
      ply.setVerbose();
      if (rank == 0)
	std::cout << ply << std::endl;
    }
    
    // Do Fit
    ply.sample(); //Also initializes
    
    if (rank == 0) {
      // Summarize results
      ply.printStatistics();

      std::vector<float> accept;
      ply.getAcceptanceFrac(accept);
      double mnacc = accept[0];
      for (unsigned int i = 1; i < nwalkers; ++i)
	mnacc += accept[i];
      std::cout << "Mean acceptance: " << mnacc / static_cast<double>(nwalkers)
		<< std::endl;

      // Write
      hdf5utils::outfiletype ftype;
      ftype = hdf5utils::getOutputFileType(outfile);
      switch(ftype) {
      case hdf5utils::TXT:
	ply.writeToFile(outfile);
	break;
      case hdf5utils::FITS:
	throw affineExcept("polyfit_example", "main",
			   "No support for FITS output");
	break;
      case hdf5utils::UNKNOWN:
      case hdf5utils::HDF5:
	ply.writeToHDF5(outfile);
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
