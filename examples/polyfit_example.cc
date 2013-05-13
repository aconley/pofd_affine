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

/*!
  \brief Fit a polynomial
*/
class polyFit : public affineEnsemble {
private :
  unsigned int ndata;
  double *x;
  double *y;
  double *iyunc;
public:
  polyFit(const std::string&, unsigned int, unsigned int, 
	  unsigned int, unsigned int, bool);
  virtual ~polyFit();
  void initChains();
  void readFile(const std::string&);
  double getLogLike(const paramSet&);
};

polyFit::polyFit(const std::string& datafile, unsigned int NWALKERS, 
		 unsigned int NPARAMS, unsigned int NSAMPLES,
		 unsigned int MIN_BURN, bool FIXED_BURN=false) :
  affineEnsemble(NWALKERS, NPARAMS, NSAMPLES, 0, MIN_BURN, FIXED_BURN) {
  ndata = 0;
  x = NULL;
  y = NULL;
  iyunc = NULL;
  readFile(datafile);
}

polyFit::~polyFit() {
  if (x != NULL) delete[] x;
  if (y != NULL) delete[] y;
  if (iyunc != NULL) delete[] iyunc;
}
  
void polyFit::readFile(const std::string& datafile) {
  std::ifstream ifs;
  ifs.open(datafile.c_str());
  if (!ifs) {
    std::stringstream errstr;
    errstr << "Couldn't open input file: " << datafile;
    throw affineExcept("polyFit","readFile",
		       errstr.str(),1);
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
    if (iyunc != NULL) delete[] iyunc;
    if (nentries > 0) {
      ndata = nentries;
      x = new double[ndata];
      y = new double[ndata];
      iyunc = new double[ndata];
    } else {
      ndata = 0;
      x = y = iyunc = NULL;
    }
  }

  //Read
  for (unsigned int i = 0; i < ndata; ++i)
    ifs >> x[i] >> y[i] >> iyunc[i];
  for (unsigned int i = 0; i < ndata; ++i)
    iyunc[i] = 1.0/iyunc[i];
  ifs.close();

}


void polyFit::initChains() {
  //Just generate random positions from -1 to 1 for
  // initial vectors

  if (rank != 0) return;

  chains.clear();
  chains.addChunk(1);

  paramSet p( getNParams() ); //nparams is the poly order + 1
  double logLike;

  unsigned int nwalk = getNWalkers();
  unsigned int npar = getNParams();
  for (unsigned int i = 0; i < nwalk; ++i) {
    for (unsigned int j = 0; j < npar; ++j)
      p[j] = 2*rangen.doub()-1;
    logLike = getLogLike(p);
    chains.addNewStep( i, p, logLike );
    naccept[i] = 1;
  }
  chains.setSkipFirst();
}

double polyFit::getLogLike(const paramSet& p) {

  //Evaluate the polynomial at each data point
  //Order starts with smallest order first
  unsigned int npar = getNParams();
  double val, xval, diffs;
  double logLike;

  logLike = 0.0;
  for (unsigned int i = 0; i < ndata; ++i) {
    xval = x[i];
    val = p[npar-1];
    for (int j=static_cast<int>(npar)-2; j >= 0; --j)
      val = val*xval+p[j];
    diffs = (y[i] - val)*iyunc[i];
    logLike -= diffs*diffs;
  }
  logLike *= 0.5;

  return logLike;
}

/////////////////////////////////////

int main(int argc, char** argv) {

  const unsigned int nterms = 4;
  const double coeffs[nterms] = {0.5, 0.04, -0.15, 0.3}; // Correct model

  unsigned int nwalkers, nsamples, min_burn;
  double scalefac;
  bool verbose, fixed_burn;

  //Defaults
  verbose     = false;
  fixed_burn  = false;
  min_burn    = 50;
  scalefac    = 2;
  
  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help",no_argument,0,'h'},
    {"fixedburn", no_argument, 0, 'f'},
    {"minburn",required_argument,0,'m'},
    {"scalefac",required_argument,0,'s'},
    {"verbose",no_argument,0,'v'},
    {"Version",no_argument,0,'V'},
    {0,0,0,0}
  };
  char optstring[] = "hfm:s:vV";

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
	std::cerr << "\tpolyfit_example -- fits a polynomial to example" 
		  << std::endl;
	std::cerr << "\t\tdata set." << std::endl;
	std::cerr << std::endl;
	std::cerr << "SYNOPSIS" << std::endl;
	std::cerr << "\tpolyfit_example nwalkers nsamples"
		  << std::endl;
	std::cerr << "DESCRIPTION:" << std::endl;
	std::cerr << "\tFits the polynomail using an affine-invariant " 
		  << "MCMC code," << std::endl;
	std::cerr << "\twith nwalkers walkers and generating (approximately)"
		  << " nsamples samples." << std::endl;
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
	std::cerr << "\t-m, --minburn MINBURN" << std::endl;
	std::cerr << "\t\tMinimum number of burn-in steps to do per "
		  << "walker (def: 50)" << std::endl;
	std::cerr << "\t-o, --outfile FILENAME\n" << std::endl;
	std::cerr << "\t\tWrite the resulting chain to this file."
		  << std::endl;
	std::cerr << "\t-s, --scalefac SCALEFAC\n" << std::endl;
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

  if (nproc < 2) {
    if (rank == 0) {
      std::cerr << "Must run on multiple processes" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  if ( optind >= argc-1 ) {
    if (rank == 0) {
      std::cerr << "Required arguments missing" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  nwalkers = atoi(argv[optind]);
  nsamples = atoi(argv[optind + 1]);

  if (nwalkers == 0) return 0;
  if (nsamples == 0) return 0;
  if (scalefac <= 0) {
    std::cerr << "Invalid (non-positive) scaling factor " << scalefac
	      << std::endl;
    return 1;
  }

  try {
    std::string infile="exampledata/polyexample.txt";
    polyFit ply(infile, nwalkers, nterms, nsamples, min_burn, fixed_burn);
    ply.setScalefac(scalefac);

    ply.setParamName(0,"c[0]");
    ply.setParamName(1,"c[1]*x");
    ply.setParamName(2,"c[2]*x^2");
    ply.setParamName(3,"c[3]*x^3");

    if (verbose) ply.setVerbose();
    
    // Do Fit
    if (verbose && rank == 0)
      std::cout << "Entering main loop" << std::endl;
    ply.sample();
    
    if (rank == 0) {
      // Summarize results
      float mn, var, low, up;
      float conflevel = 0.683;
      for (unsigned int i = 0; i < 4; ++i) {
	ply.getParamStats(i, mn, var, low, up, conflevel);
	std::cout << "Parameter: " << ply.getParamName(i) << " Mean: " 
		  << mn << " Stdev: " << sqrt(var) << " true: " << coeffs[i] 
		  << std::endl;
	std::cout << "  lower limit: " << low << " upper limit: "
		  << up << " (" << conflevel * 100.0 << "% limit)" << std::endl;
      }

      std::vector<float> accept;
      ply.getAcceptanceFrac(accept);
      double mnacc = accept[0];
      for (unsigned int i = 1; i < nwalkers; ++i)
	mnacc += accept[i];
      std::cout << "Mean acceptance: " << mnacc / static_cast<double>(nwalkers)
		<< std::endl;
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
