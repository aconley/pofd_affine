#include<iostream>
#include<string>

#include<getopt.h>

#include "../include/pofdMCMC.h"
#include "../include/pofdMCMCDouble.h"
#include "../include/affineExcept.h"

int main(int argc, char** argv) {
  bool twod;
  std::string initfile, specfile, outfile;
  unsigned int nwalkers;
  unsigned int nsamples;
  // Number of steps totally ignored
  unsigned int init_steps; 
  // Minimum number of steps before burn check after init_Steps
  unsigned int min_burn;
  // Do fixed burn-in length
  bool fixed_burn;
  // Fraction of autocorr steps to add before rechecking burn-in
  double burn_multiple; 
  // Scale fac for Z generation (controls how new steps are generated)
  double scalefac;

  //Defaults
  twod = false;
  nwalkers = 200;
  nsamples = 2000;
  init_steps = 30;
  min_burn = 30;
  fixed_burn = false;
  burn_multiple = 1.5;
  scalefac = 2.0;

  MPI_Init(&argc, &argv);

  int rank, nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  static struct option long_options[] = {
    {"burnmult",required_argument,0,'b'},
    {"help", no_argument, 0, 'h'},
    {"double",no_argument,0,'d'},
    {"fixedburn", no_argument, 0, 'f'},
    {"initsteps",required_argument,0,'i'},
    {"minburn",required_argument,0,'m'},
    {"nsamples",required_argument,0,'n'},
    {"nwalkers",required_argument,0,'N'},
    {"scalefac",required_argument,0,'s'},
    {"version",no_argument,0,'V'}, 
    {0,0,0,0}
  };
  char optstring[] = "b:hdfi:m:n:N:s:V";
  int c;
  int option_index = 0;

  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      if (rank == 0) {
	std::cerr << "NAME" << std::endl;
	std::cerr << "\tpofd_affine_mcmc -- run a MCMC using a P(D) model. "
		  << "Both" << std::endl;
	std::cerr << "\tone-dimensional and two-dimensional models are "
		  << "supported." << std::endl;
	std::cerr << std::endl;
	std::cerr << "SYNOPSIS" << std::endl;
	std::cerr << "\t pofd_affine_mcmc [options] initfile specfile outfile"
		  << std::endl;
	std::cerr << std::endl;
	std::cerr << "DESCRIPTION" << std::endl;
	std::cerr << "\tUses the P(D) formalism to run an affine invariant"
		  << "MCMC fit" << std::endl;
	std::cerr << "\tto a data set.  Broadly speaking, initfile specifies "
		  << "the model" << std::endl;
	std::cerr << "\tspecification, initial setup, and parameter limits, "
		  << "while specfile" << std::endl;
	std::cerr << "\tspecifies the data and beams as well as a number of "
		  << "items related" << std::endl;
	std::cerr << "\tto the fitting, such as priors, whether the instrument"
		  << " noise" << std::endl;
	std::cerr << "\tis marginalized over, or whether the beam is "
		  << "histogrammed.  The output is written to outfile." 
		  << std::endl;
	std::cerr << std::endl;
	std::cerr << "OPTIONS" << std::endl;
	std::cerr << "\t-b, --burnmult BURN_MULTIPLE" << std::endl;
	std::cerr << "\t\tThe number of extra steps taking between burn in"
		  << std::endl;
	std::cerr << "\t\tchecks is this times the estimated autocorrelation"
		  << std::endl;
	std::cerr << "\t\tlength (def: 1.5)" << std::endl;
	std::cerr << "\t-d, --double" << std::endl;
	std::cerr << "\t\tUse the 2D model." << std::endl;
	std::cerr << "\t-f, --fixedburn" << std::endl;
	std::cerr << "\t\tUsed a fixed burn in length before starting main"
		  << " sample," << std::endl;
	std::cerr << "\t\trather than using the autocorrelation.  With this"
		  << " set" << std::endl;
	std::cerr << "\t\tmin_burn steps will be carried out, then the"
		  << " main" << std::endl;
	std::cerr << "\t\tmain sampling will start -- in other words,"
		  << " initsteps is" << std::endl;
	std::cerr << "\t\tset to zero." << std::endl;
	std::cerr << "\t-h, --help" << std::endl;
	std::cerr << "\t\tOutput this help message and exit." << std::endl;
	std::cerr << "\t-i, --initsteps INITIAL_STEPS" << std::endl;
	std::cerr << "\t\tThis many steps per walker are taken and thrown away"
		  << std::endl;
	std::cerr << "\t\t first (def: 30)." << std::endl;
	std::cerr << "\t-m, --minburn NSTEPS" << std::endl;
	std::cerr << "\t\tThis many steps are taken per walker before the burn"
		  << std::endl;
	std::cerr << "\t\tin is checked.  This also controls how many "
		  << "additional"
		  << std::endl;
	std::cerr << "\t\tsteps get taken if burn in hasn't happened, in"
		  << " combination" << std::endl;
	std::cerr << "\t\twith bunmult (def: 30)." << std::endl;
	std::cerr << "\t-n, --nsamples NSAMPLES" << std::endl;
	std::cerr << "\t\tNumber of total samples to generate across all "
		  << "walkers." << std::endl;
	std::cerr << "\t\tThis is rounded to an integer power of the number of"
		  << std::endl;
	std::cerr << "\t\twalkers (def: 2000)." << std::endl;
	std::cerr << "\t-N, --nwalkers NWALKERS" << std::endl;
	std::cerr << "\t\tNumber of walkers to use (def: 200)." << std::endl;
	std::cerr << "\t-V, --version" << std::endl;
	std::cerr << "\t\tOutput the version number and exit" << std::endl;
	std::cerr << "\t-s, --scalefac SCALEFACTOR" << std::endl;
	std::cerr << "\t\tScale factor in stretch steps (def: 2.0)."
		  << std::endl;
      }
      MPI_Finalize();
      return 0;
      break;
    case 'b' :
      burn_multiple = atof(optarg);
      break;
    case 'd' :
      twod = true;
      break;
    case 'f':
      fixed_burn = true;
      break;
    case 'i' :
      init_steps = atoi(optarg);
      break;
    case 'm' :
      min_burn = atoi(optarg);
      break;
    case 'n' :
      nsamples = atoi(optarg);
      break;
    case 'N':
      nwalkers = atoi(optarg);
      break;
    case 's' :
      scalefac = atof(optarg);
      break;
    case 'V' :
      if (rank == 0) 
	std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
		  << std::endl;
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

  if (optind >= argc - 2) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  initfile = std::string(argv[optind]);
  specfile = std::string(argv[optind+1]);
  outfile = std::string(argv[optind+2]);

  //Test inputs
  if (nwalkers == 0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) number of walkers: "
		<< nwalkers << std::endl;
    return 1;
  }
  if (nsamples == 0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) number of samples: "
		<< nsamples << std::endl;
    return 1;
  }
  if (burn_multiple < 1.0) {
    if (rank == 0)
      std::cerr << "Burn multiple must be >= 1.0, you provided: "
		<< burn_multiple << std::endl;
    return 1;
  }
  if (min_burn == 0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) minburn: "
		<< min_burn << std::endl;
    return 1;
  }
  if (scalefac <= 0.0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) scalefac: "
		<< scalefac << std::endl;
    return 1;
  }

  if (fixed_burn) init_steps = 0;
    
  try {
    if (!twod) {
      pofdMCMC engine(initfile, specfile, nwalkers, nsamples, init_steps,
		      min_burn, fixed_burn, burn_multiple, scalefac);
      
      if (rank == 0) std::cout << "Initializing chains" << std::endl;
      engine.initChains();
      
      if (rank == 0) std::cout << "Starting main loop" << std::endl;
      engine.sample();
      
      if (rank == 0) {
	std::cout << "Done with sampling; writing to: " << outfile << std::endl;
	engine.writeToFile(outfile);
      }
      MPI_Finalize();
      
    } else {
      
      pofdMCMCDouble engine(initfile, specfile, nwalkers, nsamples, init_steps,
			    min_burn, fixed_burn, burn_multiple, scalefac);
      
      if (rank == 0) std::cout << "Initializing chains" << std::endl;
      engine.initChains();
      
      if (rank == 0) std::cout << "Starting main loop" << std::endl;
      engine.sample();
      
      if (rank == 0) {
	std::cout << "Done with sampling; writing to: " << outfile << std::endl;
	engine.writeToFile(outfile);
      }
      MPI_Finalize();
    }
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered: " << ex << std::endl;
    return 1;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Allocation error encountered: " << ba.what() 
	      << std::endl;
    return 2;
  }

  return 0;
}
