#include<iostream>
#include<string>

#include<getopt.h>

#include "../include/pofdMCMC.h"
#include "../include/pofdMCMCDouble.h"
#include "../include/affineExcept.h"

int main(int argc, char** argv) {

  bool twod, write_as_hdf;
  std::string initfile, specfile, outfile;
  unsigned int nwalkers;
  unsigned int nsamples;
  // Number of steps totally ignored
  unsigned int init_steps; 
  // Temperature during init_steps
  double init_temp;
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
  init_temp = 2.0;
  min_burn = 30;
  fixed_burn = false;
  burn_multiple = 1.5;
  scalefac = 2.0;
  write_as_hdf = false;

  MPI_Init(&argc, &argv);

  int rank, nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  static struct option long_options[] = {
    {"burnmult",required_argument,0,'b'},
    {"help", no_argument, 0, 'h'},
    {"double",no_argument,0,'d'},
    {"fixedburn", no_argument, 0, 'f'},
    {"hdf", no_argument, 0, 'H'},
    {"initsteps",required_argument,0,'i'},
    {"inittemp", required_argument, 0, 'I'},
    {"minburn",required_argument,0,'m'},
    {"nsamples",required_argument,0,'n'},
    {"nwalkers",required_argument,0,'N'},
    {"scalefac",required_argument,0,'s'},
    {"version",no_argument,0,'V'}, 
    {0,0,0,0}
  };
  char optstring[] = "b:hdfHi:I:m:n:N:s:V";
  int c;
  int option_index = 0;

  while ((c = getopt_long(argc, argv, optstring, long_options,
			  &option_index)) != -1) 
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
		  << "histogrammed.  The " << std::endl;
	std::cerr << "\toutput is written to outfile." << std::endl;
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
	std::cerr << "\t\trather than using the autocorrelation." << std::endl;
	std::cerr << "\t-h, --help" << std::endl;
	std::cerr << "\t\tOutput this help message and exit." << std::endl;
	std::cerr << "\t-H, --hdf" << std::endl;
	std::cerr << "\t\tWrite as HDF5 file rather than text file."
		  << std::endl;
	std::cerr << "\t-i, --initsteps INITIAL_STEPS" << std::endl;
	std::cerr << "\t\tThis many steps per walker are taken and thrown away"
		  << std::endl;
	std::cerr << "\t\t first.  The highest likelihood point from this set"
		  << " is then" << std::endl;
	std::cerr << "\t\tused to re-seed the initial conditions before "
		  << "burn-in" << std::endl;
	std::cerr << "\t\tis started.  If zero, no initial steps performed "
		  << "(def: 30)." << std::endl;
	std::cerr << "\t-I, --inittemp VALUE" << std::endl;
	std::cerr << "\t\tTemperature used during initial steps (def: 2.0)"
		  << std::endl;
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
    case 'H':
      write_as_hdf = true;
      break;
    case 'i' :
      init_steps = atoi(optarg);
      break;
    case 'I':
      init_temp = atof(optarg);
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

  if (optind >= argc - 2) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    MPI_Finalize();
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
    MPI_Finalize();
    return 1;
  }
  if (init_steps > 0 && init_temp <= 0.0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) initial temperature: "
		<< init_temp << std::endl;
    MPI_Finalize();
    return 1;
  }
  if (nsamples == 0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) number of samples: "
		<< nsamples << std::endl;
    MPI_Finalize();
    return 1;
  }
  if (burn_multiple < 1.0) {
    if (rank == 0)
      std::cerr << "Burn multiple must be >= 1.0, you provided: "
		<< burn_multiple << std::endl;
    MPI_Finalize();
    return 1;
  }
  if (min_burn == 0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) minburn: "
		<< min_burn << std::endl;
    MPI_Finalize();
    return 1;
  }
  if (scalefac <= 0.0) {
    if (rank == 0)
      std::cerr << "Invalid (non-positive) scalefac: "
		<< scalefac << std::endl;
    MPI_Finalize();
    return 1;
  }

  if (nproc < 2) {
    if (rank == 0) 
      std::cerr << "Must run on multiple processes" << std::endl;
    MPI_Finalize();
    return 1;
  }

  try {
    if (!twod) {
      // Single band case
      pofdMCMC engine(initfile, specfile, nwalkers, nsamples, init_steps,
		      init_temp, min_burn, fixed_burn, burn_multiple, 
		      scalefac);
      
      // Initialize
      engine.initChains();
      if (rank == 0 && (engine.getVerbosity() > 1))
	std::cout << engine << std::endl;

      // Main loop
      engine.sample();
      
      if (rank == 0) {
	if (engine.getVerbosity() > 1) {
	  std::vector<float> accept;
	  engine.getAcceptanceFrac(accept);
	  float mnval = accept[0];
	  for (unsigned int i = 1; i < accept.size(); ++i)
	    mnval += accept[i];
	  std::cout << "Mean acceptance fraction: " 
		    << mnval / (accept.size() - 1.0) << std::endl;
	}

	// Compute autocorrelation length
	if (write_as_hdf || (engine.getVerbosity() > 1)) {
	  std::vector<float> acor;
	  bool succ = engine.getAcor(acor);
	  if (engine.getVerbosity() > 1) {
	    if (succ) {
	      std::cout << "Autocorrelation length: " << acor[0];
	      for (unsigned int i = 1; i < acor.size(); ++i)
		std::cout << " " << acor[i];
	      std::cout << std::endl;
	    } else {
	      std::cout << "Failed to compute autocorrelation" 
			<< std::endl;
	    }
	  }
	}

	if (engine.getVerbosity() > 1)
	  std::cout << "Done with sampling; writing to: " << outfile 
		    << std::endl;
	if (write_as_hdf)
	  engine.writeToHDF5(outfile);
	else
	  engine.writeToFile(outfile);
      }

    } else {
      // Dual band case

      pofdMCMCDouble engine(initfile, specfile, nwalkers, nsamples, init_steps,
			    init_temp, min_burn, fixed_burn, burn_multiple, 
			    scalefac);
      
      // Initialize
      engine.initChains();
      if (rank == 0 && (engine.getVerbosity() > 1))
	std::cout << engine << std::endl;

      // Main loop
      engine.sample();
      
      if (rank == 0) {
	if (engine.getVerbosity() > 1) {
	  std::vector<float> accept;
	  engine.getAcceptanceFrac(accept);
	  float mnval = accept[0];
	  for (unsigned int i = 1; i < accept.size(); ++i)
	    mnval += accept[i];
	  std::cout << "Mean acceptance fraction: " 
		    << mnval / (accept.size() - 1.0) << std::endl;
	}

	if (write_as_hdf || (engine.getVerbosity() > 1)) {
	  std::vector<float> acor;
	  bool succ = engine.getAcor(acor);
	  if (engine.getVerbosity() > 1) {
	    if (succ) {
	      std::cout << "Autocorrelation length: " << acor[0];
	      for (unsigned int i = 1; i < acor.size(); ++i)
		std::cout << " " << acor[i];
	      std::cout << std::endl;
	    } else {
	      std::cout << "Failed to compute autocorrelation" 
			<< std::endl;
	    }
	  }
	}

	if (engine.getVerbosity() > 1)
	  std::cout << "Done with sampling; writing to: " << outfile 
		    << std::endl;
	if (write_as_hdf)
	  engine.writeToHDF5(outfile);
	else
	  engine.writeToFile(outfile);
      }
    }
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered: " << ex << std::endl;
    int jnk;
    if (rank == 0)
      for (int i = 1; i < nproc; ++i)
	MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
    else
      MPI_Send(&jnk, 1, MPI_INT, 0, mcmc_affine::ERROR, MPI_COMM_WORLD);
    MPI_Finalize();
    return 1;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Allocation error encountered: " << ba.what() 
	      << std::endl;
    int jnk;
    if (rank == 0)
      for (int i = 1; i < nproc; ++i)
	MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
    else
      MPI_Send(&jnk, 1, MPI_INT, 0, mcmc_affine::ERROR, MPI_COMM_WORLD);
    MPI_Finalize();
    return 2;
  }

  MPI_Finalize();
  return 0;
}
