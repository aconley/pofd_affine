#include<iostream>
#include<string>
#include<sstream>

#include<getopt.h>

#include "../include/errorSnakeSingle.h"
#include "../include/affineExcept.h"

//See pofd_affine_getdNdS comment to explain why I do this as a global
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"double", no_argument, 0, 'd'},
  {"version", no_argument, 0, 'V'}, //Below here not parsed in main routine
  {"nflux", required_argument, 0, 'n'},
  {"nsamples", required_argument, 0, 'N'},
  {"verbose", no_argument, 0, 'v'},
  {0, 0, 0, 0}
};

char optstring[] = "hdVn:N:v";

int getErrorSnakeSingle(int argc, char **argv) {

  unsigned int nflux, nsamples;
  bool verbose;
  
  std::string infile; // pofd_affine_mcmc output file to process
  std::string outputfile; //Ouput pofd option
  
  //Defaults
  nflux               = 1024;
  nsamples            = 1000;
  verbose             = false;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'n' :
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'N' :
      nsamples = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 1 ) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  infile = std::string(argv[optind]);
  outputfile = std::string(argv[optind + 1]);

  //Input tests
  if (nflux == 0) {
    std::cerr << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (nsamples == 0) {
    std::cerr << "Error -- nsamples is zero." << std::endl;
      return 1;
  }

  // Actual computation
  try {
    errorSnakeSingle errsnk(nflux);
    errsnk.build(infile, nsamples, false);
    errsnk.writeAsHDF5(outputfile);
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex.what() << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    return 16;
  }

  return 0;
}



///////////////////////////////////////

int main( int argc, char** argv ) {
  bool twod;

  twod = false;

  //Only interested in a) displaying help and b) figuring out
  // if this is 1D or 2D c) displaying the version number
  int c;
  int option_index = 0;
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_errorsnake -- build an error snake for a"
		<< std::endl;
      std::cerr << "\t pofd_affine_mcmc fit." << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\t pofd_affine_errorsnake infile outfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tBuilds error snakes from the results of a"
		<< "pofd_affine_mcmc" << std::endl;
      std::cerr << "\tfit." << std::endl;
      std::cerr << "\tinfile is the HDF5 output of pofd_affine_mcmc."
		<< std::endl;
      std::cerr << "\toutfile is the output error snake, also as HDF5."
		<< std::endl;
    case 'd' :
      twod = true;
      break;
    case 'V' :
      std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
		<< std::endl;
      return 0;
      break;
    }

  if (!twod)
    return getErrorSnakeSingle(argc, argv);
}