#include<iostream>
#include<string>
#include<sstream>

#include<getopt.h>

#include "../include/numberCountsKnotsSplineStats.h"
#include "../include/numberCountsDoubleLogNormalStats.h"
#include "../include/hdf5utils.h"
#include "../include/affineExcept.h"

/*! \brief Factory function for getting statsAccumulator instance */
statsAccumulator* getStatsAccumulator(const std::string& infile,
                                      unsigned int nflux) {

  // Open the input file, figure out the model type
  hid_t fileid = H5Fopen(infile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Read the model type
  hid_t groupid = H5Gopen(fileid, "LikelihoodParams/Model", H5P_DEFAULT);
  std::string mtype = hdf5utils::readAttString(groupid, "ModelType");
  H5Gclose(groupid);
  H5Fclose(fileid);

  if (mtype == "numberCountsKnotsSpline") {
    return new numberCountsKnotsSplineStats(nflux);
  } else if (mtype == "numberCountsDoubleLogNormal") {
    return new numberCountsDoubleLogNormalStats(nflux, true);
  } else {
    std::stringstream errstr;
    errstr << "Got unexpected model type: " << mtype;
    throw affineExcept("pofd_affine_getStats", "getStatsAccumulator",
                       errstr.str());
  }
}

///////////////////////////////////////

int main(int argc, char** argv) {

  struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"double", no_argument, 0, 'd'},
    {"version", no_argument, 0, 'V'},
    {"nflux", required_argument, 0, 'n'},
    {"thin", required_argument, 0, 't'},
    {0, 0, 0, 0}
  };

  char optstring[] = "hdVn:t:";

  unsigned int nflux, thin;

  std::string infile; // pofd_affine_mcmc output file to process
  std::string outputfile; //Ouput pofd option
  
  //Defaults
  nflux               = 1024;
  thin                = 5;

  int c;
  int option_index = 0;
  while ((c = getopt_long(argc, argv, optstring, long_options,
                                            &option_index)) != -1) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getStats -- get summary statistics for" 
                << std::endl;
      std::cerr << "\ta pofd_affine_mcmc fit." << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\t pofd_affine_getStats infile outfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tBuilds summary statistics from the results of a"
                << std::endl;
                  std::cerr << "\t pofd_affine_mcmc fit, 1D or 2D." << std::endl;
      std::cerr << "\tinfile is the HDF5 output of pofd_affine_mcmc."
                << std::endl;
      std::cerr << "\toutfile is the output statistics, also as HDF5."
                << std::endl;
    case 'n' :
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case 't' :
      thin = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'V' :
      std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
                << std::endl;
      return 0;
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
  if (thin == 0) {
    std::cerr << "Error -- thin is invalid (zero)." << std::endl;
      return 1;
  }

  // Actual computation
  statsAccumulator *acc = nullptr;
  try {
    acc = getStatsAccumulator(infile, nflux);
    acc->build(infile, thin, true);
    acc->writeAsHDF5(outputfile);
    delete acc;
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex.what() << std::endl;
    if (acc != nullptr) delete acc;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (acc != nullptr) delete acc;
    return 16;
  }

  return 0;
}
