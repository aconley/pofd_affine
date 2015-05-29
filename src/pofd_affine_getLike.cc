#include<iostream>
#include<string>
#include<vector>
#include<iomanip>
#include<fstream>
#include<sstream>

#include<getopt.h>

#include "../include/calcLike.h"
#include "../include/specFile.h"
#include "../include/calcLikeDouble.h"
#include "../include/specFileDouble.h"
#include "../include/affineExcept.h"

int getLikeSingle(const std::string& initfile, const std::string specfile,
                  float sigmult=1.0) {

  try {
    //Read in the initialization file knot positions, values
    initFileKnots model_info(initfile, false, false);

    //Data files, evaulation specifications
    specFile spec_info(specfile);

    //Set up model parameters
    unsigned int nknots = model_info.getNKnots();
    if (nknots == 0)
      throw affineExcept("pofd_affine_getLike", "getLikeSingle",
                         "No info read in");

    paramSet pars(nknots + 1);
    model_info.getParams(pars);
    pars[nknots] = sigmult;

    //Make sure we got some data
    if (spec_info.datafiles.size() == 0) 
      throw affineExcept("pofd_affine_getLike", "getLikeSingle",
                         "No datafiles loaded");

    calcLike likeSet(spec_info.fftsize, spec_info.ninterp, 
                     spec_info.bin_data, spec_info.nbins);

    // Set priors
    if (spec_info.has_cfirbprior)
      likeSet.setCFIRBPrior(spec_info.cfirbprior_mean,
                            spec_info.cfirbprior_stdev);
    if (spec_info.has_poissonprior)
      likeSet.setPoissonPrior(spec_info.poissonprior_mean,
                              spec_info.poissonprior_stdev);
    if (spec_info.fit_sigma && spec_info.has_sigprior)
      likeSet.setSigmaPrior(spec_info.sigprior_stdev);
    
    // Regularization
    if (spec_info.regularization_alpha > 0.0)
      likeSet.setRegularizationAlpha(spec_info.regularization_alpha);

    //Read data
    if (spec_info.verbosity >= 2)
      std::cout << "Reading in data files" << std::endl;
    likeSet.readDataFromFiles(spec_info.datafiles, spec_info.psffiles, 
                              spec_info.sigmas, spec_info.like_norm,
                              spec_info.ignore_mask, spec_info.mean_sub, 
                              spec_info.minbeamval, spec_info.beam_histogram, 
                              spec_info.nbeamhist, spec_info.exp_conf);
    
    if (spec_info.has_wisdom_file) likeSet.addWisdom(spec_info.wisdom_file);
    if (spec_info.verbosity >= 2) likeSet.setVerbose();
      
    if (spec_info.verbosity >= 1) {
      printf("  FFTsize:       %u\n", spec_info.fftsize);
      printf("  Nknots:        %u\n", model_info.getNKnots());
      if (spec_info.beam_histogram)
        printf("  Using histogramming to reduce beam size\n");  
      if (spec_info.bin_data)
        printf("  Using histogramming to reduce data size to: %u\n",
               spec_info.nbins);  
      printf("  Positions and initial values:\n");
      std::pair<double,double> pr;
      for (unsigned int i = 0; i < nknots; ++i) {
        pr = model_info.getKnot(i);
        printf("   %11.5e  %11.5e\n",pr.first,pr.second);
      }
    }
  
    //And, get that likelihood
    likeSet.setKnotPositions(model_info);
    likeSet.setRRanges(pars);
    bool pars_rejected;
    double LogLike = likeSet.getLogLike(pars, pars_rejected);
    if (pars_rejected) {
      std::cout << "Parameters rejected in likelihood computation" << std::endl;
      return 1;
    }

    std::cout << "log Likelihoood is: "
              << std::setprecision(7) << LogLike << std::endl;
  } catch (const affineExcept& ex) {
    std::cerr << ex << std::endl;
    return 4;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    return 8;
  }
  return 0;
}


/////////////////////////////////////

int getLikeDouble(const std::string& initfile, const std::string& specfile,
                  float sigmult1=1.0, float sigmult2=1.0) {

  try {
    //Read in the initialization file knot positions, values
    initFileDoubleLogNormal model_info(initfile, false, false);
    
    //data, likelihood params
    specFileDouble spec_info(specfile);
    
    //Set up model
    unsigned int ntot = model_info.getNTot();
    if (ntot == 0)
      throw affineExcept("pofd_affine_getLike", "getLikeDouble",
                         "No info read in");

    // Number of params is number of knots (ntot)
    // + 2 for the sigmas in each band -- although some may be fixed
    // + 2 for the bonus mean flux per area
    // + 2 for the bonus mean flux^2 per area
    paramSet pars(ntot + 6);
    model_info.getParams(pars);
    pars[ntot] = sigmult1;  // Sigma multipliers
    pars[ntot+1] = sigmult2;

    //Make sure we got some data
    if (spec_info.datafiles1.size() == 0) 
      throw affineExcept("pofd_affine_getLike", "getLikeDouble",
                         "No datafiles loaded");
    
    calcLikeDouble likeSet(spec_info.fftsize, spec_info.nedge, 
                           spec_info.edge_set, spec_info.bin_data, 
                           spec_info.nbins);

    //Set priors 
    if (spec_info.has_cfirbprior1)
      likeSet.setCFIRBPrior1(spec_info.cfirbprior_mean1,
                             spec_info.cfirbprior_stdev1);
    if (spec_info.has_cfirbprior2)
      likeSet.setCFIRBPrior2(spec_info.cfirbprior_mean2,
                             spec_info.cfirbprior_stdev2);
    if (spec_info.has_poissonprior1)
      likeSet.setPoissonPrior1(spec_info.poissonprior_mean1,
                               spec_info.poissonprior_stdev1);
    if (spec_info.has_poissonprior2)
      likeSet.setPoissonPrior2(spec_info.poissonprior_mean2,
                               spec_info.poissonprior_stdev2);
    if (spec_info.fit_sigma1 && spec_info.has_sigprior1)
      likeSet.setSigmaPrior1(spec_info.sigprior_stdev1);
    if (spec_info.fit_sigma2 && spec_info.has_sigprior2)
      likeSet.setSigmaPrior2(spec_info.sigprior_stdev2);

    if (spec_info.regularization_alpha > 0.0)
      likeSet.setRegularizationAlpha(spec_info.regularization_alpha);

    //Read data
    likeSet.readDataFromFiles(spec_info.datafiles1, spec_info.datafiles2, 
                              spec_info.psffiles1, spec_info.psffiles2,
                              spec_info.sigmas1, spec_info.sigmas2, 
                              spec_info.like_norm, spec_info.ignore_mask, 
                              spec_info.mean_sub, spec_info.minbeamval,
                              spec_info.beam_histogram, spec_info.nbeamhist,
                              spec_info.exp_conf1, spec_info.exp_conf2);

    if (spec_info.has_wisdom_file) likeSet.addWisdom(spec_info.wisdom_file);
    if (spec_info.verbosity >= 2) likeSet.setVerbose();
      
    if (spec_info.verbosity >= 1) {
      printf("  FFTsize:       %u\n", spec_info.fftsize);
      if (spec_info.beam_histogram)
        printf("  Using histogramming to reduce beam size\n");  
      if (spec_info.bin_data)
        printf("  Using histogramming to reduce data size to: %u by %u\n",
               spec_info.nbins, spec_info.nbins);  
      printf("  Knot Positions and initial values:\n");
      std::pair<double,double> pr;
      for (unsigned int i = 0; i < model_info.getNKnots(); ++i) {
        pr = model_info.getKnot(i);
        printf("   %11.5e  %11.5e\n", pr.first, pr.second);
      }
      printf("  Sigma Positions and initial values:\n");
      for (unsigned int i = 0; i < model_info.getNSigmas(); ++i) {
        pr = model_info.getSigma(i);
        printf("   %11.5e  %11.5e\n", pr.first, pr.second);
      }
      printf("  Offset Positions and initial values:\n");
      for (unsigned int i = 0; i < model_info.getNOffsets(); ++i) {
        pr = model_info.getOffset(i);
        printf("   %11.5e  %11.5e\n", pr.first, pr.second);
      }
    }

    //Get likelihood
    likeSet.setupModel(model_info);
    likeSet.setRRanges(pars); // Have to set parameters first

    bool pars_rejected;
    double LogLike = likeSet.getLogLike(pars, pars_rejected);
    if (pars_rejected) {
      std::cout << "Parameters rejected in likelihood computation" << std::endl;
      return 1;
    }

    std::cout << "log Likelihoood is: "
              << std::setprecision(7) << LogLike << std::endl;
  } catch (const affineExcept& ex) {
    std::cerr << ex << std::endl;
    return 4;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    return 8;
  }
  return 0;
}

/////////////////////////////////

int main(int argc, char** argv) {
  bool twod;
  float sigmult, sigmult1, sigmult2;
  std::string initfile, specfile;

  twod = false;
  sigmult = 1.0;
  sigmult1 = 1.0;
  sigmult2 = 1.0;

  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"double",no_argument, 0, 'd'},
    {"sigmult", required_argument, 0, 's'},
    {"sigmult1", required_argument, 0, '!'},
    {"sigmult2", required_argument, 0, '@'},
    {"version",no_argument, 0, 'V'},
    {0,0,0,0}
  };
  char optstring[] = "hds:!:@:V";
  int c;
  int option_index = 0;

  while ( ( c = getopt_long(argc,argv,optstring,long_options,
                            &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getLike -- evaluate the likelihood of some "
                << "data" << std::endl;
      std::cerr << "\tagainst a number counts model using the P(D) approach."
                << " Both" << std::endl;
      std::cerr << "\tone-dimensional and two-dimensional models are supported."
                << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\t pofd_mcmc_getLike [options] initfile specfile"
                << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates the likelihood of the data specified in "
                << "specfile against" << std::endl;
      std::cerr << "\tthe model in initfile. The 1D model is a log-space spline"
                << " model for" << std::endl;
      std::cerr << "\tthe number counts, and the 2D model is the 1D spline "
                << "model times" << std::endl;
      std::cerr << "\ta log-normal color function for the second band, with the"
                << std::endl;
      std::cerr << "\tlog-space variance and mean color stored as splines in"
                << " the flux" << std::endl;
      std::cerr << "\tof the first band." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tinitfile is a text file specifying the model; the exact"
                << " details" << std::endl;
      std::cerr << "\t(given below) depend on whether the 1D or 2D case is"
                << " being used." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tFor the 1D case, initfile is a text file giving the "
                << "positions" << std::endl;
      std::cerr << "\tof the spline knots and their values in the format"
                << " knotflux value." << std::endl;
      std::cerr << "\tAdditional elements on each line are ignored."
                << std::endl;
      std::cerr << "\tFor the 2D case, initfile is a text file giving the "
                << "positions" << std::endl;
      std::cerr << "\tof the knot points and their values, followed by the "
                << "sigma" << std::endl;
      std::cerr << "\tknot positions and their values, then likewise for the "
                << "colour" << std::endl;
      std::cerr << "\toffset.  The format is three numbers on the first line, "
                << "giving" << std::endl;
      std::cerr << "\tthe number of number count knots, sigma knots, and "
                << "offset knots," << std::endl;
      std::cerr << "\tfollowed by a number of lines again with the format"
                << std::endl;
      std::cerr << "\tknotpos value.  The sigmas and offsets are in log space."
                << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tspecfile is a text file containing the beam(s), "
                << "datafile(s), sigma(s)," << std::endl;
      std::cerr << "\tetc. using the same format as pofd_affine_mcmc.  It also"
                << " controls" << std::endl;
      std::cerr << "\thow the likelihood is evaluated "
                << "(beam histogramming, etc.)." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-d, --double" << std::endl;
      std::cerr << "\t\tUse the 2D model." << std::endl;
      std::cerr << "\t-h, --help" << std::endl;
      std::cerr << "\t\tOutput this help message and exit." << std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput the version number and exit" << std::endl;
      std::cerr << "ONE-D OPTIONS" << std::endl;
      std::cerr << "\t-s, --sigmult VALUE" << std::endl;
      std::cerr << "\t\tSigma multplier to use (def: 1.0)." << std::endl;
      std::cerr << "TWO-D OPTIONS" << std::endl;
      std::cerr << "\t--sigmult1 VALUE" << std::endl;
      std::cerr << "\t\tSigma multplier to use, band 1 (def: 1.0)." 
                << std::endl;
      std::cerr << "\t--sigmult2 VALUE" << std::endl;
      std::cerr << "\t\tSigma multplier to use, band 2 (def: 1.0)." 
                << std::endl;
      return 0;
      break;
    case 'd' :
      twod = true;
      break;
    case 's':
      sigmult = atof(optarg);
      break;
    case '!':
      sigmult1 = atof(optarg);
      break;
    case '@':
      sigmult2 = atof(optarg);
      break;
    case 'V' :
      std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
                << std::endl;
      return 0;
      break;
    }

  if (optind >= argc-1 ) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
              << std::endl;
    return 1;
  }
  initfile = std::string(argv[optind]);
  specfile = std::string(argv[optind+1]);

  if (! twod) {
    if (sigmult <= 0.0) {
      std::cerr << "Invalid (non-positive) user supplied sigmult: "
                << sigmult << std::endl;
      return 1;
    }
    return getLikeSingle(initfile, specfile, sigmult);
  } else {
    if (sigmult1 <= 0.0) {
      std::cerr << "Invalid (non-positive) user supplied sigmult1: "
                << sigmult1 << std::endl;
      return 1;
    }
    if (sigmult2 <= 0.0) {
      std::cerr << "Invalid (non-positive) user supplied sigmult2: "
                << sigmult2 << std::endl;
      return 1;
    }

    return getLikeDouble(initfile, specfile, sigmult1, sigmult2);
  }
}
