#include<iostream>
#include<string>
#include<iomanip>

#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/hdf5utils.h"
#include "../include/PDFactory.h"
#include "../include/PDFactoryDouble.h"
#include "../include/numberCountsKnotsSpline.h"
#include "../include/numberCountsDoubleLogNormal.h"
#include "../include/paramSet.h"
#include "../include/affineExcept.h"

//See pofd_affine_getdNdS comment to explain why I do this as a global
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"double", no_argument, 0, 'd'},
  {"version", no_argument, 0, 'V'}, //Below here not parsed in main routine
  {"histogram", no_argument, 0, 'H'},
  {"log", no_argument, 0, 'l'},
  {"nbins", required_argument, 0, '1'},
  {"nflux", required_argument, 0, 'n'},
  {"ninterp", required_argument, 0, 'N'},
  {"rfile", required_argument, 0, 'r'},
  {"sigma", required_argument, 0, 's'},
  {"verbose", no_argument, 0, 'v'},
  {"wisdom", required_argument, 0, 'w'},
  {"nedge", required_argument, 0, '5'},
  {"sigma1", required_argument, 0, '6'},
  {"sigma2", required_argument, 0, '7'},
  {0, 0, 0, 0}
};

char optstring[] = "hdVHl1:n:N:r:s:vw:5:6:7:";

int getPDSingle(int argc, char **argv) {

  double sigma_noise; //Noise
  unsigned int nflux, ninterp, nbins;
  bool has_wisdom, getLog, writeR, histogram, verbose;
  double minflux, maxflux;
  std::string wisdom_file;
  
  std::string modelfile; //Model file
  std::string outputfile; //Ouput pofd option
  std::string psffile;  //Beam file
  std::string rfile; //!< File to write R to
  
  //Defaults
  sigma_noise         = 0.0;
  has_wisdom          = false;
  nflux               = 262144;
  ninterp             = 2048;
  nbins               = 120;
  verbose             = false;
  histogram           = false;
  getLog              = false;
  writeR              = false;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'H' :
      histogram = true;
      break;
    case 'l':
      getLog = true;
      break;
    case '1':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'n' :
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'N' :
      ninterp = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'r':
      writeR = true;
      rfile = std::string(optarg);
      break;
    case 's' :
      sigma_noise = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    case 'w' :
      has_wisdom = true;
      wisdom_file = std::string(optarg);
      break;
    }

  if (optind >= argc - 4) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  psffile    = std::string(argv[optind + 1]);
  minflux    = atof(argv[optind + 2]);
  maxflux    = atof(argv[optind + 3]);
  outputfile = std::string(argv[optind + 4]);

  //Input tests
  if (nflux == 0) {
    std::cerr << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (ninterp == 0) {
    std::cerr << "Error -- ninterp is zero." << std::endl;
      return 1;
  }
  if (ninterp > nflux) {
    std::cerr << "nflux must be larger than ninterp"
	      << std::endl;
    return 1;
  }
  if (sigma_noise < 0.0) {
    std::cerr << "Invalid noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (writeR && rfile.empty()) {
    std::cerr << "Requested write of R to empty filename" << std::endl;
    return 1;
  }

  if (minflux > maxflux) std::swap(minflux, maxflux);
  
  //Actual PD computation
  try {
    initFileKnots model_info(modelfile, false, false);

    numberCountsKnotsSpline model;
    model_info.getKnotPos(model);

    PDFactory pfactory(ninterp);
    beam bm(psffile, histogram, nbins);

    paramSet pars(model_info.getNKnots());
    model_info.getParams(pars);

    model.setParams(pars);
    PD pd;

    if (has_wisdom) {
      if (verbose)
	std::cout << "Reading in wisdom file: " << wisdom_file 
		  << std::endl;
      pfactory.addWisdom(wisdom_file);
    }
    if (verbose) pfactory.setVerbose();

    if (verbose) {
      printf("  Beam area:            %0.3e [deg^2]\n",bm.getEffectiveArea());
      printf("  Pixel size:           %0.1f [arcsec]\n", bm.getPixSize());
      printf("  Flux per area:        %0.3f [Jy deg^-2]\n",
	     model.getFluxPerArea());
      printf("  Nknots:               %u\n",model_info.getNKnots());
      printf("  Inst Sigma:           %0.5f [Jy]\n",sigma_noise);
      printf("  Requesting R range:   %0.3f-%0.3f [Jy]\n",
	     minflux, maxflux);
      double modelMax = model.getMaxFlux();
      printf("  Model range est:      %0.3f-%0.3f [Jy]\n", 0., modelMax);
      if (histogram)
	printf("  Using beam histogramming to reduce beam size to: %u %u\n",
	       bm.getNPos(), bm.getNNeg());
      printf("  Interpolation length:  %u\n",ninterp);
      printf("  Positions and initial values:\n");
      std::pair<double,double> pr;
      for (unsigned int i = 0; i < model_info.getNKnots(); ++i) {
	pr = model_info.getKnot(i);
	printf("   %11.5e  %11.5e\n",pr.first,pr.second);
      }
      if (getLog)
	printf("  Getting Log_2 P(D)\n");
    }

    // Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << std::endl;
    bool succ;
    succ = pfactory.initPD(nflux, minflux, maxflux, model, bm);
    if (!succ) {
      std::cerr << "Error initializing P(D) -- parameters not valid" 
		<< std::endl;
      return 1;
    }    
    std::cerr << "Get PD" << std::endl;
    pfactory.getPD(sigma_noise, pd, getLog);
    
    // Write it
    hdf5utils::outfiletype oft = hdf5utils::getOutputFileType(outputfile);
    if (oft == hdf5utils::HDF5 || oft == hdf5utils::UNKNOWN) {
      if (verbose) std::cout << "Writing P(D) to file " << outputfile 
			     << " as HDF5" << std::endl;
      hid_t file_id, group_id;
      file_id = H5Fcreate(outputfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);

      pd.writeToHDF5Handle(file_id);

      group_id = H5Gcreate(file_id, "Model", H5P_DEFAULT, H5P_DEFAULT, 
			  H5P_DEFAULT);
      model.writeToHDF5Handle(group_id, true);
      H5Gclose(group_id);

      H5Fclose(file_id);
    } else if (oft == hdf5utils::FITS) {
      if (verbose) std::cout << "Writing P(D) to file " << outputfile 
			     << " as FITS" << std::endl;
      pd.writeToFits(outputfile);
    } else if (oft == hdf5utils::TXT) {
      if (verbose) std::cout << "Writing P(D) to file " << outputfile 
			     << " as text" << std::endl;
      std::ofstream ofs(outputfile.c_str());
      if (!ofs) {
	std::cerr << "Error opening output file: " << outputfile
		  << std::endl;
	return 64;
      }
      ofs << pd;
      ofs.close();
    } else {
      std::cerr << "Unsupported output file type for " << outputfile
		<< std::endl;
      return 128;
    }
    
    // Write R as HDF5
    if (writeR) {
      if (verbose) std::cout << "Writing R to: " << rfile 
			     << " as HDF5" << std::endl;
      pfactory.writeRToHDF5(rfile);
    }

#ifdef TIMING
    std::cout << "Timing results:" << std::endl;
    pfactory.summarizeTime(2);
#endif

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

////////////////////////
int getPDDouble(int argc, char** argv) {

  double sigma1, sigma2; //Noise
  unsigned int nflux, nedge, nbins;
  bool has_wisdom, histogram, verbose, getLog, doedge, writeR;
  double minflux1, maxflux1, minflux2, maxflux2;
  std::string wisdom_file;
  
  std::string modelfile; //Model file
  std::string outputfile; //Ouput pofd option
  std::string psffile1, psffile2;  //Beam file
  std::string rfile; //!< File to write R to

  //Defaults
  sigma1              = 0.0;
  sigma2              = 0.0;
  has_wisdom          = false;
  nflux               = 2048;
  verbose             = false;
  histogram           = false;
  nbins               = 150;
  nedge               = 256;
  getLog              = false;
  doedge              = true;
  writeR              = false;

  int c;
  int option_index = 0;
  optind = 1;
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'H':
      histogram = true;
      break;
    case 'l':
      getLog = true;
      break;
    case '1':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case '5':
      nedge = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'n':
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'r':
      writeR = true;
      rfile = std::string(optarg);
      break;
    case '6':
      sigma1 = atof(optarg);
      break;
    case '7':
      sigma2 = atof(optarg);
      break;
    case 'v':
      verbose = true;
      break;
    case 'w':
      has_wisdom = true;
      wisdom_file = std::string(optarg);
      break;
    }

  if (optind >= argc - 7) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  psffile1   = std::string(argv[optind + 1]);
  psffile2   = std::string(argv[optind + 2]);
  minflux1   = atof(argv[optind + 3]);
  maxflux1   = atof(argv[optind + 4]);
  minflux2   = atof(argv[optind + 5]);
  maxflux2   = atof(argv[optind + 6]);
  outputfile = std::string(argv[optind + 7]);

  //Input tests
  if (nflux == 0) {
    std::cerr << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (sigma1 < 0.0) {
    std::cerr << "Invalid noise level in band 1: must be >= 0.0" << std::endl;
    return 1;
  }
  if (sigma2 < 0.0) {
    std::cerr << "Invalid noise level in band 2: must be >= 0.0" << std::endl;
    return 1;
  }
  if (writeR && rfile.empty()) {
    std::cout << "Requested R be written to empty filename" << std::endl;
    return 1;
  }
  if (nedge == 0) doedge = false;
  if (minflux1 > maxflux1) std::swap(minflux1, maxflux1);
  if (minflux2 > maxflux2) std::swap(minflux2, maxflux2);

  //Actual PD computation
  try {
    initFileDoubleLogNormal model_info(modelfile, false, false);

    numberCountsDoubleLogNormal model;
    model_info.getModelPositions(model);
    paramSet pars(model_info.getNTot());
    model_info.getParams(pars);
    model.setParams(pars);

    PDFactoryDouble pfactory(nedge);
    doublebeam bm(psffile1, psffile2, histogram, nbins);

    PDDouble pd;

    if (has_wisdom) {
      if (verbose)
	std::cout << "Reading in wisdom file: " << wisdom_file 
		  << std::endl;
      pfactory.addWisdom(wisdom_file);
    }
    if (verbose) pfactory.setVerbose();

    if (verbose) {
      printf("  Beam area, band 1:       %0.3e [deg^2]\n",
	     bm.getEffectiveArea1());
      printf("  Beam area, band 2:       %0.3e [deg^2]\n",
	     bm.getEffectiveArea2());
      printf("  Pixel size:              %0.1f [arcsec]\n", bm.getPixSize());
      printf("  Flux per area, band 1:   %0.3f [Jy deg^-2]\n",
	     model.getFluxPerArea(0));

      printf("  Flux per area, band 2:   %0.3f [Jy deg^-2]\n",
	     model.getFluxPerArea(1));
      printf("  Nknots:                  %u\n", model_info.getNKnots());
      printf("  Nsigma:                  %u\n", model_info.getNSigmas());
      printf("  Noffset:                 %u\n", model_info.getNOffsets());
      printf("  Inst Sigma, band 1:      %0.5f [Jy]\n", sigma1);
      printf("  Inst Sigma, band 2:      %0.5f [Jy]\n", sigma2);
      printf("  Requesting R ranges:     %0.3f-%0.3f %0.3f-%0.3f [Jy]\n",
	     minflux1, maxflux1, minflux2, maxflux2);
      dblpair modelMax = model.getMaxFlux();
      printf("  Model range est:         %0.3f-%0.3f %0.3f-%0.3f [Jy]\n",
	     0., modelMax.first, 0., modelMax.second);
      if (histogram)
	printf("  Using beam histogramming to reduce beam size to %u %u %u %u\n",
	       bm.getNHist(0), bm.getNHist(1), bm.getNHist(2), bm.getNHist(3));
      else
	printf("  Beam size is %u %u %u %u\n",
	       bm.getNPix(0), bm.getNPix(1), bm.getNPix(2), bm.getNPix(3));

      printf("  Knot Positions and initial values:\n");
      std::pair<double,double> pr;
      for (unsigned int i = 0; i < model_info.getNKnots(); ++i) {
	pr = model_info.getKnot(i);
	printf("   %11.5e  %11.5e\n",pr.first,pr.second);
      }
      printf("  Sigma Positions and initial values:\n");
      for (unsigned int i = 0; i < model_info.getNSigmas(); ++i) {
	pr = model_info.getSigma(i);
	printf("   %11.5e  %11.5e\n",pr.first,pr.second);
      }
      printf("  Offset Positions and initial values:\n");
      for (unsigned int i = 0; i < model_info.getNOffsets(); ++i) {
	pr = model_info.getOffset(i);
	printf("   %11.5e  %11.5e\n",pr.first,pr.second);
      }
      if (getLog)
	printf("  Getting Log_2 P(D)\n");
    }

    //Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << std::endl;
    bool succ;
    succ = pfactory.initPD(nflux, minflux1, maxflux1, minflux2, maxflux2,
			   model, bm, doedge);
    if (!succ) {
      std::cerr << "Error initializing P(D) -- parameters not valid" 
		<< std::endl;
      return 1;
    }    
    pfactory.getPD(sigma1, sigma2, pd, getLog);

    // Write PD
    hdf5utils::outfiletype oft = hdf5utils::getOutputFileType(outputfile);
    if (oft == hdf5utils::HDF5 || oft == hdf5utils::UNKNOWN) {
      if (verbose) std::cout << "Writing P(D) to file " << outputfile 
			     << " as HDF5" << std::endl;
      hid_t file_id, group_id;
      file_id = H5Fcreate(outputfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);

      pd.writeToHDF5Handle(file_id);

      group_id = H5Gcreate(file_id, "Model", H5P_DEFAULT, H5P_DEFAULT, 
			  H5P_DEFAULT);
      model.writeToHDF5Handle(group_id, true);
      H5Gclose(group_id);

      H5Fclose(file_id);
    } else if (oft == hdf5utils::FITS) {
      if (verbose) std::cout << "Writing P(D) to file " << outputfile 
			     << " as FITS" << std::endl;
      pd.writeToFits(outputfile);
    } else if (oft == hdf5utils::TXT) {
      if (verbose) std::cout << "Writing P(D) to file " << outputfile 
			     << " as text" << std::endl;
      std::ofstream ofs(outputfile.c_str());
      if (!ofs) {
	std::cerr << "Error opening output file: " << outputfile
		  << std::endl;
	return 64;
      }
      ofs << pd;
      ofs.close();
    } else {
      std::cerr << "Unsupported output file type for " << outputfile
		<< std::endl;
      return 128;
    }

    // Write R as HDF5
    if (writeR) {
      if (verbose) std::cout << "Writing R to: " << rfile 
			     << " as HDF5" << std::endl;
      pfactory.writeRToHDF5(rfile);
    }

#ifdef TIMING
    std::cout << "Timing results:" << std::endl;
    pfactory.summarizeTime(2);
#endif
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
      std::cerr << "\tpofd_affine_getPD -- get the PD for a number counts "
		<< "model. Both" << std::endl;
      std::cerr << "\t one-dimensional and two-dimensional models are "
		<< "supported." << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tEither" << std::endl;
      std::cerr << std::endl;

      std::cerr << "\t pofd_affine_getPD [options] modelfile beamfile minflux"
		<< std::endl;
      std::cerr << "\t\t maxflux outputfile"
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfor the 1D case or" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\t pofd_affine_getPD -d [options] modelfile beamfile1"
		<< " beamfile2" << std::endl;
      std::cerr << "\t  minflux1 maxflux1 minflux2 maxflux2 outputfile" 
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfor the 2D case." << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates P(D) for the model in modelfile and write it"
		<< " to" << std::endl;
      std::cerr << "\toutputfile.  The 1D model is a log-space spline model for"
		<< " the" << std::endl;
      std::cerr << "\tnumber counts, and the 2D model is the 1D spline model" 
		<< " times" << std::endl;
      std::cerr << "\ta log-normal color function for the second band, with the"
		<< std::endl;
      std::cerr << "\tlog-space variance and mean color stored as splines in"
		<< " the flux" << std::endl;
      std::cerr << "\tof the first band." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tmodelfile is a text file specifying the model; the exact"
		<< " details" << std::endl;
      std::cerr << "\t(given below) depend on whether the 1D or 2D case is"
		<< " being used." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tFor the 1D case, modelfile is a text file giving the "
		<< "positions" << std::endl;
      std::cerr << "\tof the spline knots and their values in the format"
		<< " knotflux value." << std::endl;
      std::cerr << "\tAdditional elements on each line are ignored."
		<< std::endl;
      std::cerr << "\tFor the 2D case, modelfile is a text file giving the "
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
      std::cerr << "\tThe format of outputfile is controlled by the file"
		<< " extension," << std::endl;
      std::cerr << "\twith HDF5, FITS and text supported.  The default is "
		<< "HDF5." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tbeamfile gives the name of a FITS file containing the "
		<< "beam" << std::endl;
      std::cerr << "\tin the 1D case, and beamfile1, beamfile2 give the beam "
		<< "in each" << std::endl;
      std::cerr << "\t of the two bands in the 2D case." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tIn both cases the output PD is written to outputfile "
		<< "either as" << std::endl;
      std::cerr << "\ttext, FITS, or HDF5 depending on the options." 
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tThe minimum and maximum ranges of R (not the output P(D))"
		<< " are" << std::endl;
      std::cerr << "\tspecified by minflux/maxflux." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-d, --double" << std::endl;
      std::cerr << "\t\tUse the 2D model." << std::endl;
      std::cerr << "\t-H, --histogram" << std::endl;
      std::cerr << "\t\tUse beam histogramming." << std::endl;
      std::cerr << "\t-l, --log" << std::endl;
      std::cerr << "\t\tReturn log2 P(D) rather than P(D)." << std::endl;
      std::cerr << "\t--nbins VALUE" << std::endl;
      std::cerr << "\t\tNumber of bins to use in beam histogram (if "
		<< "histogramming" << std::endl;
      std::cerr << "\t\tis being used.  (def: 120 for 1D, 150 for 2D)."
		<< std::endl;
      std::cerr << "\t-n, --nflux VALUE" << std::endl;
      std::cerr << "\t\tThe number of requested fluxes, also the FFT length." 
		<< std::endl;
      std::cerr << "\t\tBest if a power of 2. For the 2D model, this is the"
		<< std::endl;
      std::cerr << "\t\tsize along each dimension. (def: 262144 for 1D, 2048"
		<< std::endl;
      std::cerr << "\t\tfor 2D)." << std::endl;
      std::cerr << "\t-r, --rfile FILENAME" << std::endl;
      std::cerr << "\t\tWrite R as an HDF5 file to this filename."
		<< std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tPrint informational messages while running."
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput the version number and exit." << std::endl;
      std::cerr << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cerr << "\t\tName of FFTW wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "ONE-D ONLY OPTIONS" << std::endl;
      std::cerr << "\t-N, --ninterp value" << std::endl;
      std::cerr << "\t\tLength of interpolation vector used in computing R"
		<< " (def: 2048)." << std::endl;
      std::cerr << "\t-s, --sigma noise" << std::endl;
      std::cerr << "\t\tThe assumed per-pixel noise (assumed Gaussian, "
		<< "def: 0)." << std::endl;
      std::cerr << "TWO-D ONLY OPTIONS" << std::endl;
      std::cerr << "\t--nedge value" << std::endl;
      std::cerr << "\t\tNumber of bins in edge integrals for R."
		<< std::endl;
      std::cerr << "\t--sigma1 noise" << std::endl;
      std::cerr << "\t\tThe assumed per-pixel noise, band 1 (assumed Gaussian,"
		<< std::endl;
      std::cerr << "\t\tdef: 0)." << std::endl;
      std::cerr << "\t--sigma2 noise" << std::endl;
      std::cerr << "\t\tThe assumed per-pixel noise, band 2 (assumed Gaussian,"
		<< std::endl;
      std::cerr << "\t\tdef: 0)." << std::endl;
      return 0;
      break;
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
    return getPDSingle(argc,argv);
  else
    return getPDDouble(argc,argv);
}
