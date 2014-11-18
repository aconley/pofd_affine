#include<iostream>

#include<getopt.h>
#include<hdf5.h>

#include "../include/global_settings.h"
#include "../include/numberCountsKnotsSpline.h"
#include "../include/numberCountsDoubleLogNormal.h"
#include "../include/paramSet.h"
#include "../include/hdf5utils.h"
#include "../include/affineExcept.h"

//All sub-parses have to have the same options to avoid
// getting warnings.  So we give them all the same long_options,
// but then only process the appropriate ones, ignoring the rest
//Yes, this is a bit complicated and error prone, but such is life
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"double", no_argument, 0, 'd'},
  {"logspace", no_argument, 0, 'l'},
  {"logcounts", no_argument, 0, 'L'},
  {"version", no_argument, 0, 'V'}, 
  {"projected", no_argument, 0, 'p'}, //After this, not parsed by main
  {"minflux", required_argument, 0, 'm'},
  {"maxflux", required_argument, 0, 'M'},
  {"nflux", required_argument, 0, 'n'},
  {"verbose", no_argument, 0, 'v'},
  {"minflux1", required_argument, 0, '1'},
  {"minflux2", required_argument, 0, '2'},
  {"maxflux1", required_argument, 0, '3'},
  {"maxflux2", required_argument, 0, '4'},
  {"nflux1", required_argument, 0, '5'},
  {"nflux2", required_argument, 0, '6'},
  {0, 0, 0, 0}
};
char optstring[] = "hdLlVpm:M:n:v1:2:3:4:5:6:"; 

//One-D version
int getNSingle(int argc, char** argv) {
  
  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  double minflux, maxflux;
  unsigned int nflux;
  bool has_user_maxflux, has_user_minflux, verbose, logspace, logcounts;

  minflux          = 0;
  maxflux          = 1; //Will always be overridden
  nflux            = 1000;
  has_user_minflux = false;
  has_user_maxflux = false;
  logspace         = false;
  logcounts        = false;
  verbose          = false;

  int c;
  int option_index = 0;
  optind = 1; //Resets parsing
  while ((c = getopt_long(argc, argv, optstring, long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'l' :
      logspace = true;
      break;
    case 'L':
      logcounts = true;
      break;
    case 'm' :
      has_user_minflux = true;
      minflux = atof(optarg);
      break;
    case 'M' :
      has_user_maxflux = true;
      maxflux = atof(optarg);
      break;
    case 'n' :
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 1) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  initfile = std::string(argv[optind]);
  outfile = std::string(argv[optind + 1]);

  //Input check
  if (has_user_maxflux && (maxflux <= minflux)) {
    std::cerr << "Maxflux must be > minflux" << std::endl;
    return 1;
  }
  if (nflux == 0) {
    std::cerr << "nflux must be positive" << std::endl;
    return 2;
  }

  double *dNdS = nullptr;
  try {
    initFileKnots model_info(initfile, false, false);

    numberCountsKnotsSpline model;
    model_info.getKnotPos(model);

    paramSet pars(model_info.getNKnots());
    model_info.getParams(pars);
    model.setParams(pars);

    //Check minflux if logspace is set
    if (logspace) {
      if (!has_user_minflux)
	minflux = model.getKnotPos(0); //Sorted inside model
      if (minflux <= 0.0)
	throw affineExcept("pofd_affine_getdNdS", "getNSingle",
			   "Minflux must be positive if using logspace");
    }

    //Get maximum flux if not set, use to compute dflux
    if (!has_user_maxflux )
      maxflux = model.getMaxFlux();

    double dflux;
    if (nflux > 1) 
      if (logspace) 
	dflux = log2(maxflux / minflux) / static_cast<double>(nflux - 1);
      else
	dflux = (maxflux - minflux)/static_cast<double>(nflux - 1);
    else
      dflux = 1.0; //Not used

    if (verbose) {
      std::cout << "Flux per sq degree: " << model.getFluxPerArea()
		<< std::endl;
      std::cout << "Number of sources per area: "
		<< model.getNS() << std::endl;
    }

    //Calculation loop
    dNdS  = new double[nflux];
    if (logspace) {
      double lmin = log2(minflux);
      for (unsigned int i = 0; i < nflux; ++i)
	dNdS[i] = 
	  model.getNumberCounts(exp2(lmin + static_cast<double>(i) * dflux));
    } else
      for (unsigned int i = 0; i < nflux; ++i)
	dNdS[i] = model.getNumberCounts(minflux + static_cast<double>(i) * 
					dflux);

    if (logcounts)
      for (unsigned int i = 0; i < nflux; ++i)
	dNdS[i] = log10(dNdS[i]);
    
    //Write out; text or HDF5
    hdf5utils::outfiletype oft = hdf5utils::getOutputFileType(outfile);
    if (oft == hdf5utils::HDF5 || oft == hdf5utils::UNKNOWN) {
      if (verbose) std::cout << "Writing to: " << outfile 
			     << " as HDF5" << std::endl;
      hid_t file_id;
      file_id = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);
      if (H5Iget_ref(file_id) < 0) {
	H5Fclose(file_id);
	throw affineExcept("pofd_affine_getdNdS", "getNSingle",
			   "Failed to open HDF5 file to write");
      }
      hsize_t adims;
      hid_t mems_id, att_id, group_id;
      hbool_t btmp;

      // Properties
      adims = 1;
      mems_id = H5Screate_simple(1, &adims, nullptr);
      btmp = static_cast<hbool_t>(logspace);
      att_id = H5Acreate2(file_id, "logspace", H5T_NATIVE_HBOOL,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_HBOOL, &btmp);
      H5Aclose(att_id);
      att_id = H5Acreate2(file_id, "dflux", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux);
      H5Aclose(att_id);
      H5Sclose(mems_id);

      // Model
      group_id = H5Gcreate(file_id, "Model", H5P_DEFAULT, H5P_DEFAULT, 
			  H5P_DEFAULT);
      model.writeToHDF5Handle(group_id, true);
      H5Gclose(group_id);

      // Fluxes
      double* fluxes;
      fluxes = new double[nflux];
      if (logspace) {
	double lmin = log2(minflux);
	for (unsigned int i = 0; i < nflux; ++i) 
	  fluxes[i] = exp2(lmin + static_cast<double>(i) * dflux);
      } else
	for (unsigned int i = 0; i < nflux; ++i) 
	  fluxes[i] = minflux + static_cast<double>(i) * dflux;
      hdf5utils::writeDataDoubles(file_id, "Flux", nflux, fluxes);
      delete[] fluxes;

      // dNdS
      if (logcounts)
	hdf5utils::writeDataDoubles(file_id, "Log10dNdS", nflux, dNdS);
      else
	hdf5utils::writeDataDoubles(file_id, "dNdS", nflux, dNdS);

      H5Fclose(file_id);
    } else if (oft == hdf5utils::TXT) {
      FILE *fp;
      fp = fopen( outfile.c_str(), "w");
      if (!fp) {
	throw affineExcept("pofd_affine_getdNdS", "getNSingle",
			   "Failed to open text file file to write");
      }
      if (logcounts)
	fprintf(fp, "%12s   %12s\n", "Flux", "Log10dNdS");
      else
	fprintf(fp, "%12s   %12s\n", "Flux", "dNdS");
      if (logspace) {
	double lmin = log2(minflux);
	for (unsigned int i = 0; i < nflux; ++i) 
	  fprintf(fp, "%12.6e   %15.9e\n",
		  exp2(lmin + static_cast<double>(i) * dflux), dNdS[i]);
      } else
	for (unsigned int i = 0; i < nflux; ++i) 
	  fprintf(fp, "%12.6e   %15.9e\n",
		  minflux + static_cast<double>(i) * dflux, dNdS[i]);
      fclose(fp);
    } else if (oft == hdf5utils::FITS)
      throw affineExcept("pofd_affine_getdNdS", "getNSingle",
			 "FITS output not supported");
    
    delete[] dNdS; 
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (dNdS != nullptr) delete[] dNdS;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (dNdS != nullptr) delete[] dNdS;
    return 16;
  }
  return 0;
}

/////////////////////////////////////////////

//Two-D projected version
int getNProjected(int argc, char** argv) {
  
  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  double minflux2, maxflux2;
  unsigned int nflux2;
  bool has_user_minflux2, has_user_maxflux2, verbose, logspace, logcounts;

  minflux2          = 0;
  maxflux2          = 1; //Will always be overridden
  nflux2            = 100;
  has_user_minflux2 = false;
  has_user_maxflux2 = false;
  logcounts         = false;
  logspace          = false;
  verbose           = false;

  int c;
  int option_index = 0;
  optind = 1; //Resets parsing
  while ((c = getopt_long(argc, argv, optstring, long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'l' :
      logspace = true;
      break;
    case 'L':
      logcounts = true;
      break;
    case 'm' :
      has_user_minflux2 = true;
      minflux2 = atof(optarg);
      break;
    case 'M' :
      has_user_maxflux2 = true;
      maxflux2 = atof(optarg);
      break;
    case 'n' :
      nflux2 = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-1) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  initfile = std::string(argv[optind]);
  outfile = std::string(argv[optind + 1]);

  //Input check
  if (has_user_minflux2 && (maxflux2 <= minflux2)) {
    std::cerr << "Maxflux must be > minflux in band 2" << std::endl;
    return 2;
  }
  if (has_user_maxflux2 && (maxflux2 <= minflux2)) {
    std::cerr << "Maxflux must be > minflux in band 2" << std::endl;
    return 2;
  }
  if (nflux2 == 0) {
    std::cerr << "nflux must be positive in band 2" << std::endl;
    return 4;
  }


  double *dNdS = nullptr;
  try {
    initFileDoubleLogNormal model_info(initfile, false, false);
    
    numberCountsDoubleLogNormal model;
    model_info.getModelPositions(model);

    paramSet pars(model_info.getNTot());
    model_info.getParams(pars);
    model.setParams(pars);

    if (logspace) {
      if (!has_user_minflux2)
	minflux2 = model.getKnotPosition(0); //Sorted inside model
      if (minflux2 <= 0.0)
	throw affineExcept("pofd_affine_getdNdS", "getNProjected",
			   "Minflux must be positive if using logspace");
    }

    //Get maximum flux if not set, use to compute dflux
    if (!has_user_maxflux2)
      maxflux2 = model.getMaxFlux().second;
    double dflux2;
    if (nflux2 > 1) {
      if (logspace)
	dflux2 = log2(maxflux2 / minflux2) / static_cast<double>(nflux2 - 1);
      else
	dflux2 = (maxflux2 - minflux2) / static_cast<double>(nflux2 - 1);
    } else
      dflux2 = 1.0; //Not used, avoids compiler warning

    if (verbose) {
      std::cout << "Flux per sq degree, band 1: "
		<< model.getFluxPerArea(0) << std::endl;
      std::cout << "Flux per sq degree, band 2: "
		<< model.getFluxPerArea(1) << std::endl;
      std::cout << "Number of sources per area: "
		<< model.getNS() << std::endl;
    }

    //Calculation loop
    dNdS  = new double[nflux2];
    double flux2;
    if (logspace) {
      double lmin = log2(minflux2);
      for (unsigned int i = 0; i < nflux2; ++i) {
	flux2 = exp2(lmin + static_cast<double>(i) * dflux2);
	dNdS[i] = model.getBand2NumberCounts(flux2);
      }
    } else {
      for (unsigned int i = 0; i < nflux2; ++i) {
	flux2 = minflux2 + static_cast<double>(i) * dflux2;
	dNdS[i] = model.getBand2NumberCounts(flux2);
      }
    }

    if (logcounts)
      for (unsigned int i = 0; i < nflux2; ++i)
	dNdS[i] = log10(dNdS[i]);
    
    //Write out
    hdf5utils::outfiletype oft = hdf5utils::getOutputFileType(outfile);
    if (oft == hdf5utils::HDF5 || oft == hdf5utils::UNKNOWN) {
      if (verbose) std::cout << "Writing to: " << outfile 
			     << " as HDF5" << std::endl;
      hid_t file_id;
      file_id = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);
      if (H5Iget_ref(file_id) < 0) {
	H5Fclose(file_id);
	throw affineExcept("pofd_affine_getdNdS", "getNProjected",
			   "Failed to open HDF5 file to write");
      }
      hsize_t adims;
      hid_t mems_id, att_id, group_id;
      hbool_t btmp;

      // Properties
      adims = 1;
      mems_id = H5Screate_simple(1, &adims, nullptr);
      btmp = static_cast<hbool_t>(logspace);
      att_id = H5Acreate2(file_id, "logspace", H5T_NATIVE_HBOOL,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_HBOOL, &btmp);
      H5Aclose(att_id);
      att_id = H5Acreate2(file_id, "dflux", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux2);
      H5Aclose(att_id);
      H5Sclose(mems_id);

      // Model
      group_id = H5Gcreate(file_id, "Model", H5P_DEFAULT, H5P_DEFAULT, 
			  H5P_DEFAULT);
      model.writeToHDF5Handle(group_id, true);
      H5Gclose(group_id);

      // Fluxes
      double* fluxes;
      fluxes = new double[nflux2];
      if (logspace) {
	double lmin = log2(minflux2);
	for (unsigned int i = 0; i < nflux2; ++i) 
	  fluxes[i] = exp2(lmin + static_cast<double>(i) * dflux2);
      } else
	for (unsigned int i = 0; i < nflux2; ++i) 
	  fluxes[i] = minflux2 + static_cast<double>(i) * dflux2;
      hdf5utils::writeDataDoubles(file_id, "Flux", nflux2, fluxes);
      delete[] fluxes;

      // dNdS
      if (logcounts)
	hdf5utils::writeDataDoubles(file_id, "Log10dNdS", nflux2, dNdS);
      else
	hdf5utils::writeDataDoubles(file_id, "dNdS", nflux2, dNdS);

      H5Fclose(file_id);
    } else if (oft == hdf5utils::TXT) {
      FILE *fp;
      fp = fopen(outfile.c_str(), "w");
      if (!fp)
	throw affineExcept("pofd_affine_getdNdS", "getNProjected",
			   "Failed to open text output file");
      if (logcounts)
	fprintf(fp, "%12s   %12s\n", "Flux", "Log10dNdS");
      else
	fprintf(fp, "%12s   %12s\n", "Flux", "dNdS");
      if (logspace) {
	double lmin = log2(minflux2);
	for (unsigned int i = 0; i < nflux2; ++i)
	  fprintf(fp, "%12.6e   %15.9e\n",
		  exp2(lmin + static_cast<double>(i) * dflux2), dNdS[i]);
      } else for (unsigned int i = 0; i < nflux2; ++i) 
	       fprintf(fp, "%12.6e   %15.9e\n",
		       minflux2 + static_cast<double>(i) * dflux2, dNdS[i]);
      fclose(fp);
    } else if (oft == hdf5utils::FITS)
	throw affineExcept("pofd_affine_getdNdS", "getNProjected",
			   "FITS output not supported");

    delete[] dNdS; 
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (dNdS != nullptr) delete[] dNdS;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (dNdS != nullptr) delete[] dNdS;
    return 16;
  }
  return 0;
}

//////////////////////////////////////////////

//Two-D version
int getNDouble(int argc, char** argv) {
  
  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  double minflux1, maxflux1, minflux2, maxflux2;
  unsigned int nflux1,nflux2;
  bool has_user_minflux1, has_user_minflux2;
  bool has_user_maxflux1, has_user_maxflux2;
  bool verbose, logspace, logcounts;

  minflux1          = 0;
  minflux2          = 0;
  maxflux1          = 1; //Will always be overridden
  maxflux2          = 1; //Will always be overridden
  nflux1            = 100;
  nflux2            = 100;
  has_user_minflux1 = false;
  has_user_minflux2 = false;
  has_user_maxflux1 = false;
  has_user_maxflux2 = false;
  logspace          = false;
  logcounts         = false;
  verbose           = false;

  int c;
  int option_index = 0;
  optind = 1; //Resets parsing
  while ((c = getopt_long(argc, argv, optstring, long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'l':
      logspace = true;
      break;
    case 'L':
      logcounts = true;
      break;
    case '1' :
      minflux1 = atof(optarg);
      has_user_minflux1 = true;
      break;
    case '2' :
      minflux2 = atof(optarg);
      has_user_minflux2 = true;
      break;
    case '3' :
      has_user_maxflux1 = true;
      maxflux1 = atof(optarg);
      break;
    case '4' :
      has_user_maxflux2 = true;
      maxflux2 = atof(optarg);
      break;
    case '5' :
      nflux1 = static_cast<unsigned int>(atoi(optarg));
      break;
    case '6' :
      nflux2 = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-1 ) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  initfile = std::string(argv[optind]);
  outfile = std::string(argv[optind + 1]);

  //Input check
  if (has_user_maxflux1 && (maxflux1 <= minflux1)) {
    std::cerr << "Maxflux must be > minflux in band 1" << std::endl;
    return 1;
  }
  if (has_user_maxflux2 && (maxflux2 <= minflux2)) {
    std::cerr << "Maxflux must be > minflux in band 2" << std::endl;
    return 2;
  }
  if (nflux1 == 0) {
    std::cerr << "nflux must be positive in band 1" << std::endl;
    return 3;
  }
  if (nflux2 == 0) {
    std::cerr << "nflux must be positive in band 2" << std::endl;
    return 4;
  }

  double *dNdS = nullptr;
  try {
    initFileDoubleLogNormal model_info(initfile, false, false);
    
    numberCountsDoubleLogNormal model;
    model_info.getModelPositions(model);

    paramSet pars(model_info.getNTot());
    model_info.getParams(pars);
    model.setParams(pars);


    //Check minflux values if logspace is set
    if (logspace) {
      if (!has_user_minflux1)
	minflux1 = model.getMinFlux().first; //Sorted inside model
      else
	if (minflux1 <= 0.0)
	  throw affineExcept("pofd_affine_getdNdS", "getNDouble",
			     "Minflux1 must be positive if using logspace");
      if (!has_user_minflux2)
	minflux2 = model.getMinFlux().second; //Sorted inside model
      else
	if (minflux2 <= 0.0)
	throw affineExcept("pofd_affine_getdNdS", "getNDouble",
			   "Minflux2 must be positive if using logspace");
    }
      
    //Get maximum flux if not set, use to compute dflux
    if (!has_user_maxflux1)
      maxflux1 = model.getMaxFlux().first;
    if (! has_user_maxflux2 )
      maxflux2 = model.getMaxFlux().second;

    double dflux1, dflux2;
    if (nflux1 > 1) {
      if (logspace)
	dflux1 = log2(maxflux1 / minflux1) / static_cast<double>(nflux1 - 1);
      else
	dflux1 = (maxflux1 - minflux1) / static_cast<double>(nflux1 - 1);
    } else
      dflux1 = 1.0; //Not used, but avoid compiler warning
    if (nflux1 > 1) {
      if (logspace)
	dflux2 = log2(maxflux2 / minflux2) / static_cast<double>(nflux2 - 1);
      else
	dflux2 = (maxflux2 - minflux2) / static_cast<double>(nflux2 - 1);
    } else
      dflux2 = 1.0;

    if (verbose) {
      std::cout << "Flux per sq degree, band 1: "
		<< model.getFluxPerArea(0) << std::endl;
      std::cout << "Flux per sq degree, band 2: "
		<< model.getFluxPerArea(1) << std::endl;
      std::cout << "Number of sources per area: "
		<< model.getNS() << std::endl;
    }

    //Calculation loop
    dNdS  = new double[nflux1 * nflux2];
    double flux1, flux2;
    if (logspace) {
      double lmin1 = log2(minflux1);
      double lmin2 = log2(minflux2);
      for (unsigned int i = 0; i < nflux1; ++i) {
	flux1 = exp2(lmin1 + static_cast<double>(i) * dflux1);
	for (unsigned int j = 0; j < nflux2; ++j) {
	  flux2 = exp2(lmin2 + static_cast<double>(j) * dflux2);
	  dNdS[i * nflux2 + j] = model.getNumberCounts(flux1, flux2);
	}
      }
    } else {
      for (unsigned int i = 0; i < nflux1; ++i) {
	flux1 = minflux1 + static_cast<double>(i) * dflux1;
	for (unsigned int j = 0; j < nflux2; ++j) {
	  flux2 = minflux2 + static_cast<double>(j) * dflux2;
	  dNdS[i * nflux2 + j] = model.getNumberCounts(flux1, flux2);
	}
      }
    }

    if (logcounts)
      for (unsigned int i = 0; i < nflux1 * nflux2; ++i)
	dNdS[i] = log10(dNdS[i]);
    
    //Write out
    hdf5utils::outfiletype oft = hdf5utils::getOutputFileType(outfile);
    if (oft == hdf5utils::HDF5 || oft == hdf5utils::UNKNOWN) {
      if (verbose) std::cout << "Writing to: " << outfile 
			     << " as HDF5" << std::endl;
      hid_t file_id;
      file_id = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);
      if (H5Iget_ref(file_id) < 0) {
	H5Fclose(file_id);
	throw affineExcept("pofd_affine_getdNdS", "getNDouble",
			   "Failed to open HDF5 file to write");
      }
      hsize_t adims;
      hid_t mems_id, att_id, group_id;

      // Properties
      adims = 1;
      mems_id = H5Screate_simple(1, &adims, nullptr);
      att_id = H5Acreate2(file_id, "dflux1", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux1);
      H5Aclose(att_id);
      att_id = H5Acreate2(file_id, "dflux2", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux2);
      H5Aclose(att_id);
      H5Sclose(mems_id);

      // Model
      group_id = H5Gcreate(file_id, "Model", H5P_DEFAULT, H5P_DEFAULT, 
			  H5P_DEFAULT);
      model.writeToHDF5Handle(group_id, true);
      H5Gclose(group_id);

      // Fluxes
      double* fluxes;
      fluxes = new double[nflux1];
      if (logspace) {
	double lmin1 = log2(minflux1);
    	for (unsigned int i = 0; i < nflux1; ++i)
	  fluxes[i] = exp2(lmin1 + static_cast<double>(i) * dflux1);
      } else
	for (unsigned int i = 0; i < nflux1; ++i) 
	  fluxes[i] = minflux1 + static_cast<double>(i) * dflux1;
      hdf5utils::writeDataDoubles(file_id, "Flux1", nflux1, fluxes);
      if (nflux1 != nflux2) {
	delete[] fluxes;
	fluxes = new double[nflux2];
      }
      if (logspace) {
	double lmin2 = log2(minflux2);
    	for (unsigned int i = 0; i < nflux2; ++i)
	  fluxes[i] = exp2(lmin2 + static_cast<double>(i) * dflux2);
      } else
	for (unsigned int i = 0; i < nflux2; ++i) 
	  fluxes[i] = minflux2 + static_cast<double>(i) * dflux2;
      hdf5utils::writeDataDoubles(file_id, "Flux2", nflux2, fluxes);
      delete[] fluxes;

      // dNdS
      if (logcounts)
	hdf5utils::writeData2DDoubles(file_id, "Log10dNdS", nflux1,
				      nflux2, dNdS);
      else
	hdf5utils::writeData2DDoubles(file_id, "dNdS", nflux1, nflux2, dNdS);

      H5Fclose(file_id);
    } else if (oft == hdf5utils::TXT) {
      FILE *fp;
      fp = fopen( outfile.c_str(),"w");
      if (!fp) 
	throw affineExcept("pofd_affine_getdNdS", "getNDouble",
			   "Failed to open text output file");
      fprintf(fp, "Log10Counts: %s\n", logcounts ? "true" : "false");
      fprintf(fp, "%s %u %12.6e %12.6e\n",logspace ? "true" : "false",
	      nflux1, minflux1, dflux1);
      fprintf(fp, "%s %u %12.6e %12.6e\n", logspace ? "true" : "false",
	      nflux2, minflux2, dflux2);
      for (unsigned int i = 0; i < nflux1; ++i) {
	for (unsigned int j = 0; j < nflux2-1; ++j)
	  fprintf(fp, "%13.7e ", dNdS[i * nflux2 + j]);
	fprintf(fp, "%13.7e\n", dNdS[i * nflux2 + nflux2 - 1]);
      }
      fclose(fp);
    } else if (oft == hdf5utils::FITS)
      throw affineExcept("pofd_affine_getdNdS", "getNDouble",
			 "FITS output not supported");

    delete[] dNdS; 
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (dNdS != nullptr) delete[] dNdS;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (dNdS != nullptr) delete[] dNdS;
    return 16;
  }
  return 0;
}

/////////////////////////////////////////////

int main( int argc, char** argv ) {
  bool twod, projected;

  twod = false;
  projected = false;

  //Only interested in a) displaying help and b) figuring out
  // if this is 1D or 2D c) displaying the version number
  int c;
  int option_index = 0;
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getdNdS -- get differential number counts "
		<< "for a model." << std::endl;
      std::cerr << "\tBoth one-dimensional and two-dimensional models are "
		<< "supported." << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\t pofd_affine_getdNdS [options] initfile outfile"
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates the differential number counts for the model "
		<< "in initfile" << std::endl;
      std::cerr << "\tand writes it to outfile.  The 1D model is a log-space "
		<< "spline" << std::endl;
      std::cerr << "\tmodel for the number counts, and the 2D model is the 1D" 
		<< " spline" << std::endl;
      std::cerr << "\tmodel times a log-normal color function for the second"
		<< " band," << std::endl;
      std::cerr << "\twith the log-space variance and mean color stored as"
		<< " splines" << std::endl;
      std::cerr << "\tin the flux of the first band." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tIn the 2D case, you can also ask for the projected"
		<< " band 2" << std::endl;
      std::cerr << "\tmodel." << std::endl;
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
      std::cerr << std::endl;
      std::cerr << "\tFor the 2D (or projected 2D) case, initfile is a text "
		<< "file" << std::endl;
      std::cerr << "\tgiving the positions of the knot points and their "
		<< "values, " << std::endl;
      std::cerr << "\tfollowed by the sigma knot positions and their values, "
		<< "then" << std::endl;
      std::cerr << "\tlikewise for the colour offset. The format is three "
		<< "numbers" << std::endl;
      std::cerr << "\ton the first line, giving the number of number count "
		<< "knots," << std::endl;
      std::cerr << "\tsigma knots, and offset knots, followed by a number of "
		<< "lines" << std::endl;
      std::cerr << "\tagain with the format knotpos value.  The sigmas and "
		<< "offsets" << std::endl;
      std::cerr << "\tare in log space." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tIn all cases the output number counts are written to"
		<< " outfile." << std::endl;
      std::cerr << "\tThe file output type is controlled by the output file"
		<< std::endl;
      std::cerr << "\textension (h5 for HDF5, txt for text)." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-d, --double" << std::endl;
      std::cerr << "\t\tUse the 2D model." << std::endl;
      std::cerr << "\t-p, --projected" << std::endl;
      std::cerr << "\t\tUse the 2D model, but project out the band 1 counts."
		<< " In" << std::endl;
      std::cerr << "\t\tother words, produce the band 2 counts from the 2-band"
		<< std::endl;
      std::cerr << "\t\tmodel." << std::endl;
      std::cerr << "\t-l, --logspace" << std::endl;
      std::cerr << "\t\tSpace the fluxes in log space rather than linearly"
		<< std::endl;
      std::cerr << "\t-L, --logcounts" << std::endl;
      std::cerr << "\t\tOutput Log10 dNdS instead of dNdS" << std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tRun in verbose mode" << std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput the version number and exit." << std::endl;
      std::cerr << std::endl;
      std::cerr << "ONE-D OPTIONS" << std::endl;
      std::cerr << "\t-M, --maxflux value" << std::endl;
      std::cerr << "\t\tThe maximum flux to calculate dN/dS for (def: highest"
		<< std::endl;
      std::cerr << "\t\tknot specified in initfile)." << std::endl;
      std::cerr << "\t-m, --minflux value" << std::endl;
      std::cerr << "\t\tThe minimum flux to calculate dN/dS for (def: 0, "
		<< "unless" << std::endl;
      std::cerr << "\t\t--logspace is set)." << std::endl;
      std::cerr << "\t-n, --nflux value" << std::endl;
      std::cerr << "\t\tNumber of flux values to output (def: 1000)" 
		<< std::endl;
      std::cerr << "PROJECTED OPTIONS" << std::endl;
      std::cerr << "\t-M, --maxflux value" << std::endl;
      std::cerr << "\t\tThe maximum band 2 flux to calculate dN/dS for"
		<< std::endl;
      std::cerr << "\t\t(def: a few sigma above mean value at highest band 1 "
		<< "knot" << std::endl;
      std::cerr << "\t\tspecified in initfile)." << std::endl;
      std::cerr << "\t-m, --minflux value" << std::endl;
      std::cerr << "\t\tThe minimum flux to calculate dN/dS for (def: 0, "
		<< "unless" << std::endl;
      std::cerr << "\t\t--logspace is set)." << std::endl;
      std::cerr << "\t-n, --nflux value" << std::endl;
      std::cerr << "\t\tNumber of flux values to output (def: 100)" 
		<< std::endl;
      std::cerr << "TWO-D OPTIONS" << std::endl;
      std::cerr << "\t--maxflux1 value" << std::endl;
      std::cerr << "\t\tThe maximum flux to calculate dN/dS for in band 1"
		<< std::endl;
      std::cerr << "\t\t(def: highest knot specified in initfile)." 
		<< std::endl;
      std::cerr << "\t--maxflux2 value" << std::endl;
      std::cerr << "\t\tThe maximum flux to calculate dN/dS for in band 2"
		<< std::endl;
      std::cerr << "\t\t(def: a few sigma above mean value at highest band 1 "
		<< "knot)." << std::endl;
      std::cerr << "\t--minflux1 value" << std::endl;
      std::cerr << "\t\tThe minimum flux to calculate dN/dS for in band 1. "
		<< "(def:" << std::endl;
      std::cerr << "\t\tlowest band 1 knot specified in initfile)" 
		<< std::endl;
      std::cerr << "\t--minflux2 value" << std::endl;
      std::cerr << "\t\tThe minimum flux to calculate dN/dS for in band 2. "
		<< "(def: 0)" << std::endl;
      std::cerr << "\t--nflux1 value" << std::endl;
      std::cerr << "\t\tNumber of flux values to output, band 1 (def: 100)" 
		<< std::endl;
      std::cerr << "\t--nflux2 value" << std::endl;
      std::cerr << "\t\tNumber of flux values to output, band 2 (def: 100)" 
		<< std::endl;
      return 0;
      break;
    case 'd' :
      twod = true;
      break;
    case 'p':
      projected = true;
      break;
    case 'V' :
      std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
		<< std::endl;
      return 0;
      break;
    }
  
  if (twod) {
    if (projected)
      std::cout << "WARNING: Both twod and projected set; just using twod"
		<< std::endl;
    return getNDouble(argc, argv);
  } else if (projected) 
    return getNProjected(argc, argv);
  else
    return getNSingle(argc, argv);

}
