#include<iostream>

#include<getopt.h>
#include<hdf5.h>

#include "../include/global_settings.h"
#include "../include/beam.h"
#include "../include/doublebeam.h"
#include "../include/numberCountsKnotsSpline.h"
#include "../include/numberCountsDoubleLogNormal.h"
#include "../include/paramSet.h"
#include "../include/utility.h"
#include "../include/affineExcept.h"

//All sub-parses have to have the same options to avoid
// getting warnings.  So we give them all the same long_options,
// but then only process the appropriate ones, ignoring the rest
//Yes, this is a bit complicated and error prone, but such is life
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"double",no_argument, 0, 'd'},
  {"verbose",no_argument, 0, 'v'},
  {"version",no_argument, 0, 'V'}, //Below here, not parsed in main routine
  {"hdf5", no_argument, 0, '1'},
  {"histogram",no_argument, 0, 'H'},
  {0, 0, 0, 0}
};
char optstring[] = "dh1HvV";

//One-D version
int getRSingle( int argc, char** argv ) {
  
  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  std::string psffile; //Beam file
  bool histogram, verbose, write_as_hdf5; 
  double minflux, maxflux;
  unsigned int nflux;

  histogram     = false;
  verbose       = false;
  write_as_hdf5 = false;

  int c;
  int option_index = 0;
  optind = 1; //!< Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case '1':
      write_as_hdf5 = true;
      break;
    case 'H' :
      histogram = true;
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 5) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  minflux = atof(argv[optind]);
  maxflux = atof(argv[optind + 1]);
  nflux = static_cast<unsigned int>(atoi(argv[optind + 2]));
  initfile = std::string(argv[optind + 3]);
  psffile = std::string(argv[optind + 4]);
  outfile = std::string(argv[optind + 5]);
  
  double dflux;
  if (nflux > 1)
    dflux = (maxflux - minflux)/ static_cast<double>(nflux - 1);
  else
    dflux = 1.0;

  double *fluxes = NULL;
  double *R = NULL;
  try {
    initFileKnots model_info(initfile, false, false);

    numberCountsKnotsSpline model;
    model_info.getKnotPos(model);

    beam bm(psffile, histogram);
    paramSet pars(model_info.getNKnots());
    model_info.getParams(pars);
    model.setParams(pars);

    if (verbose) {
      std::cout << "Flux per area: " << model.getFluxPerArea()
		<< std::endl;
      std::cout << "Beam area: " << bm.getEffectiveArea() << std::endl;
      if (bm.hasPos()) std::cout << " Beam has positive components" 
				 << std::endl;
      if (bm.hasNeg()) std::cout << " Beam has negative components" 
				 << std::endl;
    }
    fluxes = new double[nflux];
    for (unsigned int i = 0; i < nflux; ++i)
      fluxes[i] = minflux + static_cast<double>(i) * dflux;

    R = new double[nflux];
    model.getR(nflux, fluxes, bm, R);
    
    if (write_as_hdf5) {
      if (verbose) std::cout << "Writing to: " << outfile 
			     << " as HDF5" << std::endl;
      hid_t file_id;
      file_id = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);
      if (H5Iget_ref(file_id) < 0) {
	H5Fclose(file_id);
	throw affineExcept("pofd_affine_getR", "pofd_affine_getR",
			   "Failed to open HDF5 file to write");
      }
      hsize_t adims;
      hid_t mems_id, att_id, group_id, dat_id;
      
      // Properties
      adims = 1;
      mems_id = H5Screate_simple(1, &adims, NULL);
      att_id = H5Acreate2(file_id, "dflux", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux);
      H5Aclose(att_id);
      H5Sclose(mems_id);

      // Model
      group_id = H5Gcreate(file_id, "Model", H5P_DEFAULT, H5P_DEFAULT, 
			  H5P_DEFAULT);
      model.writeToHDF5Handle(group_id);
      H5Gclose(group_id);

      // Rflux
      adims = nflux;
      mems_id = H5Screate_simple(1, &adims, NULL);
      dat_id = H5Dcreate2(file_id, "Flux", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, fluxes);
      H5Dclose(dat_id);

      dat_id = H5Dcreate2(file_id, "R", H5T_NATIVE_DOUBLE, mems_id,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, R);
      H5Dclose(dat_id);
      H5Sclose(mems_id);

      H5Fclose(file_id);
    } else {
      // Text file
      if (verbose) std::cout << "Writing to: " << outfile << std::endl;
      FILE *fp;
      fp = fopen(outfile.c_str(), "w");
      if (!fp) {
	std::cerr << "Failed to open output file" << std::endl;
	return 128;
      }
      fprintf(fp, "#%-11s   %-12s\n", "Flux", "R");
      for (unsigned int i = 0; i < nflux; ++i) 
	fprintf(fp, "%12.6e   %15.9e\n", fluxes[i], R[i]);
      fclose(fp);
    }
    delete[] fluxes;
    delete[] R; 
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (fluxes != NULL) delete[] fluxes;
    if (R != NULL) delete[] R;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (fluxes != NULL) delete[] fluxes;
    if (R != NULL) delete[] R;
    return 16;
  }
  return 0;
}

int getRDouble(int argc, char** argv) {
  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  std::string psffile1, psffile2; //Beam file
  bool histogram, verbose, write_as_hdf5;
  double minflux1, maxflux1, minflux2, maxflux2;
  unsigned int nflux1, nflux2;

  histogram     = false;
  verbose       = false;
  write_as_hdf5 = false;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case '1':
      write_as_hdf5 = true;
      break;
    case 'H' :
      histogram = true;
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 9) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  minflux1 = atof(argv[optind]);
  maxflux1 = atof(argv[optind + 1]);
  nflux1   = static_cast<unsigned int>(atoi(argv[optind + 2]));
  minflux2 = atof(argv[optind + 3]);
  maxflux2 = atof(argv[optind + 4]);
  nflux2   = static_cast<unsigned int>(atoi(argv[optind + 5]));
  initfile = std::string(argv[optind + 6]);
  psffile1 = std::string(argv[optind + 7]);
  psffile2 = std::string(argv[optind + 8]);
  outfile  = std::string(argv[optind + 9]);

  if (nflux1 == 0) {
    std::cerr << "Invalid (non-positive) nflux1" << std::endl;
    return 1;
  }
  if (nflux2 == 0) {
    std::cerr << "Invalid (non-positive) nflux2" << std::endl;
    return 1;
  }
  double dflux1, dflux2;
  if (maxflux1 < minflux1) std::swap(minflux1, maxflux1);
  if (nflux1 > 1) 
    dflux1 = (maxflux1 - minflux1) / static_cast<double>(nflux1 - 1);
  else
    dflux1 = maxflux1 - minflux1;
  if (maxflux2 < minflux2) std::swap(minflux2, maxflux2);
  if (nflux2 > 1)
    dflux2 = (maxflux2 - minflux2) / static_cast<double>(nflux2 - 1);
  else
    dflux2 = maxflux2 - minflux2;
  
  //Main computation
  double *fluxes1, *fluxes2, *R;
  fluxes1 = fluxes2 = R = NULL;
  try {
    initFileDoubleLogNormal model_info(initfile, false, false);

    numberCountsDoubleLogNormal model;
    model_info.getModelPositions(model);

    doublebeam bm(psffile1, psffile2, histogram);
    paramSet pars(model_info.getNTot());
    model_info.getParams(pars);
    model.setParams(pars);

    if (verbose) {
      std::cout << "Flux per area, band 1: " << model.getFluxPerArea(0)
		<< std::endl;
      std::cout << "Beam area, band 1: " << bm.getEffectiveArea1() << std::endl;
      std::cout << "Flux per area, band 2: " << model.getFluxPerArea(1)
		<< std::endl;
      std::cout << "Beam area, band 2: " << bm.getEffectiveArea2() << std::endl;
    }
    
    std::cout << "R(0.003, 0.003): "
	      << model.getR(0.003, 0.003, bm) << std::endl;
    exit(1);

    fluxes1 = new double[nflux1];
    for (unsigned int i = 0; i < nflux1; ++i)
      fluxes1[i] = minflux1 + static_cast<double>(i) * dflux1;
    fluxes2 = new double[nflux2];
    for (unsigned int i = 0; i < nflux2; ++i)
      fluxes2[i] = minflux2 + static_cast<double>(i) * dflux2;

    R = new double[nflux1 * nflux2];
    model.getR(nflux1, fluxes1, nflux2, fluxes2, bm, R);
    
    if (write_as_hdf5) {
      if (verbose) std::cout << "Writing to: " << outfile 
			     << " as HDF5" << std::endl;
      hid_t file_id;
      file_id = H5Fcreate(outfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			  H5P_DEFAULT);
      if (H5Iget_ref(file_id) < 0) {
	H5Fclose(file_id);
	throw affineExcept("pofd_affine_getR", "pofd_affine_getR",
			   "Failed to open HDF5 file to write");
      }
      hsize_t adims;
      hid_t mems_id, att_id, dat_id, group_id;
      
      // Properties
      adims = 1;
      mems_id = H5Screate_simple(1, &adims, NULL);
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
      model.writeToHDF5Handle(group_id);
      H5Gclose(group_id);

      // Fluxes
      adims = nflux1;
      mems_id = H5Screate_simple(1, &adims, NULL);
      dat_id = H5Dcreate2(file_id, "Flux1", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, fluxes1);
      H5Dclose(dat_id);
      H5Sclose(mems_id);
      adims = nflux2;
      mems_id = H5Screate_simple(1, &adims, NULL);
      dat_id = H5Dcreate2(file_id, "Flux2", H5T_NATIVE_DOUBLE,
			  mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, fluxes2);
      H5Dclose(dat_id);
      H5Sclose(mems_id);
      
      // R, which is 2D
      hsize_t dims_steps[2] = {nflux1, nflux2};
      mems_id = H5Screate_simple(2, dims_steps, NULL);
      dat_id = H5Dcreate2(file_id, "R", H5T_NATIVE_DOUBLE, mems_id,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	       H5P_DEFAULT, R);			  
      H5Dclose(dat_id);
      H5Sclose(mems_id);

      H5Fclose(file_id);
    } else {
      // Text file
      if (verbose) std::cout << "Writing to: " << outfile << std::endl;

      FILE *fp;
      fp = fopen( outfile.c_str(),"w");
      if (!fp) {
	std::cerr << "Failed to open output file" << std::endl;
	return 128;
      }
      fprintf(fp, "#%4u %4u\n", nflux1, nflux2);
      fprintf(fp, "#minflux1: %12.6e dflux1: %12.6e\n", minflux1, dflux1);
      fprintf(fp, "#minflux2: %12.6e dflux2: %12.6e\n", minflux2, dflux2);
      for (unsigned int i = 0; i < nflux1; ++i) {
	for (unsigned int j = 0; j < nflux2-1; ++j)
	  fprintf(fp, "%13.7e ", R[nflux2*i+j]);
	fprintf(fp, "%13.7e\n", R[nflux2*i+nflux2-1]);
      }
      fclose(fp);
    }

    delete[] fluxes1;
    delete[] fluxes2;
    delete[] R; 
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (fluxes1 != NULL) delete[] fluxes1;
    if (fluxes2 != NULL) delete[] fluxes2;
    if (R != NULL) delete[] R;
    return 16;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (fluxes1 != NULL) delete[] fluxes1;
    if (fluxes2 != NULL) delete[] fluxes2;
    if (R != NULL) delete[] R;
    return 16;
  }

  return 0;

}

////////////////////////////////////////

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
      std::cerr << "\tpofd_affine_getR -- get R for a number counts model."
		<< "  Both" << std::endl;
      std::cerr << "\tone-dimensional and two-dimensional models are supported."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tEither" << std::endl;
      std::cerr << std::endl;

      std::cerr << "\t pofd_affine_getR [options] minflux maxflux nflux initfile"
		<< std::endl;
      std::cerr<< "\t  beamfile outfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfor the 1D case or" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\t pofd_affine_getR -d [options] minflux1 maxflux1 nflux1"
		<< std::endl;
      std::cerr << "\t  minflux2 manflux2 nflux2 initfile beamfile1 beamfile2"
		<< std::endl;
      std::cerr << "\t  outfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfor the 2D case." << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates R for the model in initfile using the P(D)"
		<< " formalism and" << std::endl;
      std::cerr << "\twrites it to outfile.  The 1D model is a log-space "
		<< "spline" << std::endl;
      std::cerr << "\tmodel for the number counts, and the 2D model is the 1D" 
		<< " spline" << std::endl;
      std::cerr << "\tmodel times a log-normal color function for the second"
		<< " band," << std::endl;
      std::cerr << "\twith the log-space variance and mean color stored as"
		<< " splines" << std::endl;
      std::cerr << "\tin the flux of the first band." << std::endl;
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
      std::cerr << "\tminflux, maxflux, and nflux give the minimum, maximum,"
		<< " and" << std::endl;
      std::cerr << "\tnumber of fluxes to evaluate R for in the 1D case.  For"
		<< " the 2D" << std::endl;
      std::cerr << "\tcase this is extended to the minimum, maximum, and number"
		<< " of" << std::endl;
      std::cerr << "\tfluxes in each of the two bands.  Similarly, beamfile"
		<< " gives the" << std::endl;
      std::cerr << "\tname of a FITS file containing the beam in the 1D case,"
		<< " and" << std::endl;
      std::cerr << "\tbeamfile1, beamfile2 give the beam in each of the two "
		<< "bands" << std::endl;
      std::cerr << "\tin the 2D case." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tIn both cases the output R is written to outfile as"
		<< " text." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-d, --double" << std::endl;
      std::cerr << "\t\tUse the 2D model." << std::endl;
      std::cerr << "\t-H, --histogram" << std::endl;
      std::cerr << "\t\tUse beam histogramming." << std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput the version number and exit." << std::endl;
      std::cerr << std::endl;
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

  if (! twod)
    return getRSingle(argc,argv);
  else
    return getRDouble(argc,argv);
}
