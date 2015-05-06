#include<ctime>
#include<sstream>
#include<iomanip>

#include "../include/errorSnakeSingle.h"
#include "../include/affineExcept.h"
#include "../include/paramSet.h"
#include "../include/hdf5utils.h"
#include "../include/ran.h"

/*!
  \param[in] npoints  Number of flux density points
*/
errorSnakeSingle::errorSnakeSingle(unsigned int NPOINTS):
  npoints(NPOINTS) {
  s = new float[npoints];
  snake = new float[2 * nprob * npoints];
}

errorSnakeSingle::~errorSnakeSingle() {
  delete[] s;
  delete[] snake;
}

/*!
  \param[in] infile Input file to read; must be output of pofd_affine_mcmc
  \param[in] nsamples Number of samples to use.  Should be > 1000
  \param[in] logspace Logarithmically space the interpolant points
*/
void errorSnakeSingle::build(const std::string& infile,
			     unsigned int nsamples, bool logspace) {
  if (nsamples < 1000)
    throw affineExcept("errorSnakeSingle", "build",
		       "Number of samples too small to build all snakes");

  hid_t fileid = H5Fopen(infile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Read the model
  hid_t groupid = H5Gopen(fileid, "LikelihoodParams/Model", H5P_DEFAULT);
  model.readFromHDF5Handle(groupid);
  H5Gclose(groupid);

  // Read in chains
  hid_t stepsid = H5Dopen(fileid, "Chains/Chain", H5P_DEFAULT);
  hid_t spaceid = H5Dget_space(stepsid);
  
  int ndims = H5Sget_simple_extent_ndims(spaceid);
  if (ndims != 3) {
    H5Sclose(spaceid);
    H5Dclose(stepsid);
    H5Fclose(fileid);
    throw affineExcept("errorSnakeSingle", "build",
		       "Chains of unexpected dimensionality");
  }
  hsize_t datasize[3];
  hsize_t maxdatasize[3];
  H5Sget_simple_extent_dims(spaceid, datasize, maxdatasize);
  unsigned int nwalkers = static_cast<unsigned int>(datasize[0]);
  unsigned int nsteps = static_cast<unsigned int>(datasize[1]);
  unsigned int nparams = static_cast<unsigned int>(datasize[2]);

  // Allocate and read. Flattened because HDF5 works with
  //  flattened arrays
  float *steps = new float[nwalkers * nsteps * nparams];
  H5Dread(stepsid, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
	  H5P_DEFAULT, steps);
  H5Sclose(spaceid);
  H5Dclose(stepsid);
  H5Fclose(fileid);

  // Set up flux densities
  double minflux = model.getMinFlux();
  double maxflux = model.getMaxFlux();
  maxflux -= 1e-4 * (maxflux - minflux); // maxflux counts are 0, avoid
  if (logspace) {
    double log_minflux = log(minflux);
    double dflux = (log(maxflux) - log_minflux) /
      static_cast<double>(npoints - 1);
    for (unsigned int i = 0; i < npoints; ++i)
      s[i] = exp(log_minflux + dflux * static_cast<double>(i));
  } else {
    double dflux = (maxflux - minflux) / static_cast<double>(npoints - 1);
    for (unsigned int i = 0; i < npoints; ++i)
      s[i] = minflux + dflux * static_cast<double>(i);
  }

  // Working space, npoints by nsamples
  float *working;
  working = new float[npoints * nsamples];

  // Need some random numbers
  ran r;
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(NULL));
  seed += static_cast<unsigned long long int>(clock());
  r.setSeed(seed);
  
  // Main evaluation loop
  unsigned int nknots = model.getNKnots();
  paramSet p(nknots);
  unsigned int idx1, idx2, sidx;
  for (unsigned int i = 0; i < nsamples; ++i) {
    // Grab a random set of parameters
    idx1 = r.selectFromRange(0, nwalkers);
    idx2 = r.selectFromRange(0, nsteps);
    sidx = (idx1 * nsteps + idx2) * nparams;
    p.setParamValues(nknots, &steps[sidx]);
    model.setParams(p);
    // Convert to log10 so that in same units as params
    for (unsigned int j = 0; j < npoints; ++j)
      working[j * nsamples + i] = log10(model.getNumberCounts(s[j]));    
  }

  // Build statistics at each point, linearly interpolating in
  // probability array
  float *wptr;
  float lowidx[nprob], highidx[nprob]; // Fractional index of interpolation
  for (unsigned int i = 0; i < nprob; ++i)
    lowidx[i] = 0.5 * (1.0 - plevels[i]) * nsamples;
  for (unsigned int i = 0; i < nprob; ++i)
    highidx[i] = 0.5 * (1 + plevels[i]) * nsamples;
  for (unsigned int j = 0; j < npoints; ++j) {
    wptr = working + j * nsamples;
    std::sort(wptr, wptr + nsamples);
    for (unsigned int k = 0; k < nprob; ++k) {
      sidx = 2 * (k * npoints + j);
      // Lower value
      unsigned int idxl = static_cast<unsigned int>(lowidx[k]);
      snake[sidx] = wptr[idxl] + (wptr[idxl+1] - wptr[idxl]) *
	(lowidx[k] - idxl);
      // And upper
      unsigned int idxu = static_cast<unsigned int>(highidx[k]);
      snake[sidx + 1] =
	wptr[idxu] + (wptr[idxu+1] - wptr[idxu]) * (highidx[k] - idxu);
    }
  }

  // Clean up and go home
  delete[] working;
  delete[] steps;
  
}

/*!
  \param[in] filename File to write error snake to.
*/
void errorSnakeSingle::writeAsHDF5(const std::string& filename) const {

  hid_t file_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("errorSnakeSingle", "writeToHDF5",
		       "Failed to open HDF5 file to write");
  }

  hdf5utils::writeDataFloats(file_id, "Flux", npoints, s);
  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "Snake" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(file_id, pname.str(),
				 npoints, 2, &snake[i * npoints * 2]);
  }

  H5Fclose(file_id);
}
