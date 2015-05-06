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
  nknots(0), npoints(NPOINTS) {
  knots = nullptr;
  uncertainty = nullptr;
  s = new float[npoints];
  mean = new float[npoints];
  median = new float[npoints];
  snake = new float[2 * nprob * npoints];
}

errorSnakeSingle::~errorSnakeSingle() {
  if (knots != nullptr) delete[] knots;
  if (uncertainty != nullptr) delete[] uncertainty;
  delete[] s;
  delete[] mean;
  delete[] median;
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

  nknots = model.getNKnots();
  if (knots != nullptr) delete[] knots;
  knots = new float[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = model.getKnotPos(i);
  if (uncertainty != nullptr) delete[] uncertainty;
  uncertainty = new float[2 * nprob * nknots];
  
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

  // We need some random numbers
  ran r;
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(NULL));
  seed += static_cast<unsigned long long int>(clock());
  r.setSeed(seed);

  // We will interpolate to get the ranges.  So it's convenient
  //  to pre-calculate the effective index in nsamples after we sort.
  float lowidx[nprob], highidx[nprob]; // Fractional index of interpolation
  for (unsigned int i = 0; i < nprob; ++i)
    lowidx[i] = 0.5 * (1.0 - plevels[i]) * nsamples;
  for (unsigned int i = 0; i < nprob; ++i)
    highidx[i] = 0.5 * (1 + plevels[i]) * nsamples;
  
  // Uncertainties at knots
  // Need some working space
  unsigned int idx1, idx2, sidx;
  float *working = new float[nknots * nsamples];
  for (unsigned int i = 0; i < nsamples; ++i) {
    // Grab a random set of parameters
    idx1 = r.selectFromRange(0, nwalkers);
    idx2 = r.selectFromRange(0, nsteps);
    sidx = (idx1 * nsteps + idx2) * nparams;
    for (unsigned int j = 0; j < nknots; ++j)
      working[j * nsamples + i] = steps[sidx + j];
  }

  // Do statistics
  float *wptr;
  for (unsigned int j = 0; j < nknots; ++j) {
    wptr = working + j * nsamples;
    // Sort so we can get the percentiles.  This could be slightly
    //  more efficient if we used O(n) selection techniques instead
    //  of O(n ln n) sorting, but it's easily fast enough
    std::sort(wptr, wptr + nsamples);

    for (unsigned int k = 0; k < nprob; ++k) {
      sidx = 2 * (k * nknots + j);
      // Lower value
      unsigned int idxl = static_cast<unsigned int>(lowidx[k]);
      uncertainty[sidx] = wptr[idxl] + (wptr[idxl+1] - wptr[idxl]) *
	(lowidx[k] - idxl);
      // And upper
      unsigned int idxu = static_cast<unsigned int>(highidx[k]);
      uncertainty[sidx + 1] =
	wptr[idxu] + (wptr[idxu+1] - wptr[idxu]) * (highidx[k] - idxu);
    }
  }
  delete[] working;
  
  // Now do more continuous curve.  Same thing, pretty much,
  //  except that nknots -> npoints, and we also keep the mean/median
  working = new float[npoints * nsamples];
  paramSet p(nknots);
  for (unsigned int i = 0; i < nsamples; ++i) {
    idx1 = r.selectFromRange(0, nwalkers);
    idx2 = r.selectFromRange(0, nsteps);
    sidx = (idx1 * nsteps + idx2) * nparams;
    p.setParamValues(nknots, &steps[sidx]);
    model.setParams(p);
    // Convert to log10 so that in same units as params
    for (unsigned int j = 0; j < npoints; ++j)
      working[j * nsamples + i] = log10(model.getNumberCounts(s[j]));    
  }

  double mnval;
  for (unsigned int j = 0; j < npoints; ++j) {
    // Three things to do: the mean, the median, the error snakes
    wptr = working + j * nsamples;
    std::sort(wptr, wptr + nsamples);

    // Mean
    mnval = wptr[0];
    for (unsigned int k = 0; k < nsamples; ++k)
      mnval += wptr[k];
    mean[j] = static_cast<float>(mnval / static_cast<double>(nsamples));

    // Median; should really interpolate, but meh, assume even
    // The realistic error introduced by this is pretty tiny
    median[j] = wptr[nsamples / 2];
    
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

  hdf5utils::writeDataFloats(file_id, "KnotPositions", nknots, knots);
  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "KnotValuesp=" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(file_id, pname.str(),
				 nknots, 2, &uncertainty[i * nknots * 2]);
  }
  
  hdf5utils::writeDataFloats(file_id, "Flux", npoints, s);
  hdf5utils::writeDataFloats(file_id, "MeanLog10Counts", npoints, mean);
  hdf5utils::writeDataFloats(file_id, "MedianLog10Counts", npoints, median);
  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "Log10CountsSnakep=" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(file_id, pname.str(),
				 npoints, 2, &snake[i * npoints * 2]);
  }

  H5Fclose(file_id);
}
