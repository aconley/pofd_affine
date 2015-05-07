#include<ctime>
#include<sstream>
#include<iomanip>

#include "../include/numberCountsKnotsSplineError.h"
#include "../include/affineExcept.h"
#include "../include/paramSet.h"
#include "../include/hdf5utils.h"

/*!
  \param[in] npoints  Number of flux density points
*/
numberCountsKnotsSplineError::
numberCountsKnotsSplineError(unsigned int NPOINTS):
  nparams(0), nknots(0), knotPositions(nullptr), 
  bestParams(nullptr), meanParams(nullptr), 
  medParams(nullptr), uncertaintyParams(nullptr), 
  npoints(NPOINTS) {
  s = new float[npoints];
  mean = new float[npoints];
  median = new float[npoints];
  snake = new float[2 * nprob * npoints];
}

numberCountsKnotsSplineError::~numberCountsKnotsSplineError() {
  if (knotPositions != nullptr) delete[] knotPositions;
  if (bestParams != nullptr) delete[] bestParams;
  if (meanParams != nullptr) delete[] meanParams;
  if (medParams != nullptr) delete[] medParams;
  if (uncertaintyParams != nullptr) delete[] uncertaintyParams;
  delete[] s;
  delete[] mean;
  delete[] median;
  delete[] snake;
}

/*!
  \param[in] n New number of parameters
*/
void numberCountsKnotsSplineError::setNParams(unsigned int n) {
  if (n == nparams) return;
  if (bestParams != nullptr) delete[] bestParams;
  bestParams = new float[n];
  if (meanParams != nullptr) delete[] meanParams;
  meanParams = new float[n];
  if (medParams != nullptr) delete[] medParams;
  medParams = new float[n];
  if (uncertaintyParams != nullptr) delete[] uncertaintyParams;
  uncertaintyParams = new float[2 * nprob * n];
  nparams = n;
 }

/*!
  \param[in] infile Input file to read; must be output of pofd_affine_mcmc
  \param[in] thin Thinning to use on inputs
  \param[in] logspace Logarithmically space the interpolant points
*/
void numberCountsKnotsSplineError::build(const std::string& infile,
                                         unsigned int thin, bool logspace) {

  hid_t fileid = H5Fopen(infile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Read the model
  hid_t groupid = H5Gopen(fileid, "LikelihoodParams/Model", H5P_DEFAULT);
  model.readFromHDF5Handle(groupid);
  H5Gclose(groupid);

  // Get knot positions
  nknots = model.getNKnots();
  if (knotPositions != nullptr) delete[] knotPositions;
  knotPositions = new float[nknots];
  for (unsigned int i = 0; i < nknots; ++i)
    knotPositions[i] = model.getKnotPos(i);

  // Read in chains
  hid_t stepsid = H5Dopen(fileid, "Chains/Chain", H5P_DEFAULT);
  hid_t spaceid = H5Dget_space(stepsid);
  
  int ndims = H5Sget_simple_extent_ndims(spaceid);
  if (ndims != 3) {
    H5Sclose(spaceid);
    H5Dclose(stepsid);
    H5Fclose(fileid);
    throw affineExcept("numberCountsKnotsSplineError", "build",
		                   "Chains of unexpected dimensionality");
  }
  hsize_t datasize[3];
  hsize_t maxdatasize[3];
  H5Sget_simple_extent_dims(spaceid, datasize, maxdatasize);
  unsigned int nwalkers = static_cast<unsigned int>(datasize[0]);
  unsigned int nsteps = static_cast<unsigned int>(datasize[1]);
  unsigned int fnparams = static_cast<unsigned int>(datasize[2]);

  setNParams(fnparams);

  // Get names of parameters
  groupid = H5Gopen(fileid, "ParamInfo", H5P_DEFAULT);
  paramNames = hdf5utils::readAttStrings(groupid, "ParamNames");
  H5Gclose(groupid);

  // Allocate and read steps. Flattened because HDF5 works with
  //  flattened arrays
  float *steps = new float[nwalkers * nsteps * nparams];
  H5Dread(stepsid, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
	        H5P_DEFAULT, steps);
  H5Sclose(spaceid);
  H5Dclose(stepsid);

  // And likelihood, will be nwalkers by nsteps
  hid_t likeid = H5Dopen(fileid, "Chains/Likelihood", H5P_DEFAULT);
  double *like = new double[nwalkers * nsteps];
  H5Dread(likeid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
          H5P_DEFAULT, like);
  H5Dclose(likeid);

  // Done with input
  H5Fclose(fileid);

  // Deal with parameter uncertainties
  // Get best fit point
  bestLike = like[0];
  unsigned int bestLike_loc = 0;
  for (unsigned int i = 0; i < nwalkers * nsteps; ++i)
    if (like[i] > bestLike) {
      bestLike = like[i];
      bestLike_loc = i;
    }
  for (unsigned int i = 0; i < nparams; ++i)
    bestParams[i] = steps[bestLike_loc * nparams + i];
  // Done with this
  delete[] like;

  // Uncertainties of params
  // Need some working space -- first figure out how many
  //  samples we will use  
  unsigned int nsamples = (nwalkers * nsteps + thin - 1) / thin;
  unsigned int sidx;
  float *working = new float[nparams * nsamples];
  for (unsigned int i = 0; i < nsamples; ++i) {
    sidx = thin * nparams * i;
    for (unsigned int j = 0; j < nparams; ++j)
      working[j * nsamples + i] = steps[sidx + j];
  }

  // Do statistics
  // We will interpolate to get the ranges.  So it's convenient
  //  to pre-calculate the effective index in nsamples after we sort.
  float lowidx[nprob], highidx[nprob]; // Fractional index of interpolation
  for (unsigned int i = 0; i < nprob; ++i)
    lowidx[i] = 0.5 * (1.0 - plevels[i]) * nsamples;
  for (unsigned int i = 0; i < nprob; ++i)
    highidx[i] = 0.5 * (1 + plevels[i]) * nsamples;
  
  float *wptr;
  for (unsigned int j = 0; j < nparams; ++j) {
    wptr = working + j * nsamples;
    // Sort so we can get the percentiles.  This could be slightly
    //  more efficient if we used O(n) selection techniques instead
    //  of O(n ln n) sorting, but it's easily fast enough
    std::sort(wptr, wptr + nsamples);

    // Mean
    double mnval = wptr[0];  // accumulate in double
    for (unsigned int k = 1; k < nsamples; ++k)
      mnval += wptr[k];
    meanParams[j] = 
      static_cast<float>(mnval / static_cast<double>(nsamples));

    // Median.  We should, in principle, interpolate here as well,
    //  but in practice it's hardly worth it for any realistic
    //  number of samples
    medParams[j] = wptr[nsamples / 2];

    // Error snake
    for (unsigned int k = 0; k < nprob; ++k) {
      sidx = 2 * (k * nparams + j);
      // Lower value
      unsigned int idxl = static_cast<unsigned int>(lowidx[k]);
      uncertaintyParams[sidx] = wptr[idxl] + 
        (wptr[idxl+1] - wptr[idxl]) * (lowidx[k] - idxl);
      // And upper
      unsigned int idxu = static_cast<unsigned int>(highidx[k]);
      uncertaintyParams[sidx + 1] = wptr[idxu] + 
        (wptr[idxu+1] - wptr[idxu]) * (highidx[k] - idxu);
    }
  }
  delete[] working;
  
  // Now the smoother curves -- that is, the error snake
  //  Same thing, pretty much, more points, and we have
  //  to explicitly load up the model.
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

  // Note that we only process knots here, not bonus params,
  //  since there is not such thing as 
  working = new float[npoints * nsamples];
  paramSet p(nknots);
  for (unsigned int i = 0; i < nsamples; ++i) {
    sidx = thin * nparams * i;
    p.setParamValues(nknots, &steps[sidx]);
    model.setParams(p);
    // Convert counts to log10 so that in same units as params
    for (unsigned int j = 0; j < npoints; ++j)
      working[j * nsamples + i] = log10(model.getNumberCounts(s[j]));    
  }

  for (unsigned int j = 0; j < npoints; ++j) {
    wptr = working + j * nsamples;
    std::sort(wptr, wptr + nsamples);
    double mnval = wptr[0];
    for (unsigned int k = 1; k < nsamples; ++k)
      mnval += wptr[k];
    mean[j] = 
      static_cast<float>(mnval / static_cast<double>(nsamples));
    median[j] = wptr[nsamples / 2];
    for (unsigned int k = 0; k < nprob; ++k) {
      sidx = 2 * (k * npoints + j);
      // Lower value
      unsigned int idxl = static_cast<unsigned int>(lowidx[k]);
      snake[sidx] = wptr[idxl] + 
        (wptr[idxl+1] - wptr[idxl]) * (lowidx[k] - idxl);
      // And upper
      unsigned int idxu = static_cast<unsigned int>(highidx[k]);
      snake[sidx + 1] = wptr[idxu] + 
        (wptr[idxu+1] - wptr[idxu]) * (highidx[k] - idxu);
    }
  }
  delete[] working;
  
  delete[] steps;
}

/*!
  \param[in] filename File to write error snake to.
*/
void numberCountsKnotsSplineError::
writeAsHDF5(const std::string& filename) const {

  hid_t file_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("numberCountsKnotsSplineError", "writeToHDF5",
		                   "Failed to open HDF5 file to write");
  }

  hdf5utils::writeDataFloats(file_id, "KnotPositions", 
                             nknots, knotPositions);
  hdf5utils::writeDataStrings(file_id, "ParamNames", paramNames);

  hsize_t adims;
  hid_t mems_id ,att_id;
  adims = static_cast<size_t>(1);
  mems_id = H5Screate_simple(1, &adims, nullptr);
  att_id = H5Acreate1(file_id, "BestLikelihood", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &bestLike);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  hdf5utils::writeDataFloats(file_id, "BestParams", nparams, bestParams);
  hdf5utils::writeDataFloats(file_id, "MeanParams", nparams, meanParams);
  hdf5utils::writeDataFloats(file_id, "MedianParams", nparams, medParams);

  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "Paramsp=" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(file_id, pname.str(), nparams, 2,
                                 &uncertaintyParams[i * nparams * 2]);
  }

  hdf5utils::writeDataFloats(file_id, "FluxDensity", npoints, s);
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
