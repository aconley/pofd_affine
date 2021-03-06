#include<sstream>
#include<iomanip>

#include "../include/numberCountsKnotsSplineStats.h"
#include "../include/affineExcept.h"
#include "../include/paramSet.h"
#include "../include/hdf5utils.h"

/*!
  \param[in] npoints  Number of flux density points
*/
numberCountsKnotsSplineStats::
numberCountsKnotsSplineStats(unsigned int NPOINTS):
  nparams(0), nknots(0), knotPositions(nullptr), 
  initParams(nullptr), bestParams(nullptr), meanParams(nullptr), 
  medParams(nullptr), uncertaintyParams(nullptr), 
  npoints(NPOINTS) {
  s = new float[npoints];
  mean = new float[npoints];
  median = new float[npoints];
  snake = new float[2 * nprob * npoints];
  }

numberCountsKnotsSplineStats::~numberCountsKnotsSplineStats() {
  if (knotPositions != nullptr) delete[] knotPositions;
  if (initParams != nullptr) delete[] initParams;
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
void numberCountsKnotsSplineStats::setNParams(unsigned int n) {
  if (n == nparams) return;
  if (initParams != nullptr) delete[] initParams;
  if (n > 0) initParams = new float[n]; else initParams = nullptr;
  if (bestParams != nullptr) delete[] bestParams;
  if (n > 0) bestParams = new float[n]; else bestParams = nullptr;
  if (meanParams != nullptr) delete[] meanParams;
  if (n > 0) meanParams = new float[n]; else meanParams = nullptr;
  if (medParams != nullptr) delete[] medParams;
  if (n > 0) medParams = new float[n]; else medParams = nullptr;
  if (uncertaintyParams != nullptr) delete[] uncertaintyParams;
  if (n > 0) uncertaintyParams = new float[2 * nprob * n];
  else uncertaintyParams = nullptr;
  nparams = n;
}

/*!
  \param[in] infile Input file to read; must be output of pofd_affine_mcmc
  \param[in] thin Thinning to use on inputs
  \param[in] logspace Logarithmically space the interpolant points
*/
void numberCountsKnotsSplineStats::build(const std::string& infile,
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
    throw affineExcept("numberCountsKnotsSplineStats", "build",
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
  hdf5utils::readDataFloats(groupid, "InitialPosition", 
                            nparams, initParams);
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

  accumulateStats(nparams, nsamples, working, meanParams,
                  medParams, uncertaintyParams);
  delete[] working;

  // Now the smoother curves -- that is, the error snake
  //  Same thing, pretty much, more points, and we have
  //  to explicitly load up the model.
  // Set up flux densities
  double minflux = model.getMinFlux();
  double maxflux = model.getMaxFlux();
  maxflux -= 1e-4 * (maxflux - minflux); // maxflux counts are 0, avoid
  setupFluxes(minflux, maxflux, logspace, npoints, s);
  
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

  delete[] steps;

  accumulateStats(npoints, nsamples, working, mean,
                  median, snake);
  delete[] working;
}

/*!
  \param[in] filename File to write stats to.
*/
void numberCountsKnotsSplineStats::
writeAsHDF5(const std::string& filename) const {

  hid_t file_id, group_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("numberCountsKnotsSplineStats", "writeToHDF5",
                       "Failed to open HDF5 file to write");
  }

  hdf5utils::writeDataUnsignedInt(file_id, "NDims", 1);
  hdf5utils::writeDataString(file_id, "ModelType", 
                             "numberCountsKnotsSpline");
  
  group_id = H5Gcreate(file_id, "Params", H5P_DEFAULT, H5P_DEFAULT, 
                       H5P_DEFAULT);
  hdf5utils::writeDataFloats(group_id, "KnotPositions", 
                             nknots, knotPositions);
  hdf5utils::writeDataStrings(group_id, "ParamNames", paramNames);
  hdf5utils::writeDataFloats(group_id, "ParamsInit", nparams, initParams);

  hdf5utils::writeDataDouble(group_id, "BestLikelihood", bestLike);
  hdf5utils::writeDataFloats(group_id, "ParamsBest", nparams, bestParams);

  hdf5utils::writeDataFloats(group_id, "ParamsMean", nparams, meanParams);
  hdf5utils::writeDataFloats(group_id, "ParamsMedian", nparams, medParams);

  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "Paramsp=" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(group_id, pname.str(), nparams, 2,
                                 &uncertaintyParams[i * nparams * 2]);
  }
  H5Gclose(group_id);

  group_id = H5Gcreate(file_id, "Curves", H5P_DEFAULT, H5P_DEFAULT, 
                       H5P_DEFAULT);
  hdf5utils::writeDataFloats(group_id, "FluxDensity", npoints, s);
  hdf5utils::writeDataFloats(group_id, "Log10CountsMean", npoints, mean);
  hdf5utils::writeDataFloats(group_id, "Log10CountsMedian", 
                             npoints, median);
  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "Log10Countsp=" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(group_id, pname.str(),
                                 npoints, 2, &snake[i * npoints * 2]);
  }
  H5Gclose(group_id);

  H5Fclose(file_id);
}
