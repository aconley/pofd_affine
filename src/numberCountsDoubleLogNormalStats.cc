#include<sstream>
#include<iomanip>

#include "../include/numberCountsDoubleLogNormalStats.h"
#include "../include/affineExcept.h"
#include "../include/paramSet.h"
#include "../include/hdf5utils.h"


/*! \brief Utility function for converting params to sig(f2/f1) */
inline float modelParamsToSigF1F2(float sig, float mu) {
  float s2 = sig * sig;
  return sqrt(exp(2*mu + s2) * expm1(s2));
}

/*! \brief Utility function for converting params to <f2/f1> */
inline float modelParamsToMeanF1F2(float sig, float mu) {
  return exp(mu + 0.5 * sig * sig);
}

///////////////////////////////////////////////////////////

/*!
  \param[in] npoints  Number of flux density points
  \param[in] RATOS Return the color parameter constraints as
               constraints on f2 / f1 rather than the actual 
               log normal parameters
*/
numberCountsDoubleLogNormalStats::
numberCountsDoubleLogNormalStats(unsigned int NPOINTS, 
                                 bool RATIOS):
  nparams(0), nknots(0), nsigmas(0), noffsets(0),
  knotPositions(nullptr), sigmaKnotPositions(nullptr),
  offsetKnotPositions(nullptr), bestParams(nullptr), 
  meanParams(nullptr), medParams(nullptr), uncertaintyParams(nullptr), 
  npoints(NPOINTS), areRatios(RATIOS) {
  s1D = new float[npoints];
  mean1D = new float[npoints];
  median1D = new float[npoints];
  snake1D = new float[2 * nprob * npoints];
  sSigma = new float[npoints];
  meanSigma = new float[npoints];
  medianSigma = new float[npoints];
  snakeSigma = new float[2 * nprob * npoints];
  sOffset = new float[npoints];
  meanOffset = new float[npoints];
  medianOffset = new float[npoints];
  snakeOffset = new float[2 * nprob * npoints];
}

numberCountsDoubleLogNormalStats::~numberCountsDoubleLogNormalStats() {
  if (knotPositions != nullptr) delete[] knotPositions;
  if (sigmaKnotPositions != nullptr) delete[] sigmaKnotPositions;
  if (offsetKnotPositions != nullptr) delete[] offsetKnotPositions;

  if (bestParams != nullptr) delete[] bestParams;
  if (meanParams != nullptr) delete[] meanParams;
  if (medParams != nullptr) delete[] medParams;
  if (uncertaintyParams != nullptr) delete[] uncertaintyParams;

  delete[] s1D;
  delete[] mean1D;
  delete[] median1D;
  delete[] snake1D;
  delete[] sSigma;
  delete[] meanSigma;
  delete[] medianSigma;
  delete[] snakeSigma;
  delete[] sOffset;
  delete[] meanOffset;
  delete[] medianOffset;
  delete[] snakeOffset;
}

/*!
  \param[in] n New number of parameters
*/
void numberCountsDoubleLogNormalStats::setNParams(unsigned int n) {
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
void numberCountsDoubleLogNormalStats::
  build(const std::string& infile, unsigned int thin, bool logspace) {

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
    knotPositions[i] = model.getKnotPosition(i);
  nsigmas = model.getNSigmas();
  if (sigmaKnotPositions != nullptr) delete[] sigmaKnotPositions;
  sigmaKnotPositions = new float[nsigmas];
  for (unsigned int i = 0; i < nsigmas; ++i)
    sigmaKnotPositions[i] = model.getSigmaPosition(i);
  noffsets = model.getNOffsets();
  if (offsetKnotPositions != nullptr) delete[] offsetKnotPositions;
  offsetKnotPositions = new float[noffsets];
  for (unsigned int i = 0; i < noffsets; ++i)
    offsetKnotPositions[i] = model.getOffsetPosition(i);

  // Read in chains
  hid_t stepsid = H5Dopen(fileid, "Chains/Chain", H5P_DEFAULT);
  hid_t spaceid = H5Dget_space(stepsid);
  
  int ndims = H5Sget_simple_extent_ndims(spaceid);
  if (ndims != 3) {
    H5Sclose(spaceid);
    H5Dclose(stepsid);
    H5Fclose(fileid);
    throw affineExcept("numberCountsDoubleLogNormalStats", "build",
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
  accumulateStats(nparams, nsamples, working, meanParams,
                  medParams, uncertaintyParams);
  delete[] working;
 
  // Now the smoother curves -- that is, the Stats snake
  //  Same thing, pretty much, more points, and we have
  //  to explicitly load up the model.  We also have three
  //  sets of model parameters to work with.  We do them all
  //  simultaneously to avoid model reloads.
  // Set up flux densities
  float minflux = knotPositions[0];
  float maxflux = knotPositions[nknots-1];
  maxflux -= 1e-4 * (maxflux - minflux); // maxflux counts are 0, avoid
  setupFluxes(minflux, maxflux, logspace, npoints, s1D);
  // We also have to handle the case of their being only one
  // sigma or offset knot.
  if (nsigmas > 1) {
    minflux = sigmaKnotPositions[0];
    maxflux = sigmaKnotPositions[nsigmas-1];
  } // Otherwise, leave them at the band 1 range
  setupFluxes(minflux, maxflux, logspace, npoints, sSigma);
  if (noffsets > 1) {
    minflux = offsetKnotPositions[0];
    maxflux = offsetKnotPositions[noffsets-1];
  }
  setupFluxes(minflux, maxflux, logspace, npoints, sOffset);

  // Note that we only process knots here, not bonus params
  //  We need 3 arrays this time!
  working = new float[npoints * nsamples];
  float *workingSigma = new float[npoints * nsamples];
  float *workingOffset = new float[npoints * nsamples];
  unsigned int nactualpars = model.getNParams();
  paramSet p(nactualpars);
  unsigned int widx;
  for (unsigned int i = 0; i < nsamples; ++i) {
    sidx = thin * nparams * i;
    p.setParamValues(nactualpars, &steps[sidx]);
    model.setParams(p);
    // Convert counts to log10 so that in same units as params
    for (unsigned int j = 0; j < npoints; ++j) {
      widx = i + j * nsamples;
      working[widx] = log10(model.getBand1NumberCounts(s1D[j]));
    }
    if (areRatios) {
      float mu, sig;
      for (unsigned int j = 0; j < npoints; ++j) {
        widx = i + j * nsamples;
        mu = model.getOffset(sSigma[j]);
        sig = model.getSigma(sSigma[j]);
        workingSigma[widx] = modelParamsToSigF1F2(sig, mu);
        mu = model.getOffset(sOffset[j]);
        sig = model.getSigma(sOffset[j]);
        workingOffset[widx] = modelParamsToMeanF1F2(sig, mu);
      } 
    } else {
      for (unsigned int j = 0; j < npoints; ++j) {
        widx = i + j * nsamples;
        workingSigma[widx] = model.getSigma(sSigma[j]);
        workingOffset[widx] = model.getOffset(sOffset[j]);
      }
    }
  }
  delete[] steps;

  accumulateStats(npoints, nsamples, working, mean1D,
                  median1D, snake1D);
  delete[] working;
  accumulateStats(npoints, nsamples, workingSigma, meanSigma,
                  medianSigma, snakeSigma);
  delete[] workingSigma;
  accumulateStats(npoints, nsamples, workingOffset, meanOffset,
                  medianOffset, snakeOffset);
  delete[] workingOffset;
}

/*!
  \param[in] filename File to write stats to
*/
void numberCountsDoubleLogNormalStats::
writeAsHDF5(const std::string& filename) const {

  hid_t file_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
              H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("numberCountsDoubleLogNormalStats", "writeToHDF5",
                       "Failed to open HDF5 file to write");
  }

  hdf5utils::writeDataFloats(file_id, "KnotPositions", 
                             nknots, knotPositions);
  hdf5utils::writeDataFloats(file_id, "SigmaKnotPositions", 
                             nsigmas, sigmaKnotPositions);
  hdf5utils::writeDataFloats(file_id, "offsetKnotPositions", 
                             noffsets, offsetKnotPositions);
  hdf5utils::writeDataStrings(file_id, "ParamNames", paramNames);

  hsize_t adims;
  hid_t mems_id ,att_id;
  adims = static_cast<size_t>(1);
  mems_id = H5Screate_simple(1, &adims, nullptr);
  att_id = H5Acreate1(file_id, "BestLikelihood", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &bestLike);
  H5Aclose(att_id);
  hbool_t btmp = static_cast<hbool_t>(areRatios);
  att_id = H5Acreate1(file_id, "AreF2F1Ratios", H5T_NATIVE_HBOOL,
                      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &btmp);
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

  hdf5utils::writeDataFloats(file_id, "FluxDensity1D", npoints, s1D);
  hdf5utils::writeDataFloats(file_id, "MeanLog10Counts", 
                             npoints, mean1D);
  hdf5utils::writeDataFloats(file_id, "MedianLog10Counts", 
                             npoints, median1D);
  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "Log10CountsSnakep=" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(file_id, pname.str(),
                                 npoints, 2, &snake1D[i * npoints * 2]);
  }

  hdf5utils::writeDataFloats(file_id, "FluxDensitySigma", 
                             npoints, sSigma);
  hdf5utils::writeDataFloats(file_id, "MeanSigma", 
                             npoints, meanSigma);
  hdf5utils::writeDataFloats(file_id, "MedianSigma", 
                             npoints, medianSigma);
  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "SigmaSnakep=" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(file_id, pname.str(), npoints, 2,
                                 &snakeSigma[i * npoints * 2]);
  }
  hdf5utils::writeDataFloats(file_id, "FluxDensityOffset", 
                             npoints, sOffset);
  hdf5utils::writeDataFloats(file_id, "MeanOffset", 
                             npoints, meanOffset);
  hdf5utils::writeDataFloats(file_id, "MedianOffset", 
                             npoints, medianOffset);
  for (unsigned int i = 0; i < nprob; ++i) {
    std::stringstream pname;
    pname << "OffsetSnakep=" << std::setprecision(3) << plevels[i];
    hdf5utils::writeData2DFloats(file_id, pname.str(), npoints, 2,
                                 &snakeOffset[i * npoints * 2]);
  }

  H5Fclose(file_id);
}
