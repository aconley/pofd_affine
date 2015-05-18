//numberCountsKnots.cc
#include<iostream>
#include<cmath>
#include<iomanip>
#include<limits>
#include<cstdlib>
#include<ctime>
#include<sstream>

#include "../include/global_settings.h"
#include "../include/numberCountsKnots.h"
#include "../include/utility.h"
#include "../include/hdf5utils.h"
#include "../include/affineExcept.h"

numberCountsKnots::numberCountsKnots() {
  nknots = 0;
  knots = nullptr;
  logknotvals = nullptr;
  knotvals_loaded = false;
}

/*!
  \param[in] NKNOTS Number of knots

  Knot positions are set to NaN.
*/
numberCountsKnots::numberCountsKnots(unsigned int NKNOTS) :
  nknots(NKNOTS) {
  if (nknots > 0) {
    knots = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i)
      knots[i] = std::numeric_limits<double>::quiet_NaN();
    logknotvals = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i)
      logknotvals[i] = std::numeric_limits<double>::quiet_NaN();
  } else
    knots = logknotvals = nullptr;
  knotvals_loaded = false;
}  

/*!
  \param[in] S Knot positions
*/
numberCountsKnots::numberCountsKnots(const std::vector<float>& S) {
  nknots = 0;
  knots = nullptr;
  logknotvals = nullptr;
  setKnotPositions(S);
  knotvals_loaded = false;
}

/*!
  \param[in] S Knot positions
*/
numberCountsKnots::numberCountsKnots(const std::vector<double>& S) {
  nknots = 0;
  knots = nullptr;
  logknotvals = nullptr;
  setKnotPositions(S);
  knotvals_loaded = false;
}

/*!
  \param[in] n Number of knots
  \param[in] S Knot positions
*/
numberCountsKnots::numberCountsKnots(unsigned int n, const float* const S) {
  nknots = 0;
  knots = nullptr;
  logknotvals = nullptr;
  setKnotPositions(n, S);
  knotvals_loaded = false;
}

/*!
  \param[in] n Number of knots
  \param[in] S Knot positions
*/
numberCountsKnots::numberCountsKnots(unsigned int n, const double* const S) {
  nknots = 0;
  knots = nullptr;
  logknotvals = nullptr;
  setKnotPositions(n, S);
  knotvals_loaded = false;
}

/*!
  \param[in] other Other instance to copy
*/
numberCountsKnots::numberCountsKnots(const numberCountsKnots& other) {
  if (this == &other) return; //Self-copy
  nknots = 0;
  knots = nullptr;
  logknotvals = nullptr;
  setNKnots(other.nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = other.knots[i];
  if (other.knotvals_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      knots[i] = other.logknotvals[i];
  }
  knotvals_loaded = other.knotvals_loaded;
}

numberCountsKnots::~numberCountsKnots() {
  if (knots != nullptr) delete[] knots;
  if (logknotvals != nullptr) delete[] logknotvals;
}

/*!
  \param[in] n Number of knots.

  If number of knots does not change, then previous contents are
  preserved.  Otherwise, they are cleared.
*/
void numberCountsKnots::setNKnots(unsigned int n) {
  if (nknots == n) return;
  if (knots != nullptr) delete[] knots;
  if (logknotvals != nullptr) delete[] logknotvals;
  if (n > 0) {
    knots = new double[n];
    logknotvals = new double[n];
  } else {
    knots = logknotvals = nullptr;
  }
  nknots = n;
  knotvals_loaded = false;
}

/*!
  \param[out] S Knot positions
*/
void numberCountsKnots::getKnotPositions(std::vector<double>& S) const {
  S.resize(nknots);
  if (nknots > 0)
    for (unsigned int i = 0; i < nknots; ++i)
      S[i] = knots[i];
}

/*!
  \param[in] S Input knot positions
*/
void numberCountsKnots::setKnotPositions(const std::vector<float>& S) {
  unsigned int n = S.size();
  if (n != nknots) setNKnots(n);
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsKnots", "setKnots",
                         "Negative knot positions not allowed");
  // Knot positions -must- increase, the GSL requires
  for (unsigned int i = 1; i < nknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Knot positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsKnots", "setKnotPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = static_cast<double>(S[i]);
}

/*!
  \param[in] S Input knot positions
*/
void numberCountsKnots::setKnotPositions(const std::vector<double>& S) {
  unsigned int n = S.size();
  if (n != nknots) setNKnots(n);
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsKnots", "setKnots",
                         "Negative knot positions not allowed");
  // Knot positions -must- increase, the GSL requires
  for (unsigned int i = 1; i < nknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Knot positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsKnots", "setKnotPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
}

/*!
  \param[in] n Number of knots
  \param[in] S Input knot positions
*/
void numberCountsKnots::setKnotPositions(unsigned int n, 
                                         const float* const S) {
  if (n != nknots) setNKnots(n);
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsKnots", "setKnots",
                         "Negative knot positions not allowed");
  // Knot positions -must- increase, the GSL requires
  for (unsigned int i = 1; i < nknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Knot positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsKnots", "setKnotPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = static_cast<double>(S[i]);
}

/*!
  \param[in] n Number of knots
  \param[in] S Input knot positions
*/
void numberCountsKnots::setKnotPositions(unsigned int n, 
                                         const double* const S) {
  if (n != nknots) setNKnots(n);
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsKnots", "setKnots",
                         "Negative knot positions not allowed");
  // Knot positions -must- increase, the GSL requires
  for (unsigned int i = 1; i < nknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Knot positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsKnots", "setKnotPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
}


/*!
  \param[out] F Parameters from model

  Will set the first nknot parameters, ignoring any others.
  If the knot values aren't loaded, sets them to NaN.
 */
void numberCountsKnots::getParams(paramSet& F) const {
  if (F.getNParams() < nknots)
    throw affineExcept("numberCountsKnots", "getKnots",
                       "Not enough space in output variable");
  if (! knotvals_loaded)
    for (unsigned int i = 0; i < nknots; ++i)
      F[i] = std::numeric_limits<double>::quiet_NaN();
  else
    for (unsigned int i = 0; i < nknots; ++i)
      F[i] = pofd_mcmc::ilogfac * logknotvals[i]; //Convert to base 10
}

/*!
  \param[in] F Parameters to set in model
 */
void numberCountsKnots::setParams(const paramSet& F) {
  if (F.getNParams() <= 2)
    throw affineExcept("numberCountsKnots" ,"setKnots",
                       "No knots present");
  if (nknots != (F.getNParams())) 
    throw affineExcept("numberCountsKnots", "setKnots",
                       "Number of knot values different than expected");
  // Internal storage is base 2 (and double)
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = pofd_mcmc::logfac * static_cast<double>(F[i]); 
  knotvals_loaded = true;
}

/*!
  \returns True if the model parameters are valid
 */
bool numberCountsKnots::isValid() const {
  if (nknots == 0) return false;
  if (!knotvals_loaded) return false;
  for (unsigned int i = 0; i < nknots; ++i)
    if (std::isnan(knots[i])) return false;
  if (knots[0] <= 0.0) return false;
  for (unsigned int i = 1; i < nknots; ++i)
    if (knots[i] <= knots[i-1]) return false;
  for (unsigned int i = 0; i < nknots; ++i)
    if (std::isnan(logknotvals[i])) return false;
  return true;
}

/*!
  \param[in] p1 First parameter set
  \param[in] p2 Second parameter set
  \returns Relative distance between p1, p2 over entries this model cares 
     about.  This is defined as the sqrt of the sum differences
     between the two parameter sets divided by the sqrt of the
     length of the first parameter (or, if the lenght is 0, just
     the sqrt of the sum of the squared differences).
*/
float numberCountsKnots::paramRelativeDistance(const paramSet& p1, 
                                               const paramSet& p2) 
  const throw(affineExcept) {

  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  if (p1.getNParams() < nknots)
    throw affineExcept("numberCountsKnots", "paramRelativeDistance",
                       "paramSet 1 doesn't have enough entries");
  if (p2.getNParams() < nknots)
    throw affineExcept("numberCountsKnots", "paramRelativeDistance",
                       "paramSet 2 doesn't have enough entries");
  if (&p1 == &p2) return 0.0; // Same object!
  float val = p1[0];
  float diffval = val - p2[0];
  float lensum = val * val;
  float diffsum = diffval * diffval;
  for (unsigned int i = 1; i < nknots; ++i) {
    val = p1[i];
    lensum += val * val;
    diffval = val - p2[i];
    diffsum += diffval * diffval;
  }
  if (lensum == 0)
    return sqrt(diffsum);
  return sqrt(diffsum / lensum);
}

/*!
  \returns Minimum flux density model is non-zero at
*/
double numberCountsKnots::getMinFlux() const {
  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  return knots[0];
}

/*!
  \returns Maximum flux density model is non-zero at
*/
double numberCountsKnots::getMaxFlux() const {
  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  return knots[nknots-1];
}

/*!
  \param[in] bm Beam to compute range for
  \returns A pair of min/max fluxes over which R is expected to be nonzero.

  Doesn't check for model validity.
  We also include flux=0, even though R is technically zero there.
*/
dblpair numberCountsKnots::getRRangeInternal(const beam& bm) const 
  throw(affineExcept) {

  double maxknot = knots[nknots - 1];
  
  bool haspos = bm.hasPos();
  bool hasneg = bm.hasNeg();

  double minRF, maxRF;
  if (haspos) {
    maxRF = maxknot * bm.getMinMaxPos().second;
    if (hasneg)
      minRF = - maxknot * bm.getMinMaxNeg().second;
    else
      minRF = 0.0;
  } else {
    maxRF = 0.0;
    minRF = - maxknot * bm.getMinMaxNeg().second;
  }
  return std::make_pair(minRF, maxRF);
}

/*!
  \param[in] bm Beam to compute range for
  \returns A pair of min/max fluxes over which R is expected to be nonzero.

  We also include flux=0, even though R is technically zero there.
*/
dblpair numberCountsKnots::getRRange(const beam& bm) const
  throw(affineExcept) {

  if (!isValid())
    throw affineExcept("numberCountsKnots", "getRRange",
                       "Model is not valid");
  if (!bm.hasData())
    throw affineExcept("numberCountsKnots", "getRRange",
                       "Beam is empty");
  return getRRangeInternal(bm);
}

/*!
  \param[inout] comm MPI communicator
  \param[in] dest Destination to send messages to
*/
void numberCountsKnots::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&nknots), 1, MPI_UNSIGNED, dest,
           pofd_mcmc::NCKSENDNKNOTS, comm);
  if (nknots != 0) {
    MPI_Send(knots, nknots, MPI_DOUBLE, dest, pofd_mcmc::NCKSENDKNOTS, comm);
    MPI_Send(const_cast<bool*>(&knotvals_loaded), 1, MPI::BOOL, dest,
             pofd_mcmc::NCKSENDKNOTSLOADED, comm);
    if (knotvals_loaded)
      MPI_Send(logknotvals,nknots,MPI_DOUBLE,dest,
               pofd_mcmc::NCKSENDLOGKNOTVALS, comm);
  }
}

/*!
  \param[inout] comm MPI communicator
  \param[in] src Where messages will come from
*/
void numberCountsKnots::receiveCopy(MPI_Comm comm, int src) {
  unsigned int n;
  MPI_Status Info;
  MPI_Recv(&n, 1, MPI_UNSIGNED, src, pofd_mcmc::NCKSENDNKNOTS, comm, &Info);
  if (n != 0) {
    if (n != nknots) setNKnots(n);
    MPI_Recv(knots, nknots, MPI_DOUBLE, src, pofd_mcmc::NCKSENDKNOTS,
             comm, &Info);
    MPI_Recv(&knotvals_loaded, 1, MPI::BOOL, src, 
             pofd_mcmc::NCKSENDKNOTSLOADED, comm, &Info);
    if (knotvals_loaded)
      MPI_Recv(logknotvals, nknots, MPI_DOUBLE, src,
               pofd_mcmc::NCKSENDLOGKNOTVALS, comm, &Info);
  }
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] writevals If true, write the knot values in addition to
      positions.  If the knots aren't loaded, this is ignored.
*/
void numberCountsKnots::writeToHDF5Handle(hid_t objid, bool writevals) const {
  hsize_t adims;
  hid_t mems_id, att_id;

  if (H5Iget_ref(objid) < 0)
    throw affineExcept("numberCountsKnots", "writeToHDF5Handle",
                       "Input handle is not valid");

  // Number of knots
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  att_id = H5Acreate2(objid, "NKnots", H5T_NATIVE_UINT,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nknots);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  // Knot positions
  hdf5utils::writeDataDoubles(objid, "KnotPositions", nknots, knots);

  if (writevals && knotvals_loaded) {
    // Convert to log 10 for write
    double* kv;
    kv = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i)
      kv[i] = pofd_mcmc::ilogfac * logknotvals[i];
    hdf5utils::writeDataDoubles(objid, "Log10KnotValues", nknots, kv);
    delete[] kv;
  }
}


/*!
  \param[in] objid HDF5 handle to read from

  Initializes a model from an HDF5 file, setting the knot positions
  but not values.
*/
void numberCountsKnots::readFromHDF5Handle(hid_t objid) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("numberCountsKnots", "readFromHDF5Handle",
                       "Input handle is not valid");

  unsigned int f_nknots = 
    hdf5utils::readAttUnsignedInt(objid, "NKnots");
  setNKnots(f_nknots);
  if (f_nknots > 0) {
    // We could read directly into knots, but that could be a problem
    // for subclasses, which might want to do extra work.
    double *newknots = new double[nknots];
    hdf5utils::readDataDoubles(objid, "KnotPositions", nknots, newknots);
    setKnotPositions(nknots, newknots);
    delete[] newknots;
  }
}

/*!
  \param[inout] os Stream to write to
*/
bool numberCountsKnots::writeToStream(std::ostream& os) const {
  os << "Model parameters: " << std::endl;
  if (knotvals_loaded) {
    os << " " << std::left << std::setw(13) << "#Flux knot" << "  "
       << std::setw(13) << "Knot value" << std::endl;
    //Convert to log10 for output
    for (unsigned int i = 0; i < nknots; ++i)
      os << " " << std::left << std::setw(13) << knots[i] << "  "
         << std::setw(13) << pofd_mcmc::ilogfac * logknotvals[i] << std::endl; 
  } else {
    os << "Number of knots: " << nknots << std::endl;
    if (nknots > 0) {
      os << "Knot positions: " << knots[0];
      for (unsigned int i = 1; i < nknots; ++i)
        os << " " << knots[i];
      os << std::endl;
    }
  }
  return true;
}

/*!
  \param[inout] os Stream to write to
  \param[in] b Model to write
*/
std::ostream& operator<<(std::ostream& os, const numberCountsKnots& b) {
  b.writeToStream(os);
  return os;
}
