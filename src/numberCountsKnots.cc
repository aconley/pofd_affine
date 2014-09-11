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
  knots = NULL;
  logknotvals = NULL;
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
    knots = logknotvals = NULL;
  knotvals_loaded = false;
}  

/*!
  \param[in] S Knot positions
*/
numberCountsKnots::numberCountsKnots(const std::vector<float>& S) {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setKnotPositions(S);
  knotvals_loaded = false;
}

/*!
  \param[in] S Knot positions
*/
numberCountsKnots::numberCountsKnots(const std::vector<double>& S) {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setKnotPositions(S);
  knotvals_loaded = false;
}

/*!
  \param[in] n Number of knots
  \param[in] S Knot positions
*/
numberCountsKnots::numberCountsKnots(unsigned int n, const float* const S) {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setKnotPositions(n, S);
  knotvals_loaded = false;
}

/*!
  \param[in] n Number of knots
  \param[in] S Knot positions
*/
numberCountsKnots::numberCountsKnots(unsigned int n, const double* const S) {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setKnotPositions(n, S);
  knotvals_loaded = false;
}

/*!
  \param[in] other Other instance to copy
*/
numberCountsKnots::numberCountsKnots(const numberCountsKnots& other) {
  if (this == &other) return; //Self-copy
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
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
  if (knots != NULL) delete[] knots;
  if (logknotvals != NULL) delete[] logknotvals;
}

/*!
  \param[in] n Number of knots.

  If number of knots does not change, then previous contents are
  preserved.  Otherwise, they are cleared.
*/
void numberCountsKnots::setNKnots(unsigned int n) {
  if (nknots == n) return;
  if (knots != NULL) delete[] knots;
  if (logknotvals != NULL) delete[] logknotvals;
  if (n > 0) {
    knots = new double[n];
    logknotvals = new double[n];
  } else {
    knots = logknotvals = NULL;
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
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(objid, "nknots", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nknots);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  // Knot positions
  hdf5utils::writeDataDoubles(objid, "knotpos", nknots, knots);

  if (writevals && knotvals_loaded) {
    // Convert to log 10 for write
    double* kv;
    kv = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i)
      kv[i] = pofd_mcmc::ilogfac * logknotvals[i];
    hdf5utils::writeDataDoubles(objid, "log10knotvals", nknots, kv);
    delete[] kv;
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
  } else
    os << "Number of knots: " << nknots << std::endl;
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

////////////////////////////////////////

initFileKnots::initFileKnots() : nknots(0), has_range(false), 
				 has_lower_limits(false),
				 has_upper_limits(false) {

  knotpos = knotval = range = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;

}

/*
  \param[in] flname File to read from
  \param[in] read_range      Read in (and require) knot ranges
  \param[in] read_limits     Try to read limits; this will turn on require_range
  
  see initFileKnots::readFile for more details of the file format
*/
initFileKnots::initFileKnots(const std::string& flname, 
			     bool read_range, bool read_limits) :
  nknots(0), has_range(false), has_lower_limits(false), 
  has_upper_limits(false) {
  
  knotpos = knotval = range = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;

  readFile(flname, read_range, read_limits);
}

initFileKnots::~initFileKnots() {
  if (knotpos != NULL) delete[] knotpos;
  if (knotval != NULL) delete[] knotval;
  if (range != NULL) delete[] range;
  if (has_lowlim != NULL) delete[] has_lowlim;
  if (lowlim != NULL) delete[] lowlim;
  if (has_uplim != NULL) delete[] has_uplim;
  if (uplim != NULL) delete[] uplim;
}

/*!
  \param[in] flname File to read from
  \param[in] read_range      Read in (and require) knot ranges
  \param[in] read_limits     Try to read limits; this will turn on require_sigma

  The file format is a bunch of lines of the form
  knotpos   knotval   [range [ lowlim [ uplim ]]
  So knotpos, knotval are always required
  range is optionally required if read_range is set.  If used to generate
   points, they are generated within range/2 on either side uniformly
   around knotval.
  lowlim and uplim may be present, and are looked for if read_limits is set.
    range must also be present, and the first element found is lowlim.
    If another is also found, it is interpreted as uplim
*/
void initFileKnots::readFile(const std::string& flname, 
			     bool read_range, bool read_limits) {
  if (read_limits) read_range = true;

  //Clear old data
  nknots = 0;
  if (knotpos != NULL) delete[] knotpos;
  if (knotval != NULL) delete[] knotval;
  if (range != NULL) delete[] range;
  if (has_lowlim != NULL) delete[] has_lowlim;
  if (lowlim != NULL) delete[] lowlim;
  if (has_uplim != NULL) delete[] has_uplim;
  if (uplim != NULL) delete[] uplim;
  knotpos = knotval = range = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;
  has_range = has_lower_limits = has_upper_limits = false;

  //Figure out how many elements we require
  unsigned int nreq = 2; //Pos, value
  if (read_range) nreq += 1; //Range -- read limits not required ever

  std::ifstream initfs(flname.c_str());
  if (!initfs) {
    initfs.close();
    std::stringstream errmsg;
    errmsg << "Unable to open file:" << flname << std::endl;
    throw affineExcept("initFileKnots", "readFile", errmsg.str());
  }

  //Do the read into temporary vectors, then copy
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currval;
  std::vector<double> kp, kv, ks, kl, ku;
  std::vector<bool> hl, hu;

  while (!initfs.eof()) {
    std::getline(initfs,line);
    if (line[0] == '#') continue; //Skip comments

    //Parse into words, stipping spaces
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < nreq) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currval;
    kp.push_back(currval);
    str.str(words[1]); str.clear(); str >> currval;
    kv.push_back(currval);
    if (read_range) {
      has_range = true;
      str.str(words[2]); str.clear(); str >> currval;
      ks.push_back(currval);
    }
    if (read_limits) { 
      if (words.size() > 3) {
	// ignore limits if range is zero, since they are irrelevant
	if (has_range && ks.back() == 0) {
	  hl.push_back(false);
	  kl.push_back(std::numeric_limits<double>::quiet_NaN());
	  hu.push_back(false);
	  ku.push_back(std::numeric_limits<double>::quiet_NaN());
	} else {
	  has_lower_limits = true;
	  str.str(words[3]); str.clear(); str >> currval;
	  hl.push_back(true);
	  kl.push_back(currval);
	  if (words.size() > 4) {
	    has_upper_limits = true;
	    str.str(words[4]); str.clear(); str >> currval;
	    hu.push_back(true);
	    ku.push_back(currval);
	  } else {
	    hu.push_back(false);
	    ku.push_back(std::numeric_limits<double>::quiet_NaN());
	  }
	}
      } else {
	hl.push_back(false);
	kl.push_back(std::numeric_limits<double>::quiet_NaN());
	hu.push_back(false);
	ku.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }
  }
  initfs.close();
  
  nknots = kp.size();
  if (nknots == 0) {
    std::stringstream errstr;
    errstr << "No knot positions or values found in " << flname;
    throw affineExcept("initFileKnots", "readFiles", errstr.str());
  }

  //Copy into vectors
  knotpos = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotpos[i] = kp[i];
  knotval = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotval[i] = kv[i];
  if (has_range) {
    range = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i) range[i] = ks[i];
  }
  if (has_lower_limits) {
    has_lowlim = new bool[nknots];
    lowlim = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i) has_lowlim[i] = hl[i];
    for (unsigned int i = 0; i < nknots; ++i) lowlim[i] = kl[i];
  }
  if (has_upper_limits) {
    has_uplim = new bool[nknots];
    uplim = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i) has_uplim[i] = hu[i];
    for (unsigned int i = 0; i < nknots; ++i) uplim[i] = ku[i];
  }

  //Make sure limits and range are valid
  checkLimitsDontCross();
  checkRange();

}

/*!
  Throws an exception if they aren't valid
*/
void initFileKnots::checkLimitsDontCross() const {
  if (nknots == 0) return;
  if (!(has_lower_limits && has_upper_limits)) return; //Nothing to check
  for (unsigned int i = 0; i < nknots; ++i) {
    if (!(has_uplim[i] && has_lowlim[i])) continue;
    if (uplim[i] < lowlim[i]) {
      std::stringstream errstr;
      errstr << "Lower/Upper limits cross at index: " << i << std::endl;
      errstr << " Lower limit: " << lowlim[i] 
	     << " Upper limit: " << uplim[i];
      throw affineExcept("initFileKnots", "checkLimitsDontCross", errstr.str());
    }
    if ((range[i] > 0.) && (uplim[i] == lowlim[i])) {
      std::stringstream errstr;
      errstr << "Lower/Upper limits meet at index: " << i 
	     << " but range is not zero" << std::endl;
      errstr << " Lower limit: " << lowlim[i] << " Upper limit: " << uplim[i]
	     << " range: " << range[i];
      throw affineExcept("initFileKnots", "CheckLimitsDontCross", errstr.str());
    }
  }
}

/*!
  Throws an exception if any parameter range is not valid
*/
void initFileKnots::checkRange() const {
  //Make sure that if range is 0 then the mean value falls within
  // the range of any limits

  // Quick returns
  if (!has_range) return; //Nothing to check
  if (!(has_lower_limits || has_upper_limits)) return;

  for (unsigned int i = 0; i < nknots; ++i)
    if (range[i] == 0) {
      if (has_lower_limits && has_lowlim[i] && (knotval[i] < lowlim[i])) {
	std::stringstream errstr;
	errstr << "At knot " << i << " range is zero but mean value "
	       << knotval[i] << std::endl << " lies below lower limit "
	       << lowlim[i];
	throw affineExcept("initFileKnots", "checkRange", errstr.str());
      }
      if (has_upper_limits && has_uplim[i] && (knotval[i] > uplim[i])) {
	std::stringstream errstr;
	errstr << "At knot " << i << " range is zero but mean value "
	       << knotval[i] << std::endl << " lies above upper limit "
	       << uplim[i];
	throw affineExcept("initFileKnots", "readFiles", errstr.str());
      }
    }
}

/*!
  \param[in] idx Knot index
  \returns A pair of the knot position and central value
*/
dblpair initFileKnots::getKnot(unsigned int idx) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots", "getKnot",
		       "No knot information read in");
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "getKnot", "Invalid knot index");
  return std::make_pair(knotpos[idx], knotval[idx]);
}

/*
  \param[out] kp Set to knot positions on output
*/
void initFileKnots::getKnotPos(std::vector<double>& kp) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots", "getKnotPos",
		       "No knot information read in");
  kp.resize(nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    kp[i] = knotpos[i];
}

/*
  \param[in] idx Index of knot to get
  \returns Knot position at that index
*/
double initFileKnots::getKnotPos(unsigned int idx) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots", "getKnotPos",
		       "No knot information read in");
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "getKnotPos",
		       "Invalid knot index");
  return knotpos[idx];
}

/*
  \param[inout] model Modified on output; knot positions are set

  This will change the number of knots in the model if they
  don't match.
*/
void initFileKnots::getKnotPos(numberCountsKnots& model) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots", "getKnotPos",
		       "No knot information read in");
  model.setKnotPositions(nknots, knotpos);
}

/*
  \param[out] kp Set to knot positions on output
*/
void initFileKnots::getKnotVals(std::vector<double>& kv) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots", "getKnotVals",
		       "No knot information read in");
  kv.resize(nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    kv[i] = knotval[i];
}

/*
  \param[in] idx Index of knot to get
  \returns Knot value at that index
*/
double initFileKnots::getKnotValue(unsigned int idx) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots", "getKnotValue",
		       "No knot information read in");
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "getKnotValue",
		       "Invalid knot index");
  return knotval[idx];
}


/*
  \param[out] p Parameter set to fill.  Must be pre-allocated
       to have at least nknots elements.

  This only fills the first nknots parameters
*/
void initFileKnots::getParams(paramSet& p) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots", "getParams",
		       "No information loaded");
  if (p.getNParams() < nknots)
    throw affineExcept("initFileKnots", "getParams",
		       "Not enough space in provided paramSet");
  for (unsigned int i = 0; i < nknots; ++i)
    p[i] = knotval[i];
}

/*
  \param[in] rangen Random number generator
  \param[out] pnew New parameter set generated.  Must be pre-allocated
       to have at least nknots elements.

  This only fills the first nknots parameters.  It uses the central
  values from the initFile
*/
void initFileKnots::generateRandomKnotValues(ran& rangen, 
					     paramSet& pnew) const {
  paramSet pcen(pnew.getNParams());
  getParams(pcen); //Load central values into pcen
  generateRandomKnotValues(rangen, pnew, pcen); //Get new values
}

/*
  \param[in] rangen Random number generator
  \param[out] pnew New parameter set generated.  Must be pre-allocated
     to have at least nknots elements
  \param[in] pcen  Central parameter values

  This only fills the first nknots parameters.  This version
  allows the caller to use different central values than the ones
  from the initial file, but keep the ranges, limits, etc.
*/
void initFileKnots::generateRandomKnotValues(ran& rangen, paramSet& pnew, 
					     const paramSet& pcen) const {
  const unsigned int maxiters = 1000; //Maximum number of generation attempts

  if (nknots == 0)
    throw affineExcept("initFileKnots", "generateRandomKnotValues",
		       "No knot information read in");
    
  //Make sure p is big enough; don't resize, complain
  if (pnew.getNParams() < nknots)
    throw affineExcept("initFileKnots", "generateRandomKnotValues",
		       "Not enough space in provided new paramSet");
  if (pcen.getNParams() < nknots)
    throw affineExcept("initFileKnots", "generateRandomKnotValues",
		       "Not enough values in provided central values");

  //Deal with simple case -- everything fixed
  //So just return pcen
  if (!has_range) {
    //Check pcen is within limits
    if (has_lower_limits)
      for (unsigned int i = 0; i < nknots; ++i)
	if (has_lowlim[i] && (pcen[i] < lowlim[i])) {
	  std::stringstream errstr;
	  errstr << "For parameter " << i << " user provided central value "
		 << pcen[i] << " is below lower limit " << lowlim[i];
	  throw affineExcept("initFileKnots", "generateRandomKnotValues",
			     errstr.str());
	}
    if (has_upper_limits)
      for (unsigned int i = 0; i < nknots; ++i)
	if (has_uplim[i] && (pcen[i] > uplim[i])) {
	  std::stringstream errstr;
	  errstr << "For parameter " << i << " user provided central value "
		 << pcen[i] << " is above upper limit " << uplim[i];
	  throw affineExcept("initFileKnots", "generateRandomKnotValues",
			     errstr.str());
	}
    
    for (unsigned int i = 0; i < nknots; ++i)
      pnew[i] = pcen[i];
    return;
  }

  //Now we have at least some ranges
  //The simple case is if there are no limits.  If there are, we will
  // have to do trials.
  if (!(has_lower_limits || has_upper_limits)) {
    for (unsigned int i = 0; i < nknots; ++i)
      if (range[i] > 0)
	pnew[i] = rangen.flt() * range[i] + (pcen[i] - 0.5 * range[i]);
      else
	pnew[i] = pcen[i];
  } else {
    //Both ranges and limits
    bool goodval;
    double trialval;
    unsigned int iters;
    for (unsigned int i = 0; i < nknots; ++i) {
      if (range[i] > 0) {
	//Some sanity checks
	if (has_lowlim[i] && (lowlim[i] > pcen[i] + range[i])) {
	  std::stringstream errstr;
	  errstr << "Lower limit is too far above central value; will not be "
		 << "able to" << std::endl << "generate value for param idx: "
		 << i << " with lowlim: " << lowlim[i] << " central: " 
		 << pcen[i] << " range: " << range[i];
	  throw affineExcept("initFileKnots", "generateRandomKnotValues",
			     errstr.str());
	}
	if (has_uplim[i] && (uplim[i] < pcen[i] - range[i])) {
	  std::stringstream errstr;
	  errstr << "Upper limit is too far below central value; will not be"
		 << " able to " << std::endl << "generate value for param idx: "
		 << i << " with uplim: " << uplim[i] << " central: " 
		 << pcen[i] << " range: " << range[i];
	  throw affineExcept("initFileKnots", "generateRandomKnotValues",
			     errstr.str());
	}
	
	if (has_lowlim[i] && has_uplim[i]) {
	  //If space between them is smaller than range, use
	  // those limits to generate rather than accept/reject
	  double rng = uplim[i] - lowlim[i];
	  if (rng < range[i]) {
	    pnew[i] = lowlim[i] + rng * rangen.flt();
	  } else {
	    //Trial
	    goodval = false;
	    iters = 0;
	    while (!goodval) {
	      if (iters >= maxiters) {
		std::stringstream errstr;
		errstr << "Failed to generate acceptable value for param "
		       << i << " after " << iters << " attempts";
		throw affineExcept("initFileKnots", "generateRandomKnotValues",
				   errstr.str());
	      }
	      trialval = rangen.flt() * range[i] + (pcen[i] - 0.5 * range[i]);
	      if ((trialval >= lowlim[i]) && (trialval <= uplim[i])) 
		goodval = true;
	      ++iters;
	    }
	    pnew[i] = trialval;
	  }
	} else if (has_lowlim[i]) {
	  //Lower limit only
	  goodval = false;
	  iters = 0;
	  while (!goodval) {
	    if (iters >= maxiters) {
	      std::stringstream errstr;
	      errstr << "Failed to generate acceptable value for param "
		     << i << " after " << iters << " attempts";
	      throw affineExcept("initFileKnots","generateRandomKnotValues",
				 errstr.str());
	    }
	    trialval = rangen.flt() * range[i] + (pcen[i] - 0.5 * range[i]);
	    if (trialval >= lowlim[i]) goodval = true;
	    ++iters;
	  }
	  pnew[i] = trialval;
	} else if (has_uplim[i]) {
	  //Upper limit only
	  goodval = false;
	  iters = 0;
	  while (!goodval) {
	    if (iters >= maxiters) {
	      std::stringstream errstr;
	      errstr << "Failed to generate acceptable value for param "
		     << i << " after " << iters << " attempts";
	      throw affineExcept("initFileKnots","generateRandomKnotValues",
				 errstr.str());
	    }
	    trialval = rangen.flt() * range[i] + (pcen[i] - 0.5 * range[i]);
	    if (trialval <= uplim[i]) goodval = true;
	    ++iters;
	  }
	  pnew[i] = trialval;
	} else {
	  //No limit, easy cakes
	  pnew[i] = rangen.flt() * range[i] + (pcen[i] - 0.5 * range[i]);
	}
      } else {
	//Range is 0.  Check to make sure this is within the limits
	if (has_lowlim[i] && (pcen[i] < lowlim[i])) {
	  std::stringstream errstr;
	  errstr << "For parameter " << i << " user provided central value "
		 << pcen[i] << " is below lower limit " << lowlim[i]
		 << " and range is 0.";
	  throw affineExcept("initFileKnots", "generateRandomKnotValues",
			     errstr.str());
	}
	if (has_uplim[i] && (pcen[i] > uplim[i])) {
	  std::stringstream errstr;
	  errstr << "For parameter " << i << " user provided central value "
		 << pcen[i] << " is above upper limit " << uplim[i]
		 << " and range is 0.";
	  throw affineExcept("initFileKnots", "generateRandomKnotValues",
			     errstr.str());
	}
	pnew[i] = pcen[i];
      }
    }
  }
}

/*!
  \param[in] idx Index of knot

  \returns Range for that knot
*/
double initFileKnots::getKnotRange(unsigned int idx) const {
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "getKnotRange", 
		       "Invalid knot index");
  if (!has_range) return std::numeric_limits<double>::quiet_NaN();
  return range[idx];
}

/*
  \param[in] idx Knot index
  \returns True if knot is fixed, otherwise false
*/
bool initFileKnots::isKnotFixed(unsigned int idx) const {
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "isKnotFixed",
		       "Invalid knot index");
  if (!has_range) return false;
  if (range[idx] == 0) return true;
  return false;
}

/*
  \param[in] idx Knot index
  \returns True if knot has lower limit, otherwise false
*/
bool initFileKnots::knotHasLowerLimit(unsigned int idx) const {
  if (!has_lower_limits) return false;
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "knotHasLowerLimit",
		       "Invalid knot index");
  return has_lowlim[idx];
}

/*
  \param[in] idx Knot index
  \returns Lower limit on knot, or NaN if none
*/
double initFileKnots::getLowerLimit(unsigned int idx) const {
  if (!has_lower_limits) return std::numeric_limits<double>::quiet_NaN();
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "getLowerLimit",
		       "Invalid knot index");
  if (!has_lowlim[idx]) return std::numeric_limits<double>::quiet_NaN();
  return lowlim[idx];
}

/*
  \param[in] idx Knot index
  \returns True if knot has upper limit, otherwise false
*/
bool initFileKnots::knotHasUpperLimit(unsigned int idx) const {
  if (!has_upper_limits) return false;
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "knotHasUpperLimit",
		       "Invalid knot index");
  return has_uplim[idx];
}


/*
  \param[in] idx Knot index
  \returns Upper limit on knot, or NaN if none
*/
double initFileKnots::getUpperLimit(unsigned int idx) const {
  if (!has_upper_limits) return std::numeric_limits<double>::quiet_NaN();
  if (idx >= nknots)
    throw affineExcept("initFileKnots", "getUpperLimit",
		       "Invalid knot index");
  if (!has_uplim[idx]) return std::numeric_limits<double>::quiet_NaN();
  return uplim[idx];
}


/*!
  \param[in] p Parameter set to check validity of
  \returns True if p is valid (within limits), otherwise false
*/
bool initFileKnots::isValid(const paramSet& p) const {
  if (!(has_lower_limits || has_upper_limits)) return true;
  if (p.getNParams() < nknots)
    throw affineExcept("initFileKnots", "isValid",
		       "Not enough params in paramSet to test validity");
  double val;
  for (unsigned int i = 0; i < nknots; ++i) {
    val = p[i];
    if (has_lowlim[i] && (val < lowlim[i])) return false;
    if (has_uplim[i] && (val > uplim[i])) return false;
  }
  return true;
}

/*!
  \param[in] objid Handle to write information to
*/
void initFileKnots::writeToHDF5Handle(hid_t objid) const {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("initFileKnots", "writeToHDF5Handle",
		       "Input handle is not valid");

  // Has range, limits, etc.
  hdf5utils::writeAttBool(objid, "has_init_param_range", has_range);
  hdf5utils::writeAttBool(objid, "has_param_lower_limits", has_lower_limits);
  hdf5utils::writeAttBool(objid, "has_param_upper_limits", has_upper_limits);
  
  // Ranges
  if ((nknots > 0) {
    if (has_range)
      hdf5utils::writeDataDoubles(objid, "param_init_range", nknots, range);

    if (has_lower_limits) {
      hdf5utils::writeDataBools(objid, "param_has_lowerlim", 
				nknots, has_lowlim);
      hdf5utils::writeDataDoubles(objid, "param_lowerlim", 
				  nknots, lowlim);
    }
    if (has_upper_limits) {
      hdf5utils::writeDataBools(objid, "param_has_upperlim", 
				nknots, has_uplim);
      hdf5utils::writeDataDoubles(objid, "param_upperlim", 
				  nknots, uplim);
    }
  }
}

/*!
  \param[inout] comm MPI communicator
  \param[in] dest Destination to send messages to
*/
void initFileKnots::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&nknots), 1, MPI_UNSIGNED, dest, 
	   pofd_mcmc::IFKSENDNKNOTS, comm);
  if (nknots != 0) {
    MPI_Send(knotpos, nknots, MPI_DOUBLE, dest, 
	     pofd_mcmc::IFKSENDKNOTPOS, comm);
    MPI_Send(knotval, nknots, MPI_DOUBLE, dest, 
	     pofd_mcmc::IFKSENDKNOTVAL, comm);
    MPI_Send(const_cast<bool*>(&has_range), 1, MPI::BOOL, dest,
	     pofd_mcmc::IFKHASRANGE, comm);
    if (has_range)
      MPI_Send(range, nknots, MPI_DOUBLE, dest, pofd_mcmc::IFKSENDRANGE, comm);
    MPI_Send(const_cast<bool*>(&has_lower_limits), 1, MPI::BOOL, dest,
	     pofd_mcmc::IFKHASLOWERLIMITS, comm);
    if (has_lower_limits) {
      MPI_Send(has_lowlim, nknots, MPI::BOOL, dest, 
	       pofd_mcmc::IFKSENDHASLOWLIM, comm);
      MPI_Send(lowlim, nknots, MPI_DOUBLE, dest, 
	       pofd_mcmc::IFKSENDLOWLIM, comm);
    }
    MPI_Send(const_cast<bool*>(&has_upper_limits), 1, MPI::BOOL, dest,
	     pofd_mcmc::IFKHASUPPERLIMITS, comm);
    if (has_upper_limits) {
      MPI_Send(has_uplim, nknots, MPI::BOOL, dest,
	       pofd_mcmc::IFKSENDHASUPLIM, comm);
      MPI_Send(uplim, nknots, MPI_DOUBLE, dest, pofd_mcmc::IFKSENDUPLIM, comm);
    }
  }
}

/*!
  \param[inout] comm MPI communicator
  \param[in] src Where messages will come from
*/
void initFileKnots::receiveCopy(MPI_Comm comm, int src) {
  //Delete everything for simplicity
  if (knotpos != NULL) delete[] knotpos;
  if (knotval != NULL) delete[] knotval;
  if (range != NULL) delete[] range;
  if (has_lowlim != NULL) delete[] has_lowlim;
  if (lowlim != NULL) delete[] lowlim;
  if (has_uplim != NULL) delete[] has_uplim;
  if (uplim != NULL) delete[] uplim;
  knotpos = knotval = range = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;

  MPI_Status Info;

  MPI_Recv(&nknots, 1, MPI_UNSIGNED, src, 
           pofd_mcmc::IFKSENDNKNOTS, comm, &Info);
  if (nknots > 0) {
    knotpos = new double[nknots];
    MPI_Recv(knotpos, nknots, MPI_DOUBLE, src, 
             pofd_mcmc::IFKSENDKNOTPOS, comm, &Info);
    knotval = new double[nknots];
    MPI_Recv(knotval, nknots, MPI_DOUBLE, src, 
             pofd_mcmc::IFKSENDKNOTVAL, comm, &Info);
    MPI_Recv(&has_range, 1, MPI::BOOL, src, pofd_mcmc::IFKHASRANGE,
             comm, &Info);
    if (has_range) {
      range = new double[nknots];
      MPI_Recv(range, nknots, MPI_DOUBLE, src,pofd_mcmc::IFKSENDRANGE,
               comm, &Info);
    }
    MPI_Recv(&has_lower_limits, 1, MPI::BOOL, src, 
             pofd_mcmc::IFKHASLOWERLIMITS, comm, &Info);
    if (has_lower_limits) {
      has_lowlim = new bool[nknots];
      MPI_Recv(has_lowlim, nknots, MPI::BOOL, src, 
               pofd_mcmc::IFKSENDHASLOWLIM, comm, &Info);
      lowlim = new double[nknots];
      MPI_Recv(lowlim, nknots, MPI_DOUBLE, src, 
               pofd_mcmc::IFKSENDLOWLIM, comm, &Info);
    }
    MPI_Recv(&has_upper_limits, 1, MPI::BOOL, src, 
             pofd_mcmc::IFKHASUPPERLIMITS, comm, &Info);
    if (has_upper_limits) {
      has_uplim = new bool[nknots];
      MPI_Recv(has_uplim, nknots, MPI::BOOL, src, 
               pofd_mcmc::IFKSENDHASUPLIM, comm, &Info);
      uplim = new double[nknots];
      MPI_Recv(uplim, nknots, MPI_DOUBLE, src, pofd_mcmc::IFKSENDUPLIM,
               comm, &Info);
    }
  }
}

