#include<sstream>
#include<limits>

#include "../include/initFileDoubleLogNormal.h"
#include "../include/global_settings.h"
#include "../include/hdf5utils.h"
#include "../include/affineExcept.h"

initFileDoubleLogNormal::initFileDoubleLogNormal() : 
  nknots(0), nsigmas(0), noffsets(0), sigmaidx(0), offsetidx(0),
  has_range(false), has_lower_limits(false), has_upper_limits(false) {

  knotpos = knotval = range = lowlim = uplim = nullptr;
  has_lowlim = has_uplim = nullptr;
}

/*!
  \param[in] flname File to read from
  \param[in] read_range      Read in (and require) knot ranges
  \param[in] read_limits     Try to read limits; this will turn on 
                             require_range

  See initFileDoubleLogNormal::readFile for details of the file format
*/
initFileDoubleLogNormal::initFileDoubleLogNormal(const std::string& flname, 
                                                 bool read_range, 
                                                 bool read_limits) : 
  nknots(0), nsigmas(0), noffsets(0), sigmaidx(0), offsetidx(0),
  has_range(false), has_lower_limits(false), has_upper_limits(false) {

  knotpos = knotval = range = lowlim = uplim = nullptr;
  has_lowlim = has_uplim = nullptr;

  readFile(flname, read_range, read_limits);
}

initFileDoubleLogNormal::~initFileDoubleLogNormal() {
  if (knotpos != nullptr) delete[] knotpos;
  if (knotval != nullptr) delete[] knotval;
  if (range != nullptr) delete[] range;
  if (has_lowlim != nullptr) delete[] has_lowlim;
  if (lowlim != nullptr) delete[] lowlim;
  if (has_uplim != nullptr) delete[] has_uplim;
  if (uplim != nullptr) delete[] uplim;
}

/*!
  \param[in] flname File to read from
  \param[in] read_range      Read in (and require) knot ranges
  \param[in] read_limits     Try to read limits; this will turn on require_range

  The file format is: 
  first, a line giving the n1 ns no
  where n1 is the number of knots in the band 1 model, ns is the number
  of knots in the color sigma model, and no is the number of knots in the
  color offset model.  This is followed by n1+ns+no lines of the form
  knotpos   knotval   [range [ lowlim [ uplim ]]
  So knotpos, knotval are always required
  range is optionally required if read_range is set.  If used to generate
   points, they are generated within range/2 on either side uniformly
   around knotval.
  lowlim and uplim may be present, and are looked for if read_limits is set.
    range must also be present, and the first element found is lowlim.
    If another is also found, it is interpreted as uplim
*/
void initFileDoubleLogNormal::readFile(const std::string& flname, 
                                       bool read_range, bool read_limits) {
  if (read_limits) read_range = true;

  //Clear old data
  nknots = 0; 
  nsigmas = 0;
  noffsets = 0;
  sigmaidx = 0;
  offsetidx = 0;
  if (knotpos != nullptr) { delete[] knotpos; knotpos = nullptr; }
  if (knotval != nullptr) { delete[] knotval; knotval = nullptr; }
  if (range != nullptr) { delete[] range; range = nullptr; }
  if (has_lowlim != nullptr) { delete[] has_lowlim; has_lowlim = nullptr; }
  if (lowlim != nullptr) { delete[] lowlim; lowlim = nullptr; }
  if (has_uplim != nullptr) { delete[] has_uplim; has_uplim = nullptr; }
  if (uplim != nullptr) { delete[] uplim; uplim = nullptr; }
  has_range = has_lower_limits = has_upper_limits = false;

  //Figure out how many elements we require
  unsigned int nreq = 2; //Pos, value
  if (read_range) nreq += 1; //Range -- read limits not required ever

  unsigned int nk, ns, no; //Number of knots, sigmas, offsets
  std::vector<double> wvec1, wvec2, wvec3, wvec4, wvec5;
  std::vector<bool> h4, h5;
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currval;

  std::ifstream initfs(flname.c_str());
  if (!initfs) {
    initfs.close();
    std::stringstream errmsg;
    errmsg << "Unable to open file:" << flname << std::endl;
    throw affineExcept("initFileDoubleLogNormal", "readFile", errmsg.str());
  }

  //Read in number
  initfs >> nk >> ns >> no;
  if (nk < 2) {
    initfs.close();
    throw affineExcept("initFileDoubleLogNormal", "readFile",
                       "Need at least 2 band 1 knots");
  }
  if (ns < 1) {
    initfs.close();
    throw affineExcept("initFileDoubleLogNormal", "readFile",
                       "Need at least one sigma color model knot");

  }
  if (no < 1) {
    initfs.close();
    throw affineExcept("initFileDoubleLogNormal", "readFile",
                       "Need at least one offset color model knot");
  }
  nknots = sigmaidx = nk; //These could be combined, but meh
  nsigmas = ns;
  noffsets = no;
  offsetidx = nk + ns;
  
  //Read in values
  while (!initfs.eof()) {
    std::getline(initfs, line);
    if (line[0] == '#') continue; //Comment
    utility::stringwords(line, words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < nreq) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currval;
    wvec1.push_back(currval);
    str.str(words[1]); str.clear(); str >> currval;
    wvec2.push_back(currval); 

    //Now, all the options
    if (read_range) {
      has_range = true;
      str.str(words[2]); str.clear(); str >> currval;
      wvec3.push_back(currval);
    }
    if (read_limits) { 
      if (words.size() > 3) {
        //Ignore limits if range is zero -- what would be the point?
        if (has_range && wvec3.back() == 0) {
          h4.push_back(false);
          wvec4.push_back(std::numeric_limits<double>::quiet_NaN());
          h5.push_back(false);
          wvec5.push_back(std::numeric_limits<double>::quiet_NaN());
        } else {
          has_lower_limits = true;
          str.str(words[3]); str.clear(); str >> currval;
          h4.push_back(true);
          wvec4.push_back(currval);
          if (words.size() > 4) {
            has_upper_limits = true;
            str.str(words[4]); str.clear(); str >> currval;
            h5.push_back(true);
            wvec5.push_back(currval);
          } else {
            h5.push_back(false);
            wvec5.push_back(std::numeric_limits<double>::quiet_NaN());
          }
        }
      } else {
        h4.push_back(false);
        wvec4.push_back(std::numeric_limits<double>::quiet_NaN());
        h5.push_back(false);
        wvec5.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }
  }
  initfs.close();

  unsigned int ntot = nk + ns + no;
  if (wvec1.size() != ntot) {
    std::stringstream errstr;
    errstr << "Expected " << ntot << " values, got: " 
           << wvec1.size();
    throw affineExcept("initFileDoubleLogNormal", "readFile", errstr.str());
  }

  //Copy into vectors
  knotpos = new double[ntot];
  for (unsigned int i = 0; i < ntot; ++i) knotpos[i] = wvec1[i];
  knotval = new double[ntot];
  for (unsigned int i = 0; i < ntot; ++i) knotval[i] = wvec2[i];
  if (has_range) {
    range = new double[ntot];
    for (unsigned int i = 0; i < ntot; ++i) range[i] = wvec3[i];
  }
  if (has_lower_limits) {
    has_lowlim = new bool[ntot];
    lowlim = new double[ntot];
    for (unsigned int i = 0; i < ntot; ++i) has_lowlim[i] = h4[i];
    for (unsigned int i = 0; i < ntot; ++i) lowlim[i] = wvec4[i];
  }
  if (has_upper_limits) {
    has_uplim = new bool[ntot];
    uplim = new double[ntot];
    for (unsigned int i = 0; i < ntot; ++i) has_uplim[i] = h5[i];
    for (unsigned int i = 0; i < ntot; ++i) uplim[i] = wvec5[i];
  }

  // Check limits and ranges
  checkLimitsDontCross();
  checkRange();
}

void initFileDoubleLogNormal::checkLimitsDontCross() const {
  //Make sure lower/upper limits don't cross
  if (!(has_lower_limits && has_upper_limits)) return;

  unsigned int ntot = nknots + nsigmas + noffsets;
  if (ntot == 0) return;

  for (unsigned int i = 0; i < ntot; ++i) {
    if (!(has_uplim[i] && has_lowlim[i])) continue;
    if (uplim[i] < lowlim[i]) {
      std::stringstream errstr;
      errstr << "Lower/Upper limits cross at index: " << i << std::endl;
      errstr << " Lower limit: " << lowlim[i] 
             << " Upper limit: " << uplim[i];
      throw affineExcept("initFileDoubleLogNormal", "checkLimitsDontCross", 
                         errstr.str());
    }
    if ((range[i] > 0.) && (uplim[i] == lowlim[i])) {
      std::stringstream errstr;
      errstr << "Lower/Upper limits meet at index: " << i 
             << " but range is not zero" << std::endl;
      errstr << " Lower limit: " << lowlim[i] << " Upper limit: " << uplim[i]
             << " range: " << range[i];
      throw affineExcept("initFileDoubleLogNormal", "checkLimitsDontCross", 
                         errstr.str());
    }
  }
}

void initFileDoubleLogNormal::checkRange() const {
  //Make sure that if ranges is 0 then the mean value falls within
  // the range of any limits
  
  // Quick returns
  if (!has_range) return;
  if (!(has_lower_limits || has_upper_limits)) return;

  unsigned int ntot = nknots + nsigmas + noffsets;
  if (ntot == 0) return;

  for (unsigned int i = 0; i < ntot; ++i)
    if (range[i] == 0) {
      if (has_lower_limits && has_lowlim[i] && (knotval[i] < lowlim[i])) {
        std::stringstream errstr;
        errstr << "At knot " << i << " range is zero but mean value "
               << knotval[i] << std::endl << " lies below lower limit "
               << lowlim[i];
        throw affineExcept("initFileDoubleLogNormal", "checkRange", 
                           errstr.str());
      }
      if (has_upper_limits && has_uplim[i] && (knotval[i] > uplim[i])) {
        std::stringstream errstr;
        errstr << "At knot " << i << " range is zero but mean value "
               << knotval[i] << std::endl << " lies above upper limit "
               << uplim[i];
        throw affineExcept("initFileDoubleLogNormal", "checkRange", 
                           errstr.str());
      }
    }
}

/*!
  \param[in] idx Knot index
  \returns Tuple of knot position, value
*/
std::pair<double, double> 
initFileDoubleLogNormal::getKnot(unsigned int idx) const {
  if (nknots == 0)
    throw affineExcept("initFileDoubleLogNormal", "getKnot",
                       "No knot information read in");
  if (idx >= nknots)
    throw affineExcept("initFileDoubleLogNormal", "getKnot",
                       "Invalid knot index");
  return std::make_pair(knotpos[idx], knotval[idx]);
}

/*!
  \param[in] idx Sigma index
  \returns Tuple of sigma knot position, value
*/
std::pair<double, double> 
initFileDoubleLogNormal::getSigma(unsigned int idx) const {
  if (nsigmas == 0)
    throw affineExcept("initFileDoubleLogNormal", "getSigma",
                       "No sigma information read in");
  if (idx >= nsigmas)
    throw affineExcept("initFileDoubleLogNormal", "getSigma",
                       "Invalid sigma index");
  return std::make_pair(knotpos[idx+sigmaidx], knotval[idx+sigmaidx]);
}

/*!
  \param[in] idx Offset index
  \returns Tuple of offset knot position, value
*/
std::pair<double, double> 
initFileDoubleLogNormal::getOffset(unsigned int idx) const {
  if (noffsets == 0)
    throw affineExcept("initFileDoubleLogNormal", "getOffset",
                       "No offset information read in");
  if (idx >= noffsets)
    throw affineExcept("initFileDoubleLogNormal", "getOffset",
                       "Invalid offset index");
  return std::make_pair(knotpos[idx+offsetidx], knotval[idx+offsetidx]);
}


/*
  \param[inout] model Modified on output; knot positions are set for
                      band 1 model and color model

  This will change the number of knots in the model if they
  don't match.
 */
void 
initFileDoubleLogNormal::getModelPositions(numberCountsDoubleLogNormal& model) 
  const {
  if (nknots + nsigmas + noffsets == 0)
    throw affineExcept("initFileDoubleLogNormal", "getModelPositions",
                       "No knot information read in");
  model.setKnotPositions(nknots, knotpos);
  model.setSigmaPositions(nsigmas, knotpos + sigmaidx);
  model.setOffsetPositions(noffsets, knotpos + offsetidx);
}

/*
  \param[out] kp Set to knot positions on output for band 1 number counts
*/
void initFileDoubleLogNormal::getKnotPos(std::vector<double>& kp) const {
  if (nknots == 0)
    throw affineExcept("initFileDoubleLogNormal", "getKnotPos",
                       "No knot information read in");
  kp.resize(nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    kp[i] = knotpos[i];

}

/*
  \param[out] kp Set to knot positions on output for band 1 number counts
*/
void initFileDoubleLogNormal::getKnotVals(std::vector<double>& kv) const {
  if (nknots == 0)
    throw affineExcept("initFileDoubleLogNormal", "getKnotVals",
                       "No knot information read in");
  kv.resize(nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    kv[i] = knotval[i];
}

/*
  \param[out] kp Set to knot positions on output for color model sigma
*/
void initFileDoubleLogNormal::getSigmaPos(std::vector<double>& kp) const {
  if (nsigmas == 0)
    throw affineExcept("initFileDoubleLogNormal", "getSigmaPos",
                       "No sigma information read in");
  kp.resize(nsigmas);
  for (unsigned int i = 0; i < nsigmas; ++i)
    kp[i] = knotpos[i+sigmaidx];
}

/*
  \param[out] kp Set to knot values on output for color model sigma
*/
void initFileDoubleLogNormal::getSigmaVals(std::vector<double>& kv) const {
  if (nsigmas == 0)
    throw affineExcept("initFileDoubleLogNormal", "getSigmaVals",
                       "No sigma information read in");
  kv.resize(nsigmas);
  for (unsigned int i = 0; i < nsigmas; ++i)
    kv[i] = knotval[i+sigmaidx];
}

/*
  \param[out] kp Set to knot positions on output for color model offset
*/
void initFileDoubleLogNormal::getOffsetPos(std::vector<double>& kp) const {
  if (noffsets == 0)
    throw affineExcept("initFileDoubleLogNormal", "getOffsetPos",
                       "No offset information read in");
  kp.resize(noffsets);
  for (unsigned int i = 0; i < noffsets; ++i)
    kp[i] = knotpos[i+offsetidx];
}

/*
  \param[out] kp Set to knot values on output for color model offset
*/
void initFileDoubleLogNormal::getOffsetVals(std::vector<double>& kv) const {
  if (noffsets == 0)
    throw affineExcept("initFileDoubleLogNormal", "getOffsetVals",
                       "No offset information read in");
  kv.resize(noffsets);
  for (unsigned int i = 0; i < noffsets; ++i)
    kv[i] = knotval[i+offsetidx];
}

/*
  \param[out] p Filled on output with the mean knot values.

  This only fills the first nknots + nsigmas + noffsets parameters
*/
void initFileDoubleLogNormal::getParams(paramSet& p) const {
  if (nknots == 0)
    throw affineExcept("initFileDoubleLogNormal", "getParams",
                       "No information loaded");
  if (p.getNParams() < nknots+noffsets+nsigmas)
    throw affineExcept("initFileDoubleLogNormal", "getParams",
                       "Not enough space in provided paramSet");
  for (unsigned int i = 0; i < nknots+nsigmas+noffsets; ++i)
    p[i] = knotval[i];
}

/*
  \param[in] rangen Random number generator
  \param[out] p Filled on output with the mean knot values.

  This only fills the first nknots+nsigmas+noffsets parameters.  It uses
  the central values from the initialization file
*/
void initFileDoubleLogNormal::generateRandomKnotValues(ran& rangen, 
                                                       paramSet& pnew) const {
  paramSet pcen(pnew.getNParams());
  getParams(pcen); //Load central values into pcen
  generateRandomKnotValues(rangen, pnew, pcen); //Get new values
}

/*
  \param[in] rangen Random number generator
  \param[out] pnew New parameter set generated
  \param[in] pcen  Central parameter values

  This only fills the first nknots parameters.  This version
  allows the caller to use different central values than the ones
  from the initial file, but keep the sigmas, limits, etc.
*/
void initFileDoubleLogNormal::
generateRandomKnotValues(ran& rangen, paramSet& pnew, 
                         const paramSet& pcen) const {
  const unsigned int maxiters = 1000; //Maximum number of generation attempts

  unsigned int ntot = nknots + noffsets + nsigmas;
  if (ntot == 0)
    throw affineExcept("initFileDoubleLogNormal", "generateRandomKnotValues",
                       "No knot information read in");
    
  //Make sure p is big enough; don't resize, complain
  // Note that setting the sigma multipliers is not handled here,
  // so we don't test for those slots
  if (pnew.getNParams() < ntot)
    throw affineExcept("initFileDoubleLogNormal", "generateRandomKnotValues",
                       "Not enough space in provided paramSet");

  if (pcen.getNParams() < ntot)
    throw affineExcept("initFileDoubleLogNormal", "generateRandomKnotValues",
                       "Not enough params in central param Set");

  //Deal with simple case -- everything fixed
  //So just return pcen -- but check it first
  if (!has_range) {
    if (has_lower_limits)
      for (unsigned int i = 0; i < ntot; ++i)
        if (has_lowlim[i] && (pcen[i] < lowlim[i])) {
          std::stringstream errstr;
          errstr << "For parameter " << i << " user provided central value "
                 << pcen[i] << " is below lower limit " << lowlim[i];
          throw affineExcept("initFileKnots", "generateRandomKnotValues",
                             errstr.str());
        }
    if (has_upper_limits)
      for (unsigned int i = 0; i < ntot; ++i)
        if (has_uplim[i] && (pcen[i] > uplim[i])) {
          std::stringstream errstr;
          errstr << "For parameter " << i << " user provided central value "
                 << pcen[i] << " is above upper limit " << uplim[i];
          throw affineExcept("initFileKnots", "generateRandomKnotValues",
                             errstr.str());
        }
    for (unsigned int i = 0; i < ntot; ++i)
      pnew[i] = pcen[i];
    return;
  }

  //Now we have at least some ranges
  //The simple case is if there are no limits.  If there are, we will
  // have to do trials.
  if (!(has_lower_limits || has_upper_limits)) {
    for (unsigned int i = 0; i < ntot; ++i)
      pnew[i] = rangen.flt() * range[i] + (pcen[i] - 0.5 * range[i]);
  } else {
    //Both sigmas and limits
    bool goodval;
    double trialval;
    unsigned int iters;
    for (unsigned int i = 0; i < ntot; ++i) {
      if (range[i] > 0) {
        //Some sanity checks
        if (has_lowlim[i] && (lowlim[i] > pcen[i] + range[i])) {
          std::stringstream errstr;
          errstr << "Lower limit is too far above central value; will not be "
                 << "able to" << std::endl << "generate value for param idx: "
                 << i << " with lowlim: " << lowlim[i] << " central: " 
                 << pcen[i] << " range: " << range[i];
          throw affineExcept("initFileDoubleLogNormal", 
                             "generateRandomKnotValues", errstr.str());
        }
        if (has_uplim[i] && (uplim[i] < pcen[i] - range[i])) {
          std::stringstream errstr;
          errstr << "Upper limit is too far below central value; will not be"
                 << " able to " << std::endl << "generate value for param idx: "
                 << i << " with uplim: " << uplim[i] << " central: " 
                 << pcen[i] << " range: " << range[i];
          throw affineExcept("initFileDoubleLogNormal",
                             "generateRandomKnotValues", errstr.str());
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
                throw affineExcept("initFileDoubleLogNormal",
                                   "generateRandomKnotValues", errstr.str());
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
              throw affineExcept("initFileDoubleLogNormal",
                                 "generateRandomKnotValues", errstr.str());
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
              throw affineExcept("initFileDoubleLogNormal",
                                 "generateRandomKnotValues", errstr.str());
            }
            trialval = rangen.flt() * range[i] + (pcen[i] - 0.5 * range[i]);
            if (trialval <= uplim[i]) goodval = true;
            ++iters;
          }
          pnew[i] = trialval;
        } else {
          //No limit, easy cakes
          trialval = rangen.flt() * range[i] + (pcen[i] - 0.5 * range[i]);
        }
      } else {
        //Range is 0.  Check to make sure this is within the limits
        if (has_lowlim[i] && (pcen[i] < lowlim[i])) {
          std::stringstream errstr;
          errstr << "For parameter " << i << " user provided central value "
                 << pcen[i] << " is below lower limit " << lowlim[i]
                 << " and range is 0.";
          throw affineExcept("initFileDoubleLogNormal", 
                             "generateRandomKnotValues", errstr.str());
        }
        if (has_uplim[i] && (pcen[i] > uplim[i])) {
          std::stringstream errstr;
          errstr << "For parameter " << i << " user provided central value "
                 << pcen[i] << " is above upper limit " << uplim[i]
                 << " and range is 0.";
          throw affineExcept("initFileDoubleLogNormal", 
                             "generateRandomKnotValues", errstr.str());
        }
        pnew[i] = pcen[i];
      }
    }
  }

}

/*!
  \param[in] idx Which parameter index to get
  \returns Knot sigma at idx
  
  Note that idx ranges over all the parameters -- including
  the color model.  The valid range is thus [0,nknots + nsigmas + noffsets).
*/
double initFileDoubleLogNormal::getKnotRange(unsigned int idx) const {
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal",
                       "getKnotSigma", "Invalid knot index");
  if (!has_range) return std::numeric_limits<double>::quiet_NaN();
  return range[idx];
}

/*!
  \param[in] idx Which paraemter index to return value for
  
  \returns True if the knot is fixed, false if not.

  Note that idx ranges over all the parameters -- including
  the color model. The valid range is thus [0,nknots + nsigmas + noffsets).
*/
bool initFileDoubleLogNormal::isKnotFixed(unsigned int idx) const {
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal", "isKnotFixed",
                       "Invalid knot index");
  if (!has_range) return false;
  if (range[idx] == 0) return true;
  return false;
}

/*
  \param[in] idx Knot index
  \returns True if knot has lower limit, otherwise false
*/
bool initFileDoubleLogNormal::knotHasLowerLimit(unsigned int idx) const {
  if (!has_lower_limits) return false;
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal", "knotHasLowerLimit",
                       "Invalid knot index");
  return has_lowlim[idx];
}

/*
  \param[in] idx Knot index
  \returns Lower limit on knot, or NaN if none
 */
double initFileDoubleLogNormal::getLowerLimit(unsigned int idx) const {
  if (!has_lower_limits) return std::numeric_limits<double>::quiet_NaN();
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal", "getLowerLimit",
                       "Invalid knot index");
  if (!has_lowlim[idx]) return std::numeric_limits<double>::quiet_NaN();
  return lowlim[idx];
}

/*
  \param[in] idx Knot index
  \returns True if knot has upper limit, otherwise false
*/
bool initFileDoubleLogNormal::knotHasUpperLimit(unsigned int idx) const {
  if (!has_upper_limits) return false;
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal", "knotHasUpperLimit",
                       "Invalid knot index");
  return has_uplim[idx];
}


/*
  \param[in] idx Knot index
  \returns Upper limit on knot, or NaN if none
*/
double initFileDoubleLogNormal::getUpperLimit(unsigned int idx) const {
  if (!has_upper_limits) return std::numeric_limits<double>::quiet_NaN();
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal", "getUpperLimit",
                       "Invalid knot index");
  if (!has_uplim[idx]) return std::numeric_limits<double>::quiet_NaN();
  return uplim[idx];
}

/*
  \returns Min Sigma value.  This is either 0.0, or 1/2 the smallest
    sigma model lower limit if all sigma knots have such a limit.
*/
double initFileDoubleLogNormal::getMinSigma() const {
  if (!has_lower_limits) return 0;
  bool all_sig_lim = has_lowlim[nknots];
  for (unsigned int i = nknots+1; i < nknots + nsigmas; ++i)
    all_sig_lim &= has_lowlim[i];
  if (!all_sig_lim) return 0;
  double min_sig = lowlim[nknots];
  for (unsigned int i = nknots+1; i < nknots + nsigmas; ++i)
    if (min_sig > lowlim[i]) min_sig = lowlim[i];
  if (min_sig <= 0) return 0;
  return 0.5 * min_sig;
}

/*!
  \param[in] p Parameters to test
  \returns True if p is valid (within limits)
*/
bool initFileDoubleLogNormal::isValid(const paramSet& p) const {
  if (! (has_lower_limits || has_upper_limits)) return true;
  unsigned int ntot = nknots + nsigmas + noffsets;
  if (p.getNParams() < ntot)
    throw affineExcept("initFileDoubleLogNormal", "isValid",
                       "Not enough params in paramSet to test validity");
  double val;
  for (unsigned int i = 0; i < ntot; ++i) {
    val = p[i];
    if (has_lowlim[i] && (val < lowlim[i])) return false;
    if (has_uplim[i] && (val > uplim[i])) return false;
  }
  return true;
}

/*!
  \param[in] objid Handle to write information to
*/
void initFileDoubleLogNormal::writeToHDF5Handle(hid_t objid) const {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("initFileDoubleLogNormal", "writeToHDF5Handle",
                       "Input handle is not valid");

  // Has range, limits, etc.
  hdf5utils::writeAttBool(objid, "HasInitParamRange", has_range);
  hdf5utils::writeAttBool(objid, "HasParamLowerLimits", has_lower_limits);
  hdf5utils::writeAttBool(objid, "HasParamUpperLimits", has_upper_limits);
  
  // Ranges
  unsigned int ntotknots = nknots + nsigmas + noffsets;
  if (ntotknots > 0) {
    if (has_range) 
      hdf5utils::writeDataDoubles(objid, "ParamInitRange", ntotknots, range);

    if (has_lower_limits) {
      hdf5utils::writeDataBools(objid, "ParamHasLowerlim", 
                                ntotknots, has_lowlim);
      hdf5utils::writeDataDoubles(objid, "ParamLowerlim", 
                                  ntotknots, lowlim);
    }

    if (has_upper_limits) {
      hdf5utils::writeDataBools(objid, "ParamHasUpperlim", 
                                ntotknots, has_uplim);
      hdf5utils::writeDataDoubles(objid, "ParamUpperlim", 
                                  ntotknots, uplim);
    }
  }
}

/*!
  \param[in] comm MPI communication handle
  \param[in] dest Destination for send operation
*/
void initFileDoubleLogNormal::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&nknots), 1, MPI_UNSIGNED, dest,
           pofd_mcmc::IFDLNSENDNKNOTS, comm);
  MPI_Send(const_cast<unsigned int*>(&nsigmas), 1, MPI_UNSIGNED, dest,
           pofd_mcmc::IFDLNSENDNSIGMAS, comm);
  MPI_Send(const_cast<unsigned int*>(&noffsets), 1, MPI_UNSIGNED, dest,
           pofd_mcmc::IFDLNSENDNOFFSETS, comm);
  unsigned int ntot = nknots+nsigmas+noffsets;
  if (ntot != 0) {
    MPI_Send(knotpos, ntot, MPI_DOUBLE, dest, 
             pofd_mcmc::IFDLNSENDKNOTPOS, comm);
    MPI_Send(knotval, ntot, MPI_DOUBLE, dest, 
             pofd_mcmc::IFDLNSENDKNOTVAL, comm);
    MPI_Send(const_cast<bool*>(&has_range), 1, MPI::BOOL, dest, 
             pofd_mcmc::IFDLNHASRANGE, comm);
    if (has_range)
      MPI_Send(range, ntot, MPI_DOUBLE, dest, pofd_mcmc::IFDLNSENDRANGE, comm);
    MPI_Send(const_cast<bool*>(&has_lower_limits), 1, MPI::BOOL, dest,
             pofd_mcmc::IFDLNHASLOWERLIMITS, comm);
    if (has_lower_limits) {
      MPI_Send(has_lowlim, ntot, MPI::BOOL, dest,
               pofd_mcmc::IFDLNSENDHASLOWLIM, comm);
      MPI_Send(lowlim, ntot, MPI_DOUBLE, dest, 
               pofd_mcmc::IFDLNSENDLOWLIM, comm);
    }
    MPI_Send(const_cast<bool*>(&has_upper_limits), 1, MPI::BOOL, dest,
             pofd_mcmc::IFDLNHASUPPERLIMITS, comm);
    if (has_upper_limits) {
      MPI_Send(has_uplim, ntot, MPI::BOOL, dest,
               pofd_mcmc::IFDLNSENDHASUPLIM, comm);
      MPI_Send(uplim, ntot, MPI_DOUBLE, dest, pofd_mcmc::IFDLNSENDUPLIM, comm);
    }
  }
}

/*!
  \param[inout] comm MPI communicator
  \param[in] src Where messages will come from
*/
void initFileDoubleLogNormal::receiveCopy(MPI_Comm comm, int src) {
  //Delete everything for simplicity
  if (knotpos != nullptr) delete[] knotpos;
  if (knotval != nullptr) delete[] knotval;
  if (range != nullptr) delete[] range;
  if (has_lowlim != nullptr) delete[] has_lowlim;
  if (lowlim != nullptr) delete[] lowlim;
  if (has_uplim != nullptr) delete[] has_uplim;
  if (uplim != nullptr) delete[] uplim;
  knotpos = knotval = range = lowlim = uplim = nullptr;
  has_lowlim = has_uplim = nullptr;

  MPI_Status Info;
  MPI_Recv(&nknots, 1, MPI_UNSIGNED, src, pofd_mcmc::IFDLNSENDNKNOTS,
           comm, &Info);
  MPI_Recv(&nsigmas, 1, MPI_UNSIGNED, src, pofd_mcmc::IFDLNSENDNSIGMAS,
           comm, &Info);
  MPI_Recv(&noffsets, 1, MPI_UNSIGNED, src, pofd_mcmc::IFDLNSENDNOFFSETS,
           comm, &Info);
  unsigned int ntot = nknots + nsigmas + noffsets;
  sigmaidx = nknots;
  offsetidx = nknots + nsigmas;
  if (ntot > 0) {
    knotpos = new double[ntot];
    MPI_Recv(knotpos, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDKNOTPOS,
             comm, &Info);
    knotval = new double[ntot];
    MPI_Recv(knotval, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDKNOTVAL,
             comm, &Info);
    MPI_Recv(&has_range, 1, MPI::BOOL, src, pofd_mcmc::IFDLNHASRANGE,
             comm, &Info);
    if (has_range) {
      range = new double[ntot];
      MPI_Recv(range, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDRANGE,
               comm, &Info);
    }
    MPI_Recv(&has_lower_limits, 1, MPI::BOOL, src,
             pofd_mcmc::IFDLNHASLOWERLIMITS, comm, &Info);
    if (has_lower_limits) {
      has_lowlim = new bool[ntot];
      MPI_Recv(has_lowlim, ntot, MPI::BOOL, src, 
               pofd_mcmc::IFDLNSENDHASLOWLIM, comm, &Info);
      lowlim = new double[ntot];
      MPI_Recv(lowlim, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDLOWLIM,
               comm, &Info);
    }
    MPI_Recv(&has_upper_limits, 1, MPI::BOOL, src, 
             pofd_mcmc::IFDLNHASUPPERLIMITS, comm, &Info);
    if (has_upper_limits) {
      has_uplim = new bool[ntot];
      MPI_Recv(has_uplim, ntot, MPI::BOOL, src, pofd_mcmc::IFDLNSENDHASUPLIM,
               comm, &Info);
      uplim = new double[ntot];
      MPI_Recv(uplim, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDUPLIM,
               comm, &Info);
    }
  }
}
