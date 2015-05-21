#ifndef __initFileDoubleLogNormal__
#define __initFileDoubleLogNormal__

#include "../include/numberCountsDoubleLogNormal.h"

/*!
  \brief A class to read in model specifications from init files, 2D case

  \ingroup Models
*/
class initFileDoubleLogNormal {
 private:
  //These are all c arrays rather than vectors for ease of compatability
  // with MPI.  We also store everything as one vector

  //The first two (knotpos, knotval) are required
  unsigned int nknots; //!< Number of knots in band 1 number counts
  unsigned int nsigmas; //!< Number of knots in color sigma
  unsigned int noffsets; //!< Number of knots in color offset
  unsigned int sigmaidx; //!< Index of first sigma value
  unsigned int offsetidx; //!< Index of first offset value
  double* knotpos; //!< Positions of knots (all, including sigma and offsets)
  double* knotval; //!< Initial value center for knot value (all, inc. sigmas and offsets)

  //These are optional
  bool has_range; //!< Has initial value range
  double* range; //!< Range for initial positions
  bool has_lower_limits; //!< Has some lower limit information
  bool* has_lowlim; //!< Knots have lower limit
  double* lowlim; //!< Value of lower limit
  bool has_upper_limits; //!< Has some upper limit information
  bool* has_uplim; //!< Knots have upper limit
  double* uplim; //!< Value of upper limit

  void checkLimitsDontCross() const; //!< Make sure limits make sense, throw if not
  void checkRange() const; //!< Make sure param ranges make sense, throw if not

 public:
  initFileDoubleLogNormal(); //!< Basic constructor
  /*! \brief Constructor with file read */
  initFileDoubleLogNormal(const std::string&, bool=false, bool=false);
  ~initFileDoubleLogNormal(); //!< Destructor

  unsigned int getNKnots() const { return nknots; } //!< Get number of knots in band 1
  unsigned int getNSigmas() const { return nsigmas; } //!< Get number of knots in color model sigma
  unsigned int getNOffsets() const { return noffsets; } //!< Get number of knots in color model offset
  unsigned int getNTot() const { return nknots+nsigmas+noffsets; } //!< Get total number of model knots

  std::pair<double,double> getKnot(unsigned int idx) const; //!< Get band 1 knot pos and value
  std::pair<double,double> getSigma(unsigned int idx) const; //!< Get color sigma pos and value
  std::pair<double,double> getOffset(unsigned int idx) const; //!< Get color offset pos and value

  void readFile(const std::string&, bool=false, bool=false); //!< Read file

  void getKnotPos(std::vector<double>&) const; //!< Gets the knot positions for band 1
  void getKnotVals(std::vector<double>&) const; //!< Gets the knot values for band 1
  void getSigmaPos(std::vector<double>&) const; //!< Gets the knot positions for color model sigma
  void getSigmaVals(std::vector<double>&) const; //!< Gets the knot values for color model sigma
  void getOffsetPos(std::vector<double>&) const; //!< Gets the knot positions for color model offset
  void getOffsetVals(std::vector<double>&) const; //!< Gets the knot values for color model offset

  void getModelPositions(numberCountsDoubleLogNormal&) const; //!< Sets knot locations in model for all model components
  void getParams(paramSet& p) const; //!< Sets param values to central values
  void generateRandomKnotValues(ran&, paramSet& pnew) const; //!< Seed knot values
  void generateRandomKnotValues(ran&, paramSet& pnew, const paramSet& pcen) const; //!< Seed knot values

  double getKnotRange(unsigned int) const; //!< Get knot range
  bool isKnotFixed(unsigned int) const; //!< Is a knot fixed?

  bool knotHasLowerLimit(unsigned int) const; //!< Does knot have a lower limit
  double getLowerLimit(unsigned int) const; //!< Get knot lower limit

  bool knotHasUpperLimit(unsigned int) const; //!< Does knot have a lower limit
  double getUpperLimit(unsigned int) const; //!< Get knot lower limit

  double getMinSigma() const; //!< Use sigma knot lower limits to get min sigma

  bool isValid(const paramSet&) const; //!< Checks if parameters are within allowed ranges

  void writeToHDF5Handle(hid_t) const; //!< Write to HDF5 handle

  void sendSelf(MPI_Comm, int dest) const; //!< Send self
  void receiveCopy(MPI_Comm, int src); //!< Receive
};

#endif
