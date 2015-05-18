#ifndef __initFileKnots__
#define __initFileKnots__

#include "../include/numberCountsKnots.h"

/*!
  \brief A class to read in model specifications from init files, 1D case

  \ingroup Models
*/
class initFileKnots {
 private:
  //These are all c arrays rather than vectors for ease of compatability
  // with MPI

  //The first two (knotpos, knotval) are required
  unsigned int nknots; //!< Number of knots
  double* knotpos; //!< Positions of knots
  double* knotval; //!< Initial value center for knot value

  //These are optional
  bool has_range; //!< Has initial value range
  double* range; //!< Range for initial positions
  bool has_lower_limits; //!< Has some lower limit information
  bool* has_lowlim; //!< Knots have lower limit
  double* lowlim; //!< Value of lower limit
  bool has_upper_limits; //!< Has some upper limit information
  bool* has_uplim; //!< Knots have upper limit
  double* uplim; //!< Value of upper limit

  void checkLimitsDontCross() const; //!< Make sure limits are valid, throw if not
  void checkRange() const; //!< Make sure any range is valid, throw if not

 public:
  initFileKnots(); //!< Basic constructor
  initFileKnots(const std::string&, bool=false, bool=false); //!< Constructor with file read
  ~initFileKnots(); //!< Destructor

  unsigned int getNKnots() const { return nknots; } //!< Get number of knots
  dblpair getKnot(unsigned int idx) const; //!< Get knot pos and value

  void readFile(const std::string&, bool=false, bool=false); //!< Read file

  void getKnotPos(std::vector<double>&) const; //!< Gets the knot positions
  void getKnotVals(std::vector<double>&) const; //!< Gets the knot values
  void getKnotPos(numberCountsKnots&) const; //!< Sets knot locations in model
  double getKnotPos(unsigned int) const; //!< Get knot position
  double getKnotValue(unsigned int) const; //!< Get knot value

  void getParams(paramSet& p) const; //!< Sets param values to central values
  void generateRandomKnotValues(ran&, paramSet& pnew) const; //!< Seed knot values
  void generateRandomKnotValues(ran&, paramSet& pnew, const paramSet& pcen) const; //!< Seed knot values

  double getKnotRange(unsigned int) const; //!< Get knot range
  bool isKnotFixed(unsigned int) const; //!< Is a knot fixed?
  
  bool knotHasLowerLimit(unsigned int) const; //!< Does knot have a lower limit
  double getLowerLimit(unsigned int) const; //!< Get knot lower limit

  bool knotHasUpperLimit(unsigned int) const; //!< Does knot have a lower limit
  double getUpperLimit(unsigned int) const; //!< Get knot lower limit

  bool isValid(const paramSet&) const; //!< Checks if parameters are within allowed ranges
  
  void writeToHDF5Handle(hid_t objid) const; //!< Writes parameter limits to HDF5 handle

  void sendSelf(MPI_Comm, int dest) const; //!< Send self
  void receiveCopy(MPI_Comm, int src); //!< Receive
};

#endif
