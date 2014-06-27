//numberCountsKnots.h

#ifndef __numberCountsKnots__
#define __numberCountsKnots__

#include<vector>
#include<utility>

#include "../include/numberCounts.h"
#include "../include/ran.h"
#include "../include/paramSet.h"

/*!
  \brief Number counts with knots abstract base class, 1D case
 
  The first nknots parameters in paramSet are the log knot
  values, anything after that are special (sigma, mean offset,
  possibly clustering parameters).

  The user inputs knot values in \f$\log_{10}\f$, but they are
  stored internally as base 2 logarithms.

  \ingroup Models
 */
class numberCountsKnots : public numberCounts {
 protected:
  unsigned int nknots; //!< Number of knots
  double* knots; //!< Location of knots

  double* logknotvals; //!< Log2 values of differential number counts at knots
  bool knotvals_loaded; //!< Are knot values loaded?

  virtual void setNKnots(unsigned int n); //!< Sets number of knots

  dblpair getRRangeInternal(const beam& bm) const throw(affineExcept); //!< Unchecked version of getRRange

 public:
  numberCountsKnots(); //!< Default constructor
  explicit numberCountsKnots(unsigned int); //!< Constructor with number of knots
  numberCountsKnots(const std::vector<float>&); //!< Vector constructor
  numberCountsKnots(const std::vector<double>&); //!< Vector constructor
  numberCountsKnots(unsigned int, const float* const); //!< C array constructor
  numberCountsKnots(unsigned int, const double* const); //!< C array constructor
  numberCountsKnots(const numberCountsKnots&); //!< Copy constructor
  ~numberCountsKnots(); //!< Destructor

  /*! \brief Load knot positions into vector */
  virtual void getKnotPositions(std::vector<double>&) const; 
  /*! \brief Set knot positions, vector version */
  virtual void setKnotPositions(const std::vector<double>&);
  /*! \brief Set knot positions, c array version */
  virtual void setKnotPositions(unsigned int, const double* const);
  /*! \brief Set knot positions, vector version */
  virtual void setKnotPositions(const std::vector<float>&);
  /*! \brief Set knot positions, c array version */
  virtual void setKnotPositions(unsigned int, const float* const);

  virtual void getParams(paramSet&) const; //!< Return current parameters
  virtual void setParams(const paramSet&); //!< Set parameters

  unsigned int getNKnots() const { return nknots; } //!< Returns number of knots
  /*! \brief Get position of specified knot */
  double getKnotPos(unsigned int i) const { return knots[i]; }

  /*! \brief Get knot position and log value */
  dblpair getLogKnot(unsigned int i) const {
    return std::pair<double,double>(knots[i], logknotvals[i]); }

  virtual bool isValid() const; //!< Are model parameters valid

  /*! \brief Distance between two parameter sets over params model cares about*/
  virtual float paramRelativeDistance(const paramSet& p1, const paramSet& p2) 
    const throw(affineExcept);

  double getMaxFlux() const; //!< Maximum flux supported by model
  double getMinFlux() const; //!< Minimum flux supported by model

  /*! \brief Get range over which R is expected to be nonzero */
  dblpair getRRange(const beam&) const throw(affineExcept);

  virtual void writeToHDF5Handle(hid_t objid) const; //!< Output to HDF5
  virtual bool writeToStream(std::ostream& os) const; //<! Output to stream

  virtual void sendSelf(MPI_Comm, int dest) const; //!< Send self
  virtual void receiveCopy(MPI_Comm, int src); //!< Receive
};

/*! \brief Write to stream operator */
std::ostream& operator<<(std::ostream& os, const numberCountsKnots& b);

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
