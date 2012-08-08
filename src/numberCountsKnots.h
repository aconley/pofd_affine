//numberCountsKnots.h

#ifndef __numberCountsKnots__
#define __numberCountsKnots__

#include<vector>
#include<utility>

#include<numberCounts.h>
#include<ran.h>
#include<paramSet.h>

/*!
  \brief Number counts with knots abstract base class
  \ingroup Models
 
  The first nknots parameters in paramSet are the log knot
  values, anything after that are special (sigma, mean offset,
  possibly clustering parameters).

  The user inputs knot values in \f$\log_{10}\f$, but they are
  stored internally as base 2 logarithms.
 */
class numberCountsKnots : public numberCounts {
 protected:
  unsigned int nknots; //!< Number of knots
  double* knots; //!< Location of knots

  double* logknotvals; //!< Log2 values of differential number counts at knots
  bool knotvals_loaded; //!< Are knot values loaded?

  virtual void setNKnots(unsigned int n); //!< Sets number of knots

 public:
  numberCountsKnots(); //!< Default constructor
  explicit numberCountsKnots(unsigned int); //!< Constructor with number of knots
  numberCountsKnots( const std::vector<double>& ); //!< Vector constructor
  numberCountsKnots( unsigned int, const double* const); //!< C array constructor
  numberCountsKnots( const numberCountsKnots& ); //!< Copy constructor
  ~numberCountsKnots(); //!< Destructor

  /*! \brief Set knot positions, vector version */
  virtual void setKnotPositions(const std::vector<double>&);
  /*! \brief Set knot positions, c array version */
  virtual void setKnotPositions(unsigned int, const double* const);

  virtual void setParams(const paramSet&); //!< Set parameters
  
  unsigned int getNKnots() const { return nknots; } //!< Returns number of knots
  /*! \brief Get position of specified knot */
  double getKnotPos( unsigned int i ) const { return knots[i]; }

  /*! \brief Get knot position and log value */
  std::pair<double,double> getLogKnot(unsigned int i) const {
    return std::pair<double,double>(knots[i],logknotvals[i]); }

  virtual bool isValid() const; //!< Are model parameters valid

  double getMaxFlux() const; //!< Maximum flux supported by model
  double getMinFlux() const; //!< Minimum flux supported by model
  
  bool writeToStream(std::ostream& os) const; //<! Output

  virtual void sendSelf(MPI::Comm&, int dest) const; //!< Send self
  virtual void recieveCopy(MPI::Comm&, int src); //!< Recieve
};

/*! \brief Write to stream operator */
std::ostream& operator<<(std::ostream& os, const numberCountsKnots& b);

/*!
  \brief A class to read in model specifications from init files
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
  bool has_sigma; //!< Has initial value sigma
  double* sigma; //!< Sigma value for initial positions
  bool has_lower_limits; //!< Has some lower limit information
  bool* has_lowlim; //!< Knots have lower limit
  double* lowlim; //!< Value of lower limit
  bool has_upper_limits; //!< Has some upper limit information
  bool* has_uplim; //!< Knots have upper limit
  double* uplim; //!< Value of upper limit

  mutable ran rangen; //!< Random number generator

 public:
  initFileKnots(); //!< Basic constructor
  initFileKnots(const std::string&, bool=false, bool=false); //!< Constructor with file read
  ~initFileKnots(); //!< Destructor

  unsigned int getNKnots() const { return nknots; } //!< Get number of knots
  std::pair<double,double> getKnot(unsigned int idx) const; //!< Get knot pos and value

  void readFile(const std::string&, bool=false, bool=false); //!< Read file

  /*! \brief Set seed of random number generator */
  void setSeed( unsigned long long int seed ) const { rangen.setSeed(seed); }

  void getKnotPos(std::vector<double>&) const; //!< Gets the knot positions
  void getKnotVals(std::vector<double>&) const; //!< Gets the knot values
  void getKnotPos(numberCountsKnots&) const; //!< Sets knot locations in model
  double getKnotPos(unsigned int) const; //!< Get knot position
  double getKnotValue(unsigned int) const; //!< Get knot value

  void getParams(paramSet& p) const; //!< Sets knot values to central values
  void generateRandomKnotValues(paramSet& p) const; //!< Seed knot values

  double getKnotSigma(unsigned int) const; //!< Get knot sigma
  bool isKnotFixed(unsigned int) const; //!< Is a knot fixed?
  
  bool knotHasLowerLimit(unsigned int) const; //!< Does knot have a lower limit
  double getLowerLimit(unsigned int) const; //!< Get knot lower limit

  bool knotHasUpperLimit(unsigned int) const; //!< Does knot have a lower limit
  double getUpperLimit(unsigned int) const; //!< Get knot lower limit

  bool isValid(const paramSet&) const; //!< Checks if parameters are within allowed ranges
  
  void sendSelf(MPI::Comm&, int dest) const; //!< Send self
  void recieveCopy(MPI::Comm&, int src); //!< Recieve
};

#endif
