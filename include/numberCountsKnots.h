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
  virtual ~numberCountsKnots(); //!< Destructor

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
  virtual void setParams(const paramSet&) override; //!< Set parameters

  unsigned int getNKnots() const noexcept { return nknots; } //!< Returns number of knots
  /*! \brief Get position of specified knot */
  double getKnotPos(unsigned int i) const noexcept { return knots[i]; }

  /*! \brief Get knot position and log value */
  dblpair getLogKnot(unsigned int i) const noexcept {
    return std::pair<double,double>(knots[i], logknotvals[i]); }

  virtual bool isValid() const noexcept override; //!< Are model parameters valid

  /*! \brief Distance between two parameter sets over params model cares about*/
  virtual float paramRelativeDistance(const paramSet& p1, const paramSet& p2) 
    const throw(affineExcept) override;

  double getMaxFlux() const noexcept override; //!< Maximum flux supported by model
  double getMinFlux() const noexcept override; //!< Minimum flux supported by model

  /*! \brief Get range over which R is expected to be nonzero */
  dblpair getRRange(const beam&) const throw(affineExcept) override;

  virtual void writeToHDF5Handle(hid_t objid, bool=false) const; //!< Output to HDF5
  virtual bool writeToStream(std::ostream& os) const; //<! Output to stream

  virtual void readFromHDF5Handle(hid_t); //!< Read from HDF5
  
  virtual void sendSelf(MPI_Comm, int dest) const; //!< Send self
  virtual void receiveCopy(MPI_Comm, int src); //!< Receive
};

/*! \brief Write to stream operator */
std::ostream& operator<<(std::ostream& os, const numberCountsKnots& b);

#endif
