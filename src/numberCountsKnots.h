//numberCountsKnots.h

#ifndef __numberCountsKnots__
#define __numberCountsKnots__

#include<vector>

#include<numberCounts.h>

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

  virtual void SendSelf(MPI::Comm&, int dest) const; //!< Send self
  virtual void RecieveCopy(MPI::Comm&, int src); //!< Recieve
};

/*! \brief Write to stream operator */
std::ostream& operator<<(std::ostream& os, const numberCountsKnots& b);

#endif
