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
  static const double ftol;
  unsigned int nknots; //!< Number of knots
  double* knots; //!< Location of knots

  //The knot values are stored as base e 
  double* logknotvals; //!< Log values of differential number counts at knots
  bool knotvals_loaded; //!< Are knot values loaded?

  virtual void setNKnots(unsigned int n); //!< Sets number of knots

 public:
  numberCountsKnots(); //!< Default constructor
  explicit numberCountsKnots(unsigned int);
  numberCountsKnots( const std::vector<double>& );
  numberCountsKnots( unsigned int, const double* const);
  numberCountsKnots( const numberCountsKnots& );
  ~numberCountsKnots(); //!< Destructor

  virtual void setKnotPositions(const std::vector<double>&);
  virtual void setKnotPositions(unsigned int, const double* const);

  virtual void setParams(const paramSet&); //!< Set parameters
  
  unsigned int getNKnots() const { return nknots; } //!< Returns number of knots
  double getKnotPos( unsigned int i ) const { return knots[i]; }

  std::pair<double,double> getLogKnot(unsigned int i) const {
    return std::pair<double,double>(knots[i],logknotvals[i]); }

  virtual bool isValid() const;

  virtual double getNS() const = 0; //!< Return the number of sources with \f$S > S_{min}\f$ per square degree

  virtual double getMeanFluxPerArea() const = 0; //!< Mean flux per unit area (sq deg)
  virtual double getMeanFluxSqPerArea() const = 0; //!< Mean flux squared per unit area (sq deg)

  double getMaxFlux() const;
  double getMinFlux() const;
  
  bool writeToStream(std::ostream& os) const; //<! Output

  virtual void SendSelf(MPI::Comm&, int dest) const; //!< Send self
  virtual void RecieveCopy(MPI::Comm&, int src); //!< Recieve
};

std::ostream& operator<<(std::ostream& os, const numberCountsKnots& b);

#endif
