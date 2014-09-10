//numberCountsDouble.h

#ifndef __numberCountsDouble__
#define __numberCountsDouble__

#include<string>
#include<fstream>

#include<hdf5.h>

#include "../include/global_settings.h"
#include "../include/doublebeam.h"
#include "../include/paramSet.h"

/*!
  \brief Number counts model abstract base class, 2D case

  It is up to each subclass to interpret the values in paramSet

  \ingroup Models
*/
//The use of const in getR really means don't change the model parameters,
// not any internal state data.
class numberCountsDouble {
 public:
  
  numberCountsDouble() {};
  virtual ~numberCountsDouble() {};
  
  virtual bool isValid() const = 0; //!< See if model params are valid

  /*! \brief Relative distance between two sets of params over model params*/
  virtual float paramRelativeDistance(const paramSet& p1, const paramSet& p2)
    const throw(affineExcept) = 0;

  /*! \brief Get Mean Flux per unit area */
  virtual double getFluxPerArea(unsigned int) const = 0;

  /*! \brief Get Mean Flux^2 per unit area */
  virtual double getFluxSqPerArea(unsigned int) const = 0;

  /*! \brief Minimum flux model is defined for */
  virtual dblpair getMinFlux() const = 0;

  /*! \brief Maxium flux model is defined for */
  virtual dblpair getMaxFlux() const = 0;

  /*! \brief Get range over which R is expected to be nonzero */
  virtual std::pair<dblpair, dblpair> 
    getRRange(const doublebeam&) const throw(affineExcept) = 0;

  virtual void setParams(const paramSet& params)=0; //!< Set parameters

  virtual unsigned int getNParams() const = 0; //!< Return number of parameters required

  /*! Evaluates number counts model */
  virtual double getNumberCounts(double, double) const = 0; 

  /*! Evaluates band 1 number counts model at specified value of S_1 */
  virtual double getBand1NumberCounts(double) const = 0;
  /*! Evaluates band 2 number counts model at specified value of S_2 */
  virtual double getBand2NumberCounts(double) const = 0;
  
  /*! \brief Get number of source responses, single value version */
  virtual double getR(double, double, const doublebeam&) const = 0;
  
  /*! \brief Get number of source responses, array version*/
  virtual void getR(unsigned int,const double* const,
		    unsigned int,const double* const,
		    const doublebeam&, double*) const = 0;

  virtual bool writeToStream(std::ostream& os) const=0; //!< Output to stream
};

/*! \brief Write to stream */
std::ostream& operator<<(std::ostream& os, const numberCountsDouble& b);

#endif
