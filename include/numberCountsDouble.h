//numberCountsDouble.h

#ifndef __numberCountsDouble__
#define __numberCountsDouble__

#include<string>
#include<fstream>

#include "hdf5.h"

#include "../include/doublebeam.h"
#include "../include/paramSet.h"

/*!
  \brief Galaxy counts model abstract base class, 2D case

  It is up to each subclass to interpret paramSet

  \ingroup Models
*/
//The use of const in getR really means don't change the model parameters,
// not any internal state data.
class numberCountsDouble {
 public:
  /*! \brief For user to request what R they want */
  enum rtype { BEAMPOS=1, BEAMNEG=2, BEAMPOSNEG=3, BEAMNEGPOS=4, 
	       BEAMALL=5 }; 
  
  numberCountsDouble() {};
  virtual ~numberCountsDouble() {};
  
  virtual bool isValid() const = 0; //!< See if model params are valid

  /*! \brief Get Mean Flux per unit area */
  virtual double getFluxPerArea(unsigned int) const = 0;

  /*! \brief Get Mean Flux^2 per unit area */
  virtual double getFluxSqPerArea(unsigned int) const = 0;

  /*! \brief Minimum flux model is defined for */
  virtual double getMinFlux(unsigned int) const = 0;

  /*! \brief Maxium flux model is defined for */
  virtual double getMaxFlux(unsigned int) const = 0;

  virtual void setParams(const paramSet& params)=0; //!< Set parameters

  virtual unsigned int getNParams() const = 0; //!< Return number of parameters required

  /*! Evaluates number counts model */
  virtual double getNumberCounts( double, double ) const = 0; 

  /*! \brief Get number of source responses, single value version */
  virtual double getR(double, double, const doublebeam&, 
		      rtype=BEAMALL) const = 0;
  
  /*! \brief Get number of source responses, array version*/
  virtual void getR(unsigned int,const double* const,
		    unsigned int,const double* const,
		    const doublebeam&, double*, rtype=BEAMALL) const = 0;

  virtual void writeToHDF5Handle(hid_t objid) const=0; //!< Output to HDF5
  virtual bool writeToStream(std::ostream& os) const=0; //!< Output to stream
};

/*! \brief Write to stream */
std::ostream& operator<<(std::ostream& os, const numberCountsDouble& b);

#endif
