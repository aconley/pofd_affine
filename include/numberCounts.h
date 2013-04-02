//numberCounts.h

#ifndef __numberCounts__
#define __numberCounts__

#include<string>
#include<fstream>

#include "../include/beam.h"
#include "../include/paramSet.h"

/*!
  \brief Galaxy counts model abstract base class
  \ingroup Models

  It is up to each subclass to interpret paramSet
*/
//The use of const in getR really means don't change the model parameters,
// not any internal state data.
class numberCounts {
 public:
  /*! \brief For user to request what R they want */
  enum rtype { BEAMPOS=1, BEAMNEG=2, BEAMBOTH=3 }; 
  
  numberCounts() {};
  virtual ~numberCounts() {};

  virtual void setParams(const paramSet& params)=0; //!< Set parameters  
  virtual bool isValid() const = 0; //!< See if model params are valid

  /*! \brief Minimum flux model is defined for */
  virtual double getMinFlux() const = 0;

  /*! \brief Maxium flux model is defined for */
  virtual double getMaxFlux() const = 0;

  /*! \brief Return the number of sources per unit area */
  virtual double getNS() const = 0; 

  /*! \brief Get flux per unit area */
  virtual double getFluxPerArea() const = 0; 
  /*! \brief Get mean flux^2 per unit area) */
  virtual double getFluxSqPerArea() const = 0;

  /*! Get differential number counts */
  virtual double getNumberCounts(double) const = 0; 

  /*! \brief Get number of source responses, single value version */
  virtual double getR(double,const beam&, rtype=BEAMBOTH) const = 0;
  
  /*! \brief Get number of source responses, array version*/
  virtual void getR(unsigned int n,const double* const,
		    const beam&,double*, rtype=BEAMBOTH) const = 0;

  virtual bool writeToStream(std::ostream& os) const=0; //!< Output
};

/*! \brief Write to stream */
std::ostream& operator<<(std::ostream& os, const numberCounts& b);

#endif
