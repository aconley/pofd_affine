//numberCounts.h

#ifndef __numberCounts__
#define __numberCounts__

#include<string>
#include<fstream>

#include<beam.h>
#include<paramSet.h>

/*!
  \brief Galaxy counts model abstract base class
  \ingroup Models

  It is up to each subclass to interpret paramSet
*/
//The use of const in getR really means don't change the model parameters,
// not any internal state data.
class numberCounts {
 public:
  enum rtype { BEAMPOS=1, BEAMNEG=2, BEAMBOTH=3 }; //!< For user to request what R they wayt
  
  numberCounts() {};
  virtual ~numberCounts() {};
  
  virtual bool isValid() const = 0; //!< See if model params are valid

  /*! \brief Get Mean Flux per unit area */
  virtual double getMeanFluxPerArea() const = 0;

  /*! \brief Minimum flux model is defined for */
  virtual double getMinFlux() const = 0;

  /*! \brief Maxium flux model is defined for */
  virtual double getMaxFlux() const = 0;

  virtual void setParams(const paramSet& params)=0; //!< Set parameters

  /*! Evaluates number counts model */
  virtual double getNumberCounts( double ) const = 0; 

  /*! \brief Get number of source responses, single value version */
  virtual double getR(double,const beam&, rtype=BEAMBOTH) const = 0;
  
  /*! \brief Get number of source responses, array version*/
  virtual void getR(unsigned int n,const double* const,
		    const beam&,double*, rtype=BEAMBOTH) const = 0;

  virtual bool writeToStream(std::ostream& os) const=0; //<! Output
};

std::ostream& operator<<(std::ostream& os, const numberCounts& b);

#endif
