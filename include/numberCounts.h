//numberCounts.h

#ifndef __numberCounts__
#define __numberCounts__

#include<string>
#include<fstream>

#include<hdf5.h>

#include "../include/global_settings.h"
#include "../include/beam.h"
#include "../include/paramSet.h"

/*!
  \brief Number counts model abstract base class, 2D case

  It is up to each subclass to interpret the values in paramSet

  \ingroup Models
*/
//The use of const in getR really means don't change the model parameters,
// not any internal state data.
class numberCounts {
 public:
  numberCounts() {}; //!< Constructor
  virtual ~numberCounts() {}; //!< Destructor

  virtual void setParams(const paramSet& params)=0; //!< Set parameters  
  virtual bool isValid() const noexcept = 0; //!< See if model params are valid

  /*! \brief Distance between two parameter sets over params model cares about*/
  virtual float paramRelativeDistance(const paramSet& p1, const paramSet& p2) 
    const throw(affineExcept) = 0;

  /*! \brief Minimum flux model is defined for */
  virtual double getMinFlux() const noexcept = 0;

  /*! \brief Maxium flux model is defined for */
  virtual double getMaxFlux() const noexcept = 0;

  /*! \brief Get range over which R is expected to be nonzero */
  virtual dblpair getRRange(const beam&) const throw(affineExcept) = 0;

  /*! \brief Return the number of sources per unit area */
  virtual double getNS() const noexcept = 0; 

  /*! \brief Get flux per unit area */
  virtual double getFluxPerArea() const noexcept = 0; 
  /*! \brief Get mean flux^2 per unit area */
  virtual double getFluxSqPerArea() const noexcept = 0;

  /*! Get differential number counts */
  virtual double getNumberCounts(double) const noexcept = 0; 

  /*! \brief Get number of source responses, single value version */
  virtual double getR(double, const beam&) const = 0;

  /*! \brief Tikhonov regularization log Likelihood penalty */
  virtual double differenceRegularize(double) const;

  /*! \brief Get number of source responses, array version*/
  virtual void getR(unsigned int n, const double* const,
                    const beam&, double*) const = 0;

  virtual bool writeToStream(std::ostream& os) const=0; //!< Output to stream
};

/*! \brief Write to stream */
std::ostream& operator<<(std::ostream& os, const numberCounts& b);

#endif
