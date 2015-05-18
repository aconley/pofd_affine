//paramset.h

#ifndef __paramSet__
#define __paramSet__

#include<vector>
#include<istream>
#include<ostream>
#include<stdexcept>

#include<mpi.h>
#include "../include/affineExcept.h"

/*!
  \brief Class for holding parameter values of the model
*/
class paramSet final {
 private:
  unsigned int nparams; //!< Number of parameters
  float* paramvals;    //!< Array holding parameter values

 public:
  paramSet(); //!< Basic constructor
  /*! \brief  Constructor with number of parameters */
  explicit paramSet(unsigned int NPARAMS);
  paramSet(const std::vector<float>&); //!< Constructor with parameter vector
  paramSet(const std::vector<double>&); //!< Constructor with parameter vector
  paramSet(unsigned int, const float* const); //!< Constructor with parameter c array
  paramSet(unsigned int, const double* const); //!< Constructor with parameter c array
  paramSet(const paramSet&); //!< Copy constructor
  paramSet(paramSet&&); //!< Move constructor
  ~paramSet(); //!< Destructor

  void setNParams(unsigned int); //!< Set number of parameters
  void clear(); //!< Clear parameters, setting number of params to zero

  unsigned int getNParams() const { return nparams; } //!< Get number of parameters
  /*! \brief Return particular parameter value */
  const float& operator[](unsigned int i) const { return paramvals[i]; }
  /*! \brief Return particular parameter value */
  float& operator[](unsigned int i) { return paramvals[i]; }
  const float& at(unsigned int) const throw(std::range_error); //!< Element access with range check
  float& at(unsigned int) throw(std::range_error); //!< Element access with range check

  /*! \brief Set parameter values from vector */
  void setParamValues(const std::vector<float>&) throw(affineExcept);
  /*! \brief Set parameter values from c array */
  void setParamValues(unsigned int, const float* const) throw(affineExcept);
  /*! \brief Set particular parameter value */
  void setParamValue(unsigned int i, float val) { paramvals[i]=val; }
  /*! \brief Copy from other paramSet */
  paramSet& operator=(const paramSet&);
  /*! \brief Copy with move semantics */
  paramSet& operator=(paramSet&&);
  /*! \brief Are parameter sets equal */
  bool operator==(const paramSet&) const;
  /*! \brief Get distance (euclidean) between two sets of parameters */
  float getDist(const paramSet&) const throw(affineExcept);

  //Input
  void readFromStream(std::istream& is);  //!< Read parameters from stream

  //Output
  bool writeToStream(std::ostream& os) const; //!< Write parameters to stream
  
  /*! \brief MPI copying, send operation */
  void sendSelf(MPI_Comm comm, int dest) const;
  /*! \brief MPI copying, recieving operation */
  void receiveCopy(MPI_Comm comm, int src);

};

/*! \brief Read from stream */
std::istream& operator>>(std::istream& is, paramSet& p);
/*! \brief Write to stream */
std::ostream& operator<<(std::ostream& os, const paramSet& p);

#endif
