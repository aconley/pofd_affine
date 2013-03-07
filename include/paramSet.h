//paramset.h

#ifndef __paramSet__
#define __paramSet__

#include<vector>
#include<istream>
#include<ostream>
#include<stdexcept>

#include<mpi.h>

/*!
  \brief Class for holding parameter values of the model
*/
class paramSet {
 private:
  unsigned int nparams; //!< Number of parameters
  double* paramvals;    //!< Array holding parameter values
 public:
  paramSet(); //!< Basic constructor
  paramSet(unsigned int NPARAMS); //!< Constructor with number of parameters
  paramSet(unsigned int, const double* const); //!< Constructor with parameter c array
  paramSet(const std::vector<double>&); //!< Constructor with parameter vector
  paramSet(const paramSet&); //!< Copy constructor
  ~paramSet(); //!< Destructor

  void setNParams(unsigned int); //!< Set number of parameters
  void clear(); //!< Clear parameters, setting number of params to zero

  unsigned int getNParams() const { return nparams; } //!< Get number of parameters
  /*! \brief Return particular parameter value */
  const double& operator[](unsigned int i) const { return paramvals[i]; }
  /*! \brief Return particular parameter value */
  double& operator[](unsigned int i) { return paramvals[i]; }
  const double& at(unsigned int) const throw(std::range_error); //!< Element access with range check
  double& at(unsigned int) throw(std::range_error); //!< Element access with range check

  /*! \brief Set parameter values from vector */
  void setParamValues(const std::vector<double>&);
  /*! \brief Set parameter values from c array */
  void setParamValues(unsigned int, const double* const);
  /*! \brief Set particular parameter value */
  void setParamValue(unsigned int i, double val) { paramvals[i]=val; }
  /*! \brief Copy from other paramSet */
  paramSet& operator=(const paramSet&);
  /*! \brief Are parameter sets equal */
  bool operator==(const paramSet&) const;
  /*! \brief Get distance (euclidean) between two sets of parameters */
  double getDist(const paramSet&) const;

  //Input
  void readFromStream(std::istream& is);  //!< Read parameters from stream

  //Output
  bool writeToStream(std::ostream& os) const; //!< Write parameters to stream
  
  /*! \brief MPI copying, send operation */
  void sendSelf(MPI::Comm& comm, int dest) const;
  /*! \brief MPI copying, recieving operation */
  void recieveCopy(MPI::Comm& comm, int src);

};

/*! \brief Read from stream */
std::istream& operator>>(std::istream& is, paramSet& p);
/*! \brief Write to stream */
std::ostream& operator<<(std::ostream& os, const paramSet& p);

#endif
