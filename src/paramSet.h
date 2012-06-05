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
  unsigned int nparams;
  double* paramvals;
 public:
  paramSet();
  paramSet(unsigned int NPARAMS);
  paramSet(unsigned int, const double* const);
  paramSet(const std::vector<double>&);
  paramSet(const paramSet&);
  ~paramSet();

  void setNParams(unsigned int);
  void clear();

  unsigned int getNParams() const { return nparams; }
  const double& operator[](unsigned int i) const { return paramvals[i]; }
  double& operator[](unsigned int i) { return paramvals[i]; }
  const double& at(unsigned int) const throw(std::range_error); //!< Element access with range check
  double& at(unsigned int) throw(std::range_error); //!< Element access with range check

  void setParamValues(const std::vector<double>&);
  void setParamValues(unsigned int, const double* const);
  void setParamValue(unsigned int i, double val) { paramvals[i]=val; }
  paramSet& operator=(const paramSet&);
  bool operator==(const paramSet&) const;
  double getDist(const paramSet&) const;

  //Input
  void readFromStream(std::istream& is);

  //Output
  bool writeToStream(std::ostream& os) const;
  
  void sendSelf(MPI::Comm& comm, int dest) const;
  void recieveCopy(MPI::Comm& comm, int src);

};

std::istream& operator>>(std::istream& is, paramSet& p);
std::ostream& operator<<(std::ostream& os, const paramSet& p);

#endif
