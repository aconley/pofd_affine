#include<cmath>

#include "../include/global_settings.h"
#include "../include/paramSet.h"
#include "../include/affineExcept.h"

paramSet::paramSet() {
  nparams = 0;
  paramvals = NULL;
}

paramSet::paramSet(unsigned int NPARAMS) {
  nparams = 0;
  paramvals = NULL;
  setNParams(NPARAMS);
}

paramSet::paramSet(const std::vector<float>& vec) {
  nparams = 0;
  paramvals = NULL;
  setNParams(vec.size());
  if (nparams > 0) for (unsigned int i = 0; i < nparams; ++i)
		     paramvals[i]=vec[i];
}

paramSet::paramSet(unsigned int N, const float* const VAL) {
  nparams = 0;
  paramvals = NULL;
  setNParams(N);
  if (nparams > 0) 
    for (unsigned int i = 0; i < nparams; ++i)
      paramvals[i]=VAL[i];
}

paramSet::paramSet(const paramSet& other) {
  nparams = 0;
  paramvals = NULL;
  setNParams(other.nparams);
  for (unsigned int i = 0; i < other.nparams; ++i)
    paramvals[i] = other.paramvals[i];
}

paramSet::~paramSet() {
  if (paramvals != NULL) delete[] paramvals;
}

void paramSet::clear() {
  if (paramvals != NULL) delete[] paramvals;
  paramvals = NULL;
  nparams = 0;
}

/*!  Doesn't preserve old params 
 */
void paramSet::setNParams(unsigned int npar) {
  if (npar == nparams) return;
  if (paramvals != NULL) delete[] paramvals;
  if (npar > 0) paramvals = new float[npar]; else paramvals=NULL;
  nparams=npar;
}

paramSet& paramSet::operator=(const paramSet& other) {
  if (this == &other) return *this;
  setNParams(other.nparams);
  for (unsigned int i = 0; i < other.nparams; ++i)
    paramvals[i] = other.paramvals[i];
  return *this;
}

/*!
  Floating point equals, always dangerous.
 */
bool paramSet::operator==(const paramSet& other) const {
  if (nparams != other.nparams) return false;
  if (nparams == 0) return true; //ambiguous...
  bool retval = (paramvals[0] == other.paramvals[0]);
  for (unsigned int i = 1; i < nparams; ++i)
    retval &= (paramvals[i] == other.paramvals[i]);
  return retval;
}

/*!
  \returns Square root of the sum of the differences
 */
float paramSet::getDist(const paramSet& other) const {
  if (other.nparams != nparams)
    throw affineExcept("paramSet","getDist",
		       "Input paramSets don't have the same size",1);
  if (nparams == 0) return 0.0; //ambigouous case
  float val, distsq;
  val = paramvals[0] - other.paramvals[0];
  distsq = val*val;
  for (unsigned int i = 1; i < nparams; ++i) {
    val = paramvals[i] - other.paramvals[i];
    distsq += val*val;
  }
  return sqrt(val);
}

const float& paramSet::at(unsigned int i) const throw(std::range_error) {
  if (i >= nparams) throw std::range_error("paramSet at out of range access");
  return paramvals[i];
}

float& paramSet::at(unsigned int i) throw(std::range_error) {
  if (i >= nparams) throw std::range_error("paramSet at out of range access");
  return paramvals[i];
}

/*!
  \param[in] vec Input parameter vector
  Doesn't allow for resizing, doesn't change noise values
 */
void paramSet::setParamValues(const std::vector<float>& vec) {
  if (vec.size() != nparams)
    throw affineExcept("paramSet","setParamValues",
		     "Input vector wrong length",1);
  for (unsigned int i = 0; i < nparams; ++i)
    paramvals[i]=vec[i];
}

/*!
  Doesn't allow for resizing
 */
void paramSet::setParamValues(unsigned int N, const float* const VAL) {
  if (N != nparams)
    throw affineExcept("paramSet","setParamValues",
		       "Input array wrong length",1);
  for (unsigned int i = 0; i < nparams; ++i)
    paramvals[i]=VAL[i];
}

/*!
  Doesn't allow resizing
*/
void paramSet::readFromStream(std::istream& is) {
  for (unsigned int i = 0; i < nparams; ++i)
    is >> paramvals[i];
}

bool paramSet::writeToStream(std::ostream& os) const {
  for (unsigned int i = 0; i < nparams; ++i)
    os << "   " << paramvals[i];
  return true;
}

void paramSet::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&nparams), 1, MPI_UNSIGNED, 
	   dest, mcmc_affine::PSSENDNPARS, comm);
  MPI_Send(paramvals, nparams, MPI_FLOAT, dest, mcmc_affine::PSSENDPVALS,
	   comm);
}

void paramSet::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  unsigned int newpars;
  MPI_Recv(&newpars, 1, MPI_UNSIGNED, src, mcmc_affine::PSSENDNPARS, 
	   comm, &Info);
  setNParams(newpars);
  MPI_Recv(paramvals, newpars, MPI_FLOAT, src, mcmc_affine::PSSENDPVALS,
	   comm, &Info);
}


std::istream& operator>>(std::istream& is, paramSet& p) {
  p.readFromStream(is);
  return is;
}

std::ostream& operator<<(std::ostream& os, const paramSet& p) {
  p.writeToStream(os);
  return os;
}

