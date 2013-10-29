#include<cmath>

#include "../include/global_settings.h"
#include "../include/paramSet.h"
#include "../include/affineExcept.h"

paramSet::paramSet() {
  nparams = 0;
  paramvals = NULL;
}

/*!
  \param[in] NPARAMS Number of parameters
*/
paramSet::paramSet(unsigned int NPARAMS) {
  nparams = 0;
  paramvals = NULL;
  setNParams(NPARAMS);
}

/*!
  \param[in] vec Vector of parameters
*/
paramSet::paramSet(const std::vector<float>& vec) {
  nparams = 0;
  paramvals = NULL;
  setNParams(vec.size());
  if (nparams > 0) 
    for (unsigned int i = 0; i < nparams; ++i)
      paramvals[i] = vec[i];
}

/*!
  \param[in] vec Vector of parameters

  Converted to float internally
*/
paramSet::paramSet(const std::vector<double>& vec) {
  nparams = 0;
  paramvals = NULL;
  setNParams(vec.size());
  if (nparams > 0) 
    for (unsigned int i = 0; i < nparams; ++i)
      paramvals[i] = static_cast<float>(vec[i]);
}

/*!
  \param[in] N Number of parameters
  \param[in] VAL Array of paramters
*/
paramSet::paramSet(unsigned int N, const float* const VAL) {
  nparams = 0;
  paramvals = NULL;
  setNParams(N);
  if (nparams > 0) 
    for (unsigned int i = 0; i < nparams; ++i)
      paramvals[i] = VAL[i];
}

/*!
  \param[in] N Number of parameters
  \param[in] VAL Array of paramters, converted to float internally
*/
paramSet::paramSet(unsigned int N, const double* const VAL) {
  nparams = 0;
  paramvals = NULL;
  setNParams(N);
  if (nparams > 0) 
    for (unsigned int i = 0; i < nparams; ++i)
      paramvals[i] = static_cast<float>(VAL[i]);
}

/*!
  \param[in] other Parameter set to copy
*/
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

/*!  
  \param[in] npar Number of new parameters

  Doesn't preserve old params 
 */
void paramSet::setNParams(unsigned int npar) {
  if (npar == nparams) return;
  if (paramvals != NULL) delete[] paramvals;
  if (npar > 0) paramvals = new float[npar]; else paramvals=NULL;
  nparams=npar;
}

/*!
  \param[in] other paramSet to copy

  Will resize as needed
*/
paramSet& paramSet::operator=(const paramSet& other) {
  if (this == &other) return *this;
  setNParams(other.nparams);
  for (unsigned int i = 0; i < other.nparams; ++i)
    paramvals[i] = other.paramvals[i];
  return *this;
}

/*!
  \param[in] other paramSet to check equality against
  \returns True if they have the same number of parameters and same values

  Floating point equals, always dangerous.  Returns true if no
  parameters.
 */
bool paramSet::operator==(const paramSet& other) const {
  if (this == &other) return true;
  if (nparams != other.nparams) return false;
  if (nparams == 0) return true; //ambiguous...
  bool retval = (paramvals[0] == other.paramvals[0]);
  for (unsigned int i = 1; i < nparams; ++i)
    retval &= (paramvals[i] == other.paramvals[i]);
  return retval;
}

/*!
  \param[in] other paramSet to compute distance with respect to
  \returns Square root of the sum of the differences
*/
float paramSet::getDist(const paramSet& other) const {
  if (this == &other) return 0.0;
  if (other.nparams != nparams)
    throw affineExcept("paramSet","getDist",
		       "Input paramSets don't have the same size",1);
  if (nparams == 0) return 0.0; //ambigouous case
  float val, distsq;
  val = paramvals[0] - other.paramvals[0];
  distsq = val * val;
  for (unsigned int i = 1; i < nparams; ++i) {
    val = paramvals[i] - other.paramvals[i];
    distsq += val * val;
  }
  return sqrt(distsq);
}

/*!
  \param[in] i Index to access
  \returns A constant reference to the value at index.

  Throws a range error if the index is invalid
*/
const float& paramSet::at(unsigned int i) const throw(std::range_error) {
  if (i >= nparams) throw std::range_error("paramSet at out of range access");
  return paramvals[i];
}

/*!
  \param[in] i Index to access
  \returns A reference to the value at index.

  Throws a range error if the index is invalid
*/
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
  \param[in] N Number of input values -- must be the same as the current
     number of parameters
  \param[in] VAL Array of values

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
  \param[inout] ifs Input stream

  Doesn't allow resizing
*/
void paramSet::readFromStream(std::istream& ifs) {
  for (unsigned int i = 0; i < nparams; ++i)
    ifs >> paramvals[i];
}

/*!
  \param[inout] os Output stream
  \returns True
*/
bool paramSet::writeToStream(std::ostream& os) const {
  for (unsigned int i = 0; i < nparams; ++i)
    os << "   " << paramvals[i];
  return true;
}

/*!
  \param[in] comm Communicator
  \param[in] dest Destination of messages
*/
void paramSet::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&nparams), 1, MPI_UNSIGNED, 
	   dest, mcmc_affine::PSSENDNPARS, comm);
  MPI_Send(paramvals, nparams, MPI_FLOAT, dest, mcmc_affine::PSSENDPVALS,
	   comm);
}

/*!
  \param[in] comm Communicator
  \param[in] src Source of messages
*/
void paramSet::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  unsigned int newpars;
  MPI_Recv(&newpars, 1, MPI_UNSIGNED, src, mcmc_affine::PSSENDNPARS, 
	   comm, &Info);
  setNParams(newpars);
  MPI_Recv(paramvals, newpars, MPI_FLOAT, src, mcmc_affine::PSSENDPVALS,
	   comm, &Info);
}

/*!
  \param[inout] is Input stream
  \param[inout] p paramSet to load
*/
std::istream& operator>>(std::istream& is, paramSet& p) {
  p.readFromStream(is);
  return is;
}

/*!
  \param[inout] os Output stream
  \param[in] p paramSet to write.
*/
std::ostream& operator<<(std::ostream& os, const paramSet& p) {
  p.writeToStream(os);
  return os;
}

