//numberCountsKnots.cc
#include<iostream>
#include<cmath>
#include<iomanip>
#include<limits>
#include<cstdlib>

#include<global_settings.h>
#include<numberCountsKnots.h>

#include<utility.h>
#include<affineExcept.h>

numberCountsKnots::numberCountsKnots() {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  knotvals_loaded = false;
}

numberCountsKnots::numberCountsKnots( unsigned int NKNOTS ) :
  nknots(NKNOTS) {
  if (nknots > 0) {
    knots = new double[nknots];
    logknotvals = new double[nknots];
  } else
    knots = logknotvals = NULL;
  knotvals_loaded = false;
}  

numberCountsKnots::numberCountsKnots( const std::vector<double>& S ) {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setKnotPositions(S);
  knotvals_loaded = false;
}

numberCountsKnots::numberCountsKnots( unsigned int n, const double* const S) {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setKnotPositions(n,S);
  knotvals_loaded = false;
}

numberCountsKnots::numberCountsKnots( const numberCountsKnots& other ) {
  if ( this == &other ) return; //Self-copy
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setNKnots(other.nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = other.knots[i];
  if (other.knotvals_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      knots[i] = other.logknotvals[i];
  }
  knotvals_loaded = other.knotvals_loaded;
}

numberCountsKnots::~numberCountsKnots() {
  if (knots != NULL) delete[] knots;
  if (logknotvals != NULL) delete[] logknotvals;
}

void numberCountsKnots::setNKnots(unsigned int n) {
  if ( nknots == n ) return;
  if ( knots != NULL ) delete[] knots;
  if ( logknotvals != NULL ) delete[] logknotvals;
  if ( n > 0 ) {
    knots = new double[n];
    logknotvals = new double[n];
  } else {
    knots = logknotvals = NULL;
  }
  nknots = n;
  knotvals_loaded = false;
}

/*!
  \param[in] S Input knot positions
*/
void numberCountsKnots::setKnotPositions(const std::vector<double>& S) {
  unsigned int n = S.size();
  if (n != nknots) setNKnots(n);
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsKnots","setKnots",
			 "Negative knot positions not allowed",1);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
}

/*!
  \param[in] n Number of knots
  \param[in] S Input knot positions
*/
void numberCountsKnots::setKnotPositions(unsigned int n, 
					 const double* const S) {
  if (n != nknots) setNKnots(n);
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsKnots","setKnots",
			 "Negative knot positions not allowed",1);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
}

/*!
  \param[in] F Parameters to set in model
 */
void numberCountsKnots::setParams(const paramSet& F) {
  if (F.getNParams() <= 2)
    throw affineExcept("numberCountsKnots","setKnots",
		       "No knots present",1);
  if (nknots != (F.getNParams())) 
    throw affineExcept("numberCountsKnots","setKnots",
		       "Number of knot values different than expected",2);
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = pofd_mcmc::logfac*F[i]; //convert to base 2
  knotvals_loaded = true;
}

/*!
  \returns True if the model parameters are valid
 */
bool numberCountsKnots::isValid() const {
  if (nknots == 0) return false;
  if (!knotvals_loaded) return false;
  for (unsigned int i = 0; i < nknots; ++i)
    if ( std::isnan(knots[i]) ) return false;
  if ( knots[0] <= 0.0 ) return false;
  for (unsigned int i = 1; i < nknots; ++i)
    if (knots[i] <= knots[i-1] ) return false;
  for (unsigned int i = 0; i < nknots; ++i)
    if ( std::isnan(logknotvals[i]) ) return false;
  return true;
}

double numberCountsKnots::getMinFlux() const {
  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  return knots[0];
}

double numberCountsKnots::getMaxFlux() const {
  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  return knots[nknots-1];
}

void numberCountsKnots::SendSelf(MPI::Comm& comm, int dest) const {
  comm.Send(&nknots,1,MPI::UNSIGNED,dest,pofd_mcmc::NCKSENDNKNOTS);
  if (nknots != 0) {
    comm.Send(knots,nknots,MPI::DOUBLE,dest,pofd_mcmc::NCKSENDKNOTS);
    comm.Send(&knotvals_loaded,1,MPI::BOOL,dest,pofd_mcmc::NCKSENDKNOTSLOADED);
    if (knotvals_loaded)
      comm.Send(logknotvals,nknots,MPI::DOUBLE,dest,
		pofd_mcmc::NCKSENDLOGKNOTVALS);
  }
}

void numberCountsKnots::RecieveCopy(MPI::Comm& comm, int src) {
  unsigned int n;
  comm.Recv(&n,1,MPI::UNSIGNED,src,pofd_mcmc::NCKSENDNKNOTS);
  if (n != 0) {
    if (n != nknots) setNKnots(n);
    comm.Recv(knots,nknots,MPI::DOUBLE,src,pofd_mcmc::NCKSENDKNOTS);
    comm.Recv(&knotvals_loaded,1,MPI::BOOL,src,pofd_mcmc::NCKSENDKNOTSLOADED);
    if (knotvals_loaded)
      comm.Recv(logknotvals,nknots,MPI::DOUBLE,src,
		pofd_mcmc::NCKSENDLOGKNOTVALS);
  }
}

bool numberCountsKnots::writeToStream(std::ostream& os) const {
  os << "Model parameters: " << std::endl;
  if (knotvals_loaded) {
    os << " " << std::left << std::setw(13) << "#Flux knot" << "  "
       << std::setw(13) << "Knot value" << std::endl;
    //Convert to log10 for output
    for (unsigned int i = 0; i < nknots; ++i)
      os << " " << std::left << std::setw(13) << knots[i] << "  "
	 << std::setw(13) << pofd_mcmc::ilogfac * logknotvals[i] << std::endl; 
  } else
    os << "Number of knots: " << nknots << std::endl;
  return true;
}

std::ostream& operator<<(std::ostream& os, const numberCountsKnots& b) {
  b.writeToStream(os);
  return os;
}
