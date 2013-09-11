#ifndef __global_settings__
#define __global_settings__

#include<cmath>
#include<unistd.h>

/*!
  \brief Global convenience variables
*/
namespace mcmc_affine {
  const char version[] = "0.2.4"; //!< Version number

  const double pi = 3.141592653589793238462643383279502884197; //!< \f$\pi\f$

 //Messages we understand
  //These are hardwired to particular values to ease debugging if
  // we get a message we don't expect.
  /*! \brief MPI message numbers */
  enum affine_messages { ERROR=1, STOP=2, REQUESTPOINT=3, SENDINGPOINT=4,
			 SENDINGRESULT=5, SENDINGPARREJ=6,
			 PSSENDNPARS=100, PSSENDPVALS=101,
			 PSTSENDIDX=110, PSTSENDOLIKE=111, 
			 PSTSENDNLIKE=112, PSTSENDZ=113};
			 
  /*! \brief Keeps track of special states for parameters */
  enum { FIXED=1, ACIGNORE=2 };
  
  const useconds_t usleeplen = 1; //!< Sleep time if no message ready in usec
}
#endif
