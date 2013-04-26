#ifndef __global_settings__
#define __global_settings__

#include<cmath>

/*!
  \brief Global convenience variables
*/
namespace mcmc_affine {
  const char version[] = "0.2.0"; //!< Version number

  const double pi = 3.141592653589793238462643383279502884197; //!< \f$\pi\f$

 //Messages we understand
  //These are hardwired to particular values to ease debugging if
  // we get a message we don't expect.
  /*! \brief MPI message numbers */
  enum affine_messages { ERROR=1, STOP=2, REQUESTPOINT=3, SENDINGPOINT=4,
			 SENDINGRESULT=5,
			 PSSENDNPARS=1000, PSSENDPVALS=1002,
			 PSTSENDIDX=1100, PSTSENDOLIKE=1101, 
			 PSTSENDNLIKE=1102, PSTSENDZ=1103};
			 

}
#endif
