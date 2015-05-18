#ifndef __global_settings__
#define __global_settings__

#include<cmath>
#include<unistd.h>

/*!
  \brief Global convenience variables
*/
namespace mcmc_affine {
  const char version[] = "0.3.0"; //!< Version number

  // Just in case this cmath doesn't have PI (not all do)
  const double pi = 3.141592653589793238462643383279502884197; //!< \f$\pi\f$

  // Messages we understand
  // These are hardwired to particular values to ease debugging if
  //  we get a message we don't expect.
  // Note that we don't use enum class because MPI excpects integers as
  //  tags.
  /*! \brief MPI message numbers */
  enum affine_messages { ERROR=1, STOP=2, REQUESTPOINT=3, SENDINGPOINT=4,
                         SENDINGRESULT=5, SENDINGPARREJ=6,
                         PSSENDNPARS=100, PSSENDPVALS=101,
                         PSTSENDIDX=110, PSTSENDOLIKE=111, 
                         PSTSENDNLIKE=112, PSTSENDZ=113 };
                         
  /*! \brief Keeps track of special states for parameters 
    FIXED means the parameter is fixed, so all steps should have
          the same value.  Additionally, this parameter should not be used when
          checking convergence.
    ACIGNORE means that the parameter should be generated as usual,
          but should not be included when computing the autocorrelation
          (and hence plays no role in computing the convergence).
    BONUS means that this is not a parameter of the model, but is a
          quantity calculated from the model and stored as a bonus
          piece of information.  It can't be generated as usual, but
          must be supplied by the likelihood model on computation.
          As an example, a number counts model might not formally have the mean
          flux per area as a parameter, but it may be useful to compute
          and store that information in the chain as a bonus parameter.
          Such a parameter will also always be ignored in the autocorrelation.
  */
  // Note that we don't use enum class because we do bitwise operations here
  enum { FIXED=1, ACIGNORE=2, BONUS=4 };
  
  const useconds_t msleeplen = 1000; //!< Master sleep time if no message ready in usec
  const useconds_t ssleeplen = 50; //!< Slave sleep time if no message ready in usec
}
#endif
