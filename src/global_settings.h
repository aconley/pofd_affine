#ifndef __global_settings__
#define __global_settings__

#include<cmath>

/*!
  \brief Global convenience variables
*/
namespace mcmc_affine {
  const char version[] = "0.1.0"; //Version number of MCMC affine library

  const double pi = 3.141592653589793238462643383279502884197; //!< \f$\pi\f$
  const double two_pi = 2.0*pi; //!< \f$2 \pi \f$

 //Messages we understand
  //These are hardwired to particular values to ease debugging if
  // we get a message we don't expect.
  /*! \brief MPI message numbers */
  enum affine_messages { ERROR=1, STOP=2, REQUESTPOINT=3, SENDINGPOINT=4,
			 SENDINGRESULT=5,
			 PSSENDNPARS=1000, PSSENDPVALS=1002,
			 PSTSENDIDX=1100, PSTSENDOLIKE=1101, 
			 PSTSENDNLIKE=1102 };
}

namespace pofd_mcmc {
  const char version[] = "0.1.0"; //Version number of P(D) library
  const double n_sigma_pad = 12.0; //!< Noise padding size in sigma
  const double n_sigma_pad2d = 5.0; //!< Noise padding size in sigma, 2D
  const double logfac = log2(10.0); //!< Conversion base 2, base 10
  const double ilogfac = 1.0/logfac; //!< Inverse conversion factor
  const double smalllogval=-100; //!< Log of small number

  /*! \brief Maximum transform size.  Make sure it's a power of 2*/
  const unsigned int nmax = 16777216;

  //Powers of 2
  const int npow2 = 24; //!< Number of powers of 2 in pow2
  const int pow2[npow2+1] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 
			     2048, 4096, 8192, 16384, 32768, 65536, 131072, 
			     262144, 524288, 1048576, 2097152, 4194304,
			     8388608, 16777216}; //!< Powers of 2


  //Messages specific to P(D) stuff
  enum pofd_messages { BEAMSENDPIXSIZE=1000, BEAMSENDNPOS=1001,
		       BEAMSENDHASPOSWEIGHTS=1002, BEAMSENDPOSPIXARR=1003,
		       BEAMSENDINVPOSPIXARR=1004, BEAMSENDPOSWEIGHTS=1005,
		       BEAMSENDTOTPOS=1006, BEAMSENDTOTSQPOS=1007,
		       BEAMSENDNNEG=1008,BEAMSENDHASNEGWEIGHTS=1009,
		       BEAMSENDNEGPIXARR=1010,BEAMSENDINVNEGPIXARR=1011,
		       BEAMSENDNEGWEIGHTS=1012,BEAMSENDTOTNEG=1013,
		       BEAMSENDTOTSQNEG=1014,
		       NCKSENDNKNOTS=2000, NCKSENDKNOTS=2001,
		       NCKSENDKNOTSLOADED=2002, NCKSENDLOGKNOTVALS=2003,
		       NCKSSENDLOGKNOTS=2100};

}

#endif
