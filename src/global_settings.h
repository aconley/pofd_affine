#ifndef __global_settings__
#define __global_settings__

#include<cmath>

/*!
\mainpage pofd_affine

This is the source documentation for pofd_affine, a package which
implements one and two-dimensional P(D) fitting using the affine
invariant MCMC method described in Foreman-Mackey et al. 2012,
arXiv1202.3665F.

*/

/*!
  \brief Global convenience variables
*/
namespace mcmc_affine {
  const char version[] = "0.1.0"; //!< Version number of MCMC affine library

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

/*!
  \brief P(D) convenience variables
*/
namespace pofd_mcmc {
  const char version[] = "0.1.1"; //!<Version number of P(D) library
  const double n_sigma_shift = 8.0; //!< Shift amount
  const double n_sigma_pad = 10.0; //!< Noise padding size in sigma
  const double n_sigma_shift2d = 4.0; //!< Shift amount
  const double n_sigma_pad2d = 6.0; //!< Noise padding size in sigma, 2D
  const double logfac = log2(10.0); //!< Conversion to base 2 from base 10
  const double ilogfac = 1.0/logfac; //!< Inverse conversion factor
  const double smallval=exp2(-100); //!< Log of small number (base 2)
  const double log2toe = log(2); //!< Multiply by this to go from log2 to ln

  const double smalllogval=-100; //!< Log of small number

  /*! \brief Maximum transform size.  Make sure it's a power of 2*/
  const unsigned int nmax = 16777216;

  //Powers of 2
  const int npow2 = 24; //!< Number of powers of 2 in pow2
  /*! \brief Powers of 2 */
  const int pow2[npow2+1] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 
			     2048, 4096, 8192, 16384, 32768, 65536, 131072, 
			     262144, 524288, 1048576, 2097152, 4194304,
			     8388608, 16777216}; 


  /*! \brief MPI message codes specific to P(D) routines */
  enum pofd_messages { BEAMSENDPIXSIZE=1000, BEAMSENDNPOS=1001,
		       BEAMSENDHASPOSWEIGHTS=1002, BEAMSENDPOSPIXARR=1003,
		       BEAMSENDINVPOSPIXARR=1004, BEAMSENDPOSWEIGHTS=1005,
		       BEAMSENDTOTPOS=1006, BEAMSENDTOTSQPOS=1007,
		       BEAMSENDNNEG=1008,BEAMSENDHASNEGWEIGHTS=1009,
		       BEAMSENDNEGPIXARR=1010,BEAMSENDINVNEGPIXARR=1011,
		       BEAMSENDNEGWEIGHTS=1012,BEAMSENDTOTNEG=1013,
		       BEAMSENDTOTSQNEG=1014,
		       DOUBLEBEAMSENDPIXSIZE=1500,
		       DOUBLEBEAMSENDN=1501, DOUBLEBEAMSENDHASWEIGHTS=1502,
		       DOUBLEBEAMSENDPIXARR1=1503,DOUBLEBEAMSENDPIXARR2=1504, 
		       DOUBLEBEAMSENDIPIXARR1=1505,DOUBLEBEAMSENDIPIXARR2=1506,
		       DOUBLEBEAMSENDWEIGHTS=1507,
		       DOUBLEBEAMSENDTOT1=1508,DOUBLEBEAMSENDTOT2=1509,
		       DOUBLEBEAMSENDTOTSQ1=1510,DOUBLEBEAMSENDTOTSQ2=1511,
		       DOUBLEBEAMSENDTOTSM1=1512,DOUBLEBEAMSENDTOTSM2=1513,
		       NCKSENDNKNOTS=2000, NCKSENDKNOTS=2001,
		       NCKSENDKNOTSLOADED=2002, NCKSENDLOGKNOTVALS=2003,
		       NCKSSENDLOGKNOTS=2100,
		       NCDCSENDNKNOTS=2500, NCDCSENDKNOTS=2501,
		       NCDCSENDLOGKNOTS=2502,NCDCSENDKNOTSLOADED=2503,
		       NCDCSENDLOGKNOTVALS=2504,NCDCSENDNSIGMAKNOTS=2505,
		       NCDCSENDSIGMAKNOTS=2506, NCDCSENDSIGMAVALS=2507, 
		       NCDCSENDNOFFSETKNOTS=2508, NCDCSENDOFFSETKNOTS=2509,
		       NCDCSENDOFFSETVALS=2510,
		       PDFSENDPLANSTYLE=3000, PDFHASWISDOM=3001,
		       PDFWISLEN=3002, PDFWISNAME=3003,
		       PDFVERBOSE=3004, PDFNINTERP=3004,
		       PDFDSENDPLANSTYLE=3100, PDFDHASWISDOM=3101,
		       PDFDWISLEN=3102, PDFDWISNAME=3103,
		       PDFDVERBOSE=3104, PDFDNEDGE=3104,
		       FDSENDN=4000, FDSENDDATA=4001, FDSENDISBINNED=4002,
		       FDSENDNBINS=4003, FDSENDBINVAL=4004, 
		       FDSENDBINCENT0=4005, FDSENDBINDELTA=4006,
		       FDDSENDN=4100, FDDSENDDATA1=4101, FDDSENDDATA2=4102,
		       FDDSENDISBINNED=4103, FDDSENDNBINS1=4104,
		       FDDSENDNBINS2=4105, FDDSENDBINVAL=4106,
		       FDDSENDBINCENT01=4107, FDDSENDBINDELTA1=4108,
		       FDDSENDBINCENT02=4109, FDDSENDBINDELTA2=4110,
		       CLSENDFFTSIZE=5000, CLSENDEDGEFIX=5001,
		       CLSENDNDATA=5002, CLSENDDATAREAD=5003,
		       CLSENDSIGMABASE=5004, CLSENDLIKEOFFSET=5005,
		       CLSENDLIKENORM=5006, CLSENDMAXFLUX=5007,
		       CLSENDMAXSIGMABASE=5008, CLSENDHASBEAM=5009,
		       CLSENDHASCFIRBPRIOR=5010, CLSENDCFIRBPRIORMEAN=5011,
		       CLSENDCFIRBPRIORSIGMA=5012, CLSENDHASSIGMAPRIOR=5013,
		       CLSENDSIGMAPRIORWIDTH=5014, CLSENDNINTERP=5015,
		       CLSENDNBEAM=5016, CLSENDSETNUM=5017, 
		       CLSENDBINDATA=5018, CLSENDNBINS=5019,
		       CLDSENDNDATA=5100, CLDSENDDATAREAD=5101,
		       CLDSENDSIGMABASE1=5102, CLDSENDSIGMABASE2=5103,
		       CLDSENDMAXSIGMABASE1=5104, CLDSENDMAXSIGMABASE2=5105,
		       CLDSENDLIKEOFFSET=5106, CLDSENDLIKENORM=5107,
		       CLDSENDMAXFLUX1=5108, CLDSENDMAXFLUX2=5109,
		       CLDSENDHASBEAM=5110, CLDSENDFFTSIZE=5120,
		       CLDSENDEDGEFIX=5121, CLDSENDEDGEINTEG=5122,
		       CLDSENDNEDGE=5123, CLDSENDNBEAM=5124,
		       CLDSENDSETNUM=5125, CLDSENDBINDATA=5126,
		       CLDSENDNBINS=5127, CLDSENDHASCFIRBPRIOR1=5128,
		       CLDSENDCFIRBPRIORMEAN1=5129, 
		       CLDSENDCFRIBPRIORSIGMA1=5130, CLDSENDHASCFIRBPRIOR2=5131,
		       CLDSENDCFIRBPRIORMEAN2=5132, 
		       CLDSENDCFRIBPRIORSIGMA2=5133, CLDSENDHASSIGMAPRIOR1=5134,
		       CLDSENDSIGMAPRIORWIDTH1=5135, CLDSENDHASSIGMAPRIOR2=5136,
		       CLDSENDSIGMAPRIORWIDTH2=5137,
		       IFKSENDNKNOTS=6000, IFKSENDKNOTPOS=6001,
		       IFKSENDKNOTVAL=6002, IFKHASSIGMA=6003,
		       IFKSENDSIGMA=6004, IFKHASLOWERLIMITS=6005,
		       IFKSENDHASLOWLIM=6006, IFKSENDLOWLIM=6007,
		       IFKHASUPPERLIMITS=6008, IFKSENDHASUPLIM=6009,
		       IFKSENDUPLIM=6010, IFDLNSENDNKNOTS=6500, 
		       IFDLNSENDNSIGMAS=6501, IFDLNSENDNOFFSETS=6502,
		       IFDLNSENDKNOTPOS=6503, IFDLNSENDKNOTVAL=6504, 
		       IFDLNHASSIGMA=6505, IFDLNSENDSIGMA=6506, 
		       IFDLNHASLOWERLIMITS=6507, IFDLNSENDHASLOWLIM=6508, 
		       IFDLNSENDLOWLIM=6509, IFDLNHASUPPERLIMITS=6510, 
		       IFDLNSENDHASUPLIM=6511, IFDLNSENDUPLIM=6512,
		       PMCMCSENDINIT=10000, PMCMCSENDINGINIT=10001,
		       PMCMCSENDNPAR=10002, PMCMCISREADY=10003};

}

#endif
