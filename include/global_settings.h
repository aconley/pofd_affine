#ifndef __global_settings__
#define __global_settings__

#include<cmath>
#include<utility>
#include<unistd.h>

/*!
\mainpage pofd_affine

This is the source documentation for pofd_affine, a package which
implements one and two-dimensional P(D) fitting using the affine
invariant MCMC method described in Foreman-Mackey et al. 2013,
PASP 125, 925 (arXiv1202.3665F).

*/

/*!
  \brief Global convenience variables
*/
namespace mcmc_affine {
  const char version[] = "0.2.9"; //!< Version number

  // Just in case this cmath doesn't have PI (not all do)
  const double pi = 3.141592653589793238462643383279502884197; //!< \f$\pi\f$
  const double two_pi = 2.0*pi; //!< \f$2 \pi \f$

 //Messages we understand
  //These are hardwired to particular values to ease debugging if
  // we get a message we don't expect.  These shouldn't overlap with
  // pofd_messages, so are all < 1000
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
  enum { FIXED=1, ACIGNORE=2, BONUS=4 };
  
  const useconds_t msleeplen = 1000; //!< Master sleep time if no message ready in usec
  const useconds_t ssleeplen = 50; //!< Slave sleep time if no message ready in usec
}

/*!
  \brief P(D) convenience variables
*/
namespace pofd_mcmc {
  const char version[] = "0.4.1"; //!<Version number of P(D) library
  const double n_zero_pad = 7.5; //!< Zero padding size in sigma
  const double logfac = log2(10.0); //!< Conversion to base 2 from base 10
  const double ilogfac = 1.0 / logfac; //!< Inverse conversion factor
  const double smalllogval = -30.; //!< Log2 of small number
  const double smallval = exp2(smalllogval); //!< Log of small number (base 2)
  const double log2toe = log(2.); //!< Multiply by this to go from log2 to ln
  const double ilog2toe = 1.0 / log(2.); //!< exponent conversion factor

  // All 1000 or more to avoid overlap with affine_messages
  // Why do this?  Purely for debugging purposes, so if there is a problem
  // I can trace back which messages are flying around.
  /*! \brief MPI message codes specific to P(D) routines */
  enum pofd_messages { BEAMSENDPIXSIZE=1000, BEAMSENDMINVAL=1001,
		       BEAMSENDNBINS=1002,
		       BEAMSENDNPOS=1003, BEAMSENDPOSPIXARR=1004,
		       BEAMSENDINVPOSPIXARR=1005, BEAMSENDISPOSHIST=1006,
		       BEAMSENDPOSNBINS=1007, BEAMSENDPOSWEIGHTS=1008, 
		       BEAMSENDPOSHISTVAL=1009,
		       BEAMSENDNNEG=1010, BEAMSENDNEGPIXARR=1011,
		       BEAMSENDINVNEGPIXARR=1012, BEAMSENDISNEGHIST=1013,
		       BEAMSENDNEGNBINS=1014, BEAMSENDNEGWEIGHTS=1015, 
		       BEAMSENDNEGHISTVAL=1016,
		       BEAMSENDTOTPOS=1017, BEAMSENDTOTSQPOS=1018,
		       BEAMSENDTOTNEG=1019, BEAMSENDTOTSQNEG=1020,
		       DOUBLEBEAMSENDPIXSIZE=1500, DOUBLEBEAMSENDMINVAL=1501,
		       DOUBLEBEAMSENDHASSIGN=1502, DOUBLEBEAMSENDMIN1=1503,
		       DOUBLEBEAMSENDMAX1=1504, DOUBLEBEAMSENDMIN2=1505,
		       DOUBLEBEAMSENDMAX2=1506, DOUBLEBEAMSENDNPIX=1507,
		       DOUBLEBEAMSENDPIXARR1=1508,DOUBLEBEAMSENDPIXARR2=1509,
		       DOUBLEBEAMSENDIPIXARR1=1510,DOUBLEBEAMSENDIPIXARR2=1511,
		       DOUBLEBEAMSENDNBINS=1512, 
		       DOUBLEBEAMSENDISHISTOGRAMMED=1513,
		       DOUBLEBEAMSENDNHIST=1514, DOUBLEBEAMSENDBINWEIGHTS=1515,
		       DOUBLEBEAMSENDBINVALS1=1516, DOUBLEBEAMSENDBINVALS2=1517,
		       DOUBLEBEAMSENDTOT1=1518,DOUBLEBEAMSENDTOT2=1519,
		       DOUBLEBEAMSENDTOTSQ1=1520,DOUBLEBEAMSENDTOTSQ2=1521,
		       DOUBLEBEAMSENDHASLOGRATIO=1522,
		       DOUBLEBEAMSENDLOGRATIO=1523,
		       DOUBLEBEAMSENDHASBINLOGRATIO=1524,
		       DOUBLEBEAMSENDBINLOGRATIO=1524,
		       NCKSENDNKNOTS=2000, NCKSENDKNOTS=2001,
		       NCKSENDKNOTSLOADED=2002, NCKSENDLOGKNOTVALS=2003,
		       NCKSSENDLOGKNOTS=2100, NCDCSENDKPLOADED=2500,
		       NCDCSENDNKNOTS=2501, NCDCSENDKNOTS=2502,
		       NCDCSENDLOGKNOTS=2503, NCDCSENDKVLOADED=2504,
		       NCDCSENDLOGKNOTVALS=2505, NCDCSENDSPLOADED=2506,
		       NCDCSENDNSIGMAKNOTS=2507, NCDCSENDSIGMAKNOTS=2508, 
		       NCDCSENDSVLOADED=2509, NCDCSENDSIGMAVALS=2510, 
		       NCDCSENDOPLOADED=2511, NCDCSENDNOFFSETKNOTS=2512,
		       NCDCSENDOFFSETKNOTS=2513, NCDCSENDOVLOADED=2514,
		       NCDCSENDOFFSETVALS=2515,
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
		       CLSENDLIKENORM=5006, CLSENDMINDATAFLUX=5007,
		       CLSENDMAXDATAFLUX=5008, CLSENDMINRFLUX=5009,
		       CLSENDMAXRFLUX=5010, CLSENDMAXSIGMABASE=5011, 
		       CLSENDEXPCONF=5012, CLSENDHASBEAM=5013, 
		       CLSENDHASCFIRBPRIOR=5014, CLSENDCFIRBPRIORMEAN=5015, 
		       CLSENDCFIRBPRIORSIGMA=5016, 
		       CLSENDHASPOISSONPRIOR=5017, CLSENDPOISSONPRIORMEAN=5018,
		       CLSENDPOISSONPRIORSIGMA=5019, CLSENDHASSIGMAPRIOR=5020, 
		       CLSENDSIGMAPRIORWIDTH=5021, CLSENDNINTERP=5022, 
		       CLSENDNBEAM=5023, CLSENDSETNUM=5024, 
		       CLSENDBINDATA=5025, CLSENDNBINS=5026,
		       CLSENDREGULARIZATIONALPHA=5027,
		       CLDSENDNDATA=5100, CLDSENDDATAREAD=5101,
		       CLDSENDSIGMABASE1=5102, CLDSENDSIGMABASE2=5103,
		       CLDSENDMAXSIGMABASE1=5104, CLDSENDMAXSIGMABASE2=5105,
		       CLDSENDEXPCONF1=5106, CLDSENDEXPCONF2=5107,
		       CLDSENDLIKEOFFSET=5108, CLDSENDLIKENORM=5109,
		       CLDSENDMINDATAFLUX1=5110, CLDSENDMAXDATAFLUX1=5111,
		       CLDSENDMINDATAFLUX2=5112, CLDSENDMAXDATAFLUX2=5113,
		       CLDSENDMINRFLUX1=5114, CLDSENDMAXRFLUX1=5115,
		       CLDSENDMINRFLUX2=5116, CLDSENDMAXRFLUX2=5117,
		       CLDSENDHASBEAM=5118, CLDSENDFFTSIZE=5119,
		       CLDSENDEDGEFIX=5120, CLDSENDEDGEINTEG=5121,
		       CLDSENDNEDGE=5122, CLDSENDNBEAM=5123,
		       CLDSENDSETNUM=5124, CLDSENDBINDATA=5125,
		       CLDSENDNBINS=5126, 
		       CLDSENDHASCFIRBPRIOR1=5127,
		       CLDSENDCFIRBPRIORMEAN1=5128, 
		       CLDSENDCFRIBPRIORSIGMA1=5129, 
		       CLDSENDHASCFIRBPRIOR2=5130,
		       CLDSENDCFIRBPRIORMEAN2=5131, 
		       CLDSENDCFRIBPRIORSIGMA2=5132,
		       CLDSENDHASPOISSONPRIOR1=5133,
		       CLDSENDPOISSONPRIORMEAN1=5134, 
		       CLDSENDPOISSONPRIORSIGMA1=5135, 
		       CLDSENDHASPOISSONPRIOR2=5136,
		       CLDSENDPOISSONPRIORMEAN2=5137, 
		       CLDSENDPOISSONPRIORSIGMA2=5138,
		       CLDSENDHASSIGMAPRIOR1=5139,
		       CLDSENDSIGMAPRIORWIDTH1=5140, 
		       CLDSENDHASSIGMAPRIOR2=5141,
		       CLDSENDSIGMAPRIORWIDTH2=5142,
		       CLDSENDREGULARIZATIONALPHA=5143,
		       IFKSENDNKNOTS=6000, IFKSENDKNOTPOS=6001,
		       IFKSENDKNOTVAL=6002, IFKHASRANGE=6003,
		       IFKSENDRANGE=6004, IFKHASLOWERLIMITS=6005,
		       IFKSENDHASLOWLIM=6006, IFKSENDLOWLIM=6007,
		       IFKHASUPPERLIMITS=6008, IFKSENDHASUPLIM=6009,
		       IFKSENDUPLIM=6010, IFDLNSENDNKNOTS=6500, 
		       IFDLNSENDNSIGMAS=6501, IFDLNSENDNOFFSETS=6502,
		       IFDLNSENDKNOTPOS=6503, IFDLNSENDKNOTVAL=6504, 
		       IFDLNHASRANGE=6505, IFDLNSENDRANGE=6506, 
		       IFDLNHASLOWERLIMITS=6507, IFDLNSENDHASLOWLIM=6508, 
		       IFDLNSENDLOWLIM=6509, IFDLNHASUPPERLIMITS=6510, 
		       IFDLNSENDHASUPLIM=6511, IFDLNSENDUPLIM=6512,
		       PMCMCSENDINIT=10000, PMCMCSENDINGINIT=10001,
		       PMCMCSENDNPAR=10002, PMCMCISREADY=10003};
}

typedef std::pair<double, double> dblpair;

#endif
