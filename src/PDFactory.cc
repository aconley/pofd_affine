#include<sstream>
#include<cmath>
#include<cstring>
#include<limits>

#include<global_settings.h>
#include<PDFactory.h>
#include<affineExcept.h>

const double PDFactory::subedgemult = 1e-5;

/*!
 */
PDFactory::PDFactory(unsigned int NINTERP) { 
  init(NINTERP);
}

/*
  \param[in] wisfile Wisdom file filename
 */
PDFactory::PDFactory(const std::string& wisfile, unsigned int NINTERP ) {
  init(NINTERP);
  addWisdom(wisfile);
}

PDFactory::~PDFactory() {
  if (RinterpFlux != NULL) delete[] RinterpFlux;
  if (RinterpVals != NULL) delete[] RinterpVals;
  if (acc != NULL)  gsl_interp_accel_free(acc);
  if (spline != NULL) gsl_spline_free(spline);

  if (rvals != NULL) fftw_free(rvals);
  if (rtrans != NULL) fftw_free(rtrans);
  if (pofd != NULL) fftw_free(pofd);
  if (pval != NULL) fftw_free(pval);

  if (plan != NULL) fftw_destroy_plan(plan); 
  if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
}

void PDFactory::init(unsigned int NINTERP) {
  lastfftlen = 0;
  currsize = 0;
  ninterp = NINTERP;

#ifdef TIMING
  resetTime();
#endif

  acc = gsl_interp_accel_alloc();
  spline = NULL;

  interpvars_allocated = false;
  RinterpFlux = NULL;
  RinterpVals = NULL;

  plan = plan_inv = NULL;
  rvars_allocated = false;
  rvals = NULL;
  rtrans = NULL;
  pval = NULL;
  pofd = NULL;

  dflux = 0.0;

  verbose = false;
  has_wisdom = false;
  fftw_plan_style = FFTW_ESTIMATE;

  max_sigma = 0.0;
  initialized = false;

}

#ifdef TIMING
void PDFactory::resetTime() {
  RTime = p0Time = fftTime = posTime = copyTime = normTime = 0;
  edgeTime = meanTime = logTime = 0;
}

void PDFactory::summarizeTime(unsigned int nindent) const {
  std::string prestring(nindent,' ');
  prestring = "  ";
    
  std::cout << "R time: " << prestring 
	    << 1.0*RTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "p0 time: " << prestring 
	    << 1.0*p0Time/CLOCKS_PER_SEC << std::endl;
  std::cout << "fft time: " << prestring 
	    << 1.0*fftTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "pos time: " << prestring 
	    << 1.0*posTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "copy time: " << prestring 
	    << 1.0*copyTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "norm time: " << prestring 
	    << 1.0*normTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "edge time: " << prestring 
	    << 1.0*edgeTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "mean time: " << prestring 
	    << 1.0*meanTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "log time: " << prestring 
	    << 1.0*logTime/CLOCKS_PER_SEC << std::endl;
}
#endif


/*
  Doesn't force resize if there is already enough room available

  \param[in] NSIZE new size
  \returns True if a resize was needed
 */
bool PDFactory::resize(unsigned int NSIZE) {
  if (NSIZE > currsize) {
    strict_resize(NSIZE);
    return true;
  } else return false;
}

/*!
  Forces a resize

  \param[in] NSIZE new size (must be > 0)
 */
//I don't check for 0 NSIZE because that should never happen
// in practice, and it isn't worth the idiot-proofing
// rtrans always gets nulled in here, then refilled when you
//  call initPD
void PDFactory::strict_resize(unsigned int NSIZE) {
  if (NSIZE == currsize) return;
  freeRvars();
  currsize = NSIZE;
  initialized = false;
}

void PDFactory::allocateRvars() {
  if (rvars_allocated) return;
  if (currsize == 0)
    throw affineExcept("PDFactory","allocate_rvars",
		       "Invalid (0) currsize",1);
  rvals = (double*) fftw_malloc(sizeof(double)*currsize);
  pofd = (double*) fftw_malloc(sizeof(double)*currsize);
  unsigned int fsize = currsize/2+1;
  rtrans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fsize);
  pval = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fsize);
  rvars_allocated = true;
}

void PDFactory::freeRvars() {
  if (rvals != NULL) { fftw_free(rvals); rvals=NULL; }
  if (rtrans != NULL) {fftw_free(rtrans); rtrans=NULL; }
  if (pval != NULL) { fftw_free(pval); pval = NULL; }
  if (pofd != NULL) { fftw_free(pofd); pofd = NULL; }
  rvars_allocated = false;
  initialized = false;
}

/*!
  \param[in] NINTERP New interpolation length
 */
void PDFactory::setNInterp(unsigned int NINTERP) {
  if (NINTERP == ninterp) return;
  ninterp = NINTERP;
  if (interpvars_allocated) {
    freeInterp();
    allocateInterp();
  }
}

void PDFactory::allocateInterp() {
  if (interpvars_allocated) return;
  if (ninterp == 0)
    throw affineExcept("PDFactory","allocateInterp",
		       "Invalid (0) ninterp",1);
  
  //Note that acc is always allocated
  spline = gsl_spline_alloc( gsl_interp_cspline, 
			     static_cast<size_t>(ninterp));
  RinterpFlux = new double[ninterp];
  RinterpVals = new double[ninterp];
  interpvars_allocated = true;
}


void PDFactory::freeInterp() {
  if (RinterpFlux != NULL) { delete[] RinterpFlux; RinterpFlux = NULL; }
  if (RinterpVals != NULL) { delete[] RinterpVals; RinterpVals = NULL; }
  interpvars_allocated = false;
  initialized = false;
}

/*!
  Frees all internal memory
*/
void PDFactory::free() {
  freeRvars();
  freeInterp();
}

/*!
  \param[in] filename Name of wisdom file
*/
bool PDFactory::addWisdom(const std::string& filename) {
  FILE *fp = NULL;
  fp = fopen( filename.c_str(), "r" );
  if (!fp) {
    std::stringstream str;
    str << "Error opening wisdom file: " << filename;
    throw affineExcept("PDFactory","addWisdom",str.str(),1);
  }
  if (fftw_import_wisdom_from_file(fp) == 0) {
    std::stringstream str;
    str << "Error reading wisdom file: " << wisdom_file;
    throw affineExcept("PDFactory","addWisdom",str.str(),1);
  }
  fclose(fp);
  fftw_plan_style = FFTW_WISDOM_ONLY || FFTW_ESTIMATE;
  has_wisdom = true;
  wisdom_file = filename;
  if (plan != NULL) {
    fftw_destroy_plan(plan); 
    plan = NULL;
  }
  if (plan_inv != NULL) {
    fftw_destroy_plan(plan_inv);
    plan_inv = NULL;
  }

  initialized = false;

  return true;
}


/*!
  Gets ready for P(D) computation by preparing R, including forward
  transforming it.  The idea is to compute everything that doesn't
  require knowing sigma here so we can call this multiple times on
  maps with the same beam but different noise values.
 
  \param[in] n       Size of transform 
  \param[in] sigma   Maximum allowed sigma
  \param[in] maxflux Maximum flux generated in R
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] bm      Beam 

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux will end up being -less- than maxflux in practice
  by about the mean flux + 10 sigma.
 */
void PDFactory::initPD(unsigned int n, double sigma,
		       double maxflux, const numberCounts& model,
		       const beam& bm ) {

  //Make sure we have enough room
  bool did_resize = resize(n);

  if (!interpvars_allocated) allocateInterp();

  //Set whether we have beam
  bool has_pos = bm.hasPos();
  bool has_neg = bm.hasNeg();
  if ( ! ( has_pos || has_neg ) )
    throw affineExcept("PDFactory","initPD",
		       "Beam has neither positive nor negative bits",1);

  //Set min/max interpolation is filled in for; the lower
  // limit is tricky.  For now controlled by user
  //The max must be slightly less than the top value because R is zero there
  //The max value is the max model flux times the max pixel
  // so we have to have different versions for pos/neg
  double ininterpm1 = 1.0/static_cast<double>(ninterp-1);
  modelmin = model.getMinFlux();
  modelmax = model.getMaxFlux();

  double inm1 = 1.0 / static_cast<double>(n-1);
  dflux = maxflux * inm1;

  //Compute R.  We do this via interpolation -- compute R
  // for ninterp positions, then fill that into rvals.
  // We do pos and neg seperately because the interpolation works
  // better on each component, rather than interpolating on
  // the sum of R.  At least, that's the theory.

  //Allocate memory if needed; this is a way of not allocating
  // these until we run into a beam that needs them
  if (!rvars_allocated) allocateRvars();

  if (has_pos) {
    double mininterpflux = modelmin * subedgemult;
    double maxinterpflux = modelmax * bm.getMaxPos();
    //R[maxinterpflux] = 0.0, which we don't want to include in our
    // log/log interpolation, so step it back slightly
    double dinterp = 0.001*(maxinterpflux-mininterpflux)*ininterpm1;
    maxinterpflux -= dinterp;

    //Note that the interpolation is in log space
    double dinterpflux = (log2(maxinterpflux/mininterpflux))*ininterpm1;
    for (unsigned int i = 0; i < ninterp; ++i)
      RinterpFlux[i] = mininterpflux*exp2(static_cast<double>(i)*dinterpflux);

#ifdef TIMING
    std::clock_t starttime = std::clock();
#endif
    model.getR(ninterp,RinterpFlux,bm,RinterpVals,numberCounts::BEAMPOS);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

    //Load the values into the spline.  Note we take the log --
    // we are interpolating in log space 
    double val;
    for (unsigned int i = 0; i < ninterp; ++i) {
      val = RinterpVals[i];
      if (val > 0) RinterpVals[i] = log2(val); //log2 is faster than ln
      else RinterpVals[i] = pofd_mcmc::smalllogval;
    }
    //Load the Spline
    gsl_spline_init( spline, RinterpFlux, RinterpVals, 
                     static_cast<size_t>(ninterp) );

    //Now interpolate out; note r still needs to be multiplied by dflux
    // for use later.  We figure out which bins are covered by the
    // interpolation for efficiency
    int st = static_cast<int>( mininterpflux/dflux + 0.9999999999999999 );
    unsigned int minitidx = (st < 0) ? 0 : static_cast<unsigned int>(st);
    unsigned int maxitidx = static_cast<unsigned int>(maxinterpflux/dflux);
    if (maxitidx >= n) maxitidx=n-1;

    //Now interpolate, setting to zero outside the range
    for (unsigned int i = 0; i < minitidx; ++i)
      rvals[i] = 0.0;
    double cflux, splval;
    for (unsigned int i = minitidx; i <= maxitidx; ++i) {
      cflux = static_cast<double>(i)*dflux; //Min is always 0 in R
      splval = gsl_spline_eval( spline, cflux, acc );
      rvals[i] = exp2(splval);
    }
    for (unsigned int i = maxitidx+1; i < n; ++i)
      rvals[i] = 0.0;
  }
  if (has_neg) {
    //And now, we do basically the same thing for the negative beam
    double mininterpflux = modelmin * subedgemult;
    double maxinterpflux = modelmax*bm.getMinAbsNeg();
    double dinterp = 0.001*(maxinterpflux-mininterpflux)*ininterpm1;
    maxinterpflux -= dinterp;
    double dinterpflux = (log2(maxinterpflux/mininterpflux))*ininterpm1;
    for (unsigned int i = 0; i < ninterp; ++i)
      RinterpFlux[i] = 
        mininterpflux*exp2(static_cast<double>(i)*dinterpflux);

#ifdef TIMING
    std::clock_t starttime = std::clock();
#endif
    model.getR(ninterp,RinterpFlux,bm,RinterpVals,numberCounts::BEAMNEG);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

    double val;
    for (unsigned int i = 0; i < ninterp; ++i) {
      val = RinterpVals[i];
      if (val > 0) RinterpVals[i] = log2(val); //log2 is faster
      else RinterpVals[i] = pofd_mcmc::smalllogval;
    }
    gsl_spline_init( spline, RinterpFlux, RinterpVals, 
                     static_cast<size_t>(ninterp) );

    //Now interpolate out, same as before
    int st = static_cast<int>( mininterpflux/dflux + 0.9999999999999999 );
    unsigned int minitidx = (st < 0) ? 0 : static_cast<unsigned int>(st);
    unsigned int maxitidx = static_cast<unsigned int>(maxinterpflux/dflux);
    if (maxitidx >= n) maxitidx=n-1;

    //Interpolate out
    if (has_pos) {
      //In this case we are adding onto previous r values.  We can ignore
      // the ends, since we would just be adding 0 anyways.
      double cflux, splval;
      for (unsigned int i = minitidx; i <= maxitidx; ++i) {
	cflux = static_cast<double>(i)*dflux; //Min is always 0 in R
	splval = gsl_spline_eval( spline, cflux, acc );
	rvals[i] += exp2(splval);
      }
    } else {
      for (unsigned int i = 0; i < minitidx; ++i)
	rvals[i] = 0.0;
      double cflux, splval;
      for (unsigned int i = minitidx; i <= maxitidx; ++i) {
	cflux = static_cast<double>(i)*dflux; //Min flux is always 0 in R
	splval = gsl_spline_eval( spline, cflux, acc );
	rvals[i] = exp2(splval);
      }
      for (unsigned int i = maxitidx+1; i < n; ++i)
	rvals[i] = 0.0;
    }
  }

  //Now that we have R, use it to compute the mean and variance
  // mn = \int x R dx.
  mn = rvals[1]; //Noting that RFlux[0] = 0
  for (unsigned int i = 2; i < n-1; ++i)
    mn += rvals[i]*static_cast<double>(i);
  mn += 0.5*rvals[n-1]*static_cast<double>(n-1);
  mn *= dflux*dflux;  //Once for x, once for the bin size in the integral

  // var = \int x^2 R dx
  var_noi= rvals[1]; //Again, using the fact that RFlux[0] = 0
  double idbl;
  for (unsigned int i = 2; i < n-1; ++i) {
    idbl = static_cast<double>(i);
    var_noi += rvals[i]*idbl*idbl;
  }
  idbl = static_cast<double>(n-1);
  var_noi += 0.5*rvals[n-1]*idbl*idbl;
  var_noi *= dflux*dflux*dflux;

  //Now, compute the sigma for the maximum instrumental sigma
  // supported by this call to init
  sg = sqrt( var_noi + sigma*sigma );

  //Multiply R by dflux factor to represent the actual
  // number of sources in each bin
  for (unsigned int i = 0; i < n; ++i)
    rvals[i] *= dflux;

  //Make the plans, or keep the old ones if possible
  //The transform plan is that the forward transform
  // takes rvals to rtrans.  We then do things to rtrans
  // to turn it into pval, including shifting, adding noise,
  // etc, and then transform from pval to pofd.
  //Only the rvals->rtrans transformation happens in this
  // routine.  The other (pval->pofd), which depends on the value of
  // sigma for each map, happens in getPD.  But we can make
  // both plans here...

  //If we resized, we must make the new plans because the
  // addresses changed

  int intn = static_cast<int>(n);
  if ( did_resize || (lastfftlen != n) || (plan == NULL) ) {
    if (plan != NULL) fftw_destroy_plan(plan);
    plan = fftw_plan_dft_r2c_1d(intn, rvals, rtrans,
				fftw_plan_style);
  }
  if (plan == NULL) {
    std::stringstream str;
    str << "Plan creation failed for forward transform of size: " << 
      n << std::endl;
    throw affineExcept("PDFactory","initPD",str.str(),32);
  }

  if ( did_resize || (lastfftlen != n) || (plan_inv == NULL) ) {
    if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
    plan_inv = fftw_plan_dft_c2r_1d(intn, pval, pofd,
				    fftw_plan_style);
  }
  if (plan_inv == NULL) {
    std::stringstream str;
    str << "Plan creation failed for inverse transform of size: " << 
      n << std::endl;
    throw affineExcept("PDFactory","initPD",str.str(),64);
  }

  //Decide if we will shift and pad, and if so by how much
  //Only do shift if the noise is larger than one actual step size
  // Otherwise we can't represent it well.
  bool dopad = (sigma > dflux);
  doshift = ( dopad && ( mn < pofd_mcmc::n_sigma_shift*sg) );
  if ( doshift ) shift = pofd_mcmc::n_sigma_shift*sg - mn; else
    shift=0.0;

  if (verbose) {
    std::cout << " Initial mean estimate: " << mn << std::endl;
    std::cout << " Initial stdev estimate: " << sg << std::endl;
    if (doshift)
      std::cout << " Additional shift applied: " << shift << std::endl;
  }

  //Make sure that maxflux is large enough that we don't get
  // bad aliasing wrap from the top around into the lower P(D) values.
  if (maxflux <= pofd_mcmc::n_sigma_pad*sg)
    throw affineExcept("PDFactory","initPD","Top wrap problem",
		       128);

  //The other side of the equation is that we want to zero-pad the
  // top, and later discard that stuff.
  // The idea is as follows:
  // the 'target mean' of the calculated P(D) will lie at mn+shift.
  // We assume that anything within n_sigma_pad2d*sg
  // is 'contaminated'.  That means that, if n_sigma_pad*sg >
  // mn+shift, there will be some wrapping around the bottom of the P(D)
  // to contaminate the top by an amount n_sigma_pad*sg - (mn+shift).
  // We therefore zero pad and discard anything above
  // maxflux - (n_sigma_pad*sg - (mn+shift))
  if (dopad) {
    double contam = pofd_mcmc::n_sigma_pad*sg - (mn+shift);
    if (contam < 0) maxidx = n; else {
      double topflux = maxflux - contam;
      if (topflux < 0)
	throw affineExcept("PDFactory","initPD","Padding problem",
			   256);
      maxidx = static_cast< unsigned int>( topflux/dflux );
      if (maxidx > n)
	throw affineExcept("PDFactory","initPD",
			 "Padding problem",
			 512);
      //Actual padding
      for (unsigned int i = maxidx; i < n; ++i)
	rvals[i] = 0.0;
    }
  } else maxidx = n;


  //Compute forward transform of this r value, store in rtrans
  //Have to use argument version, since the address of rtrans can move
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute_dft_r2c(plan,rvals,rtrans); 
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif
  
  lastfftlen = n;
  max_sigma = sigma;
  initialized = true;

}

/*!
  Computes the P(D) for a specific value of the instrument noise,
  assuming that initPD has been called first.
 
  \param[in] sigma Instrument noise sigma
  \param[out] pd Holds P(D) on output
  \param[in] setLog If true, pd is log(P(D) on output; convenient
              for likelihood evaluation.
  \param[in] edgeFix  Apply a fix to the lower edges to minimize wrapping
                      effects using a Gaussian to each row/col

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux often won't quite match the target values.

  You must call initPD first, or bad things will probably happen.
*/
void PDFactory::getPD( double sigma, PD& pd, bool setLog, 
		       bool edgeFix) {

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise and all that fun stuff, filling pd
  // for output

  if (! initialized )
    throw affineExcept("PDFactory","getPD",
		       "Must call initPD first",1);
  if (sigma > max_sigma) {
    std::stringstream errstr("");
    errstr << "Sigma value " << sigma
	   << " larger than maximum prepared value " << max_sigma
	   << std::endl;
    errstr << "initPD should have been called with at least " << sigma;
    throw affineExcept("PDFactory","getPD",errstr.str(),2);
  }

  //Output array from 2D FFT is n/2+1
  unsigned int n = lastfftlen;
  unsigned int ncplx = n/2 + 1;
      
  //Calculate p(omega) = exp( r(omega) - r(0) ),
  // convolving together all the bits into pval, which
  // is what we will transform back into pofd.
  // The forward transform of R is stored in rtrans 
  // There are some complications because of shifts and all that.
  // The frequencies are:
  //  f = i/dflux*n  
  // We work in w instead of f (2 pi f)
  // and actually compute 
  //  exp( r(omega) - r(0) - i*shift*omega - 1/2 sigma^2 omega^2 )
  if (verbose) std::cout << "  Computing p(w)" << std::endl;

#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif

  double r0, expfac, rval, ival;
  r0 = rtrans[0][0]; //r[0] is pure real
  double iflux = mcmc_affine::two_pi / (n * dflux);

  if (doshift) {
    //This should be the most common case,
    // and corresponds to having some noise
    double sigfac = 0.5*sigma*sigma;
    double w;
    for (unsigned int idx = 1; idx < ncplx; ++idx) {
      w    = iflux * static_cast<double>(idx);
      rval = rtrans[idx][0] - r0 - sigfac*w*w;
      ival = rtrans[idx][1] - shift*w;
      expfac = exp( rval );
      pval[idx][0] = expfac*cos(ival);
      pval[idx][1] = expfac*sin(ival);
    } 
  } else {
    //No shift, sigma must be zero
    for (unsigned int idx = 1; idx < ncplx; ++idx) {
      expfac = exp(rtrans[idx][0]-r0);
      ival = rtrans[idx][1];
      pval[idx][0] = expfac*cos(ival);
      pval[idx][1] = expfac*sin(ival);
    }
  }
  //p(0) is special
  pval[0][0] = 1.0;
  pval[0][1] = 0.0;

#ifdef TIMING
  p0Time += std::clock() - starttime;
#endif

  //Transform back
  if (verbose) std::cout << " Reverse transform" << std::endl;
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute(plan_inv); //overwrites pofd with reverse transform of pval
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif

  //Enforce positivity
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int idx = 0; idx < n; ++idx)
    if (pofd[idx] < 0) pofd[idx] = 0.0;
#ifdef TIMING
  posTime += std::clock() - starttime;
#endif

  //Copy into output variable
#ifdef TIMING
  starttime = std::clock();
#endif
  pd.resize(maxidx);
  for (unsigned int i = 0; i < maxidx; ++i)
    pd.pd_[i] = pofd[i];
  pd.logflat = false;
  pd.minflux = 0.0; pd.dflux = dflux;
#ifdef TIMING
  copyTime += std::clock() - starttime;
#endif

  //Normalize
#ifdef TIMING
  starttime = std::clock();
#endif
  pd.normalize();
#ifdef TIMING
  normTime += std::clock() - starttime;
#endif

  //Fix up the edge.  Must be done after normalization
#ifdef TIMING
  starttime = std::clock();
#endif
  if (edgeFix) pd.edgeFix();
#ifdef TIMING
  edgeTime += std::clock() - starttime;
#endif

  //Now mean subtract flux axis
#ifdef TIMING
  starttime = std::clock();
#endif
  double tmn; //True mean
  pd.getMean(tmn,false);
  if ( std::isinf(tmn) || std::isnan(tmn) ) {
    std::stringstream str;
    str << "Un-shift amounts not finite: " << tmn << " " << std::endl;
    str << "At length: " << n << " with noise: " << sigma;
    throw affineExcept("PDFactory","getPD",str.str(),8);
  }
  if (verbose) {
    std::cerr << " Expected mean: " << shift+mn 
	      << " Realized mean: " << tmn << std::endl;
  }
  pd.minflux = -tmn;
#ifdef TIMING
  meanTime += std::clock() - starttime;
#endif

  //Turn PD to log for more efficient log computation of likelihood
#ifdef TIMING
  starttime = std::clock();
#endif
  if (setLog) pd.applyLog(false);
#ifdef TIMING
  logTime += std::clock() - starttime;
#endif
}
 
void PDFactory::SendSelf(MPI::Comm& comm, int dest) const {
  comm.Send(&fftw_plan_style,1,MPI::UNSIGNED,dest,pofd_mcmc::PDFSENDPLANSTYLE);
  comm.Send(&has_wisdom,1,MPI::UNSIGNED,dest,pofd_mcmc::PDFHASWISDOM);
  if (has_wisdom) {
    //Send wisdom file name
    unsigned int nstr = wisdom_file.size()+1;
    char *cstr = new char[nstr];
    std::strncpy( cstr, wisdom_file.c_str(), nstr );
    comm.Send(&nstr,1,MPI::UNSIGNED,dest,pofd_mcmc::PDFWISLEN);
    comm.Send(cstr, nstr, MPI::CHAR, dest,pofd_mcmc::PDFWISNAME);
    delete[] cstr;
  }
  comm.Send(&verbose,1,MPI::BOOL,dest,pofd_mcmc::PDFVERBOSE);
  comm.Send(&ninterp,1,MPI::UNSIGNED,dest,pofd_mcmc::PDFNINTERP);
}

//Note this doesn't copy over interal variables
void PDFactory::RecieveCopy(MPI::Comm& comm, int src) {
  comm.Recv(&fftw_plan_style,1,MPI::UNSIGNED,src,pofd_mcmc::PDFSENDPLANSTYLE);
  comm.Recv(&has_wisdom,1,MPI::UNSIGNED,src,pofd_mcmc::PDFHASWISDOM);
  if (has_wisdom) {
    //Recieve wisdom file name
    unsigned int nstr;
    comm.Recv(&nstr,1,MPI::UNSIGNED,src,pofd_mcmc::PDFWISLEN);
    char *cstr = new char[nstr];
    comm.Recv(cstr, nstr, MPI::CHAR, src,pofd_mcmc::PDFWISNAME);    
    wisdom_file = std::string(cstr);
    delete[] cstr;
    addWisdom(wisdom_file);
  }
  comm.Recv(&verbose,1,MPI::BOOL,src,pofd_mcmc::PDFVERBOSE);
  unsigned int nnewinterp;
  comm.Recv(&nnewinterp,1,MPI::UNSIGNED,src,pofd_mcmc::PDFNINTERP);
  if (nnewinterp != ninterp) {
    freeInterp();
    ninterp = nnewinterp;
  }
  initialized = false;
}
