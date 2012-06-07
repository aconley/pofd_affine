#ifndef __pdfactory__
#define __pdfactory__

#include<string>

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<fftw3.h>

#include<numberCounts.h>
#include<PD.h>
#include<paramSet.h>



/*!
  Always call initPD before using getPD for a given model.
 */

//This uses an interpolative approach -- fill R into an interpolation
// array, then get the full R from that, which is much faster than
// a full r fill
class PDFactory {
 private :
  static const double subedgemult; //!< Controls interpolation lower edge

  bool initialized; //!< some R transform is filled

  unsigned int currsize; //!< Current memory allocation
  unsigned int lastfftlen; //!< FFT length of last transform
  double max_sigma; //!< Current supported max sigma
  double mn; //!< Expected mean
  double var_noi; //!< Expected variance without instrument noise
  double sg; //!< Expected sigma, inc instrument noise

  //Working variables for transformation.
  fftw_plan plan, plan_inv; //!< Hold plans
  bool rvars_allocated; //!< Are R variables allocated (rest of this block)
  double* rvals; //!< Working space for R computation 
  fftw_complex *rtrans; //!< Holds FFTed rvals 
  fftw_complex* pval; //!< Working variable holding p = exp( stuff ) 
  double* pofd; //!< Internal P(D) variable.  
  void allocateRvars(); //!< Allocates R variables if needed
  void freeRvars();

  //Flux variables
  double dflux; //!< Flux size step of last computation
  unsigned int maxidx; //!< Max non-zero index in Rs
  bool doshift; //!< Apply shifting
  double shift; //!< Shift amount

  unsigned int fftw_plan_style; //!< FFTW plan flags
  bool has_wisdom; //!< Has FFTW wisdom 
  std::string wisdom_file; //!< Name of wisdom file

  bool verbose; //!< Outputs information about each iter

  //GSL interpolation stuff
  unsigned int ninterp; //!< Interpolation length
  double modelmin; //!< Min flux in model
  double modelmax; //!< Max flux in model
  gsl_interp_accel *acc; //!< Spline lookup accelerator
  bool interpvars_allocated; //!< Are interpolation variables allocated
  gsl_spline *spline; //!< Spline holding R to interpolate
  double *RinterpFlux; //!< Flux values to interpolate with
  double *RinterpVals; //!< R values 
  void allocateInterp(); //!< Allocates interpolation variables if needed
  void freeInterp();

  void init(unsigned int); //!< Initializes memory
  bool resize(unsigned int); //!< Sets transform size arrays
  void strict_resize(unsigned int); //!< Sets transform size arrays
  
#ifdef TIMING
  std::clock_t RTime, p0Time, fftTime, posTime, copyTime, normTime, edgeTime;
  std::clock_t meanTime, logTime;
#endif

 public :
  PDFactory(unsigned int NINTERP=1024); //!< Default constructor
  PDFactory(const std::string&,unsigned int NINTERP=1024); //!< Constructor with wisdom file
  ~PDFactory(); //!< Destructor

  void setVerbose() { verbose = true; } //!< Sets verbose mode
  void unsetVerbose() { verbose = false; } //!< Unset verbose mode

  unsigned int getLastFFTLen() const { return lastfftlen; }

  unsigned int getNInterp() const { return ninterp; }

  /*! \brief Adds wisdom file*/
  bool addWisdom(const std::string& filename);

  /*! \brief Initializes P(D) by computing R */
  void initPD(unsigned int, double, double, 
	      numberCounts&, const beam&);

  /*! \brief Gets P(D) of specified transform size */
  void getPD(double, PD&, bool setLog=true, 
	     bool edgeFix=false);

  void SendSelf(MPI::Comm&, int dest) const;
  void RecieveCopy(MPI::Comm&, int dest);

#ifdef TIMING
  void resetTime();
  void summarizeTime(unsigned int=0) const;
#endif

};

#endif
