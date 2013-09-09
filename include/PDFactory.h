#ifndef __pdfactory__
#define __pdfactory__

#include<string>

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<fftw3.h>

#include "../include/numberCounts.h"
#include "../include/PD.h"
#include "../include/paramSet.h"

/*!
  \brief Class for computing P(D) values from a set of parameters,
  1D case

  Always call initPD before using getPD for a given model.
  Uses and interpolative approach -- fills R into an interpolation array,
  then get the full R from that, which is much faster.
 */
class PDFactory {
 private :
  static const double subedgemult; //!< Controls interpolation lower edge

  bool initialized; //!< some R transform is filled

  unsigned int currsize; //!< Current memory allocation for R, rtrans, etc.
  unsigned int lastfftlen; //!< FFT length of last transform
  double max_sigma; //!< Current supported max sigma
  double mn; //!< Expected mean
  double var_noi; //!< Expected variance without instrument noise
  double sg; //!< Expected sigma, inc instrument noise

  //Working variables for transformation.
  fftw_plan plan;     //!< Holds forward transform plan
  fftw_plan plan_inv; //!< Holds inverse transform plan
  bool rvars_allocated; //!< Are R variables allocated (rest of this block)
  double* rvals; //!< Working space for R computation 
  fftw_complex *rtrans; //!< Holds FFTed rvals 
  fftw_complex* pval; //!< Working variable holding p = exp( stuff ) 
  double* pofd; //!< Internal P(D) variable.  
  void allocateRvars(); //!< Allocates R variables if needed
  void freeRvars(); //!< Free R variables

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
  void freeInterp(); //!< Free interpolation variables

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

  void free(); //!< Frees memory

  void setVerbose() { verbose = true; } //!< Sets verbose mode
  void unsetVerbose() { verbose = false; } //!< Unset verbose mode

  /*! \brief Get fft size of last transformation */
  unsigned int getLastFFTLen() const { return lastfftlen; }

  /*! \brief Get number of interpolation elements for R */
  unsigned int getNInterp() const { return ninterp; }
  /*! \brief Set number of interpolation elements for R */
  void setNInterp(unsigned int);

  /*! \brief Adds FFTW wisdom file*/
  bool addWisdom(const std::string& filename);

  /*! \brief Initializes P(D) by computing R and forward transforming it*/
  bool initPD(unsigned int, double, double, 
	      const numberCounts&, const beam&);

  /*! \brief Gets P(D) of with specified noise level */
  void getPD(double, PD&, bool setLog=true, 
	     bool edgeFix=false);

  void sendSelf(MPI_Comm, int dest) const; //!< MPI copy send operation
  void recieveCopy(MPI_Comm, int dest); //!< MPI copy recieve operation

#ifdef TIMING
  void resetTime(); //!< Reset timing information
  void summarizeTime(unsigned int=0) const; //!< Output timing information
#endif

};

#endif
