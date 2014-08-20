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

  bool rinitialized; //!< R initialized
  bool initialized; //!< R forward transform is filled

  unsigned int currsize; //!< Current memory allocation for R, rtrans, etc.
  double mn; //!< Expected mean
  double var_noi; //!< Expected variance without instrument noise
  double sg; //!< Expected sigma, inc instrument noise, or most recent P(D)

  // Plans
  unsigned int fftw_plan_style; //!< FFTW plan flags
  bool has_wisdom; //!< Has FFTW wisdom 
  std::string wisdom_file; //!< Name of wisdom file
  fftw_plan plan;     //!< Holds forward transform plan
  fftw_plan plan_inv; //!< Holds inverse transform plan

  // Variables setting up the flux values for R.  For efficiency,
  // we don't store these as an array, just the information we
  // need to compute.  If wrapRidx is currsize-1, there are no neg values
  double minflux_R; //!< Minimum value of R held
  unsigned int wrapRidx; //!< Maximum index of positive R values; above this neg

  //Working variables for transformation.
  bool rvars_allocated; //!< Are R variables allocated (rest of this block)
  double* rvals; //!< Working space for R computation 
  bool rdflux; //!< Has R been multiplied by dflux
  fftw_complex *rtrans; //!< Holds FFTed rvals 
  fftw_complex* pval; //!< Working variable holding p = exp( stuff ) 
  double* pofd; //!< Internal P(D) variable.  
  void allocateRvars(); //!< Allocates R variables if needed
  void freeRvars(); //!< Free R variables

  //Flux variables
  double dflux; //!< Flux size step of last computation
  bool doshift; //!< Apply shifting
  double shift; //!< Shift amount

  bool verbose; //!< Outputs information about each iter

  //GSL interpolation stuff.  We need two sets, one for pos, one for
  //  neg
  unsigned int ninterp; //!< Interpolation length
  gsl_interp_accel *acc_pos; //!< Spline lookup accelerator
  bool interpvars_allocated_pos; //!< Are interpolation variables allocated, pos
  gsl_spline *spline_pos; //!< Spline holding R to interpolate
  double *RinterpFlux_pos; //!< Flux values to interpolate with
  double *RinterpVals_pos; //!< R values in interpolation
  void allocateInterpPos(); //!< Allocates Pos interpolation variables if needed
  void freeInterpPos(); //!< Free interpolation variables
  gsl_interp_accel *acc_neg; //!< Spline lookup accelerator
  bool interpvars_allocated_neg; //!< Are interpolation variables allocated, neg
  gsl_spline *spline_neg; //!< Spline holding R to interpolate
  double *RinterpFlux_neg; //!< Flux values to interpolate with
  double *RinterpVals_neg; //!< R values in interpolation
  void allocateInterpNeg(); //!< Allocates Pos interpolation variables if needed
  void freeInterpNeg(); //!< Free interpolation variables

  void init(unsigned int); //!< Initializes memory
  bool resize(unsigned int); //!< Sets transform size arrays
  
  void setupTransforms(unsigned int); //!< Sets up FFTW plans

  /*! \brief Fills in R interpolation variables */
  void initRInterp(const numberCounts&, const beam&); 

  /*! \brief Sets RFlux, with wrapping */
  void initRFlux(unsigned int n, double minflux, double maxflux);

  /*! \brief Fills in rvals */
  void initR(unsigned int n, double minflux, double maxflux,
	     const numberCounts&, const beam&, bool muldr=false); 

  void getMeanVarFromR(); //!< Compute mean and variance from R

  /*! \brief Helper function for unwrapAndNormalizePD */
  unsigned int findSplitPoint() const;

  /*! \brief Moves P(D) over to output variable inside getPD */
  void unwrapAndNormalizePD(PD& pd) const;

#ifdef TIMING
  mutable std::clock_t RTime, p0Time, fftTime, posTime, copyTime;
  mutable std::clock_t normTime, edgeTime, meanTime, logTime, starttime;
#endif

 public :
  PDFactory(unsigned int NINTERP=2048); //!< Default constructor
  PDFactory(const std::string&,unsigned int NINTERP=2048); //!< Constructor with wisdom file
  ~PDFactory(); //!< Destructor

  void free(); //!< Frees memory

  void setVerbose() { verbose = true; } //!< Sets verbose mode
  void unsetVerbose() { verbose = false; } //!< Unset verbose mode

  /*! \brief Get number of interpolation elements for R */
  unsigned int getNInterp() const { return ninterp; }
  /*! \brief Set number of interpolation elements for R */
  void setNInterp(unsigned int);

  /*! \brief Adds FFTW wisdom file*/
  void addWisdom(const std::string& filename);

  /*! \brief Initializes P(D) by computing R and forward transforming it*/
  bool initPD(unsigned int n, double minflux, double maxflux,
	      const numberCounts&, const beam&);

  /*! \brief Gets P(D) of with specified noise level */
  void getPD(double, PD&, bool setLog=true);

  /*! \brief Write computed R out to HDF5 file */
  void writeRToHDF5(const std::string& filename) const;

  void sendSelf(MPI_Comm, int dest) const; //!< MPI copy send operation
  void receiveCopy(MPI_Comm, int dest); //!< MPI copy receive operation

#ifdef TIMING
  void resetTime(); //!< Reset timing information
  void summarizeTime(unsigned int=0) const; //!< Output timing information
#endif

};

#endif
