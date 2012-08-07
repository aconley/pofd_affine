#ifndef __pdfactorydouble__
#define __pdfactorydouble__

#include<string>

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_errno.h>
#include<fftw3.h>

#include<numberCountsDouble.h>
#include<PDDouble.h>
#include<paramSet.h>

/*!
  \brief Class for computing P(D) values from a set of parameters,
  two band version

  Only square FFTs are supported.
  
  Always call initPD before using getPD for a given model.
  Unlike the one band version, doesn't interpolate on R.  Two dimensional
  interpolation was either inaccurate or not faster than directly computing
  R in tests.
 */

class PDFactoryDouble {
 private :
  //Constants controlling edge integrations
  static const double lowEdgeRMult; //!< How far down to go from model limits in edge integrals
  static const bool use_edge_log_x; //!< Integrate Rx in x log space
  static const bool use_edge_log_y; //!< Integrate Ry in y log space

  bool initialized; //!< forward transformed R is filled

  unsigned int currsize; //!< Current memory allocation
  unsigned int lastfftlen; //!< FFT length of last transform
  double max_sigma1; //!< Current max supported sigma, band 1
  double max_sigma2; //!< Current max supported sigma, band 2
  double mn1; //!< Expected mean, band 1
  double mn2; //!< Expected mean, band 2
  double var_noi1; //!< Expected variance without instrumental noise, band 1
  double var_noi2; //!< Expected variance without instrumental noise, band 2
  double sg1; //!< Expected sigma (inc instrument noise), band 1
  double sg2; //!< Expected sigma (inc instrument noise), band 2

  fftw_plan plan;     //!< Holds forward transformation plan
  fftw_plan plan_inv; //!< Holds inverse transformation plan

  //Working variables for transformation
  bool rvars_allocated; //!< Are R variables (rest of this block) allocated
  double* rvals; //!< Working space for R computation, row major order
  fftw_complex *rtrans; //!< Holds FFTed rvals 
  fftw_complex* pval; //!< Working variable holding p = exp( stuff )
  double* pofd; //!< Internal P(D) variable.  
  double *RFlux1; //!< Holds R flux values for fill
  double *RFlux2; //!< Holds R flux values for fill

  void allocateRvars(); //!< Allocates R variables if needed
  void freeRvars(); //!< Free R variables

  //Edge variables
  bool edgevars_allocated; //!< Are Edge variables (this block) allocated
  unsigned int nedge; //!< Number of edge integral steps
  double* REdgeFlux1; //!< Holds flux for R edge integration
  double* REdgeFlux2; //!< Holds flux for R edge integration
  double* REdgeWork; //!< Holds R in edge integration
    
  void allocateEdgevars(); //!< Allocate Edge variables
  void freeEdgevars(); //!< Free edge variables

  double dflux1; //!< Flux size step of last computation, band 1
  double dflux2; //!< Flux size step of last computation, band 2
  unsigned int maxidx1; //!< Max non-zero index in R, band 2
  unsigned int maxidx2; //!< Max non-zero index in R, band 2
  bool doshift1; //!< Apply shifting, band 2
  bool doshift2; //!< Apply shifting, band 2
  double shift1; //!< Shift amount, band 2
  double shift2; //!< Shift amount, band 2

  unsigned int fftw_plan_style; //!< FFTW plan flags
  bool has_wisdom; //!< Has FFTW wisdom 
  std::string wisdom_file; //!< Name of wisdom file

  bool verbose; //!< Outputs information about each iter

  void init(unsigned int); //!< Initializes memory
  bool resize(unsigned int); //!< Sets transform size arrays
  void strict_resize(unsigned int); //!< Sets transform size arrays
  
  /*! \brief Get P(D) statistics from R computation (so without inst. noise) */
  void getRStats(unsigned int n, double& mn1, double& mn2,
		 double& var1, double& var2) const;

#ifdef TIMING
  std::clock_t RTime, RStatsTime, p0Time, fftTime, posTime, copyTime;
  std::clock_t normTime, edgeTime, meanTime, logTime;
#endif

 public :

  PDFactoryDouble(unsigned int nedge=256); //!< Default constructor
  PDFactoryDouble(const std::string&, unsigned int nedge=256 ); //!< Constructor with wisdom file
  ~PDFactoryDouble(); //!< Destructor

  void free(); //!< Frees memory

  void setVerbose() { verbose = true; } //!< Sets verbose mode
  void unsetVerbose() { verbose = false; } //!< Unset verbose mode

  /*! \brief Get size of last FFT */
  unsigned int getLastFFTLen() const { return lastfftlen; }

  /*! \brief Returns edge integration length */
  unsigned int getNEdge() const { return nedge; }
  void setNEdge(unsigned int); //!< Set nedge

  /*! \brief Adds FFTW wisdom file*/
  bool addWisdom(const std::string& filename);

  /*! \brief Initializes P(D) by computing R and forward transforming it*/
  void initPD(unsigned int, double, double, double, double, 
	      const numberCountsDouble&, const doublebeam&,
	      bool setEdge=true);

  /*! \brief Gets P(D) with specified noise levels */
  void getPD(double, double, PDDouble&, bool setLog=true, 
	     bool edgeFix=false);

  void SendSelf(MPI::Comm&, int dest) const; //!< MPI copy send operation
  void RecieveCopy(MPI::Comm&, int dest); //!< MPI copy recieve operation


#ifdef TIMING
  void resetTime(); //!< Reset timing information
  void summarizeTime(unsigned int=0) const; //!< Summarize timing information
#endif
};

#endif
