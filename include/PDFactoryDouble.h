#ifndef __pdfactorydouble__
#define __pdfactorydouble__

#include<string>

#ifdef TIMING
#include<ctime>
#endif

#include<gsl/gsl_errno.h>
#include<fftw3.h>

#include "../include/numberCountsDouble.h"
#include "../include/PDDouble.h"
#include "../include/paramSet.h"

/*!
  \brief Class for computing P(D) values from a set of parameters,
  2D case.

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
  
  bool rinitialized; //!< R, RFluxes are filled
  bool initialized; //!< forward transformed R is filled

  unsigned int currsize; //!< Current memory allocation
  double mn1; //!< Expected mean, band 1
  double mn2; //!< Expected mean, band 2
  double var_noi1; //!< Expected variance without instrumental noise, band 1
  double var_noi2; //!< Expected variance without instrumental noise, band 2
  double sg1; //!< Expected sigma (inc instrument noise), band 1 of most recent
  double sg2; //!< Expected sigma (inc instrument noise), band 2 of most recent

  unsigned int fftw_plan_style; //!< FFTW plan flags
  bool has_wisdom; //!< Has FFTW wisdom 
  std::string wisdom_file; //!< Name of wisdom file
  fftw_plan plan;     //!< Holds forward transformation plan
  fftw_plan plan_inv; //!< Holds inverse transformation plan

  //Working variables for transformation
  bool rvars_allocated; //!< Are R variables (rest of this block) allocated
  double* rvals; //!< Working space for R computation, row major order
  double* rsum; //!< Sum of rvals along one dimension
  bool rdflux; //!< Has rvals been multiplied by dflux1 * dflux2?
  fftw_complex *rtrans; //!< Holds FFTed rvals 
  fftw_complex* pval; //!< Working variable holding p = exp(stuff)
  double* pofd; //!< Internal P(D) variable.  
  double *RFlux1; //!< Holds R flux values for fill
  unsigned int RwrapIdx1; //!< Pos/neg wrap index in RFlux1
  double *RFlux2; //!< Holds R flux values for fill
  unsigned int RwrapIdx2; //!< Pos/neg wrap index in RFlux2

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
  bool doshift1; //!< Apply shifting, band 2
  bool doshift2; //!< Apply shifting, band 2
  double shift1; //!< Shift amount, band 2
  double shift2; //!< Shift amount, band 2

  bool verbose; //!< Outputs information about each iter

  void init(unsigned int); //!< Initializes memory
  void resize(unsigned int); //!< Sets transform size arrays

  /*! \brief Sets RFlux1, Rflux2, with wrapping */
  void initRFlux(unsigned int n, double minflux1, double maxflux1,
		 double minflux2, double maxflux2);

  /*! \brief Fills in rvals */
  void initR(unsigned int n, double minflux1, double maxflux1,
	     double minflux2, double maxflux2, 
	     const numberCountsDouble&, const doublebeam&, 
	     bool setEdge=true, bool muldr=false); 

  /*! \brief Get P(D) statistics from R computation (so without inst. noise) */
  void getRStats();

  /*! \brief Helper function for unwrapAndNormalizePD */
  unsigned int findSplitPoint(const double* const, double,
			      double) const;

  void unwrapAndNormalizePD(PDDouble& pd) const; //!< Moves computed P(D) to output var

  void setupTransforms(unsigned int); //!< Sets up FFTW plans and resizes
#ifdef TIMING
  mutable std::clock_t RTime, p0Time, fftTime, posTime, copyTime, splitTime;
  mutable std::clock_t normTime, edgeTime, meanTime, logTime, starttime;
#endif

 public :

  PDFactoryDouble(unsigned int nedge=256); //!< Default constructor
  PDFactoryDouble(const std::string&, unsigned int nedge=256 ); //!< Constructor with wisdom file
  ~PDFactoryDouble(); //!< Destructor

  void free(); //!< Frees memory

  void setVerbose() { verbose = true; } //!< Sets verbose mode
  void unsetVerbose() { verbose = false; } //!< Unset verbose mode

  /*! \brief Returns edge integration length */
  unsigned int getNEdge() const { return nedge; }
  void setNEdge(unsigned int); //!< Set nedge

  /*! \brief Adds FFTW wisdom file*/
  void addWisdom(const std::string& filename);

  /*! \brief Initializes P(D) by computing R and forward transforming it*/
  bool initPD(unsigned int n, double minflux1, double maxflux1, 
	      double minflux2, double maxflux2, 
	      const numberCountsDouble&, const doublebeam& bm,
	      bool setEdge=true);

  /*! \brief Gets P(D) with specified noise levels */
  void getPD(double, double, PDDouble&, bool setLog=true);

  /*! \brief Write out current R as an HDF5 file */
  void writeRToHDF5(const std::string&) const;
  
  void sendSelf(MPI_Comm, int dest) const; //!< MPI copy send operation
  void receiveCopy(MPI_Comm, int dest); //!< MPI copy receive operation

#ifdef TIMING
  void resetTime(); //!< Reset timing information
  void summarizeTime(unsigned int=0) const; //!< Summarize timing information
#endif
};

#endif
