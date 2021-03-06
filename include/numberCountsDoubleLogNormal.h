//numberCountsDoubleLogNormal.h

#ifndef __numberCountsDoubleLogNormal__
#define __numberCountsDoubleLogNormal__

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_roots.h>

#include "../include/global_settings.h"
#include "../include/numberCountsDouble.h"
#include "../include/ran.h"

/*!
  \brief Spline number counts model for 2D case.

  The counts are modeled as a spline at one frequency times a Log-Normal 
  colour (flux ratio) model.  

  The full expression is
  \f[
  \frac{dN}{dS_1\, dS_2} = \frac{n_1\left(S_1 \right)}{S_1} 
  \mathrm{L} \left( \frac{S_2}{S_1}; \mu\left(S_1\right),
  \sigma\left(S_1\right) \right)
  \f] 
  where \f$n_1\f$ is just the spline as in the one-dimensional case
  as implemented in numberCountsKnotsSpline.
  The division by \f$S_1\f$ may seem odd, but assures that the
  distribution marginalized over \f$S_2\f$ is just \f$ n_1 \f$.
  \f$L\f$ is just the Log-Normal distribution:
  \f[
  \mathrm{L}\left( x ; \mu, \sigma \right) =
  \frac{1}{\sqrt{2 \pi} \sigma} \frac{1}{x}
  \exp\left[ \frac{ - \left(\log x - \mu\right)^2 }{ 2 \sigma^2 } \right] .
  \f]
  \f$\sigma\f$ and \f$\mu\f$ are stored as splines as functions of
  \f$S_1\f$.  Note that \f$\sigma\f$ and \f$\mu\f$ (the actual parameters 
  of the model)
  are not the mean and square root of the variance of \f$S_2 / S_1\f$,
  but rather the mean and square root of the variance of the 
  (natural) log quantities.  These are related by
  \f[
  \left< S_2/S_1 \right> = \exp \left[ \mu\left(S_1\right) +
  \frac{1}{2}\sigma^2\left(S_1\right) \right]
  \f]
  and
  \f[
  Var[ S_2/S_1 ] = \exp \left( 2 \mu \left( S_1 \right) +
  \sigma^2\left(S_1\right) \right) 
  \left[ \exp \left( \sigma^2\left(S_1\right) \right) - 1 \right] .
  \f]
  Note that the number counts require that \f$S_2/S_1 > 0\f$.
  We also explicitly require that \f$S_1 > 0\f$.  Both are justified
  on physical grounds.
    
  \ingroup Models
*/
class numberCountsDoubleLogNormal : public numberCountsDouble {
 protected :
  
  //Number counts in band one
  unsigned int nknots; //!< Number of knot positions for counts in band one
  double *knots; //!< Locations of knots
  double *logknots; //!< Log of knot positions (stored as log 2)
  double *logknotvals; //!< Log 2 values of differential number counts at knots, band 1
  gsl_interp_accel *acc; //!< Spline lookup accelerator
  gsl_spline *splinelog; //!< Spline in log2/log2 space of values

  //LogNormal for band 2
  unsigned int nsigmaknots; //!< Number of knot positions in sigma
  double *sigmaknots; //!< Positions of sigma knots
  double *sigmavals; //!< Sigma values at positions of knots
  double min_sigma; //!< Minimum allowable sigma value
  gsl_interp_accel *accsigma; //!< Sigma accelerator
  gsl_interp *sigmainterp; //!< Sigma interpolator
  unsigned int noffsetknots; //!< Number of knot positions in offset
  double *offsetknots; //!< Positions of offset knots
  double *offsetvals; //!< Offset values at positions of knots
  gsl_interp_accel *accoffset; //!< Offset accelerator
  gsl_interp *offsetinterp; //!< Offset interpolator

  //Model validity stuff
  mutable bool knots_valid; //!< Knots ready
  mutable bool sigmas_valid; //!< Sigmas ready
  mutable bool offsets_valid; //!< Offsets ready
  bool knotpos_loaded; //!< Knot positions loaded (not validity checked)
  bool sigmapos_loaded; //!< Sigma positions loaded (not validity checked)
  bool offsetpos_loaded; //!< Offset positions loaded (not validity checked)
  bool knotvals_loaded; //!< Knot values loaded (not validity checked)
  bool sigmavals_loaded; //!< Sigma values loaded (not validity checked)
  bool offsetvals_loaded; //!< Offset values loaded (not validity checked)
  void checkKnotsValid() const noexcept; //!< Set knots_valid by checking
  void checkSigmasValid() const noexcept; //!< Set sigmas_valid by checking
  void checkOffsetsValid() const noexcept; //!< Set offsets_valid by checking

  // Workspace.  These are temporary arrays used for efficiency.
  // They are useful if we are going to compute R many times (which we
  // often are).  The sign convention for f2work_sgn may seem rather odd,
  // but it has to do with the way things are stored in RWork, where it
  // also operates like a sign offset.
  gsl_integration_workspace *gsl_work; //!< Integration workspace for QAG
  mutable gsl_root_fsolver *fslv; //!< Root solver for use in finding model max flux
  mutable unsigned int nRWork; //!< Number of elements in R working arrays
  mutable double* RWork; //!< R working array; for pre-computing R bits
  void setRWorkSize(unsigned int) const; //!< Controls R working arrays
  mutable unsigned int nf2work; //!< Number of elements in f2work array
  mutable unsigned int* f2work_sgn; //!< Sign of f2; 0 for +, 1 for -, 2 for 0
  mutable double* f2work_inv; //!< 1 / |f2|
  mutable double* f2work_log; //!< log(|f2|)
  void setf2WorkSize(unsigned int) const; //!< Controls f2work arrays

  //Convenience functions for spline computations without input checking
  double getSigmaInner(double) const noexcept; //!< Inner sigma computation
  double getOffsetInner(double) const noexcept; //!< Inner offset computation
  double getNumberCountsInner(double,double) const noexcept; //!< Inner number counts computation

  // Finds value where Log Normal value goes to specified value down from peak
  double logNormalSolver(double mu, double sig, double nsig) const;

  /*! \brief Get range over which R is expected to be nonzero, unchecked version */
  std::pair<dblpair, dblpair> getRRangeInternal(const doublebeam&) const 
    throw(affineExcept);

  /*! \brief Integrate powers of fluxes over number counts */
  double splineInt(double alpha, double beta) const noexcept;

  /*! \brief Number of elements in varr (lots of model params!) */
  static constexpr unsigned int nvarr = 18; 
  void **varr; //!< Internal evil casting array for splineInt

  bool isValidLoaded() const noexcept; //!< Are loaded values valid?

 public :
  numberCountsDoubleLogNormal(); //!< Default
  explicit numberCountsDoubleLogNormal(unsigned int, unsigned int, 
                                       unsigned int); //!< Constructor

  /*! \brief Constructor */
  numberCountsDoubleLogNormal(const std::vector<float>&, 
                              const std::vector<float>&,
                              const std::vector<float>&); //!< Constructor
  /*! \brief Constructor */
  numberCountsDoubleLogNormal(const std::vector<double>&, 
                              const std::vector<double>&,
                              const std::vector<double>&); //!< Constructor
  /*! \brief Constructor */
  numberCountsDoubleLogNormal(unsigned int, const float* const,
                              unsigned int, const float* const,
                              unsigned int, const float* const);
  /*! \brief Constructor */
  numberCountsDoubleLogNormal(unsigned int, const double* const,
                              unsigned int, const double* const,
                              unsigned int, const double* const);
  /*! \brief Copy constructor */
  numberCountsDoubleLogNormal(const numberCountsDoubleLogNormal& other);
  virtual ~numberCountsDoubleLogNormal(); //!< Destructor

  /*! \brief Copy operator */
  numberCountsDoubleLogNormal& operator=(const numberCountsDoubleLogNormal&);

  bool isValid() const noexcept override; //!< Are model settings loaded and valid?

  /*! \brief Relative distance between two sets of params over model params*/
  float paramRelativeDistance(const paramSet& p1, const paramSet& p2)
    const throw(affineExcept) override;

  void setNKnots(unsigned int n); //!< Sets number of knots in band 1
  unsigned int getNKnots() const noexcept { return nknots; } //!< Number of knots in band 1 spline
  /*! \brief Sets knot positions, band 1, vector version */
  void setKnotPositions(const std::vector<float>&);
  /*! \brief Sets knot positions, band 1, vector version */
  void setKnotPositions(const std::vector<double>&);
  /*! \brief Sets knot positions, band 1, c array version */
  void setKnotPositions(unsigned int, const float* const);
  /*! \brief Sets knot positions, band 1, c array version */
  void setKnotPositions(unsigned int, const double* const);
  /*! \brief Sets number of sigma positions for color model */
  void setNSigmas(unsigned int);
  /*! \brief Gets number of sigma knots in band 2 color model */
  unsigned int getNSigmas() const noexcept { return nsigmaknots; }
  /*! \brief Sets knot positions for band 2 sigma model, vector version */
  void setSigmaPositions(const std::vector<float>&);
  /*! \brief Sets knot positions for band 2 sigma model, vector version */
  void setSigmaPositions(const std::vector<double>&);
  /*! \brief Sets knot positions for band 2 sigma model, c array version */
  void setSigmaPositions(unsigned int, const float* const);
  /*! \brief Sets knot positions for band 2 sigma model, c array version */
  void setSigmaPositions(unsigned int, const double* const);
  /*! \brief Set minimum sigma value */
  void setMinSigma(double);
  /*! \brief Get minimum sigma value */
  double getMinSigma() const noexcept { return min_sigma; }
  /*! \brief Sets number of offset positions for color model */
  void setNOffsets(unsigned int);
  /*! \brief Gets number of offset knots in band 2 color model */
  unsigned int getNOffsets() const noexcept { return noffsetknots; }
  /*! \brief Sets knot positions for band 2 offset model, vector version */
  void setOffsetPositions(const std::vector<float>&);
  /*! \brief Sets knot positions for band 2 offset model, vector version */
  void setOffsetPositions(const std::vector<double>&);
  /*! \brief Sets knot positions for band 2 offset model, c array version */
  void setOffsetPositions(unsigned int, const float* const);
  /*! \brief Sets knot positions for band 2 offset model, c array version */
  void setOffsetPositions(unsigned int, const double* const);

  //Get individual positions
  double getKnotPosition(unsigned int i) const; //!< Get band 1 model knot position (unchecked)
  double getSigmaPosition(unsigned int i) const; //!< Get sigma knot position (unchecked)
  double getOffsetPosition(unsigned int i) const; //!< Get offset knot position (unchecked)

  /*! \brief Gets total number of parameters */
  unsigned int getNParams() const noexcept { return nknots + nsigmaknots + noffsetknots; }
  /*! \brief Gets positions of all knots */
  void getPositions(std::vector<double>&, std::vector<double>&,
                    std::vector<double>&) const;
  /*! \brief Sets positions of all knots */
  void setPositions(const std::vector<float>&,const std::vector<float>&,
                    const std::vector<float>&); 
  /*! \brief Sets positions of all knots */
  void setPositions(const std::vector<double>&,const std::vector<double>&,
                    const std::vector<double>&);
  
  void getParams(paramSet&) const override; //!< Get parameters
  void setParams(const paramSet&) override; //!< Set parameters
 
  /*! \brief Evaluates Sigma at specified value of S_1 */
  double getSigma(double) const noexcept;
  /*! \brief Evaluates Offset at specified value of S_1 */
  double getOffset(double) const noexcept;
  
  /*! \brief Evaluates number counts model */
  double getNumberCounts(double, double) const noexcept override;

  /*! Evaluates band 1 number counts model at specified value of S_1 */
  double getBand1NumberCounts(double) const noexcept override;
  /*! Evaluates band 2 number counts model at specified value of S_2 */
  double getBand2NumberCounts(double) const noexcept override;

  /*! \brief Minimum flux model is defined for */
  dblpair getMinFlux() const noexcept override;
  
  /*! \brief Maxium flux model is defined for */
  dblpair getMaxFlux() const noexcept override;
  
  /*! \brief Get range over which R is expected to be nonzero */
  std::pair<dblpair, dblpair> getRRange(const doublebeam&) const 
    throw(affineExcept) override;

  double getNS() const noexcept override; //!< Total number of sources per area
  double getFluxPerArea(unsigned int) const noexcept override; //!< Flux per unit area (sq deg)
  double getFluxSqPerArea(unsigned int) const noexcept override; //!< Flux^2 per unit area (sq deg)
  double getFluxPowerPerArea(double p1, double p2) const noexcept; //!< Integral of s1 to some power and s2 to some power over number counts
  
  //Routines for getting R
  /*! \brief Get R value at single flux1,flux2 value */
  double getR(double, double, const doublebeam&) const override;
  /*! \brief Gets R values, array version */
  void getR(unsigned int n1, const double* const, unsigned int n2,
            const double* const, const doublebeam&, double*) const override;
  
  double differenceRegularize(double) const; //!< Tikhonov regularization log Likelihood penalty

  void writeToHDF5Handle(hid_t, bool=false) const override; //!< Write to HDF5 handle
  bool writeToStream(std::ostream& os) const override; //<! Output to stream
  void readFromHDF5Handle(hid_t) override; //!< Read from HDF5 Handle
  
  void sendSelf(MPI_Comm, int dest) const; //!< Send self
  void receiveCopy(MPI_Comm, int src); //!< Receive
};

#endif
