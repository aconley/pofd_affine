//binneddata

#ifndef __PD__
#define __PD__

#include<string>
#include<ostream>

#include<hdf5.h>

#include "../include/fitsData.h"

/*!
  \brief Class to hold P(D), 1D case.  Supports interpolation.

  Only monotonic flux value grids are allowed.

  By default, the log of the P(D) is stored (and interpolated on),
  since this is what you want for the likelihood calculation.
  The user can specify the real-space P(D).
  The memory usage is something like std::vector -- each object
  has a capacity and a current size, with size < capacity.  This
  tries to avoid resizing if the memory is already allocated.
*/
class PD {
private:
  static constexpr double lowsigval = 2.5; //!< Constant used in edgeFix

  unsigned int n; //!< Current size
  unsigned int capacity; //!< Current capacity

  double getLogLikeBinned(const fitsData&) const;
  double getLogLikeUnbinned(const fitsData&) const;

public:
  /* \brief Constructor */
  explicit PD(unsigned int N=0, double MINFLUX=0.0, double DFLUX=1.0);
  PD(const PD&)=delete;
  PD(PD&&); //!< Move constructor
  ~PD(); //!< Destructor

  //Public for efficient filling -- bad form, but speed matters here
  bool logflat; //!< True if log( P(D) ) is stored instead of P(D)
  double minflux; //!< Minimum flux
  double dflux; //!< Flux step along axis
  double* pd_; //!< Actual P(D)

  void shrink(); //!< Shrink memory requirements to user size
  void strict_resize(unsigned int); //!< Resize data arrays, forcing actual resizing in all cases
  void resize(unsigned int); //!< Resize data arrays

  double getTotal() const; //!< Get sum of entries
  double getIntegral() const; //!< Get integral of entries
  void normalize(); //!< Normalize integral of model
  
  bool isLog() const { return logflat; } //!< Is log(P(D)) stored?

  void applyLog(bool=false); //!< Logify if not already
  void deLog(); //!< de-logify if not already

  void edgeFix(bool donorm=true); //!< Apply Gaussian replacement to bottom edges

  double getMean(bool donorm=true) const; //!< Get mean 

  dblpair getMeanAndVar(bool donorm=true) const; //!< Get mean and variance 
  
  PD& operator=(const PD&); //!< Copy
  PD& operator=(PD&&); //!< Move assignment

  /*! \brief Fill contents from array*/
  void fill(unsigned int N, double MIN, double DF,
            const double* const DATA, bool LOG=true); 

  /*! \brief Get flux value for specified index */
  double getFluxVal(unsigned int i) const { return minflux+static_cast<double>(i)*dflux; }
  /*! \brief Get PD value for specified index */
  double getPDVal(unsigned int i) const { return pd_[i]; }
  /*! \brief Get PD value for specified flux value (with interpolation) */
  double getPDVal(double) const;

  /*! \brief Element access */
  double operator[](unsigned int i) const { return pd_[i]; }
  unsigned int getDim() const { return n; } //!< Get number of elements in PD

  /*! \brief Get Log likelihood of data set*/
  double getLogLike(const fitsData&) const;

  std::ostream& writeToStream(std::ostream& os) const; //!< Write summary

  int writeToFits(const std::string& file) const; //!< Write as fits file
  void writeToHDF5Handle(hid_t) const; //!< Write as HDF5 handle
  void writeToHDF5(const std::string& file) const; //!< Write as HDF5 file
};

/*!\brief Write to stream */
std::ostream& operator<<(std::ostream& os, const PD&);

#endif
