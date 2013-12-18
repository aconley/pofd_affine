//beam.h

//PSF (pixelated) with positive and negative parts

#ifndef __beam__
#define __beam__

#include<mpi.h>

#include <string>
#include <algorithm>
#include <limits>

#include "../include/utility.h"

/*!
  \brief Represents beam in 1D case.

  Note that zero beam values are not stored

  Can also store the histogrammed inverse beam.

  \ingroup Beams
*/
class beam {
 private :
  static const unsigned int histothresh; //!< Minimum number of pix to histogram
  unsigned int npos; //!< Number of positive pixels
  unsigned int nneg; //!< Number of negative pixels

  // Raw beam
  bool haspos; //!< Has positive pixels
  bool hasneg; //!< Has negative pixels
  double *pospixarr; //!< Holds positive pixel values (sorted)
  double *negpixarr; //!< Holds absolute value of negative pixel values (sorted)
  double *posinvpixarr; //!< Inverse positive pixels
  double *neginvpixarr; //!< Inverse negative pixels (abs value)

  // Histogrammed beam, if present
  unsigned int nbins; //!< Number of histogram bins
  bool is_pos_histogrammed; //!< Is the positive inverse beam histogrammed?
  bool is_neg_histogrammed; //!< Is the negative inverse beam histogrammed?
  unsigned int posnbins; //!< Number of non-zero weights in pos hist beam
  unsigned int negnbins; //!< Number of non-zero weights in neg hist beam
  double *posweights; //!< Positive weights; double to avoid casting
  double *negweights; //!< Negative weights; double to avoid casting
  double *poshistval; //!< Positive inverse histogrammed beam values
  double *neghistval; //!< Negative inverse histogrammed beam values (abs value)

  // Descriptive parameters
  double pixsize; //!< Size of pixel in arcsec (if oversampled, not real pix)
  double totpos; //!< Sum of all positive elements
  double totpossq; //!< Sum of positive elements squared
  double totneg; //!< Sum of all negative elements
  double totnegsq; //!< Sum of negative elements squared

  double minval; //!< Smallest value to keep (based on abs for neg beam)

  void cleanup(); //!< Frees internal structures
  bool revSort(const double&, const double&) const; //!< Reverse sort comparison

 public :
  beam(); //!< Default constructor
  beam(const std::string& filename, bool histogram=false, 
       unsigned int NBINS=120, double MINVAL=1e-5); //!< Reads beam from file
  beam(const beam&); //!< Copy constructor
  ~beam() { cleanup(); } //!< Destructor

  void free() { cleanup(); } //!< Free all memory

  void readFile(const std::string& filename, double MINVAL=1e-5); //!< Read in file
  void makeHistogram(unsigned int NBINS=120); //!< Prepare the histogram

  double getEffectiveArea() const; //!< Get effective area of beam in sq deg
  double getEffectiveAreaPos() const; //!< Get effective area of positive beam in sq deg
  double getEffectiveAreaNeg() const; //!< Get effective area of negative beam in sq deg
  double getEffectiveAreaPix() const; //!< Get effective area of beam in pixels
  double getEffectiveAreaSq() const; //!< Get effective area of squared beam in sq deg

  beam& operator=(const beam&); //!< Copy

  unsigned int getNPos() const { return npos; } //!< Number of positive beam pix
  unsigned int getNNeg() const { return nneg; } //!< Number of negative beam pix
  bool hasPos() const { return haspos; } //!< Beam has positive pixels
  bool hasNeg() const { return hasneg; } //!< Beam has negative pixels
  bool isPosHist() const { return is_pos_histogrammed; } //!< Inverse positive hstogrammed beam is available
  bool isNegHist() const { return is_pos_histogrammed; } //!< Inverse negative hstogrammed beam is available
  unsigned int getNHistPos() const { return posnbins; } //!< Number of positive hist bins
  unsigned int getNHistNeg() const { return negnbins; } //!< Number of negative hist bins

  double getPixSize() const { return pixsize; } //!< Pixel size (1d)
  double getPos(unsigned int i) const { return pospixarr[i]; } //!< Get positive pixel
  double& getPos(unsigned int i) { return pospixarr[i]; } //!< Get positive pixel
  double getNeg(unsigned int i) const { return negpixarr[i]; } //!< Get negative pixel
  double& getNeg(unsigned int i) { return negpixarr[i]; } //!< Get negative pixel
  
  // Access to un-histogrammed beam
  /*! \brief Access inverse positive pixel array */
  const double* const getPosInvPixArr() const { return posinvpixarr; }
  /*! \brief Access inverse negative pixel array */
  const double* const getNegInvPixArr() const { return neginvpixarr; }

  // Access to histogrammed beam
  /*! \brief Access positive pixel weight array */
  const double* const getPosHistWeights() const { return posweights; }
  /*! \brief Access negative pixel weight array */
  const double* const getNegHistWeights() const { return negweights; }
  const double* const getPosHist() const { return poshistval; }
  const double* const getNegHist() const { return neghistval; }
  
  double getMinPos() const; //!< Minimum positive pixel
  double getMaxPos() const; //!< Maximum positive pixel
  double getMinAbsNeg() const; //!< Minimum absolute value negative pixel
  double getMaxAbsNeg() const; //!< Maximum absolute value negative pixel
  
  std::pair<double,double> getRangePos() const; //!< Get range of positive pixels
  std::pair<double,double> getRangeNeg() const; //!< Get range of negative pixels

  /*! \brief Returns positive pixel values raised to some power */
  void powerPos(double, double*) const;
  /*! \brief Returns negative pixel values raised to some power */
  void powerNeg(double, double*) const;

  /*! \brief Returns first index greater than specified value in pos pixel map,
   or npos if there isn't one*/
  unsigned int idxFirstGtValPos(double val) const {
    if (!haspos) return std::numeric_limits<double>::quiet_NaN();
    return utility::binary_search_gt(val, pospixarr, npos);
  }
  /*! \brief Returns first index with abs greater than specified value 
    in neg pixel map, or nneg if there isn't one*/
  unsigned int idxFirstGtValNeg(double val) const {
    if (! hasneg) return std::numeric_limits<double>::quiet_NaN();
    return utility::binary_search_gt(val, negpixarr, nneg);
  }

  /*! \brief Returns last index less than specified value in pos pixel map,
   or npos if there isn't one*/
  unsigned int idxLastLtValPos(double val) const {
    if (!haspos) return std::numeric_limits<double>::quiet_NaN();
    return utility::binary_search_lt(val, pospixarr, npos);
  }
  /*! \brief Returns last index with abs less than specified value in 
    negative pixel map, or nneg if there isn't one*/
  unsigned int idxLastLtValNeg(double val) const {
    if (! hasneg) return std::numeric_limits<double>::quiet_NaN();
    return utility::binary_search_lt(val, negpixarr, nneg);
  }

  /*! \brief Returns last index less than or equal to specified value in 
    positive pixel map, or npos if there isn't one*/
  unsigned int idxLastLteValPos(double val) const {
    if (!haspos) return std::numeric_limits<double>::quiet_NaN();
    return utility::binary_search_lte(val, pospixarr, npos);
  }

  /*! \brief Returns last index with abs less than or equal to specified 
    value in negative pixel map, or nneg if there isn't one*/
  unsigned int idxLastLteValNeg(double val) const {
    if (! hasneg) return std::numeric_limits<double>::quiet_NaN();
    return utility::binary_search_lte(val, negpixarr, nneg);
  }

  /*! \brief MPI copy send operation */
  void sendSelf(MPI_Comm, int dest) const;
  /*! \brief MPI copy recieve operation */
  void recieveCopy(MPI_Comm, int dest);
};

#endif
