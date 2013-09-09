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
  \brief Represents PSF parameters for beam.

  Note that zero beam values are not stored

  \ingroup Beams
*/
class beam {
 private :
  static const unsigned int histothresh; //!< Minimum number of pix to histogram
  unsigned int npos; //!< Number of positive pixels
  unsigned int nneg; //!< Number of negative pixels
  bool haspos; //!< Has positive pixels
  bool hasneg; //!< Has negative pixels
  bool hasposweights; //!< Has weights for positive pixels (histogrammed)
  bool hasnegweights; //!< Has weights for negative pixels (histogrammed)
  double *posweights; //!< Positive weights; double to avoid casting
  double *negweights; //!< Negative weights; double to avoid casting
  double *pospixarr; //!< Holds positive pixel values (sorted)
  double *negpixarr; //!< Holds absolute value of negative pixel values (sorted)
  double *posinvpixarr; //!< Inverse positive pixels
  double *neginvpixarr; //!< Inverse negative pixels
  double pixsize; //!< Size of pixel in arcsec (if oversampled, not real pix)
  double totpos; //!< Sum of all positive elements
  double totpossq; //!< Sum of positive elements squared
  double totneg; //!< Sum of all negative elements
  double totnegsq; //!< Sum of negative elements squared

  void cleanup(); //!< Frees internal structures
  bool revSort(const double&, const double&) const; //!< Reverse sort comparison

 public :
  beam(); //!< Default constructor
  beam(const std::string& filename, bool histogram=false, 
	double histogramlogstep=0.2); //!< Reads beam from file
  beam(const beam&); //!< Copy constructor
  ~beam() { cleanup(); } //!< Destructor

  void free() { cleanup(); } //!< Free all memory

  void readFile(const std::string& filename, bool histogram=false, 
		double histogramlogstep=0.2); //!< Read in file

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
  bool hasPosWeights() const { return hasposweights; } //!< Has positive pixel weights
  bool hasNegWeights() const { return hasnegweights; } //!< Has negative pixel weights

  double getPixSize() const { return pixsize; } //!< Pixel size (1d)
  double getPos(unsigned int i) const { return pospixarr[i]; } //!< Get positive pixel
  double& getPos(unsigned int i) { return pospixarr[i]; } //!< Get positive pixel
  double getNeg(unsigned int i) const { return negpixarr[i]; } //!< Get negative pixel
  double& getNeg(unsigned int i) { return negpixarr[i]; } //!< Get negative pixel

  /*! \brief Access inverse positive pixel array */
  const double* const getPosInvPixArr() const { return posinvpixarr; }
  /*! \brief Access inverse negative pixel array */
  const double* const getNegInvPixArr() const { return neginvpixarr; }
  /*! \brief Access positive pixel weight array */
  const double* const getPosWeights() const { return posweights; }
  /*! \brief Access negative pixel weight array */
  const double* const getNegWeights() const { return negweights; }

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
