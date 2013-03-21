//doublebeam.h

//2-band PSF (pixelated) with positive and negative parts

#ifndef __doublebeam__
#define __doublebeam__

#include<mpi.h>

#include <string>
#include <algorithm>
#include <limits>

#include "../include/utility.h"

/*!
  \brief Represents PSF parameters for 2 beams with positive
          and negative components.  Note that zero beam values
	  are not stored.  The two beams must have the same
	  pixel size and extent, and are assumed co-aligned.

  \ingroup Beams
*/
class doublebeam {
 private :
  double pixsize; //!< Size of pixel in arcsec (if oversampled, not real pix)

  //There are 4 combinations to worry about:
  // positive+positive, positive+negative, negative+positive, negative+negative
  //They will always be stored in that order
  bool hassign[4]; //!< Do we have pp,pn,np,nn pixels
  unsigned int npix[4]; //!< Number of pp, pn, np, nn pixels

  static const unsigned int histothresh; //!< Don't bother histogramming for this many or fewer

  bool has_weights[4]; //!< Using histogrammed weights? (pp,pn,np,nn)
  double *weights[4]; //!< Histogram weights, if has weights (pp,pn,np,nn); double to avoid casting
  double *pixarr1[4]; //!< Array of pixels, beam 1 (pp,pn,np,nn)
  double *pixarr2[4]; //!< Array of pixels, beam 2 (pp,pn,np,nn)
  double *ipixarr1[4]; //!< Array of inverse pixels, beam 1 (pp,pn,np,nn)
  double *ipixarr2[4]; //!< Array of inverse pixels, beam 2 (pp,pn,np,nn)

  double tot1[4]; //!< Sum of pp,pn,np,nn elements, beam 1
  double tot2[4]; //!< Sum of pp,pn,np,nn elements, beam 2
  double totsq1[4]; //!< Sum of pp,pn,np,nn squared elements, beam 1
  double totsq2[4]; //!< Sum of pp,pn,np,nn squared elements, beam 2
  double totsm1; //!< Sum of tot1
  double totsm2; //!< Sum of tot2

  void cleanup(); //!< Frees internal structures

 public :
  doublebeam(); //!< Default constructor
  doublebeam(const std::string&, const std::string&,
	     bool histogram=false, double histogramlogstep=0.2); //!< Reads beam from files
  doublebeam(const doublebeam&); //!< Copy constructor
  ~doublebeam() { cleanup(); } //!< Destructor

  void free() { cleanup(); } //!< Free all memory

  /*! \brief Read in files */
  void readFiles(const std::string&, const std::string& filename2,
		 bool histogram=false, double histogramlogstep=0.2); 
  
  void setBeams(unsigned int, const double* const,
		const double* const, double, bool histogram=false,
		double histogramlogstep=0.2 ); //!< Set beams by hand

  double getEffectiveArea1() const; //!< Get effective area of beam1 in sq deg
  double getEffectiveAreaSign1(unsigned int) const; //!< Get effective area of either pp,pn,np, or nn beam 1, in pixels
  double getEffectiveAreaSqSign1(unsigned int) const; //!< Get effective squared area of either pp,pn,np, or nn, beam 1
  double getEffectiveArea2() const; //!< Get effective area of beam2 in sq deg
  double getEffectiveAreaSign2(unsigned int) const; //!< Get effective area of either pp,pn,np, or nn beam 2, in pixels
  double getEffectiveAreaSqSign2(unsigned int) const; //!< Get effective squared area of either pp,pn,np, or nn, beam 2

  double getEffectiveAreaPixGeoMean() const; //!< Get the geometric mean of the effective areas in pixels
  
  double getMinAreaPix() const; //!< Get area of minimum area beam

  doublebeam& operator=(const doublebeam&); //!< Copy

  unsigned int getNPix(unsigned int i) const { return npix[i]; } //!< Number of pixels in pp,pn,np,nn beams
  unsigned int getTotalNPix() const; //!< Total number of pix
  bool hasSign(unsigned int i) const { return hassign[i]; } //!< Do beams have certain sign combinations? (pp,pn,np,nn)
  bool hasWeights(unsigned int i) const { return has_weights[i]; } //!< Has beam weights in sign combination (pp,pn,np,nn)

  /*! \brief Get max values for pp,pn,np,nn pieces of beam 1*/
  double getMax1(unsigned int) const;
  /*! \brief Get max values for pp,pn,np,nn pieces of beam 2*/
  double getMax2(unsigned int) const;

  /*! \brief Get min/max values for pp,pn,np,nn pieces of beam 1*/
  void getMinMax1(unsigned int, double&, double&) const;
  /*! \brief Get min/max values for pp,pn,np,nn pieces of beam 2*/
  void getMinMax2(unsigned int, double&, double&) const;

  /*! \brief Get pixel array element, band 1 */
  const double* const getPixArr1(unsigned int idx) const { return pixarr1[idx]; }
  /*! \brief Get pixel array element, band 2 */
  const double* const getPixArr2(unsigned int idx) const { return pixarr2[idx]; }

  /*! \brief Get inverse pixel array element, band 1 */
  const double* const getInvPixArr1(unsigned int idx) const { return ipixarr1[idx]; }
  /*! \brief Get inverse pixel array element, band 2 */
  const double* const getInvPixArr2(unsigned int idx) const { return ipixarr2[idx]; }

  /*! \brief Get weights element */
  const double* const getWeights(unsigned int idx) const { return weights[idx]; }

  double getPixSize() const { return pixsize; } //!< Get pixel size (1d)

  /*! \brief MPI copy send operation */
  void sendSelf(MPI::Comm&, int dest) const;
  /*! \brief MPI copy recieve operation */
  void recieveCopy(MPI::Comm&, int dest);
};

/*! \brief Structure for beam histogramming */
struct bmhist {
public:
  //For some reason averaging over 1/beam works badly
  unsigned int cnt; //!< Number of elements in bin
  double tot1; //!< Total in band 1
  double tot2; //!< Total in band 2
  bmhist() { cnt = 0; tot1=tot2=0.0; } //!< Constructor
  /*! \brief Copy constructor */
  bmhist(const bmhist& other) {
    cnt = other.cnt; tot1=other.tot1; tot2=other.tot2;
  }
  /*! \brief Copy operator */
  bmhist& operator=(const bmhist& other) {
    if (this == &other) return *this;
    cnt = other.cnt; tot1=other.tot1; tot2=other.tot2;
    return *this;
  }
    
};

#endif
