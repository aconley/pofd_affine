//doublebeam.h

//2-band PSF (pixelated) with positive and negative parts

#ifndef __doublebeam__
#define __doublebeam__

#include<mpi.h>

#include<string>
#include<utility>

#include "../include/utility.h"
#include "../include/global_settings.h"

/*!
  \brief Represents PSF parameters for 2 beams.  

  Note that zero beam values are not stored.  The two beams must have the same
  pixel size and extent, and are assumed co-aligned.

  This can also store the histogram of the inverse beam.

  There are four sign components to keep track of.  These are always
  in the order pos-pos, pos-neg, neg-pos, neg-neg

  It will also store the logarithm of the beam ratio of the two
  bands.  This is created the first time it is asked for, so is only
  computed if the user needs it.

  \ingroup Beams
*/
class doublebeam {
 private :
  static const unsigned int histothresh; //!< Don't bother histogramming for this many or fewer

  double pixsize; //!< Size of pixel in arcsec

  // Raw beam
  bool hassign[4]; //!< Do we have any elements of given sign component?
  unsigned int npix[4];  //!< Number of pp, pn, np, nn pixels
  double *pixarr1[4]; //!< Array of pixels, beam 1.  Length npix
  double *pixarr2[4]; //!< Array of pixels, beam 2.  Length npix
  double *invpixarr1[4]; //!< Array of inverse pixels, beam 1.  Length npix
  double *invpixarr2[4]; //!< Array of inverse pixels, beam 2.  Length npix
  
  // Log of beam ratio, unbinned
  mutable bool haslogratio[4]; //!< Is the log of the beam ratios computed?
  mutable double *logratio[4]; //!< Array of log(beam1/beam2), length npix

  // Histogrammed inverse beam
  // We build the nbins x nbins histogram, but only keep the non-zero entries
  unsigned int nbins; //!< Number of histogram bins
  bool ishistogrammed[4]; //!< Is given histogram sign available?
  unsigned int nhist[4]; //!< Number of bins actually stored for each histogram
  double *binweights[4]; //!< Weights (double to avoid casting expense)
  double *binvals1[4]; //!< Bin values, inverse beam, band 1
  double *binvals2[4]; //!< Bin values, inverse beam, band 2

  // Histogrammed log beam ratio
  mutable bool hasbinlogratio[4]; //!< Is the log of the histogrammed beam ratios there?
  mutable double *binlogratio[4]; //!< Array of log(beam1/beam2), histogrammed

  // Descriptive parameters. 
  double tot1[4]; //!< Sum of pp,pn,np,nn elements, beam 1
  double tot2[4]; //!< Sum of pp,pn,np,nn elements, beam 2
  double totsq1[4]; //!< Sum of pp,pn,np,nn squared elements, beam 1
  double totsq2[4]; //!< Sum of pp,pn,np,nn squared elements, beam 2
  double minbm1[4]; //!< Minimum value of beam1
  double maxbm1[4]; //!< Maximum value of beam1
  double minbm2[4]; //!< Minimum value of beam2
  double maxbm2[4]; //!< Maximum value of beam2

  double minval; //!< Minimum absolute value kept for all components

  void cleanup(); //!< Frees internal structures

 public :
  doublebeam(); //!< Default constructor
  doublebeam(const std::string&, const std::string&, bool histogram=false, 
	     unsigned int NBINS=150, double MINVAL=1e-6); //!< Reads beam from files
  ~doublebeam() { cleanup(); } //!< Destructor

  void free() { cleanup(); } //!< Free all memory

  /*! \brief Read in files */
  void readFiles(const std::string& filename1, const std::string& filename2,
		 double MINVAL=1e-6); 

  /*! \brief Sets beams from arrays */
  void setBeams(unsigned int n, const double* const beam1,
		const double* const beam2, double PIXSIZE, double MINVAL=1e-6);

  /*! \brief Build histograms */
  void makeHistogram(unsigned int NBINS=150); //!< Prepare the histogram

  bool hasData() const; //!< Is there any data set?

  double getEffectiveArea1() const; //!< Get effective area of beam1 in sq deg
  double getEffectiveAreaSign1(unsigned int) const; //!< Get effective area of either pp, pn, np, or nn beam 1, in pixels
  double getEffectiveAreaSqSign1(unsigned int) const; //!< Get effective squared area of either pp, pn, np, or nn, beam 1
  double getEffectiveArea2() const; //!< Get effective area of beam2 in sq deg
  double getEffectiveAreaSign2(unsigned int) const; //!< Get effective area of either pp, pn, np, or nn beam 2, in pixels
  double getEffectiveAreaSqSign2(unsigned int) const; //!< Get effective squared area of either pp, pn, np, or nn, beam 2
  
  bool hasSign(unsigned int i) const { return hassign[i]; } //!< Does beam have given sign component?
  unsigned int getNPix(unsigned int i) const { return npix[i]; } //!< Number of pixels in pp, pn, np, nn beams
  unsigned int getMaxNPix() const; //!< Largest number of pix
  unsigned int getTotalNPix() const; //!< Total number of pix
  double getMinval() const { return minval; } //!< Get minimum abs value
  double getPixSize() const { return pixsize; } //!< Return pixel size

  /*! \brief Get max values for pp,pn,np,nn pieces of beam 1*/
  dblpair getMinMax1(unsigned int) const;
  /*! \brief Get max values for pp,pn,np,nn pieces of beam 2*/
  dblpair getMinMax2(unsigned int) const;

  // Direct access to inverse pixel array.  Bad style, but important
  // for efficiency reasons
  /*! \brief Get inverse pixel array element, band 1 */
  const double* const getInvPixArr1(unsigned int idx) const { return invpixarr1[idx]; }
  /*! \brief Get inverse pixel array element, band 2 */
  const double* const getInvPixArr2(unsigned int idx) const { return invpixarr2[idx]; }

  // Access to log ratio of beams.  Computes it if needed
  const double* const getLogRatio(unsigned int idx) const;

  // Histogram information
  unsigned int getNBins() const { return nbins; } //!< Get number of bins
  /*! \brief Is the beam histogrammed? */
  bool isHistogrammed(unsigned int i) const { return ishistogrammed[i]; }
  /*! \brief Get number of histogram bins per sign component */
  unsigned int getNHist(unsigned int i) const { return nhist[i]; }
  
  // More direct acess 
  /*! \brief Access to histogram bin weights */
  const double* const getBinWeights(unsigned int i) const { return binweights[i]; }
  /*! \brief Access to histogram bin values, band 1 */
  const double* const getBinVals1(unsigned int i) const { return binvals1[i]; }
  /*! \brief Access to histogram bin values, band 2 */
  const double* const getBinVals2(unsigned int i) const { return binvals2[i]; }
  
  // Access to log ratio of histogrammed beams.  Computes it if needed
  const double* const getBinLogRatio(unsigned int idx) const;

  /*! \brief MPI copy send operation */
  void sendSelf(MPI_Comm, int dest) const;
  /*! \brief MPI copy receive operation */
  void receiveCopy(MPI_Comm, int dest);
};

#endif
