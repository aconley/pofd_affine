//fitsData.h

#ifndef __fitsData__
#define __fitsData__

#include <string>
#include <mpi.h>

#include "../include/global_settings.h"

/*!
  \brief Data class for holding 1-band data from a fits file.
  Supports binning.
*/

class fitsData {
 private:
  unsigned int n; //!< Number of data points
  double *data; //!< Data

  //Binning variables
  bool is_binned; //!< Has data been binned
  unsigned int nbins; //!< Number of bins
  double bincent0; //!< Center of bin 0
  double bindelta; //!< Delta of bins
  unsigned int* binval; //!< Number of elements in each bin

 public:
  fitsData(); //!< Default constructor
  fitsData(const std::string&, bool ignore_mask=false, 
	   bool meansub=false); //!< Constructor from file
  ~fitsData();

  void readData(const std::string&, bool ignore_mask=false, 
		bool domeansub=false); //!< Read data from file

  bool isBinned() const { return is_binned; } //!< Is data binned?
  void applyBinning(unsigned int); //!< Takes an unbinned image and bins it
  void removeBinning(); //!< Removes binning
  unsigned int getNBins() const { return nbins; } //!< Gets number of bins
  double getBinCent0() const { return bincent0; } //!< Gets center of 0th bin
  double getBinDelta() const { return bindelta; } //!< Gets size of bins

  bool hasData() const { return n != 0; } //!< Has data been read
  unsigned int getN() const { return n; } //!< Number of data points
  const double* const getData() const { return data; } //!< Access to data
  /*! \brief Access to binned data */
  const unsigned int* const getBinnedData() const { return binval; }

  double getData(unsigned int i) const { return data[i]; } //!< Unchecked data access
  unsigned int getBinnedData(unsigned int i) const { return binval[i]; } //!< Unckecked binned data access
  
  double meanSubtract(); //!< Subract the mean from the data. May require rebinning

  double getMax() const; //!< Get maximum flux
  double getMin() const; //!< Get minimum flux
  dblpair getMinMax() const; //!< Get minima and maxima
  double getMean() const; //!< Get mean flux

  void sendSelf(MPI_Comm, int dest) const; //!< Send self
  void recieveCopy(MPI_Comm, int src); //!< Recieve
};

#endif
