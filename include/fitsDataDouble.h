//fitsDataDouble.h

#ifndef __fitsDataDouble__
#define __fitsDataDouble__

#include<string>
#include<utility>

#include <mpi.h>

#include "../include/global_settings.h"

/*!
  \brief Data class for holding 2-band data from a fits files.
  Supports binning of data.
*/

class fitsDataDouble {
 private:
  unsigned int n; //!< Number of data points.  Must be the same in both
  double *data1; //!< Data, band 1
  double *data2; //!< Data, band 2

  //Binning variables
  bool is_binned; //!< Has data been binned
  unsigned int nbins1; //!< Number of bins, band 1
  double bincent01; //!< Center of bin 0, band 1
  double bindelta1; //!< Delta of bins, band 1
  unsigned int nbins2; //!< Number of bins, band 2
  double bincent02; //!< Center of bin 0, band 2
  double bindelta2; //!< Delta of bins, band 2
  unsigned int* binval; //!< Number of elements in each bin, row major array

  // We want to keep the filenames and extensions we read data from
  std::string file1; //!< File we read data 1 from
  unsigned int dataext1; //!< File extension for data 1
  bool has_mask1; //!< Did we have a mask from file 1?
  unsigned int maskext1; //!< Mask 1 extension
  std::string file2; //!< File we read data 2 from
  unsigned int dataext2; //!< File extension for data 2
  bool has_mask2; //!< Did we have a mask from file 2?
  unsigned int maskext2; //!< Mask 2 extension
  
  /*! \brief Reads data from single file*/
  bool readFile(const std::string&, long&, double*&, unsigned int&,
		unsigned int*&, unsigned int&, bool=false);

 public:
  fitsDataDouble(); //!< Default constructor
  fitsDataDouble(const std::string&, const std::string&, 
		 bool ignore_mask=false, 
		 bool meansub=false); //!< Constructor from file
  ~fitsDataDouble();

  void readData(const std::string&, const std::string&, bool=false, 
		bool=false); //!< Read data from file

  bool isBinned() const { return is_binned; } //!< Is data binned
  void applyBinning(unsigned int, unsigned int); //!< Takes an unbinned image and bins it
  void removeBinning(); //!< Removes binning
  std::pair<unsigned int, unsigned int> getNBins() const; //!< Get number of bins along each dimension
  dblpair getBinCent0() const; //!< Gets center of 0th bin along each dimension
  dblpair getBinDelta() const; //!< Gets size of bins along each dimension
  const unsigned int* getBinnedData() const { return binval; } //!< Access to binned elements
  unsigned int getBinnedData(unsigned int i, unsigned int j) const 
  { return binval[i*nbins2+j]; } //!< Unchecked binned data access
  unsigned int getBinnedData(unsigned int i) const { return binval[i]; } //!< Unckecked binned data access

  bool hasData() const { return n != 0; } //!< Has data been read?
  unsigned int getN() const { return n; } //!< Get number of data points
  const double* getData1() const { return data1; } //!< Direct data access, band 1
  const double* getData2() const { return data2; } //!< Direct data access, band 2
  double getData1(unsigned int i) const { return data1[i]; } //!< Unchecked
  double getData2(unsigned int i) const { return data2[i]; } //!< Unchecked
  
  dblpair meanSubtract(); //!< Subract the mean from the data. May require rebinning

  dblpair getMin() const; //!< Gets minimum in each band
  dblpair getMax() const; //!< Gets maximum in each band
  std::pair<dblpair, dblpair> getMinMax() const; //!< Get minima and maxima for both bands
  dblpair getMean() const; //!< Gets mean value in each band

  /*! \brief Get the filenames of the data */
  std::pair<std::string, std::string> getFiles() const 
    { return std::make_pair(file1, file2); }
  /*! \brief Get the data extensions */
  std::pair<unsigned int, unsigned int> getDataExt() const 
    { return std::make_pair(dataext1, dataext2); }
  /*! \brief Get whether masks were present */
  std::pair<bool, bool> hasMask() const 
    { return std::make_pair(has_mask1, has_mask2); }
  /*! \brief Get mask extensions */
  std::pair<unsigned int, unsigned int> getMaskExt() const 
    { return std::make_pair(maskext1, maskext2); }
  
  void sendSelf(MPI_Comm, int dest) const; //!< Send self
  void receiveCopy(MPI_Comm, int src); //!< Receive
};

#endif
