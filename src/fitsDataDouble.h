//fitsDataDouble.h

#ifndef __fitsDataDouble__
#define __fitsDataDouble__

#include<string>
#include<utility>

#include <mpi.h>

/*!
  \brief Data class for holding 2-band data from a fits files.
  Supports binning.
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

  /*! \brief Reads data from single file*/
  bool readFile(const std::string&,unsigned int&,double*&, 
		int*&,bool=false);


 public:
  fitsDataDouble(); //!< Default constructor
  fitsDataDouble(const std::string&, const std::string&, 
		 bool ignore_mask=false, 
		 bool meansub=false); //!< Constructor from file
  ~fitsDataDouble();

  void readData(const std::string&, const std::string&, bool=false, 
		bool=false); //!< Read data from file

  bool isBinned() const { return is_binned; }
  void applyBinning(unsigned int, unsigned int); //!< Takes an unbinned image and bins it
  std::pair<unsigned int, unsigned int> getNBins() const;
  std::pair<double,double> getBinCent0() const;
  std::pair<double,double> getBinDelta() const;
  const unsigned int* const getBinnedData() const { return binval; }
  unsigned int getBinnedData(unsigned int i, unsigned int j) const 
  { return binval[i*nbins2+j]; } //!< Unckecked
  unsigned int getBinnedData(unsigned int i) const { return binval[i]; } //!< Unckecked

  unsigned int getN() const { return n; } //!< Get number of data points
  const double* const getData1() const { return data1; }
  const double* const getData2() const { return data2; }
  double getData1(unsigned int i) const { return data1[i]; } //!< Unchecked
  double getData2(unsigned int i) const { return data2[i]; } //!< Unchecked
  
  std::pair<double,double> meanSubtract(); //!< Subract the mean from the data. May require rebinning

  std::pair<double,double> getMin() const;
  std::pair<double,double> getMax() const;
  void getMinMax(double&, double&, double&, double&) const; //!< Get minima and maxima
  std::pair<double,double> getMean() const;

  void SendSelf(MPI::Comm&, int dest) const; //!< Send self
  void RecieveCopy(MPI::Comm&, int src); //!< Recieve
};

#endif
