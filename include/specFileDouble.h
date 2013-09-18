//specFileDouble.h

#ifndef __specFileDouble__
#define __specFileDouble__

#include<string>
#include<vector>

#include "hdf5.h"

/*!
  \brief Structure to read in spec files for fits, 2D case

  Lines in the spec file control the likelihood calculations and fits.
  The format is: identifier= values
  The appropriate entries for values depends upon identifier=

  The valid values for identifier= and their associated values are:
   dataset=      datafile1 datafile2 inst_sigma1 inst_sigma2 psffile1 \
                   psffile2 [like_norm]
                    datafile1 and datafile2 are the names of fits data files,
	              in band 1 and band 2, respectively
	            inst_sigma1, inst_sigma2 are the instrumental noise 
		      (Gaussian) in the same units as the data in datafile 
		      (e.g., Jy) for datafile1 and datafile2
	            psffile1, psffile2 are FITS files giving the beam in each 
                      band
	            like_norm is an optional input controlling how the
	              likelihood is normalized for the beam area. The likelihood
		      is divided by this number (def: 1.0).
   bin_data=    nbins
                  Turns on binning of data files with specified number 
                  of bins.
   mean_sub=    bool
                  Do mean subtraction on input data if value is true
   ignore_mask= bool
                  Ignore mask information in data files
   fftsize=     nfft
                  Sets fft transform length along each dimension (def: 4096)
   edge_set=    bool
                  Turn on R edge integration if true (on by default, so this
                  is only useful for turning it off)
   nedge=       value
                  Set size of edge interpolation (def: 256).
   beam_histogram= bool
                  Activates beam histogramming.  On by default.
   hist_logstep= value
                  If beam histogramming, log bin size of histogram (def: 0.2)  
   edge_fix=    bool
                  Apply edge fix to P(D). On by default.
   fit_sigma1=  value
                  If value is true, this turns on fitting for sigma in band 1
   fit_sigma2=  value
                  If value is true, this turns on fitting for sigma in band 2
   sigmaprior1=  stdev
                  This activates the Gaussian prior on the instrumental noise
		  multiplier in band 1.  It has mean one and standard 
		  deviation stdev.  Will automatically set fit_sigma1
   sigmaprior2=  stdev
                  Same as sigmaprior1, but in band 2
   exp_conf1=    value
                  Sets the expected confusion noise value in Jy in band 1; 
		  used to set the expected likelihood values.
   exp_conf2=    value
                  Same as exp_conf1, but in the other band
   cfirbprior1=  mean stdev
                  This activates a Gaussian prior on the integrated CFIRB 
		  luminoisty in band 1. It has a mean of mean and a standard 
                  deviation stdev expressed in the same flux and area units as 
		  the model.  So if your model dN/dS is in Jy^-1 deg^-2,
		  then your CFIRB prior should be in Jy deg^-2
   cfirbprior2= mean stdev
                  Same as cfirbprior1 but in band 2
   wisdom_file= filename
                  Uses filename as a FFTW wisdom file		  
   verbose=     bool
                  Turns on verbose mode (verbosity = 1)
   ultraverbose= bool
                  Turns on ultraverbose mode (verbosity = 2)
   verbosity=    value		  
                  Sets verbosity to this level.
   seed=         value
                  Sets random number generator seed
  Having multiple lines with dataset= is additive.
  Having multiple lines of the others results in only
  the values from the final such line being used.
*/
struct specFileDouble {
  //dataset
  std::vector<std::string> datafiles1; //!< Data files, band 1
  std::vector<std::string> datafiles2; //!< Data files, band 2
  std::vector<double> sigmas1; //!< Instrumental noise values, band 1
  std::vector<double> sigmas2; //!< Instrumental noise values, band 2
  std::vector<std::string> psffiles1; //!< Beam file, band 1
  std::vector<std::string> psffiles2; //!< Beam file, band 2
  std::vector<double> like_norm; //!< Likelihood normalization

  //Various fit options
  bool bin_data; //!< Bin input data
  unsigned int nbins; //!< Number of data bins
  bool ignore_mask; //!< Ignore mask information in data files
  bool mean_sub; //!< Mean subtract data
  unsigned int fftsize; //!< Length of fft
  bool edge_set; //!< Do edge setting integration
  unsigned int nedge; //!< Number of edge integration points
  bool edge_fix; //!< Apply edge fix to P(D)
  bool beam_histogram; //!< Do beam histogramming
  double hist_logstep; //!< Size of Beam histogram log step
  bool fit_sigma1; //!< Do fit to sigma in band 1
  bool fit_sigma2; //!< Do fit to sigma in band 2
  bool has_sigprior1; //!< Is the sigma prior on, band 1
  double sigprior_stdev1; //!< Stdev of sigma prior, band 1
  bool has_sigprior2; //!< Is the sigma prior on, band 2
  double sigprior_stdev2; //!< Stdev of sigma prior, band 2
  double exp_conf1; //!< Expected confusion noise, band 1
  double exp_conf2; //!< Expected confusion noise, band 2
  bool has_cfirbprior1; //!< Is cfirb prior on, band 1
  double cfirbprior_mean1; //!< Mean value of prior, band 1
  double cfirbprior_stdev1; //!< Stdev of prior, band 1
  bool has_cfirbprior2; //!< Is cfirb prior on, band 2
  double cfirbprior_mean2; //!< Mean value of prior, band 2
  double cfirbprior_stdev2; //!< Stdev of prior, band 2
  bool has_wisdom_file; //!< Has FFTW wisdom file
  std::string wisdom_file; //!< Name of wisdom file
  unsigned int verbosity; //!< Verbosity level
  bool has_user_seed; //!< User set seed
  unsigned long long int seed; //!< RNG seed

  specFileDouble(); //!< Default constructor
  specFileDouble(const std::string&); //!< Constructor with file read

  void init(); //!< Reset values

  void writeToHDF5Handle(hid_t) const; //!< Write to HDF5 Handle

  void readFile(const std::string&); //!< Read in file
};

#endif
