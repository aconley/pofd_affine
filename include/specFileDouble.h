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

  The valid values for identifier= and their associated values are:\n\n
   -dataset=      datafile1 datafile2 inst_sigma1 inst_sigma2 psffile1 \
                   psffile2 [like_norm]\n
                    datafile1 and datafile2 are the names of fits data files,
	              in band 1 and band 2, respectively
	            inst_sigma1, inst_sigma2 are the instrumental noise 
		      (Gaussian) in the same units as the data in datafile 
		      (e.g., Jy) for datafile1 and datafile2
	            psffile1, psffile2 are FITS files giving the beam in each 
                      band
	            like_norm is an optional input controlling how the
	              likelihood is normalized for the beam area. The likelihood
		      is divided by this number (def: 1.0).\n
   -bin_data=    nbins (int)\n
                  Turns on binning of data files with specified number 
                  of bins.\n
   -mean_sub=    bool\n
                  Do mean subtraction on input data if value is true\n
   -ignore_mask= bool\n
                  Ignore mask information in data files\n
   -fftsize=     nfft (int)\n
                  Sets fft transform length along each dimension (def: 4096)\n
   -edge_set=    bool\n
                  Turn on R edge integration if true (on by default, so this
                  is only useful for turning it off)\n
   -nedge=       value (int)\n
                  Set size of edge interpolation (def: 256).\n
   -minbeamval=  value (float)\n
                  Minimum beam value used.  Everything with an absolute
		  value <= than this is ignored.  This can be zero if you
		  are not histogramming, but with histogramming a non-zero
		  value is advisable.  You should look at your
		  beam to decide a reasonable cutoff (def: 1e-6)\n
   -beam_histogram= bool\n
                  Activates beam histogramming.  On by default.\n
   -nbeamhist= value (int)\n
                  If beam histogramming, number of bins in histogram (def: 150)
   -fit_sigma1=  bool\n
                  If value is true, this turns on fitting for sigma in band 1
   -fit_sigma2=  bool\n
                  If value is true, this turns on fitting for sigma in band 2
   -sigmaprior1=  stdev (float)\n
                  This activates the Gaussian prior on the instrumental noise
		  multiplier in band 1.  It has mean one and standard 
		  deviation stdev.  Will automatically set fit_sigma1\n
   -sigmaprior2=  stdev (float)\n
                  Same as sigmaprior1, but in band 2\n
   -exp_conf1=    value (float)\n
                  Sets the expected confusion noise value in Jy in band 1; 
		  used to set the expected likelihood values.\n
   -exp_conf2=    value (float)\n
                  Same as exp_conf1, but in the other band\n
   -cfirbprior1=  mean stdev (floats)\n
                  This activates a Gaussian prior on the integrated CFIRB 
		  luminoisty in band 1. It has a mean of mean and a standard 
                  deviation stdev expressed in the same flux and area units as 
		  the model.  So if your model dN/dS is in Jy^-1 deg^-2,
		  then your CFIRB prior should be in Jy deg^-2\n
   -cfirbprior2= mean stdev (floats)\n
                  Same as cfirbprior1 but in band 2\n
   -poissonprior1=  mean stdev (floats)\n
                  This activates a Gaussian prior on the Poisson noise in
		  band 1. It has a mean of mean and a standard deviation 
		  stdev expressed in the same flux and area units as 
		  the model.  So if your model dN/dS is in Jy^-1 deg^-2,
		  then your Poisson prior should be in Jy^2 deg^-2\n
   -poissonprior2= mean stdev (floats)\n
                  Same as poissonprior1 but in band 2\n
   -regularize_alpha= value (float)\n
                  Sets regularization multiplier for difference operator
		  penalty in log likelihood for band one model (def: 0).\n
   -wisdom_file= filename (string)\n
                  Uses filename as a FFTW wisdom file\n
   -verbose=     bool\n
                  Turns on verbose mode (verbosity = 1)\n
   ultraverbose= bool\n
                  Turns on ultraverbose mode (verbosity = 2)\n
   -verbosity=    value	(int)\n
                  Sets verbosity to this level.\n
   -seed=         value (int)\n
                  Sets random number generator seed\n

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
  double minbeamval; //!< Minimum beam value used
  bool beam_histogram; //!< Do beam histogramming
  unsigned int nbeamhist; //!< Number of beam histogram bins
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
  bool has_poissonprior1; //!< Is poisson prior on, band 1
  double poissonprior_mean1; //!< Mean value of prior, band 1
  double poissonprior_stdev1; //!< Stdev of prior, band 1
  bool has_poissonprior2; //!< Is poisson prior on, band 2
  double poissonprior_mean2; //!< Mean value of prior, band 2
  double poissonprior_stdev2; //!< Stdev of prior, band 2
  double regularization_alpha; //!< Regularization multiplier
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
