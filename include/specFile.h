//specFile.h

#ifndef __specFile__
#define __specFile__

#include<string>
#include<vector>

#include "hdf5.h"

/*!
  \brief Structure to read in spec files for fits, 1D case

  Lines in the spec file control the likelihood calculations and fits.\n
  The format is: identifier= values\n
  The appropriate entries for values depends upon identifier=\n
  Where values is bool the option is turned on if value is either
  true or yes, and is not activated for any other value.\n


  The valid values for identifier= and their associated values are:\n\n

  -dataset=     datafile inst_sigma psffile [like_norm]\n
  datafile is the name of a fits data file
  inst_sigma is the instrumental noise (Gaussian) in the same
  units as the data in datafile (e.g., Jy)
  psffile is a FITS file giving the beam
  like_norm is an optional input controlling how the
  likelihood is normalized for the beam area.  The likelihood
  is divided by this number (def: 1.0).\n
  -bin_data=    nbins (int)\n
  Turns on binning of data files with specified number 
  of bins\n
  -mean_sub=    bool\n
  Do mean subtraction on input data if value is true\n
  -ignore_mask= bool\n
  Ignore mask information in input files\n
  -fftsize=     nfft (int)\n
  Sets fft transform length (def: 131072)\n
  -ninterp=     value (int)\n
  Number of R interpolation points (def: 2048)\n
  -minbeamval=  value (float)\n
  Minimum beam value used.  Everything with an absolute
  value <= than this is ignored.  This can be zero if you
  are not histogramming, but with histogramming a non-zero
  value is advisable.  You should look at your
  beam to decide a reasonable cutoff (def: 1e-5)\n
  -beam_histogram= bool\n
  Activates beam histogramming.  On by default.\n
  -nbeamhist= value (int)\n
  If beam histogramming, number of bins in histogram 
  (def: 120)\n
  -fit_sigma=   bool\n
  If value is true, this turns on fitting for sigma
  -sigmaprior=  stdev (float)\n
  This activates the Gaussian prior on the instrumental noise
  multiplier.  It has mean one and standard deviation stdev.
  Setting this will automatically turn on sigma fitting.\n
  -exp_conf=    value (float)\n
  Sets the expected confusion noise value in Jy; used to
  set the expected likelihood values.\n
  -cfirbprior=  mean stdev (floats)\n
  This activates a Gaussian prior on the integrated CFIRB 
  luminoisty. It has a mean of mean and a standard deviation 
  stdev expressed in the same flux and area units as the model.
  So if your model dN/dS is in Jy^-1 deg^-2,
  then your CFIRB prior should be in Jy deg^-2\n
  -poissonprior=  mean stdev (floats)\n
  This activates a Gaussian prior on the Poisson noise.
  It has a mean of mean and a standard deviation 
  stdev expressed in the same flux and area units as the model.
  So if your model dN/dS is in Jy^-1 deg^-2,
  then your Poisson prior should be in Jy^2 deg^-2\n
  -regularize_alpha= alpha (float)\n
  Sets regularization multiplier for difference operator
  penalty in log likelihood (def: 0).\n
  -wisdom_file= filename (string)\n
  Uses filename as a FFTW wisdom file\n
  -verbose=     bool\n
  Turns on verbose mode (verbosity = 1)\n
  -ultraverbose= bool\n
  Turns on ultraverbose mode (verbosity = 2)\n
  -verbosity=   value (int)\n
  Sets verbosity to this level.\n
  -seed=        value (int)\n
  Sets random number generator seed\n

  Having multiple lines with dataset= is additive.
  Having multiple lines with any of the others results in only
  the value(s) from the final such line being used.

*/
struct specFile {
  //dataset
  std::vector<std::string> datafiles; //!< Data files
  std::vector<double> sigmas; //!< Instrumental noise values
  std::vector<std::string> psffiles; //!< Beam file
  std::vector<double> like_norm; //!< Likelihood normalization
  
  //Various fit options
  bool bin_data; //!< Bin input data
  unsigned int nbins; //!< Number of data bins
  bool mean_sub; //!< Mean subtract data
  bool ignore_mask; //!< Ignore mask info in input files
  unsigned int fftsize; //!< Length of fft
  unsigned int ninterp; //!< Size of R interpolation
  double minbeamval; //!< Minimum beam value used
  bool beam_histogram; //!< Do beam histogramming
  unsigned int nbeamhist; //!< Number of beam histogram bins
  bool fit_sigma; //!< Are we fitting for sigma
  bool has_sigprior; //!< Is the sigma prior on
  double sigprior_stdev; //!< Stdev of sigma prior
  double exp_conf; //!< Expected confusion noise in Jy
  bool has_cfirbprior; //!< Is cfirb prior on
  double cfirbprior_mean; //!< Mean value of prior
  double cfirbprior_stdev; //!< Stdev of prior
  bool has_poissonprior; //!< Is poisson prior on
  double poissonprior_mean; //!< Mean value of prior
  double poissonprior_stdev; //!< Stdev of prior
  double regularization_alpha; //!< Regularization multiplier
  bool has_wisdom_file; //!< Has FFTW wisdom file
  std::string wisdom_file; //!< Name of wisdom file
  unsigned int verbosity; //!< Verbosity level
  bool has_user_seed; //!< User set seed
  unsigned long long int seed; //!< RNG seed
  
  specFile(); //!< Default constructor
  specFile(const std::string&); //!< Constructor with file read

  void init(); //!< Clear all values

  void writeToHDF5Handle(hid_t) const; //!< Write filenames to HDF5 handle

  void readFile(const std::string&); //!< Read in file
};

#endif
