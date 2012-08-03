//specFile.h

#ifndef __specFile__
#define __specFile__

#include<string>
#include<vector>

/*!
  \brief Structure to read in spec files for fits

  Lines in the spec file control the likelihood calculations and fits.
  The format is: identifier= values
  The appropriate entries for values depends upon identifier=
  Where values is bool the option is turned on if value is either
    true or yes, and is not activated for any other value.

  The valid values for identifier= and their associated values are:
   dataset=     datafile inst_sigma psffile [like_norm]
                  datafile is the name of a fits data file
		  inst_sigma is the instrumental noise (Gaussian) in the same
		   units as the data in datafile (e.g., Jy)
		  psffile is a FITS file giving the beam
		  like_norm is an optional input controlling how the
		   likelihood is normalized for the beam area (def: 1)
   bin_data=    nbins
                  Turns on binning of data files with specified number 
                  of bins
   meansub=     bool
                  Do mean subtraction on input data if value is true
   fftsize=     nfft
                  Sets fft transform length (def: 131072)
   beam_histogram= bool
                  Activates beam histogramming
   fit_sigma=   value
                  If value is true, this turns on fitting for sigma
   sigmaprior=  stdev
                  This activates the Gaussian prior on the instrumental noise
		  multiplier.  It has mean one and standard deviation stdev.
                  Setting this will automatically turn on sigma fitting.
   cfirbprior=  mean stdev
                  This activates a Gaussian prior on the integrated CFIRB 
		  luminoisty. It has a mean of mean and a standard deviation 
		  stdev expressed in the same flux and area units as the model.
		  So if your model dN/dS is in Jy^-1 deg^-2,
		  then your CFIRB prior should be in Jy deg^-2
   wisdom_file= filename
                  Uses filename as a FFTW wisdom file		  
   verbose=     bool
                  Turns on verbose mode
   ultraverbose= bool
                  Turns on ultraverbose mode

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
  bool meansub; //!< Mean subtract data
  unsigned int fftsize; //!< Length of fft
  bool beam_histogram; //!< Do beam histogramming
  bool fit_sigma; //!< Are we fitting for sigma
  bool has_sigprior; //!< Is the sigma prior on
  double sigprior_stdev; //!< Stdev of sigma prior
  bool has_cfirbprior; //!< Is cfirb prior on
  double cfirbprior_mean; //!< Mean value of prior
  double cfirbprior_stdev; //!< Stdev of prior
  bool has_wisdom_file; //!< Has FFTW wisdom file
  std::string wisfile; //!< Name of wisdom file
  bool verbose; //!< Run in verbose mode
  bool ultraverbose; //!< Run in ultraverbose mode

  specFile(); //!< Default constructor
  specFile(const std::string&); //!< Constructor with file read

  void init(); //!< Clear all values

  void readFile(const std::string&); //!< Read in file
};

#endif
