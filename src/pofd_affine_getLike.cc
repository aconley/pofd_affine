#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<iomanip>
#include<sstream>

#include<getopt.h>

#include<utility.h>
#include<calcLike.h>
#include<beam.h>
#include<numberCountsKnotsSpline.h>
#include<paramSet.h>
#include<affineExcept.h>

int main(int argc, char **argv) {
  double sigma_mult;
  unsigned int ninterp, fftsize, nbins;

  bool meansub, histogram, bindata;
  bool verbose, ultraverbose, has_wisdom;
  std::string wisdom_file;
  
  std::string specfile; //Input file
  std::string initfile; //Init file (having model we want)

  //Defaults
  sigma_mult          = 1.0;
  ninterp             = 1024;
  fftsize             = 131072;
  nbins               = 1000;
  meansub             = false;
  histogram           = false;
  bindata             = false;
  verbose             = false;
  ultraverbose        = false;
  has_wisdom          = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"bindata",no_argument, 0, 'b'},
    {"histogram",no_argument,0,'H'},
    {"meansub",no_argument,0,'$'},
    {"fftsize",required_argument,0,'f'},
    {"ninterp",required_argument,0,'n'},
    {"nbins",required_argument,0,'N'},
    {"sigma_mult",required_argument,0,'s'},
    {"ultraverbose",no_argument,0,'u'},
    {"verbose",no_argument,0,'v'},
    {"version",no_argument,0,'V'},
    {"wisdom",required_argument,0,'w'},
    {0,0,0,0}
  };
  
  char optstring[] = "hbH$f:n:N:s:uvVw:";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getLike -- Evaluate the likelihood of "
		<< "a multiply broken" << std::endl;
      std::cerr << "\tpower law model" << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tpofd_affine_getLike [options] specfile initfile" 
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates the likelihood of the data specified by"
		<< " specfile" << std::endl;
      std::cerr << "\tfor the model in in initfile using a multiply broken "
		<< "power" << std::endl;
      std::cerr << "\tlaw model and the P(D) formalism." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tspecfile is a text file containing the beams, datafiles,"
		<< " sigma," << std::endl;
      std::cerr << "\tetc. using the same format as pofd_affine."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tinitfile is a text file giving the positions of the"
		<< "knot points," << std::endl;
      std::cerr << "\ttheir initial values, and their estimated errors in"
		<< std::endl;
      std::cerr << "\tthe format knotflux value error.  The error is ignored."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-b, --bindata" << std::endl;
      std::cerr << "\t\tBin the data when computing the likelihood."
		<< std::endl;
      std::cerr << "\t-f, --fftsize minsize" << std::endl;
      std::cerr << "\t\tFFT size Must be a power of 2. (def: 131072)"
		<< std::endl;
      std::cerr << "\t-h --help" << std::endl;
      std::cerr << "\t\tPrint this message and exit." << std::endl;
      std::cerr << "\t-H, --histogram" << std::endl;
      std::cerr << "\t\tUse beam histogramming" << std::endl;
      std::cerr << "\t--meansub" << std::endl;
      std::cerr << "\t\tMean-subtract the input data."
		<< std::endl;
      std::cerr << "\t-n, --ninterp ninterp" << std::endl;
      std::cerr << "\t\tNumber of interpolation points in R calculation" 
		<< " (def: 1024)" << std::endl;
      std::cerr << "\t-N, --nbins VAL" << std::endl;
      std::cerr << "\t\tNumber of bins if binning data (def: 1000)."
		<< std::endl;
      std::cerr << "\t-s, --sigma_mult value" << std::endl;
      std::cerr << "\t\tAmount to multiply noise in spec file by (def: 1.0)."
		<< std::endl;
      std::cerr << "\t-u, --ultraverbose" << std::endl;
      std::cerr << "\t\tLike verbose, but more so." << std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tTurns on verbose mode." << std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput version number and exit" << std::endl;
      std::cerr << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cerr << "\t\tName of wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
      return 0;
      break;
    case 'b' :
      bindata = true;
      break;
    case 'f' :
      fftsize = static_cast< unsigned int >( atoi(optarg) );
      break;
    case 'H' :
      histogram = true;
      break;
    case 'n' :
      ninterp = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'N' :
      nbins =  static_cast<unsigned int>(atoi(optarg));
      break;
    case '$' :
      meansub = true;
      break;
    case 's' :
      sigma_mult = atof(optarg);
      break;
    case 'u' :
      verbose=ultraverbose=true;
      break;
    case 'v' :
      verbose = true;
      break;
    case 'V' :
      std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
		<< std::endl;
      return 0;
      break;
    case 'w' :
      has_wisdom = true;
      wisdom_file = std::string( optarg );
      break;
    }

  if (optind >= argc-1 ) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  specfile = std::string( argv[optind] );
  initfile = std::string( argv[optind+1] );

  //Input tests
  if (fftsize < 32768) {
    std::cerr << "Inadvisably small minimum fft size" << std::endl;
    return 1;
  }
  if (fftsize & (fftsize-1)) {
    std::cerr << "fftsize must be power of 2" << std::endl;
    return 1;
  }
  if (ninterp == 0) {
    std::cerr << "Invalid ninterp -- must be positive" << std::endl;
    return 1;
  }
  if (sigma_mult < 0) {
    std::cerr << "Invalid (non-positive) sigma mult" << std::endl;
    return 1;
  }
  if (bindata && (nbins == 0)) {
    std::cerr << "Invalid number of bins (non-positive)"
	      << std::endl;
  }

  //Read in the initialization file knot positions, values, errors
  if (verbose || ultraverbose)
    std::cout << "Reading initialization file: " << initfile << std::endl;
  std::ifstream initfs( initfile.c_str() );
  if (!initfs) {
    std::cerr << "Error readining in initialization file: "
	      << initfile << std::endl;
    return 1;
  }
  std::vector<double> knotpos, knotval;
  unsigned int nknots;
  std::string line;
  std::stringstream str;
  std::vector<std::string> words;
  double currpos, currval;
  while (!initfs.eof()) {
    std::getline(initfs,line);
    if (line[0] == '#') continue;

    //Parse into words, stipping spaces
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < 2) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currpos;
    str.str(words[1]); str.clear(); str >> currval;
    knotpos.push_back(currpos);
    knotval.push_back(currval);
  }
  initfs.close();
  if (knotpos.size() == 0) {
    std::cerr << "No knot positions or values found in " << initfile
	      << std::endl;
    return 1;
  }
  nknots = knotpos.size();

  //Process the spec file
  if (ultraverbose)
    std::cout << "Reading spec file: " << specfile << std::endl;
  std::ifstream ifs( specfile.c_str() );
  if ( ! ifs ) {
    std::cerr << "Error reading spec file: " << specfile << std::endl;
    return 1;
  }
  std::vector< std::string > datafiles;
  std::vector< double > sigmas, like_norm;
  std::vector< std::string > psffiles;
  double lnorm;
  while (! ifs.eof() ) {
    std::getline(ifs,line);
    if (line[0] == '#') continue;
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < 3) continue; //Has wrong number of entries
    
    double curr_sigma;
    datafiles.push_back( words[0] );
    str.str(words[1]); str.clear(); str >> curr_sigma;
    if (curr_sigma < 0.0) {
      std::cerr << "Invalid (negative) sigma " << curr_sigma 
		<< " from line: " << line << std::endl;
      return 1;
    }
    sigmas.push_back( curr_sigma );
    

    psffiles.push_back( words[2] );
    if (words.size() >= 4) {
      str.str(words[3]); str.clear(); str >> lnorm;
      like_norm.push_back(lnorm);
    } else like_norm.push_back(1.0);
  }
  ifs.close();
  unsigned int ndata = datafiles.size();
  if (ndata == 0) {
    std::cerr << "No datafiles loaded" << std::endl;
    return 1;
  }

  try {
    calcLike likeSet( fftsize, ninterp, true, bindata, nbins );

    //Read data
    if (verbose || ultraverbose)
      std::cout << "Reading in data files" << std::endl;
    likeSet.readDataFromFiles( datafiles, psffiles, sigmas, like_norm,
			       false, meansub, histogram );
    
    if (has_wisdom) likeSet.addWisdom(wisdom_file);

    if (ultraverbose) likeSet.setVerbose();
      
    if (ultraverbose && histogram)
      std::cout << "  Using beam histogramming to reduce beam size"
		<< std::endl;

    if ( ultraverbose ) {
      printf("  FFTsize:       %u\n",fftsize);
      printf("  Nknots:        %u\n",nknots);
      if (histogram)
	printf("  Using histogramming to reduce beam size\n");  
      if (bindata)
	printf("  Using histogramming to reduce data size to: %u\n",nbins);  
      printf("  Positions and initial values:\n");
      for (unsigned int i = 0; i < nknots; ++i)
	printf("  %11.5e  %11.5e\n",knotpos[i],knotval[i]);
    }
  
    //And, get that likelihood
    likeSet.setKnotPositions( knotpos );
    paramSet params(nknots+1);
    for (unsigned int i = 0; i < nknots; ++i)
      params[i] = knotval[i];
    params[nknots] = sigma_mult;

    double LogLike = likeSet.getLogLike(params);

    std::cout << "log Likelihoood is: "
	      << std::setprecision(10) << LogLike << std::endl;
  } catch (const affineExcept& ex) {
    std::cerr << ex << std::endl;
    return 4;
  }
  return 0;
}
