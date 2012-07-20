#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<iomanip>
#include<sstream>

#include<getopt.h>

#include<calcLike.h>
#include<calcLikeDouble.h>
#include<affineExcept.h>

int getLikeSingle(int argc, char **argv) {
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
    {"bindata",no_argument, 0, 'b'},
    {"histogram",no_argument,0,'H'},
    {"meansub",no_argument,0,'$'},
    {"fftsize",required_argument,0,'f'},
    {"ninterp",required_argument,0,'n'},
    {"nbins",required_argument,0,'N'},
    {"sigma_mult",required_argument,0,'s'},
    {"ultraverbose",no_argument,0,'u'},
    {"verbose",no_argument,0,'v'},
    {"wisdom",required_argument,0,'w'},
    {0,0,0,0}
  };
  
  char optstring[] = "bH$f:n:N:s:uvw:";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
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


/////////////////////////////////////

int getLikeDouble(int argc, char** argv) {

  double sigma_mult;
  unsigned int nedge, fftsize, nbins;

  bool meansub, histogram, bindata, edgeInterp;
  bool verbose, ultraverbose, has_wisdom;
  std::string wisdom_file;
  
  std::string specfile; //Input file
  std::string initfile; //Init file (having model we want)

  //Defaults
  sigma_mult          = 1.0;
  nedge               = 256;
  fftsize             = 2048;
  nbins               = 100;
  meansub             = false;
  histogram           = false;
  bindata             = false;
  verbose             = false;
  ultraverbose        = false;
  has_wisdom          = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"bindata",no_argument, 0, 'b'},
    {"histogram",no_argument,0,'H'},
    {"meansub",no_argument,0,'$'},
    {"fftsize",required_argument,0,'f'},
    {"nedge",required_argument,0,'n'},
    {"nbins",required_argument,0,'N'},
    {"sigma_mult",required_argument,0,'s'},
    {"ultraverbose",no_argument,0,'u'},
    {"verbose",no_argument,0,'v'},
    {"wisdom",required_argument,0,'w'},
    {0,0,0,0}
  };
  
  char optstring[] = "bH$f:n:N:s:uvw:";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
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
      nedge = static_cast<unsigned int>(atoi(optarg));
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
  if (fftsize < 512) {
    std::cerr << "Inadvisably small minimum fft size " << fftsize << std::endl;
    return 1;
  }
  if (fftsize & (fftsize-1)) {
    std::cerr << "fftsize must be power of 2" << std::endl;
    return 1;
  }
  if (nedge == 0) edgeInterp = false;

  if (sigma_mult < 0) {
    std::cerr << "Invalid (non-positive) sigma mult" << std::endl;
    return 1;
  }
  if (bindata && (nbins == 0)) {
    std::cerr << "Invalid number of bins (non-positive)"
	      << std::endl;
  }

  //////////////////////////////
  //Read in the initialization file knot positions, values, errors
  if (verbose || ultraverbose)
    std::cout << "Reading initialization file: " << initfile << std::endl;
  std::ifstream initfs( initfile.c_str() );
  if (!initfs) {
    std::cerr << "Error readining in initialization file: "
	      << initfile << std::endl;
    return 2;
  }
  unsigned int nk, ns, no; //Number of knots, sigmas, offsets
  std::vector<double> wvec1, wvec2;
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currpos, currval;

  //Read in number
  initfs >> nk >> ns >> no;
  if ( nk < 2 ) {
    initfs.close();
    std::cerr << "Need at least 2 knots!" << std::endl;
    return 4;
  }
  if ( ns < 1 ) {
    initfs.close();
    std::cerr << "Need at least 1 sigma knots!" << std::endl;
    return 4;
  }
  if ( no < 1 ) {
    initfs.close();
    std::cerr << "Need at least 1 offset knots!" << std::endl;
    return 4;
  }
  //Read in values
  while (!initfs.eof()) {
    std::getline(initfs,line);
    if (line[0] == '#') continue; //Comment
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < 2) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currpos;
    str.str(words[1]); str.clear(); str >> currval;
    wvec1.push_back( currpos );
    wvec2.push_back( currval ); 
  }
  initfs.close();
  if (wvec1.size() != nk+ns+no) {
   std::cerr << "Expected " << nk+ns+no << " values, got: " 
              << wvec1.size() << std::endl;
    return 8;
  }

  //Parse out into knot positions, etc.
  //Keep values in wvec2
  std::vector<double> knotpos(nk), sigmapos(ns), offsetpos(no);
  for (unsigned int i = 0; i < nk; ++i)
    knotpos[i] = wvec1[i];
  for (unsigned int i = nk; i < ns+nk; ++i)
    sigmapos[i-nk] = wvec1[i];
  for (unsigned int i = nk+ns; i < no+ns+nk; ++i)
    offsetpos[i-nk-ns] = wvec1[i];

  ////////////////////////
  //Process the spec file
  if (ultraverbose)
    std::cout << "Reading spec file: " << specfile << std::endl;
  std::ifstream ifs( specfile.c_str() );
  if ( ! ifs ) {
    std::cerr << "Error reading spec file: " << specfile << std::endl;
    return 1;
  }
  std::vector< std::string > datafiles1, datafiles2;
  std::vector< double > sigmas1, sigmas2, like_norm;
  std::vector< std::string > psffiles1, psffiles2;
  double lnorm;
  while (! ifs.eof() ) {
    std::getline(ifs,line);
    if (line[0] == '#') continue;
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < 6) continue; //Has wrong number of entries

    datafiles1.push_back( words[0] );
    datafiles2.push_back( words[1] );
    
    double curr_sigma;
    str.str(words[2]); str.clear(); str >> curr_sigma;
    if (curr_sigma < 0.0) {
      std::cerr << "Invalid (negative) sigma1 " << curr_sigma 
		<< " from line: " << line << std::endl;
      return 1;
    }
    sigmas1.push_back( curr_sigma );
    str.str(words[3]); str.clear(); str >> curr_sigma;
    if (curr_sigma < 0.0) {
      std::cerr << "Invalid (negative) sigma2 " << curr_sigma 
		<< " from line: " << line << std::endl;
      return 1;
    }
    sigmas2.push_back( curr_sigma );
    
    psffiles1.push_back( words[4] );
    psffiles2.push_back( words[5] );

    //like_norm (optional)
    if (words.size() >= 7) {
      str.str(words[6]); str.clear(); str >> lnorm;
      like_norm.push_back(lnorm);
    } else like_norm.push_back(1.0);
  }
  ifs.close();
  unsigned int ndata = datafiles1.size();
  if (ndata == 0) {
    std::cerr << "No datafiles loaded" << std::endl;
    return 1;
  }

  /////////////////////////////////////
  // Actual calculation
  try {
    calcLikeDouble likeSet( fftsize, nedge, true, edgeInterp,
			    bindata, nbins );

    //Read data
    if (verbose || ultraverbose)
      std::cout << "Reading in data files" << std::endl;
    likeSet.readDataFromFiles( datafiles1, datafiles2, psffiles1, psffiles2,
			       sigmas1, sigmas2, like_norm,
			       false, meansub, histogram );
    
    if (has_wisdom) likeSet.addWisdom(wisdom_file);

    if (ultraverbose) likeSet.setVerbose();
      
    if ( ultraverbose ) {
      printf("  FFTsize:       %u\n",fftsize);
      if (histogram)
	printf("  Using histogramming to reduce beam size\n");  
      if (bindata)
	printf("  Using histogramming to reduce data size to: %u by %u\n",
	       nbins,nbins);  
      printf("  Knot Positions and initial values:\n");
      for (unsigned int i = 0; i < nk; ++i)
	printf("  %11.5e  %11.5e\n",knotpos[i],wvec2[i]);
      printf("  Sigma Positions and initial values:\n");
      for (unsigned int i = 0; i < ns; ++i)
	printf("  %11.5e  %11.5e\n",sigmapos[i],wvec2[i+nk]);
      printf("  Offset Positions and initial values:\n");
      for (unsigned int i = 0; i < no; ++i)
	printf("  %11.5e  %11.5e\n",offsetpos[i],wvec2[i+nk+ns]);
    }
  
    //And, get that likelihood
    likeSet.setPositions( knotpos, offsetpos, sigmapos );
    paramSet params(wvec2);

    double LogLike = likeSet.getLogLike(params);

    std::cout << "log Likelihoood is: "
	      << std::setprecision(10) << LogLike << std::endl;
  } catch (const affineExcept& ex) {
    std::cerr << ex << std::endl;
    return 4;
  }
  return 0;
}

/////////////////////////////////

int main( int argc, char** argv ) {
  bool twod;

  twod = false;

  //Only interested in a) displaying help and b) figuring out
  // if this is 1D or 2D c) displaying the version number
  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"double",no_argument,0,'d'},
    {"version",no_argument,0,'V'},
    {0,0,0,0}
  };
  char optstring[] = "hdv";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getLike -- evaluate the likelihood of some "
		<< "data" << std::endl;
      std::cerr << "\tagainst a number counts model using the P(D) approach."
		<< " Both" << std::endl;
      std::cerr << "\tone-dimensional and two-dimensional models are supported."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\t pofd_mcmc_getLike [options] initfile specfile outputfile"
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates the likelihood of the data specified in "
		<< "specfile against" << std::endl;
      std::cerr << "\tthe model in initfile. The 1D model is a log-space spline"
		<< " model for" << std::endl;
      std::cerr << "\tthe number counts, and the 2D model is the 1D spline "
		<< "model times" << std::endl;
      std::cerr << "\ta log-normal color function for the second band, with the"
		<< std::endl;
      std::cerr << "\tlog-space variance and mean color stored as splines in"
		<< " the flux" << std::endl;
      std::cerr << "\tof the first band." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tinitfile is a text file specifying the model; the exact"
		<< " details" << std::endl;
      std::cerr << "\t(given below) depend on whether the 1D or 2D case is"
		<< " being used." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tFor the 1D case, initfile is a text file giving the "
		<< "positions" << std::endl;
      std::cerr << "\tof the spline knots and their values in the format"
		<< " knotflux value." << std::endl;
      std::cerr << "\tAdditional elements on each line are ignored."
		<< std::endl;
      std::cerr << "\tFor the 2D case, initfile is a text file giving the "
		<< "positions" << std::endl;
      std::cerr << "\tof the knot points and their values, followed by the "
		<< "sigma" << std::endl;
      std::cerr << "\tknot positions and their values, then likewise for the "
		<< "colour" << std::endl;
      std::cerr << "\toffset.  The format is three numbers on the first line, "
		<< "giving" << std::endl;
      std::cerr << "\tthe number of number count knots, sigma knots, and "
		<< "offset knots," << std::endl;
      std::cerr << "\tfollowed by a number of lines again with the format"
		<< std::endl;
      std::cerr << "\tknotpos value.  The sigmas and offsets are in log space."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tspecfile is a text file containing the beam(s), "
		<< "datafile(s), sigma(s)," << std::endl;
      std::cerr << "\tetc. using the same format as pofd_affine."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-b, --bindata" << std::endl;
      std::cerr << "\t\tBin the data when computing the likelihood."
		<< std::endl;
      std::cerr << "\t-d, --double" << std::endl;
      std::cerr << "\t\tUse the 2D model." << std::endl;
      std::cerr << "\t-f, --fftsize minsize" << std::endl;
      std::cerr << "\t\tFFT size Must be a power of 2. The size along each"
		<< " dimension" << std::endl;
      std::cerr << "\t\tfor the 2D model. (def: 131072 for 1D, 2048 for 2D)"
		<< std::endl;
      std::cerr << "\t-H, --histogram" << std::endl;
      std::cerr << "\t\tUse beam histogramming." << std::endl;
      std::cerr << "\t--meansub" << std::endl;
      std::cerr << "\t\tMean-subtract the input data."
		<< std::endl;
      std::cerr << "\t-N, --nbins VAL" << std::endl;
      std::cerr << "\t\tNumber of bins if binning data, along each dimension in"
		<< " the" << std::endl;
      std::cerr << "\t\t2D model case. (def: 1000 in 1D, 100 in 2D)."
		<< std::endl;
      std::cerr << "\t-s, --sigma_mult value" << std::endl;
      std::cerr << "\t\tAmount to multiply noise in spec file by (def: 1.0)."
		<< std::endl;
      std::cerr << "\t-u, --ultraverbose" << std::endl;
      std::cerr << "\t\tLike verbose, but more so." << std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput the version number and exit." << std::endl;
      std::cerr << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cerr << "\t\tName of FFTW wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "ONE-D ONLY OPTIONS" << std::endl;
      std::cerr << "\t-n, --ninterp ninterp" << std::endl;
      std::cerr << "\t\tNumber of interpolation points in R calculation" 
		<< " (def: 1024)" << std::endl;
      std::cerr << "TWO-D ONLY OPTIONS" << std::endl;
      std::cerr << "\t-n, --nedge nedge" << std::endl;
      std::cerr << "\t\tNumber of edge integral points in R calculation" 
		<< " (def: 256)" << std::endl;
      return 0;
      break;
    case 'd' :
      twod = true;
      break;
    case 'V' :
      std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
		<< std::endl;
      return 0;
      break;
    }

  if (! twod)
    return getLikeSingle(argc,argv);
  else
    return getLikeDouble(argc,argv);
}

