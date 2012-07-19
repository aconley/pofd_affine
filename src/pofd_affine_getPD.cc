#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<sstream>

#include<getopt.h>

#include<global_settings.h>
#include<utility.h>
#include<PDFactory.h>
#include<PDFactoryDouble.h>
#include<numberCountsKnotsSpline.h>
#include<numberCountsDoubleLogNormal.h>
#include<paramSet.h>
#include<affineExcept.h>

int getPDSingle(int argc, char **argv) {

  double sigma_noise; //Noise
  unsigned int nflux, ninterp;
  bool has_wisdom, has_user_maxflux;
  bool histogram, verbose, write_to_fits;
  double maxflux;
  std::string wisdom_file;
  
  std::string initfile; //Init file (having model we want)
  std::string outputfile; //Ouput pofd option
  std::string psffile;  //Beam file

  //Defaults
  sigma_noise         = 0.002;
  has_wisdom          = false;
  nflux               = 262144;
  ninterp             = 1024;
  maxflux             = std::numeric_limits<double>::quiet_NaN();
  has_user_maxflux    = false;
  verbose             = false;
  histogram           = false;
  write_to_fits       = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"fits", no_argument, 0, 'f'},
    {"histogram", no_argument, 0, 'H'},
    {"maxflux",required_argument,0,'2'},
    {"nflux",required_argument,0,'n'},
    {"ninterp",required_argument,0,'N'},
    {"sigmanoise",required_argument,0,'s'},
    {"verbose",no_argument,0,'v'},
    {"wisdom",required_argument,0,'w'},
    {0,0,0,0}
  };
  
  char optstring[] = "fH2:n:N:s:vw:";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'f' :
      write_to_fits = true;
      break;
    case 'H' :
      histogram = true;
      break;
    case '2' :
      maxflux = atof(optarg);
      has_user_maxflux = true;
      break;
    case 'n' :
      nflux = static_cast<unsigned int>( atoi(optarg) );
      break;
    case 'N' :
      ninterp = static_cast<unsigned int>( atoi(optarg) );
      break;
    case 's' :
      sigma_noise = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    case 'w' :
      has_wisdom = true;
      wisdom_file = std::string( optarg );
      break;
    }

  if (optind >= argc-2 ) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  initfile   = std::string( argv[optind] );
  psffile    = std::string( argv[optind+1] );
  outputfile = std::string( argv[optind+2] );

  //Input tests
  if (nflux == 0) {
    std::cerr << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (nflux & (nflux-1)) {
    std::cerr << "nflux must be power of 2" << std::endl;
    std::cerr << " yours is: " << nflux << std::endl;
    return 1;
  }
  if (ninterp == 0) {
    std::cerr << "Error -- ninterp is zero." << std::endl;
      return 1;
  }
  if (ninterp > nflux) {
    std::cerr << "nflux must be larger than ninterp"
	      << std::endl;
    return 1;
  }
  if (sigma_noise < 0.0) {
    std::cerr << "Invalid noise level: must be >= 0.0" << std::endl;
    return 1;
  }

  //Read in the initialization file knot positions, values, errors
  std::ifstream initfs( initfile.c_str() );
  if (!initfs) {
    std::cerr << "Error readining in initialization file: "
	      << initfile << std::endl;
    return 2;
  }
  std::vector<double> knotpos, knotval, knoterr;
  unsigned int nknots;
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
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
    knotval.push_back( currval );
  }
  initfs.close();
  if (knotpos.size() == 0) {
    std::cerr << "No knot positions or values found in " << initfile
	      << std::endl;
    return 2;
  }
  nknots = knotpos.size();

  //Actual PD computation
  try {
    PDFactory pfactory(ninterp);
    beam bm( psffile, histogram );
    numberCountsKnotsSpline model(knotpos);
    paramSet params(knotval);
    model.setParams(params);
    PD pd;

    bool success;
    if (has_wisdom) {
      std::cout << "Reading in wisdom file: " << wisdom_file 
		<< std::endl;
      success = pfactory.addWisdom(wisdom_file);
      if (!success) {
	std::cerr << "Error reading wisdom file: " << wisdom_file << std::endl;
	return 4;
      }
    }
    if (verbose) pfactory.setVerbose();

    if (! has_user_maxflux) {
      //Need to decide on some estimate
      double meanFluxPerBeam = model.getMeanFluxPerArea() * 
	bm.getEffectiveArea();
      maxflux = model.getMaxFlux() - meanFluxPerBeam;
    }

    if (verbose) {
      printf("  Beam area: %f\n",bm.getEffectiveArea());
      printf("  Mean flux per area: %0.8f\n",
	     model.getMeanFluxPerArea());
      printf("  Nknots: %u\n",nknots);
      printf("  Sigma:  %0.5f\n",sigma_noise);
      if (histogram)
	printf("  Using beam histogramming to reduce beam size to: %u %u\n",
	       bm.getNPos(),bm.getNNeg());
      printf("  Interpolation length:  %u\n",ninterp);
      printf("  Positions and initial values:\n");
      for (unsigned int i = 0; i < nknots; ++i)
	printf("  %11.5e  %11.5e\n",knotpos[i],knotval[i]);
    }

    //Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << std::endl;
    pfactory.initPD(nflux,sigma_noise,maxflux,model,bm);
    pfactory.getPD( sigma_noise, pd, false, true );
    
    //Write it
    if (verbose) std::cout << "Writing P(D) to file " << outputfile 
			   << std::endl;
    if (write_to_fits) {
      pd.writeToFits(outputfile);
    } else {
      std::ofstream ofs(outputfile.c_str());
      if (!ofs) {
	std::cerr << "Error opening output file: " << outputfile
		  << std::endl;
	return 64;
      }
      ofs << pd;
      ofs.close();
    }

  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    return 8;
  }

  return 0;
}

////////////////////////
int getPDDouble(int argc, char** argv) {

  double sigma1, sigma2; //Noise
  unsigned int nflux, nedge;
  bool has_wisdom, has_user_maxflux1, has_user_maxflux2;
  bool histogram, verbose, write_to_fits, doedge;
  double maxflux1, maxflux2;
  std::string wisdom_file;
  
  std::string initfile; //Init file (having model we want)
  std::string outputfile; //Ouput pofd option
  std::string psffile1, psffile2;  //Beam file

  //Defaults
  sigma1              = 0.002;
  sigma2              = 0.002;
  has_wisdom          = false;
  nflux               = 2048;
  maxflux1            = std::numeric_limits<double>::quiet_NaN();
  maxflux2            = std::numeric_limits<double>::quiet_NaN();
  has_user_maxflux1   = false;
  has_user_maxflux2   = false;
  verbose             = false;
  histogram           = false;
  write_to_fits       = false;
  nedge               = 256;
  doedge              = true;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"fits", no_argument, 0, 'f'},
    {"histogram", no_argument, 0, 'H'},
    {"maxflux1",required_argument,0,'2'},
    {"maxflux2",required_argument,0,'3'},
    {"nflux",required_argument,0,'n'},
    {"nedge",required_argument,0,'N'},
    {"sigma1",required_argument,0,'s'},
    {"sigma2",required_argument,0,'S'},
    {"verbose",no_argument,0,'v'},
    {"wisdom",required_argument,0,'w'},
    {0,0,0,0}
  };
  
  char optstring[] = "fH2:3:n:N:s:S:vw:";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'f' :
      write_to_fits = true;
      break;
    case 'H' :
      histogram = true;
      break;
    case '2' :
      maxflux1 = atof(optarg);
      has_user_maxflux1 = true;
      break;
    case '3' :
      maxflux2 = atof(optarg);
      has_user_maxflux2 = true;
      break;
    case 'n' :
      nflux = static_cast<unsigned int>( atoi(optarg) );
      break;
    case 'N' :
      nedge = static_cast<unsigned int>( atoi(optarg) );
      break;
    case 's' :
      sigma1 = atof(optarg);
      break;
    case 'S' :
      sigma2 = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    case 'w' :
      has_wisdom = true;
      wisdom_file = std::string( optarg );
      break;
    }

  if (optind >= argc-3 ) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  initfile   = std::string( argv[optind] );
  psffile1   = std::string( argv[optind+1] );
  psffile2   = std::string( argv[optind+2] );
  outputfile = std::string( argv[optind+3] );

  //Input tests
  if (nflux == 0) {
    std::cerr << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (nflux & (nflux-1)) {
    std::cerr << "nflux must be power of 2" << std::endl;
    std::cerr << " yours is: " << nflux << std::endl;
    return 1;
  }
  if (sigma1 < 0.0) {
    std::cerr << "Invalid noise level in band 1: must be >= 0.0" << std::endl;
    return 1;
  }
  if (sigma2 < 0.0) {
    std::cerr << "Invalid noise level in band 2: must be >= 0.0" << std::endl;
    return 1;
  }
  if (nedge == 0) doedge = false;

  //Read in the initialization file knot positions, values, errors
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

  //Actual PD computation
  try {
    PDFactoryDouble pfactory(nedge);
    doublebeam bm( psffile1, psffile2, histogram );
    numberCountsDoubleLogNormal model(knotpos,sigmapos,offsetpos);
    paramSet params(wvec2);
    model.setParams(params);
    PDDouble pd;

    bool success;
    if (has_wisdom) {
      std::cout << "Reading in wisdom file: " << wisdom_file 
		<< std::endl;
      success = pfactory.addWisdom(wisdom_file);
      if (!success) {
	std::cerr << "Error reading wisdom file: " << wisdom_file << std::endl;
	return 4;
      }
    }
    if (verbose) pfactory.setVerbose();

    if (! has_user_maxflux1) {
      //Need to decide on some estimate
      double meanFluxPerBeam = model.getMeanFluxPerArea(0) * 
	bm.getEffectiveArea1();
      maxflux1 = model.getMaxFlux(0) - meanFluxPerBeam;
    }
    if (! has_user_maxflux2) {
      double meanFluxPerBeam = model.getMeanFluxPerArea(1) * 
	bm.getEffectiveArea2();
      maxflux2 = model.getMaxFlux(1) + 5*sigmapos[ns-1] - meanFluxPerBeam;
    }

    if (verbose) {
      printf("  Beam area, band 1: %f\n",bm.getEffectiveArea1());
      printf("  Mean flux per area, band 1: %0.8f\n",
	     model.getMeanFluxPerArea(0));
      printf("  Beam area, band 2: %f\n",bm.getEffectiveArea2());
      printf("  Mean flux per area, band 2: %0.8f\n",
	     model.getMeanFluxPerArea(1));
      printf("  Nknots: %u\n",nk);
      printf("  Nsigma: %u\n",ns);
      printf("  Noffset: %u\n",no);
      printf("  Sigma, band 1:  %0.5f\n",sigma1);
      printf("  Sigma, band 2:  %0.5f\n",sigma2);
      if (histogram)
	printf("  Using beam histogramming to reduce beam size\n");
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

    //Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << std::endl;
    pfactory.initPD(nflux,sigma1,sigma2,maxflux1,maxflux2,model,bm,doedge);
    pfactory.getPD( sigma1, sigma2, pd, false, true );
    
    //Write it
    if (verbose) std::cout << "Writing P(D) to file " << outputfile 
			   << std::endl;
    if (write_to_fits) {
      pd.writeToFits(outputfile);
    } else {
      std::ofstream ofs(outputfile.c_str());
      if (!ofs) {
	std::cerr << "Error opening output file: " << outputfile
		  << std::endl;
	return 64;
      }
      ofs << pd;
      ofs.close();
    }

  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    return 8;
  }

  return 0;
}



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
      std::cerr << "\tpofd_affine_getPD -- get the PD for a number counts "
		<< "model. Both" << std::endl;
      std::cerr << "\tone-dimensional and two-dimensional models are supported."
		<< std::endl;
      std::cerr << "\tmodel" << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tEither" << std::endl;
      std::cerr << std::endl;

      std::cerr << "\t pofd_mcmc_getPD [options] initfile beamfile outputfile"
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfor the 1D case or" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\t pofd_mcmc_getPD -d [options] initfile beamfile1"
		<< " beamfile2" << std::endl;
      std::cerr << "\t  outfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfor the 2D case." << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates P(D) for the model in initfile and write it"
		<< " to" << std::endl;
      std::cerr << "\toutfile.  The 1D model is a log-space spline model for"
		<< " the" << std::endl;
      std::cerr << "\tnumber counts, and the 2D model is the 1D spline model" 
		<< " times" << std::endl;
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
      std::cerr << "\tbeamfile gives the name of a FITS file containing the "
		<< "beam" << std::endl;
      std::cerr << "\tin the 1D case, and beamfile1, beamfile2 give the beam "
		<< "in each" << std::endl;
      std::cerr << "\t of the two bands in the 2D case." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tIn both cases the output PD is written to outfile either"
		<< " as" << std::endl;
      std::cerr << "\ttext or as a FITS file." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-d, --double" << std::endl;
      std::cerr << "\t\tUse the 2D model." << std::endl;
      std::cerr << "\t-f, --fits" << std::endl;
      std::cerr << "\t\tWrite the output as a FITS file." << std::endl;
      std::cerr << "\t-H, --histogram" << std::endl;
      std::cerr << "\t\tUse beam histogramming." << std::endl;
      std::cerr << "\t-n, --nflux value" << std::endl;
      std::cerr << "\t\tThe number of requested fluxes, also the FFT length." 
		<< std::endl;
      std::cerr << "\t\tShould be a power of 2. For the 2D model, this is the"
		<< std::endl;
      std::cerr << "\t\tsize along each dimension. (def: 262144 for 1D, 2048"
		<< std::endl;
      std::cerr << "\t\tfor 2D)" << std::endl;
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
      std::cerr << "\t--maxflux value" << std::endl;
      std::cerr << "\t\tMaximum flux value. (def: Largest knot minus the"
		<< std::endl;
      std::cerr << "\t\tmean flux from the model)" << std::endl;
      std::cerr << "\t-N, --ninterp value" << std::endl;
      std::cerr << "\t\tLength of interpolation vector used in computing R"
		<< std::endl;
      std::cerr << "\t-s, --sigmanoise noise" << std::endl;
      std::cerr << "\t\tThe assumed per-pixel noise (assumed Gaussian, "
		<< "def: 0.002)" << std::endl;
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
    return getPDSingle(argc,argv);
  else
    return getPDDouble(argc,argv);
}
