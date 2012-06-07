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
#include<numberCountsKnotsSpline.h>
#include<paramSet.h>
#include<affineExcept.h>

int main(int argc, char **argv) {

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
    {"help", no_argument, 0, 'h'},
    {"fits", no_argument, 0, 'f'},
    {"histogram", no_argument, 0, 'H'},
    {"maxflux",required_argument,0,'2'},
    {"nflux",required_argument,0,'n'},
    {"ninterp",required_argument,0,'N'},
    {"sigmanoise",required_argument,0,'s'},
    {"verbose",no_argument,0,'v'},
    {"version",no_argument,0,'V'},
    {"wisdom",required_argument,0,'w'},
    {0,0,0,0}
  };
  
  char optstring[] = "hfH2:n:N:s:vVw:";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getPD -- get the P(D) for a" << std::endl;
      std::cerr << "\tpower law model" << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tpofd_affine_getPD [options] initfile beamfile outputfile"
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tWrites out the P(D) for the model in initfile using"
		<< std::endl;
      std::cerr << "\ta multiply broken power law model."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tinitfile is a text file giving the positions of the"
		<< " knot points," << std::endl;
      std::cerr << "\ttheir initial values, and their estimated errors in"
		<< std::endl;
      std::cerr << "\tthe format knotflux value error."
		<< std::endl;
      std::cerr << "\tBeamfile is the beam in FITS format." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-f, --fits" << std::endl;
      std::cerr << "\t\tWrite the output as a FITS file." << std::endl;
      std::cerr << "\t-h, --help" << std::endl;
      std::cerr << "\t\tPrint this message and exit." << std::endl;
      std::cerr << "\t-H, --histogram" << std::endl;
      std::cerr << "\t\tUse beam histogramming" << std::endl;
      std::cerr << "\t--maxflux value" << std::endl;
      std::cerr << "\t\tMaximum flux value. (def: Largest knot minus the"
		<< std::endl;
      std::cerr << "\t\tmean flux from the model)" << std::endl;
      std::cerr << "\t-n, --nflux value" << std::endl;
      std::cerr << "\t\tThe number of requested fluxes, also the FFT length." 
		<< std::endl;
      std::cerr << "\t\tShould be a power of 2. (def: 262144)"
		<< std::endl;
      std::cerr << "\t-N, --ninterp value" << std::endl;
      std::cerr << "\t\tLength of interpolation vector used in computing R"
		<< std::endl;
      std::cerr << "\t-s, --sigmanoise noise" << std::endl;
      std::cerr << "\t\tThe assumed per-pixel noise (assumed Gaussian, "
		<< "def: 0.002)" << std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput version number and exit" << std::endl;
      std::cerr << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cerr << "\t\tName of wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
      return 0;
      break;
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
