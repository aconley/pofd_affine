#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

#include<getopt.h>

#include<global_settings.h>
#include<numberCountsKnotsSpline.h>
#include<numberCountsDoubleLogNormal.h>
#include<paramSet.h>
#include<utility.h>
#include<affineExcept.h>

//One-D version
int getNSingle( int argc, char** argv ) {
  
  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  double minflux, maxflux;
  unsigned int nflux;
  bool has_user_maxflux, has_user_minflux, verbose, logspace;

  minflux          = 0;
  maxflux          = 1; //Will always be overridden
  nflux            = 1000;
  has_user_minflux = false;
  has_user_maxflux = false;
  logspace         = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"logspace",no_argument,0,'l'},
    {"minflux",required_argument,0,'m'},
    {"maxflux",required_argument,0,'M'},
    {"nflux",required_argument,0,'n'},
    {"verbose",no_argument,0,'v'},
    {0,0,0,0}
  };

  //Have to have global options here, even if ignored
  char optstring[] = "hdlm:M:n:vV"; 
  optind = 1; //Resets parsing
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'l' :
      logspace = true;
      break;
    case 'm' :
      has_user_minflux = true;
      minflux = atof(optarg);
      break;
    case 'M' :
      has_user_maxflux = true;
      maxflux = atof(optarg);
      break;
    case 'n' :
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-1 ) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  initfile = std::string( argv[optind] );
  outfile = std::string( argv[optind+1] );

  //Input check
  if (has_user_maxflux && (maxflux <= minflux)) {
    std::cerr << "Maxflux must be > minflux" << std::endl;
    return 1;
  }
  if (nflux == 0) {
    std::cerr << "nflux must be positive" << std::endl;
    return 2;
  }

  //Read in the model params
  std::ifstream initfs( initfile.c_str() );
  if (!initfs) {
    std::cerr << "Error readining in initialization file: "
	      << initfile << std::endl;
    return 1;
  }
  std::vector<double> knotpos, knotval;
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
    knotpos.push_back( currpos );
    knotval.push_back( currval ); 
  }
  initfs.close();
  if (knotpos.size() == 0) {
    std::cerr << "No knot positions or values found in " << initfile
	      << std::endl;
    return 1;
  }

  double *dNdS = NULL;
  try {
    numberCountsKnotsSpline model(knotpos);
    paramSet pars( knotval );
    model.setParams( pars );

    //Check minflux if logspace is set
    if (logspace) {
      if (! has_user_minflux)
	minflux = model.getKnotPos(0); //Sorted inside model
      if (minflux <= 0.0)
	throw affineExcept("pofd_affine_getdNdS","getNSingle",
			   "Minflux must be positive if using logspace",1);
    }

    //Get maximum flux if not set, use to compute dflux
    if (! has_user_maxflux )
      maxflux = model.getMaxFlux();

    double dflux;
    if (nflux > 1) 
      if (logspace) 
	dflux = log2(maxflux/minflux)/static_cast<double>(nflux-1);
      else
	dflux = (maxflux-minflux)/static_cast<double>(nflux-1);
    else
      dflux = 1.0; //Not used

    if (verbose) {
      std::cout << "Mean flux per sq degree: " << model.getMeanFluxPerArea()
		<< std::endl;
      std::cout << "Number of sources per area: "
		<< model.getNS() << std::endl;
    }

    //Calculation loop
    dNdS  = new double[nflux];
    if (logspace) {
      double lmin = log2(minflux);
      for (unsigned int i = 0; i < nflux; ++i)
	dNdS[i] = 
	  model.getNumberCounts(exp2(lmin+static_cast<double>(i)*dflux));
    } else
      for (unsigned int i = 0; i < nflux; ++i)
	dNdS[i] = model.getNumberCounts(minflux+static_cast<double>(i)*dflux);
    
    //Write out
    FILE *fp;
    fp = fopen( outfile.c_str(),"w");
    if (!fp) {
      std::cerr << "Failed to open output file" << std::endl;
      return 128;
    }
    fprintf(fp,"%12s   %12s\n","Flux","dNdS");
    if (logspace) {
      double lmin = log2(minflux);
      for (unsigned int i = 0; i < nflux; ++i) 
	fprintf(fp,"%12.6e   %15.9e\n",
		exp2(lmin+static_cast<double>(i)*dflux),dNdS[i]);
    } else for (unsigned int i = 0; i < nflux; ++i) 
	     fprintf(fp,"%12.6e   %15.9e\n",
		     minflux+static_cast<double>(i)*dflux,dNdS[i]);
    fclose(fp);
    
    delete[] dNdS; 
  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (dNdS != NULL) delete[] dNdS;
    return 8;
  }
  return 0;
}

/////////////////////////////////////////////

//Two-D version
int getNDouble( int argc, char** argv ) {
  
  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  double minflux1, maxflux1, minflux2, maxflux2;
  unsigned int nflux1,nflux2;
  bool has_user_maxflux1, has_user_maxflux2, verbose;

  minflux1          = 0;
  minflux2          = 0;
  maxflux1          = 1; //Will always be overridden
  maxflux2          = 1; //Will always be overridden
  nflux1            = 100;
  nflux2            = 100;
  has_user_maxflux1 = false;
  has_user_maxflux2 = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"minflux1",required_argument,0,'m'},
    {"minflux2",required_argument,0,'1'},
    {"maxflux1",required_argument,0,'M'},
    {"maxflux1",required_argument,0,'2'},
    {"nflux1",required_argument,0,'n'},
    {"nflux2",required_argument,0,'N'},
    {"verbose",no_argument,0,'v'},
    {0,0,0,0}
  };

  char optstring[] = "hdm:1:M:2:n:N:vV";
  optind = 1; //Resets parsing
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'm' :
      minflux1 = atof(optarg);
      break;
    case '1' :
      minflux2 = atof(optarg);
      break;
    case 'M' :
      has_user_maxflux1 = true;
      maxflux1 = atof(optarg);
      break;
    case '2' :
      has_user_maxflux2 = true;
      maxflux2 = atof(optarg);
      break;
    case 'n' :
      nflux1 = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'N' :
      nflux2 = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-1 ) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  initfile = std::string( argv[optind] );
  outfile = std::string( argv[optind+1] );

  //Input check
  if (has_user_maxflux1 && (maxflux1 <= minflux1)) {
    std::cerr << "Maxflux must be > minflux in band 1" << std::endl;
    return 1;
  }
  if (has_user_maxflux2 && (maxflux2 <= minflux2)) {
    std::cerr << "Maxflux must be > minflux in band 2" << std::endl;
    return 2;
  }
  if (nflux1 == 0) {
    std::cerr << "nflux must be positive in band 1" << std::endl;
    return 3;
  }
  if (nflux2 == 0) {
    std::cerr << "nflux must be positive in band 2" << std::endl;
    return 4;
  }
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

  double *dNdS = NULL;
  try {
    numberCountsDoubleLogNormal model(knotpos,sigmapos,offsetpos);
    paramSet pars( wvec2 );
    model.setParams( pars );

    //Get maximum flux if not set, use to compute dflux
    if (! has_user_maxflux1 )
      maxflux1 = model.getMaxFlux(0);
    if (! has_user_maxflux2 )
      maxflux2 = model.getMaxFlux(1);
    double dflux1, dflux2;
    if (nflux1 > 1) 
      dflux1 = (maxflux1-minflux1)/static_cast<double>(nflux1-1);
    else
      dflux1 = 1.0; //Not used
    if (nflux2 > 1) 
      dflux2 = (maxflux2-minflux2)/static_cast<double>(nflux2-1);
    else
      dflux2 = 1.0; //Not used

    if (verbose) {
      std::cout << "Mean flux per sq degree, band 1: "
		<< model.getMeanFluxPerArea(0) << std::endl;
      std::cout << "Mean flux per sq degree, band 2: "
		<< model.getMeanFluxPerArea(1) << std::endl;
      std::cout << "Number of sources per area: "
		<< model.getNS() << std::endl;
    }

    //Calculation loop
    dNdS  = new double[nflux1*nflux2];
    double flux1, flux2;
    for (unsigned int i = 0; i < nflux1; ++i) {
      flux1 = minflux1+static_cast<double>(i)*dflux1;
      for (unsigned int j = 0; j < nflux2; ++j) {
	flux2 = minflux2+static_cast<double>(j)*dflux2;
	dNdS[i*nflux2+j] = model.getNumberCounts(flux1,flux2);
      }
    }
    
    //Write out
    FILE *fp;
    fp = fopen( outfile.c_str(),"w");
    if (!fp) {
      std::cerr << "Failed to open output file" << std::endl;
      return 128;
    }
    fprintf(fp,"%u %12.6e %12.6e\n",nflux1,minflux1,dflux1);
    fprintf(fp,"%u %12.6e %12.6e\n",nflux2,minflux2,dflux2);
    for (unsigned int i = 0; i < nflux1; ++i) {
      for (unsigned int j = 0; j < nflux2; ++j)
	fprintf(fp,"%13.7e",dNdS[i*nflux2+j]);
      fprintf(fp,"%13.7e\n",dNdS[i*nflux2+nflux2-1]);
    }
    fclose(fp);
    
    delete[] dNdS; 
  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (dNdS != NULL) delete[] dNdS;
    return 8;
  }
  return 0;
}



/////////////////////////////////////////////

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
  //Also have to include options for 1D/2D in this string, although
  // we ignore them
  char optstring[] = "hdlm:1:M:2:n:N:v";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getdNdS -- get differential number counts for"
		<< " a model." << std::endl;
      std::cerr << "\tBoth one-dimensional and two-dimensional models are "
		<< "supported." << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tEither" << std::endl;
      std::cerr << std::endl;

      std::cerr << "\t pofd_mcmc_getdNdS [options] initfile outfile"
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates the differential number counts for the model "
		<< "in initfile" << std::endl;
      std::cerr << "\tand writes it to outfile.  The 1D model is a log-space "
		<< "spline" << std::endl;
      std::cerr << "\tmodel for the number counts, and the 2D model is the 1D" 
		<< " spline" << std::endl;
      std::cerr << "\tmodel times a log-normal color function for the second"
		<< " band," << std::endl;
      std::cerr << "\twith the log-space variance and mean color stored as"
		<< " splines" << std::endl;
      std::cerr << "\tin the flux of the first band." << std::endl;
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
      std::cerr << "\tIn both cases the output R is written to outfile as"
		<< " text." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-d, --double" << std::endl;
      std::cerr << "\t\tUse the 2D model." << std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tRun in verbose mode" << std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput the version number and exit." << std::endl;
      std::cerr << std::endl;
      std::cerr << "ONE-D ONLY OPTIONS" << std::endl;
      std::cerr << "\t-l, --logspace" << std::endl;
      std::cerr << "\t\tSpace the fluxes in log space rather than linearly"
		<< std::endl;
      std::cerr << "\t-M, --maxflux value" << std::endl;
      std::cerr << "\t\tThe maximum flux to calculate dN/dS for (def: highest"
		<< std::endl;
      std::cerr << "\t\tknot specified in initfile)." << std::endl;
      std::cerr << "\t-m, --minflux value" << std::endl;
      std::cerr << "\t\tThe minimum flux to calculate dN/dS for (def: 0, "
		<< "unless" << std::endl;
      std::cerr << "\t\t--logspace is set)." << std::endl;
      std::cerr << "\t-n, --nflux value" << std::endl;
      std::cerr << "\t\tNumber of flux values to output (def: 1000)" 
		<< std::endl;
      std::cerr << "TWO-D ONLY OPTIONS" << std::endl;
      std::cerr << "\t--maxflux1 value" << std::endl;
      std::cerr << "\t\tThe maximum flux to calculate dN/dS for in band 1"
		<< std::endl;
      std::cerr << "\t\t(def: highest knot specified in initfile)." 
		<< std::endl;
      std::cerr << "\t--maxflux2 value" << std::endl;
      std::cerr << "\t\tThe maximum flux to calculate dN/dS for in band 2"
		<< std::endl;
      std::cerr << "\t\t(def: a few sigma above mean value at highest band 1 "
		<< "knot)." << std::endl;
      std::cerr << "\t--minflux2 value" << std::endl;
      std::cerr << "\t\tThe minimum flux to calculate dN/dS for in band 2. "
		<< "(def: 0)" << std::endl;
      std::cerr << "\t--nflux2 value" << std::endl;
      std::cerr << "\t\tNumber of flux values to output, band 2 (def: 100)" 
		<< std::endl;
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
    return getNSingle(argc,argv);
  else
    return getNDouble(argc,argv);
}