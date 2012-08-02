#include<iostream>

#include<getopt.h>

#include<global_settings.h>
#include<numberCountsKnotsSpline.h>
#include<numberCountsDoubleLogNormal.h>
#include<paramSet.h>
#include<utility.h>
#include<affineExcept.h>

//All sub-parses have to have the same options to avoid
// getting warnings.  So we give them all the same long_options,
// but then only process the appropriate ones, ignoring the rest
//Yes, this is a bit complicated and error prone, but such is life
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"double",no_argument,0,'d'},
  {"version",no_argument,0,'V'}, //After this, not parsed here
  {"logspace",no_argument,0,'l'},
  {"minflux",required_argument,0,'m'},
  {"maxflux",required_argument,0,'M'},
  {"nflux",required_argument,0,'n'},
  {"verbose",no_argument,0,'v'},
  {"minflux1",required_argument,0,'1'},
  {"minflux2",required_argument,0,'2'},
  {"maxflux1",required_argument,0,'3'},
  {"maxflux2",required_argument,0,'4'},
  {"nflux1",required_argument,0,'5'},
  {"nflux2",required_argument,0,'6'},
  {0,0,0,0}
};
char optstring[] = "hdVlm:M:n:v1:2:3:4:5:6:"; 

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

  double *dNdS = NULL;
  try {
    initFileKnots model_info;
    model_info.readFile(initfile, false, false);

    numberCountsKnotsSpline model;
    model_info.getKnotPos(model);

    paramSet pars(model_info.getNKnots());
    model_info.getParams(pars);
    model.setParams(pars);

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
  optind = 1; //Resets parsing
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case '1' :
      minflux1 = atof(optarg);
      break;
    case '2' :
      minflux2 = atof(optarg);
      break;
    case '3' :
      has_user_maxflux1 = true;
      maxflux1 = atof(optarg);
      break;
    case '4' :
      has_user_maxflux2 = true;
      maxflux2 = atof(optarg);
      break;
    case '5' :
      nflux1 = static_cast<unsigned int>(atoi(optarg));
      break;
    case '6' :
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


  double *dNdS = NULL;
  try {
    initFileDoubleLogNormal model_info;
    model_info.readFile(initfile, false, false);
    
    numberCountsDoubleLogNormal model;
    model_info.getModelPositions(model);

    paramSet pars(model_info.getNTot());
    model_info.getParams(pars);
    model.setParams(pars);

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
      for (unsigned int j = 0; j < nflux2-1; ++j)
	fprintf(fp,"%13.7e ",dNdS[i*nflux2+j]);
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
