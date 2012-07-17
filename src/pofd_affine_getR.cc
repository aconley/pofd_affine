#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

#include<getopt.h>

#include<global_settings.h>
#include<beam.h>
#include<numberCountsKnotsSpline.h>
#include<paramSet.h>
#include<utility.h>
#include<affineExcept.h>

int main( int argc, char** argv ) {

  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  std::string psffile; //Beam file
  bool histogram, posonly, negonly; //Histogram beam
  double minflux, maxflux;
  unsigned int nflux;

  histogram = false;
  posonly   = false;
  negonly   = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"histogram",no_argument,0,'H'},
    {"negonly",no_argument,0,'n'},
    {"posonly",no_argument,0,'p'},
    {"version",no_argument,0,'V'},
    {0,0,0,0}
  };

  char optstring[] = "hHnpV";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getR -- get R for a number counts" 
		<< std::endl;
      std::cerr << "\tmodel" << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tpofd_mcmc_getR [options] minflux maxflux nflux initfile"
		<< std::endl;
      std::cerr<< "\t beamfile outfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates R for the model in initfile using the P(D)"
		<< std::endl;
      std::cerr << "\tformalism and writes it to outfile." 
		     << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tinitfile is a text file giving the positions of the"
		<< "knot points," << std::endl;
      std::cerr << "\ttheir initial values, and their estimated errors in"
		<< std::endl;
      std::cerr << "\tthe format knotflux value error."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tWrites positive and negative beam parts." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tminflux, maxflux, and nflux describe the flux values"
		<< std::endl;
      std::cerr << "\tto calculate R for." << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-H, --histogram" << std::endl;
      std::cerr << "\t\tUse beam histogramming" << std::endl;
      std::cerr << "\t-n, --negonly" << std::endl;
      std::cerr << "\t\tOnly use the negative parts of the beam (if they exist)"
		<< std::endl;
      std::cerr << "\t-p, --posonly" << std::endl;
      std::cerr << "\t\tOnly use the positive parts of the beam (if they exist)"
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput version number and exit" << std::endl;
      return 0;
      break;
    case 'H' :
      histogram = true;
      break;
    case 'n' :
      negonly = true;
      break;
    case 'p' :
      posonly = true;
      break;
    case 'V' :
      std::cerr << "pofd_mcmc version number: " << pofd_mcmc::version 
		<< std::endl;
      return 0;
      break;
    }

  if (optind >= argc-5 ) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  minflux = atof( argv[optind] );
  maxflux = atof( argv[optind+1] );
  nflux = static_cast<unsigned int>( atoi(argv[optind+2]) );
  initfile = std::string( argv[optind+3] );
  psffile = std::string( argv[optind+4] );
  outfile = std::string( argv[optind+5] );

  if (posonly && negonly) {
    std::cerr << "Can't set both posonly and negonly" << std::endl;
    return 1;
  }

  double dflux = (maxflux-minflux)/static_cast<double>(nflux-1);

  //Read in the params
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

  double *fluxes = NULL;
  double *R = NULL;
  try {
    numberCountsKnotsSpline model(knotpos);
    beam bm( psffile, histogram );
    paramSet pars( knotval );
    model.setParams( pars );

    numberCounts::rtype rt = numberCounts::BEAMBOTH;
    if (posonly)
      rt = numberCounts::BEAMPOS;
    else if (negonly) 
      rt = numberCounts::BEAMNEG;
    
    std::cout << "Mean flux per sq degree: " << model.getMeanFluxPerArea()
	      << std::endl;
    std::cout << "Beam area: " << bm.getEffectiveArea() << std::endl;
    if (bm.hasPos()) std::cout << " Beam has positive components" << std::endl;
    if (bm.hasNeg()) std::cout << " Beam has negative components" << std::endl;
    fluxes = new double[nflux];
    for (unsigned int i = 0; i < nflux; ++i)
      fluxes[i] = minflux + static_cast<double>(i)*dflux;

    R = new double[nflux];
    model.getR(nflux,fluxes,bm,R,rt);
    
    FILE *fp;
    fp = fopen( outfile.c_str(),"w");
    if (!fp) {
      std::cerr << "Failed to open output file" << std::endl;
      return 128;
    }
    fprintf(fp,"%12s   %12s\n","Flux","R");
    for (unsigned int i = 0; i < nflux; ++i) 
      fprintf(fp,"%12.6e   %15.9e\n",fluxes[i],R[i]);
    fclose(fp);
    
    delete[] fluxes;
    delete[] R; 
  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (fluxes != NULL) delete[] fluxes;
    if (R != NULL) delete[] R;
    return 8;
  }

  return 0;

}
