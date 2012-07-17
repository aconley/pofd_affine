#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

#include<getopt.h>

#include<global_settings.h>
#include<doublebeam.h>
#include<numberCountsDoubleLogNormal.h>
#include<paramSet.h>
#include<utility.h>
#include<affineExcept.h>

int main( int argc, char** argv ) {

  std::string initfile; //Init file (having model we want)
  std::string outfile; //File to write to
  std::string psffile1, psffile2; //Beam file
  bool histogram, posonly; //Histogram beam, use only positive part
  double minflux1, maxflux1, minflux2, maxflux2;
  unsigned int nflux1, nflux2;

  histogram = false;
  posonly   = false;

  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"histogram",no_argument,0,'H'},
    {"posonly",no_argument,0,'p'},
    {"version",no_argument,0,'V'},
    {0,0,0,0}
  };

  char optstring[] = "hHpV";
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_affine_getDoubleR -- get R for a number counts" 
		<< std::endl;
      std::cerr << "\tmodel consisting of a spline in the first dimension" 
		<< std::endl;
      std::cerr << "\tand a log-normal distribution in the second."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tpofd_mcmc_getDoubleR [options] minflux1 maxflux1 nflux1"
		<< std::endl;
      std::cerr << "\t minflux2 minflux2 initfile beamfile1 beamfile2 outfile"
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates R for the model in initfile using the P(D1,D2)"
		<< std::endl;
      std::cerr << "\tformalism and writes it to outfile." 
		     << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tminflux[12], maxflux[12], and nflux[12] describe the "
		<< "flux values" << std::endl;
      std::cerr << "\tto calculate R for in band 1 and band 2, respectively." 
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tinitfile is a text file giving the positions of the"
                << std::endl;
      std::cerr << "\tknot points and their initial values, followed by the "
		<< "sigma" << std::endl;
      std::cerr << "\tknot points and their values, then likewise for the "
		<< "colour" << std::endl;
      std::cerr << "\toffset.  The format is three numbers on the first line, "
		<< "giving" << std::endl;
      std::cerr << "\tthe number of number count knots, sigma knots, and "
		<< "offset knots," << std::endl;
      std::cerr << "\tfollowed by a number of lines with the format"
                << std::endl;
      std::cerr << "\tknotpos value." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tbeamfile1 and beamfile2 are FITS files containing the"
		<< " beam in" << std::endl;
      std::cerr << "\teach band.  They must have the same pixel size."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-H, --histogram" << std::endl;
      std::cerr << "\t\tUse beam histogramming" << std::endl;
      std::cerr << "\t-p, --posonly" << std::endl;
      std::cerr << "\t\tOnly use the positive-positive parts of the beams "
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput version number and exit" << std::endl;
      return 0;
      break;
    case 'H' :
      histogram = true;
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

  if (optind >= argc-9 ) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  minflux1 = atof( argv[optind] );
  maxflux1 = atof( argv[optind+1] );
  nflux1   = static_cast<unsigned int>( atoi(argv[optind+2]) );
  minflux2 = atof( argv[optind+3] );
  maxflux2 = atof( argv[optind+4] );
  nflux2   = static_cast<unsigned int>( atoi(argv[optind+5]) );
  initfile = std::string( argv[optind+6] );
  psffile1 = std::string( argv[optind+7] );
  psffile2 = std::string( argv[optind+8] );
  outfile  = std::string( argv[optind+9] );

  if (nflux1 == 0) {
    std::cerr << "Invalid (non-positive) nflux1" << std::endl;
    return 1;
  }
  if (nflux2 == 0) {
    std::cerr << "Invalid (non-positive) nflux2" << std::endl;
    return 1;
  }
  double dflux1, dflux2;
  if (maxflux1 < minflux1) {
    double tmp = minflux1;
    minflux1 = maxflux1; 
    maxflux1=tmp;
  }
  if (nflux1 > 1) 
    dflux1 = (maxflux1-minflux1)/static_cast<double>(nflux1-1);
  else
    dflux1 = maxflux1-minflux1;
  if (maxflux2 < minflux2) {
    double tmp = minflux2;
    minflux2 = maxflux2; 
    maxflux2=tmp;
  }
  if (nflux2 > 1) 
    dflux2 = (maxflux2-minflux2)/static_cast<double>(nflux2-1);
  else
    dflux2 = maxflux2-minflux2;
  
  //Read in the params
  std::ifstream initfs( initfile.c_str() );
  if (!initfs) {
    std::cerr << "Error readining in initialization file: "
	      << initfile << std::endl;
    return 1;
  }
  unsigned int nk, ns, no;
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

  //Main computation
  double *fluxes1, *fluxes2, *R;
  fluxes1 = fluxes2 = R = NULL;
  try {
    numberCountsDoubleLogNormal model(knotpos,sigmapos,offsetpos);
    doublebeam bm( psffile1, psffile2, histogram );
    paramSet pars( wvec2 );
    model.setParams( pars );

    numberCountsDouble::rtype rt; 
    if (posonly)
      rt = numberCountsDouble::BEAMPOS;
    else
      rt = numberCountsDouble::BEAMALL;
    
    std::cout << "Mean flux per sq degree, band 1: " 
	      << model.getMeanFluxPerArea(0)
	      << std::endl;
    std::cout << "Beam area, band 1: " << bm.getEffectiveArea1() << std::endl;
    std::cout << "Mean flux per sq degree, band 2: " 
	      << model.getMeanFluxPerArea(1)
	      << std::endl;
    std::cout << "Beam area, band 2: " << bm.getEffectiveArea2() << std::endl;
    fluxes1 = new double[nflux1];
    for (unsigned int i = 0; i < nflux1; ++i)
      fluxes1[i] = minflux1 + static_cast<double>(i)*dflux1;
    fluxes2 = new double[nflux2];
    for (unsigned int i = 0; i < nflux2; ++i)
      fluxes2[i] = minflux2 + static_cast<double>(i)*dflux2;

    R = new double[nflux1*nflux2];
    model.getR(nflux1,fluxes1,nflux2,fluxes2,bm,R,rt);
    
    FILE *fp;
    fp = fopen( outfile.c_str(),"w");
    if (!fp) {
      std::cerr << "Failed to open output file" << std::endl;
      return 128;
    }
    fprintf(fp,"%4u %4u\n",nflux1,nflux2);
    fprintf(fp,"minflux1: %12.6e dflux1: %12.6e\n",minflux1,dflux1);
    fprintf(fp,"minflux2: %12.6e dflux2: %12.6e\n",minflux2,dflux2);
    for (unsigned int i = 0; i < nflux1; ++i) {
      for (unsigned int j = 0; j < nflux2-1; ++j)
        fprintf(fp,"%13.7e ",R[nflux2*i+j]);
      fprintf(fp,"%13.7e\n",R[nflux2*i+nflux2-1]);
    }
    fclose(fp);
    
    delete[] fluxes1;
    delete[] fluxes2;
    delete[] R; 
  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (fluxes1 != NULL) delete[] fluxes1;
    if (fluxes2 != NULL) delete[] fluxes2;
    if (R != NULL) delete[] R;
    return 16;
  }

  return 0;

}
