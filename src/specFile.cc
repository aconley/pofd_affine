#include<fstream>
#include<algorithm>
#include<sstream>

#include<specFile.h>
#include<utility.h>
#include<affineExcept.h>

specFile::specFile() { init(); }

/*!
  \param[in] flname Name of file to read
 */
specFile::specFile(const std::string& flname) {
  readFile(flname); //Calls init as well
}

void specFile::init() {
  datafiles.clear();
  sigmas.clear();
  psffiles.clear();
  like_norm.clear();
  bin_data = false;
  nbins = 0;
  meansub = false;
  fftsize = 131072;
  beam_histogram = false;
  fit_sigma = false;
  has_sigprior = false;
  sigprior_stdev = 0.0;
  has_cfirbprior = false;
  cfirbprior_mean = 0.0;
  cfirbprior_stdev = 0.0;
  has_wisdom_file = false;
  wisfile.clear();
  verbose = false;
  ultraverbose = false;
}

/*!
  \param[in] flname Name of file to read
 */
void specFile::readFile(const std::string& flname) {
  init();

  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double lnorm;

  std::ifstream ifs( flname.c_str() );
  if ( ! ifs ) {
    ifs.close();
    std::stringstream errstr;
    errstr << "Error reading spec file: " << flname;
    throw affineExcept("specFile","readFile",errstr.str(),1);
  }

  //Do read
  double dblval;
  int ival;
  std::stringstream errstr;
  while (! ifs.eof() ) {

    std::getline(ifs,line);
    
    //Comment line
    if (line[0] == '#') continue;

    //Break up line
    utility::stringwords_eq(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line, leading spaces removed

    //Figure out what type of line
    if (words[0] == "dataset") {
      if (words.size() < 4) {
	errstr << "dataset line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 2);
      }

      //Datafile
      datafiles.push_back( words[1] );

      //Sigma
      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	errstr << "Invalid (negative) sigma " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 4);
      }
      sigmas.push_back( dblval );
    
      //Psffile
      psffiles.push_back( words[3] );

      //Optional like norm
      if (words.size() >= 5) {
	str.str(words[4]); str.clear(); str >> lnorm;
	like_norm.push_back(lnorm);
      } else like_norm.push_back(1.0);
    } else if (words[0] == "bin_data") {
      if (words.size() < 2) {
	errstr << "bin_data line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 2048);
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "bin_data value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFile", "readFile", errstr.str(), 4096);
      }
      
      bin_data = true;
      nbins = static_cast<unsigned int>(ival);
    } else if (words[0] == "meansub") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "meansub line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 8192);
      }

      meansub = utility::string_true(words[1]);
    } else if (words[0] == "fftsize") {
      if (words.size() < 2) {
	errstr << "fftsize line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 16384);
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "fftsize value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFile", "readFile", errstr.str(), 32768);
      }
      
      //Note we don't require this to be a power of 2
      fftsize = static_cast<unsigned int>(ival);

    } else if (words[0] == "beam_histogram") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "beam_histogram line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 65536);
      }

      beam_histogram = utility::string_true(words[1]);

    } else if (words[0] == "fit_sigma") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "fit_sigma line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 512);
      }

      fit_sigma = utility::string_true(words[1]);

    } else if (words[0] == "sigmaprior") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "sigmaprior line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 8);
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) sigma prior stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 16);
      }

      has_sigprior = true;
      fit_sigma = true;
      sigprior_stdev = dblval;

    } else if (words[0] == "cfirbprior") {

      if (words.size() < 3) {
	std::stringstream errstr;
	errstr << "cfirbprior line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 32);
      }

      has_cfirbprior = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) cfirb prior mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 64);
      }
      cfirbprior_mean = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) cfirb prior stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 128);
      }
      cfirbprior_stdev = dblval;

    } else if (words[0] == "wisdom_file") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "wisdom_file line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 131072);
      }

      has_wisdom_file = true;
      wisfile = words[1];

    } else if (words[0] == "verbose") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "verbose line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 262144);
      }

      verbose = utility::string_true(words[1]);

    } else if (words[0] == "ultraverbose") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "ultraverbose line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 524288);
      }

      ultraverbose = utility::string_true(words[1]);

    } else {
      std::stringstream errstr;
      errstr << "Couldn't determine line type for: " << line;
      throw affineExcept("specFile","readFile",errstr.str(),256);
    }
  }

  ifs.close();
}
