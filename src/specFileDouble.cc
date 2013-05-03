#include<fstream>
#include<algorithm>
#include<sstream>

#include "../include/specFileDouble.h"
#include "../include/utility.h"
#include "../include/affineExcept.h"

specFileDouble::specFileDouble() {
  init();
}

/*!
  \param[in] flname Name of file to read
 */
specFileDouble::specFileDouble(const std::string& flname) {
  readFile(flname); //Also calls init
}

void specFileDouble::init() {
  datafiles1.clear();
  datafiles2.clear();
  sigmas1.clear();
  sigmas2.clear();
  psffiles1.clear();
  psffiles2.clear();
  like_norm.clear();
  bin_data = false;
  ignore_mask = false;
  nbins = 0;
  mean_sub = false;
  fftsize = 4096;
  edge_set = true;
  nedge = 256;
  edge_fix = true;
  beam_histogram = true;
  fit_sigma1 = false;
  has_sigprior1 = false;
  sigprior_stdev1 = 0.0;
  fit_sigma2 = false;
  has_sigprior2 = false;
  sigprior_stdev2 = 0.0;
  has_cfirbprior1 = false;
  cfirbprior_mean1 = 0.0;
  cfirbprior_stdev1 = 0.0;
  has_cfirbprior2 = false;
  cfirbprior_mean2 = 0.0;
  cfirbprior_stdev2 = 0.0;
  has_wisdom_file = false;
  wisdom_file.clear();
  verbose = false;
  ultraverbose = false;
}


/*!
  \param[in] flname Name of file to read
 */
void specFileDouble::readFile(const std::string& flname) {
  std::string line;
  std::vector<std::string> words;
  std::stringstream str, errstr;
  double lnorm;

  init();

  std::ifstream ifs( flname.c_str() );
  if ( ! ifs ) {
    ifs.close();
    errstr << "Error reading spec file: " << flname;
    throw affineExcept("specFileDouble","readFile",errstr.str(),1);
  }

  //Do read
  double dblval;
  int ival;
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

      if (words.size() < 7) {
	errstr << "dataset line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 2);
      }

      //Datafiles
      datafiles1.push_back( words[1] );
      datafiles2.push_back( words[2] );

      //Sigmas
      str.str(words[3]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	errstr << "Invalid (negative) sigma1 " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 3);
      }
      sigmas1.push_back( dblval );
      str.str(words[4]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	errstr << "Invalid (negative) sigma2 " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 4);
      }
      sigmas2.push_back( dblval );
    
      //Psffiles
      psffiles1.push_back( words[5] );
      psffiles2.push_back( words[6] );

      //Optional like norm
      if (words.size() >= 8) {
	str.str(words[7]); str.clear(); str >> lnorm;
	like_norm.push_back(lnorm);
      } else like_norm.push_back(1.0);
    } else if (words[0] == "fit_sigma1") {
      if (words.size() < 2) {
	errstr << "fit_sigma1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 5);
      }

      fit_sigma1 = utility::string_true(words[1]);
    } else if (words[0] == "fit_sigma2") {
      if (words.size() < 2) {
	errstr << "fit_sigma2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 6);
      }

      fit_sigma2 = utility::string_true(words[1]);

    } else if (words[0] == "sigmaprior1") {
      if (words.size() < 2) {
	errstr << "sigmaprior1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 7);
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) sigma prior1 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 8);
      }

      fit_sigma1 = true;
      has_sigprior1 = true;
      sigprior_stdev1 = dblval;
    } else if (words[0] == "sigmaprior2") {
      if (words.size() < 2) {
	errstr << "sigmaprior2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 9);
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) sigma prior2 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 10);
      }

      fit_sigma2 = true;
      has_sigprior2 = true;
      sigprior_stdev2 = dblval;

    } else if (words[0] == "cfirbprior1") {

      if (words.size() < 3) {
	errstr << "cfirbprior1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 11);
      }

      has_cfirbprior1 = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior1 mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 12);
      }
      cfirbprior_mean1 = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior1 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 13);
      }
      cfirbprior_stdev1 = dblval;

    } else if (words[0] == "cfirbprior2") {

      if (words.size() < 3) {
	errstr << "cfirbprior2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 14);
      }

      has_cfirbprior2 = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior2 mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 15);
      }
      cfirbprior_mean2 = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior2 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 16);
      }
      cfirbprior_stdev2 = dblval;

    } else if (words[0] == "bin_data") {
      if (words.size() < 2) {
	errstr << "bin_data line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 17);
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "bin_data value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFile", "readFile", errstr.str(), 18);
      }
      
      bin_data = true;
      nbins = static_cast<unsigned int>(ival);
    } else if (words[0] == "mean_sub") {
      if (words.size() < 2) {
	errstr << "mean_sub line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 19);
      }

      mean_sub = utility::string_true(words[1]);

    } else if (words[0] == "ignore_mask") {
      if (words.size() < 2) {
	errstr << "ignore_mask line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 8);
      }

      ignore_mask = utility::string_true(words[1]);


    } else if (words[0] == "fftsize") {
      if (words.size() < 2) {
	errstr << "fftsize line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 20);
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "fftsize value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFile", "readFile", errstr.str(), 21);
      }
      
      //Note we don't require this to be a power of 2
      fftsize = static_cast<unsigned int>(ival);

    } else if (words[0] == "edge_set") {
      if (words.size() < 2) {
	errstr << "edge_set line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 22);
      }

      edge_set = utility::string_true(words[1]);

    } else if (words[0] == "nedge") {
      if (words.size() < 2) {
	errstr << "nedge line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 23);
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "nedge value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFile", "readFile", errstr.str(), 24);
      }
      
      nedge = static_cast<unsigned int>(ival);

    } else if (words[0] == "edge_fix") {
      if (words.size() < 2) {
	errstr << "edge_fix line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 25);
      }

      edge_fix = utility::string_true(words[1]);

    } else if (words[0] == "beam_histogram") {
      if (words.size() < 2) {
	errstr << "beam_histogram line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 26);
      }

      beam_histogram = utility::string_true(words[1]);

    } else if (words[0] == "wisdom_file") {
      if (words.size() < 2) {
	errstr << "wisdom_file line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 27);
      }

      has_wisdom_file = true;
      wisdom_file = words[1];

    } else if (words[0] == "verbose") {
      if (words.size() < 2) {
	errstr << "verbose line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 28);
      }

      verbose = utility::string_true(words[1]);

    } else if (words[0] == "ultraverbose") {
      if (words.size() < 2) {
	errstr << "ultraverbose line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 29);
      }

      ultraverbose = utility::string_true(words[1]);

    } else {
      errstr << "Couldn't determine line type for: " << line;
      throw affineExcept("specFileDouble", "readFile", errstr.str(), 30);
    }
  }

  ifs.close();
}
