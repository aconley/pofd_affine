#include<fstream>
#include<algorithm>
#include<sstream>

#include<specFileDouble.h>
#include<utility.h>
#include<affineExcept.h>

specFileDouble::specFileDouble() : 
  fit_sigma1(false), fit_sigma2(false),
  has_sigprior1(false), sigprior_stdev1(0.0),
  has_sigprior2(false), sigprior_stdev2(0.0),
  has_cfirbprior1(false), cfirbprior_mean1(0.0), cfirbprior_stdev1(0.0),
  has_cfirbprior2(false), cfirbprior_mean2(0.0), cfirbprior_stdev2(0.0) {}

/*!
  \param[in] flname Name of file to read
 */
specFileDouble::specFileDouble(const std::string& flname) : 
  fit_sigma1(false), fit_sigma2(false),
  has_sigprior1(false), sigprior_stdev1(0.0),
  has_sigprior2(false), sigprior_stdev2(0.0),
  has_cfirbprior1(false), cfirbprior_mean1(0.0), cfirbprior_stdev1(0.0),
  has_cfirbprior2(false), cfirbprior_mean2(0.0), cfirbprior_stdev2(0.0) {
  
  readFile(flname);
}

/*!
  \param[in] flname Name of file to read
 */
void specFileDouble::readFile(const std::string& flname) {
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double lnorm;

  std::ifstream ifs( flname.c_str() );
  if ( ! ifs ) {
    ifs.close();
    std::stringstream errstr;
    errstr << "Error reading spec file: " << flname;
    throw affineExcept("specFileDouble","readFile",errstr.str(),1);
  }

  //Do read
  double dblval;
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
	std::stringstream errstr;
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
	std::stringstream errstr;
	errstr << "Invalid (negative) sigma1 " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 4);
      }
      sigmas1.push_back( dblval );
      str.str(words[4]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (negative) sigma2 " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 8);
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
	std::stringstream errstr;
	errstr << "fit_sigma1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 16834);
      }

      fit_sigma1 = utility::string_true(words[1]);
    } else if (words[0] == "fit_sigma2") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "fit_sigma2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str(), 32768);
      }

      fit_sigma2 = utility::string_true(words[1]);

    } else if (words[0] == "sigmaprior1") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "sigmaprior1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 16);
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) sigma prior1 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 32);
      }

      fit_sigma1 = true;
      has_sigprior1 = true;
      sigprior_stdev1 = dblval;
    } else if (words[0] == "sigmaprior2") {
      if (words.size() < 2) {
	std::stringstream errstr;
	errstr << "sigmaprior2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 32);
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) sigma prior2 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 64);
      }

      fit_sigma2 = true;
      has_sigprior2 = true;
      sigprior_stdev2 = dblval;

    } else if (words[0] == "cfirbprior1") {

      if (words.size() < 3) {
	std::stringstream errstr;
	errstr << "cfirbprior1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 128);
      }

      has_cfirbprior1 = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) cfirb prior1 mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 256);
      }
      cfirbprior_mean1 = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) cfirb prior1 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 512);
      }
      cfirbprior_stdev1 = dblval;

    } else if (words[0] == "cfirbprior2") {

      if (words.size() < 3) {
	std::stringstream errstr;
	errstr << "cfirbprior2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 1024);
      }

      has_cfirbprior2 = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) cfirb prior2 mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 2048);
      }
      cfirbprior_mean2 = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	std::stringstream errstr;
	errstr << "Invalid (non positive) cfirb prior2 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str(), 4096);
      }
      cfirbprior_stdev2 = dblval;


    } else {
      std::stringstream errstr;
      errstr << "Couldn't determine line type for: " << line;
      throw affineExcept("specFileDouble","readFile",errstr.str(),8192);
    }
  }

  ifs.close();
}
