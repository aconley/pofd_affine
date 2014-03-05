#include<fstream>
#include<algorithm>
#include<sstream>

#include "../include/specFile.h"
#include "../include/utility.h"
#include "../include/affineExcept.h"

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
  mean_sub = false;
  ignore_mask = false;
  fftsize = 131072;
  ninterp = 2048;
  minbeamval = 1e-5;
  beam_histogram = true;
  nbeamhist = 120;
  fit_sigma = false;
  has_sigprior = false;
  sigprior_stdev = 0.0;
  exp_conf = 0.0;
  has_cfirbprior = false;
  cfirbprior_mean = 0.0;
  cfirbprior_stdev = 0.0;
  has_poissonprior = false;
  poissonprior_mean = 0.0;
  poissonprior_stdev = 0.0;
  has_wisdom_file = false;
  wisdom_file.clear();
  verbosity = 0;
  seed = 0LL;
  has_user_seed = false;
}

/*!
  \param[in] flname Name of file to read
*/
void specFile::readFile(const std::string& flname) {
  init();

  std::string line;
  std::vector<std::string> words;
  std::stringstream str, errstr;
  double lnorm;

  std::ifstream ifs(flname.c_str());
  if (!ifs) {
    ifs.close();
    errstr << "Error reading spec file: " << flname;
    throw affineExcept("specFile", "readFile", errstr.str());
  }

  //Do read
  double dblval;
  int ival;
  while (!ifs.eof()) {

    std::getline(ifs, line);
    
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
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      //Datafile
      datafiles.push_back(words[1]);

      //Sigma
      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	errstr << "Invalid (negative) sigma " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
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
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "bin_data value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      
      bin_data = true;
      nbins = static_cast<unsigned int>(ival);

    } else if (words[0] == "mean_sub") {
      if (words.size() < 2) {
	errstr << "mean_sub line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      mean_sub = utility::string_true(words[1]);

    } else if (words[0] == "ignore_mask") {
      if (words.size() < 2) {
	errstr << "ignore_mask line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      ignore_mask = utility::string_true(words[1]);

    } else if (words[0] == "fftsize") {
      if (words.size() < 2) {
	errstr << "fftsize line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "fftsize value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      
      //Note we don't require this to be a power of 2
      fftsize = static_cast<unsigned int>(ival);

    } else if (words[0] == "ninterp") {
      if (words.size() < 2) {
	errstr << "ninterp line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "ninterp value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      
      ninterp = static_cast<unsigned int>(ival);

    } else if (words[0] == "minbeamval") {
      if (words.size() < 2) {
	errstr << "minbeamval line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	errstr << "Invalid (non-negative) minbeamval value " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      minbeamval = dblval;

    } else if (words[0] == "beam_histogram") {
      if (words.size() < 2) {
	errstr << "beam_histogram line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      beam_histogram = utility::string_true(words[1]);

    } else if (words[0] == "nbeamhist") {
      if (words.size() < 2) {
	errstr << "nbeamhist line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> ival;
      if (ival <= 0) {
	errstr << "Invalid (non positive) nbeamhist value " << ival
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      nbeamhist = static_cast<unsigned int>(ival);

    } else if (words[0] == "fit_sigma") {
      if (words.size() < 2) {
	errstr << "fit_sigma line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      fit_sigma = utility::string_true(words[1]);

    } else if (words[0] == "sigmaprior") {
      if (words.size() < 2) {
	errstr << "sigmaprior line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) sigma prior stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      has_sigprior = true;
      fit_sigma = true;
      sigprior_stdev = dblval;

    } else if (words[0] == "exp_conf") {
      if (words.size() < 2) {
	errstr << "exp_conf line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) expected confusion noise " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      exp_conf = dblval;

    } else if (words[0] == "cfirbprior") {

      if (words.size() < 3) {
	errstr << "cfirbprior line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      has_cfirbprior = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      cfirbprior_mean = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      cfirbprior_stdev = dblval;

    } else if (words[0] == "poissonprior") {

      if (words.size() < 3) {
	errstr << "poissonprior line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      has_poissonprior = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) poisson prior mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      poissonprior_mean = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) poisson prior stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }
      poissonprior_stdev = dblval;

    } else if (words[0] == "wisdom_file") {
      if (words.size() < 2) {
	errstr << "wisdom_file line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      has_wisdom_file = true;
      wisdom_file = words[1];

    } else if (words[0] == "verbose") {
      if (words.size() < 2) {
	errstr << "verbose line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      if (utility::string_true(words[1])) verbosity = 1;

    } else if (words[0] == "ultraverbose") {
      if (words.size() < 2) {
	errstr << "ultraverbose line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      if (utility::string_true(words[1])) verbosity = 2;

    } else if (words[0] == "verbosity") {
      if (words.size() < 2) {
	errstr << "verbosity line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> ival;
      if (ival <= 0) {
	errstr << "Invalid (non positive) verbosity " << ival
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      verbosity = static_cast<unsigned int>(ival);

    } else if (words[0] == "seed") {
      if (words.size() < 2) {
	errstr << "seed line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> ival;
      if (ival <= 0) {
	errstr << "Invalid (non positive) seed " << ival
	       << " from line: " << line;
	throw affineExcept("specFile", "readFile", errstr.str());
      }

      seed = static_cast<unsigned long long int>(ival);
      has_user_seed = true;

    } else {
      errstr << "Couldn't determine line type for: " << line;
      throw affineExcept("specFile", "readFile", errstr.str());
    }
  }

  ifs.close();
}

/*!
  \param[in] objid HDF5 group to write to

  Writes the name of input files (datafiles, beamfiles),
  plus instrumental sigmas and likelihood normalizations
*/
void specFile::writeToHDF5Handle(hid_t objid) const {

  if (H5Iget_ref(objid) < 0)
    throw affineExcept("specFile", "writeToHDF5Handle",
		       "Input handle is not valid");

  hsize_t adims;
  hid_t mems_id, att_id;

  // Set up string writing
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  const char ** catmp;

  if (datafiles.size() > 0) {
    catmp = new const char*[datafiles.size()];
    for (unsigned int i = 0; i < datafiles.size(); ++i)
      catmp[i] = datafiles[i].c_str();
    adims = static_cast<hsize_t>(datafiles.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "datafiles", datatype,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, datatype, catmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] catmp;
  }

  if (psffiles.size() > 0) {
    catmp = new const char*[psffiles.size()];
    for (unsigned int i = 0; i < psffiles.size(); ++i)
      catmp[i] = psffiles[i].c_str();
    adims = static_cast<hsize_t>(psffiles.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "psffiles", datatype,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, datatype, catmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] catmp;
  }

  // Numerical data
  double *dtmp;
  if (sigmas.size() > 0) {
    dtmp = new double[sigmas.size()];
    for (unsigned int i = 0; i < sigmas.size(); ++i)
      dtmp[i] = sigmas[i];
    adims = static_cast<hsize_t>(sigmas.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "inst_sigma", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, dtmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] dtmp;
  }

  if (like_norm.size() > 0) {
    dtmp = new double[like_norm.size()];
    for (unsigned int i = 0; i < like_norm.size(); ++i)
      dtmp[i] = like_norm[i];
    adims = static_cast<hsize_t>(like_norm.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "like_norm", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, dtmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] dtmp;
  }
}
