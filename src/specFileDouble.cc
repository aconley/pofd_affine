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
  minbeamval = 1e-6;
  beam_histogram = true;
  nbeamhist = 150;
  fit_sigma1 = false;
  has_sigprior1 = false;
  sigprior_stdev1 = 0.0;
  fit_sigma2 = false;
  has_sigprior2 = false;
  sigprior_stdev2 = 0.0;
  exp_conf1 = 0.0;
  exp_conf2 = 0.0;
  has_cfirbprior1 = false;
  cfirbprior_mean1 = 0.0;
  cfirbprior_stdev1 = 0.0;
  has_cfirbprior2 = false;
  cfirbprior_mean2 = 0.0;
  cfirbprior_stdev2 = 0.0;
  has_wisdom_file = false;
  wisdom_file.clear();
  verbosity = 0;
  has_user_seed = false;
  seed = 0LL;
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
    throw affineExcept("specFileDouble", "readFile", errstr.str());
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
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      //Datafiles
      datafiles1.push_back( words[1] );
      datafiles2.push_back( words[2] );

      //Sigmas
      str.str(words[3]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	errstr << "Invalid (negative) sigma1 " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      sigmas1.push_back( dblval );
      str.str(words[4]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	errstr << "Invalid (negative) sigma2 " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
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
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      fit_sigma1 = utility::string_true(words[1]);

    } else if (words[0] == "fit_sigma2") {
      if (words.size() < 2) {
	errstr << "fit_sigma2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      fit_sigma2 = utility::string_true(words[1]);

    } else if (words[0] == "sigmaprior1") {
      if (words.size() < 2) {
	errstr << "sigmaprior1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) sigma prior1 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      fit_sigma1 = true;
      has_sigprior1 = true;
      sigprior_stdev1 = dblval;

    } else if (words[0] == "sigmaprior2") {
      if (words.size() < 2) {
	errstr << "sigmaprior2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) sigma prior2 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      fit_sigma2 = true;
      has_sigprior2 = true;
      sigprior_stdev2 = dblval;

    } else if (words[0] == "exp_conf1") {
      if (words.size() < 2) {
	errstr << "exp_conf1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) exp_conf1 " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      exp_conf1 = dblval;

    } else if (words[0] == "exp_conf2") {
      if (words.size() < 2) {
	errstr << "exp_conf2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) exp_conf2 " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      exp_conf2 = dblval;

    } else if (words[0] == "cfirbprior1") {

      if (words.size() < 3) {
	errstr << "cfirbprior1 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      has_cfirbprior1 = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior1 mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      cfirbprior_mean1 = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior1 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      cfirbprior_stdev1 = dblval;

    } else if (words[0] == "cfirbprior2") {

      if (words.size() < 3) {
	errstr << "cfirbprior2 line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      has_cfirbprior2 = true;

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior2 mean " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      cfirbprior_mean2 = dblval;

      str.str(words[2]); str.clear(); str >> dblval;
      if (dblval <= 0.0) {
	errstr << "Invalid (non positive) cfirb prior2 stdev " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      cfirbprior_stdev2 = dblval;

    } else if (words[0] == "bin_data") {
      if (words.size() < 2) {
	errstr << "bin_data line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "bin_data value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      
      bin_data = true;
      nbins = static_cast<unsigned int>(ival);

    } else if (words[0] == "mean_sub") {
      if (words.size() < 2) {
	errstr << "mean_sub line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      mean_sub = utility::string_true(words[1]);

    } else if (words[0] == "ignore_mask") {
      if (words.size() < 2) {
	errstr << "ignore_mask line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      ignore_mask = utility::string_true(words[1]);


    } else if (words[0] == "fftsize") {
      if (words.size() < 2) {
	errstr << "fftsize line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "fftsize value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      
      //Note we don't require this to be a power of 2
      fftsize = static_cast<unsigned int>(ival);

    } else if (words[0] == "edge_set") {
      if (words.size() < 2) {
	errstr << "edge_set line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      edge_set = utility::string_true(words[1]);

    } else if (words[0] == "nedge") {
      if (words.size() < 2) {
	errstr << "nedge line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      
      str.str(words[1]); str.clear(); str >> ival;
      
      if (ival <= 0) {
	errstr << "nedge value " << ival << " is invalid (non-positive)";
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }
      
      nedge = static_cast<unsigned int>(ival);

    } else if (words[0] == "minbeamval") {
      if (words.size() < 2) {
	errstr << "minbeamval line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> dblval;
      if (dblval < 0.0) {
	errstr << "Invalid (non-negative) minbeamval value " << dblval
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      minbeamval = dblval;

    } else if (words[0] == "beam_histogram") {
      if (words.size() < 2) {
	errstr << "beam_histogram line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      beam_histogram = utility::string_true(words[1]);

    } else if (words[0] == "nbeamhist") {
      if (words.size() < 2) {
	errstr << "nbeamhist line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> ival;
      if (ival <= 0) {
	errstr << "Invalid (non positive) nbeamhist value " << ival
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      nbeamhist = static_cast<unsigned int>(ival);

    } else if (words[0] == "wisdom_file") {
      if (words.size() < 2) {
	errstr << "wisdom_file line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      has_wisdom_file = true;
      wisdom_file = words[1];

    } else if (words[0] == "verbose") {
      if (words.size() < 2) {
	errstr << "verbose line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      if (utility::string_true(words[1])) verbosity = 1;

    } else if (words[0] == "ultraverbose") {
      if (words.size() < 2) {
	errstr << "ultraverbose line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      if (utility::string_true(words[1])) verbosity = 2;

    } else if (words[0] == "verbosity") {
      if (words.size() < 2) {
	errstr << "verbosity line doesn't have right number of entries: "
	       << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
      }

      str.str(words[1]); str.clear(); str >> ival;
      if (ival <= 0.0) {
	errstr << "Invalid (non positive) verbosity " << ival
	       << " from line: " << line;
	throw affineExcept("specFileDouble", "readFile", errstr.str());
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
      throw affineExcept("specFileDouble", "readFile", errstr.str());
    }
  }

  ifs.close();
}

/*!
  \param[in] objid HDF5 group to write to

  Writes the name of input files (datafiles, beamfiles),
  plus instrumental sigmas and likelihood normalizations.
*/
void specFileDouble::writeToHDF5Handle(hid_t objid) const {

  if (H5Iget_ref(objid) < 0)
    throw affineExcept("specFileDouble", "writeToHDF5Handle",
		       "Input handle is not valid");
  hsize_t adims;
  hid_t mems_id, att_id;

  // Set up string writing
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  const char ** catmp;

  if (datafiles1.size() > 0) {
    catmp = new const char*[datafiles1.size()];
    for (unsigned int i = 0; i < datafiles1.size(); ++i)
      catmp[i] = datafiles1[i].c_str();
    adims = static_cast<hsize_t>(datafiles1.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "datafiles1", datatype,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, datatype, catmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] catmp;
  }

  if (datafiles2.size() > 0) {
    catmp = new const char*[datafiles2.size()];
    for (unsigned int i = 0; i < datafiles2.size(); ++i)
      catmp[i] = datafiles2[i].c_str();
    adims = static_cast<hsize_t>(datafiles2.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "datafiles2", datatype,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, datatype, catmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] catmp;
  }

  if (psffiles1.size() > 0) {
    catmp = new const char*[psffiles1.size()];
    for (unsigned int i = 0; i < psffiles1.size(); ++i)
      catmp[i] = psffiles1[i].c_str();
    adims = static_cast<hsize_t>(psffiles1.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "psffiles1", datatype,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, datatype, catmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] catmp;
  }

  if (psffiles2.size() > 0) {
    catmp = new const char*[psffiles2.size()];
    for (unsigned int i = 0; i < psffiles2.size(); ++i)
      catmp[i] = psffiles2[i].c_str();
    adims = static_cast<hsize_t>(psffiles2.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "psffiles2", datatype,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, datatype, catmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] catmp;
  }

  // Numerical data
  double *dtmp;
  if (sigmas1.size() > 0) {
    dtmp = new double[sigmas1.size()];
    for (unsigned int i = 0; i < sigmas1.size(); ++i)
      dtmp[i] = sigmas1[i];
    adims = static_cast<hsize_t>(sigmas1.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "inst_sigma1", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, dtmp);
    H5Aclose(att_id);
    H5Sclose(mems_id);
    delete[] dtmp;
  }
  if (sigmas2.size() > 0) {
    dtmp = new double[sigmas2.size()];
    for (unsigned int i = 0; i < sigmas2.size(); ++i)
      dtmp[i] = sigmas2[i];
    adims = static_cast<hsize_t>(sigmas2.size());
    mems_id = H5Screate_simple(1, &adims, NULL);
    att_id = H5Acreate1(objid, "inst_sigma2", H5T_NATIVE_DOUBLE,
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
