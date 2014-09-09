//hdf5utils.h

#ifndef __hdf5utils__
#define __hdf5utils__

#include<string>
#include<vector>

#include "hdf5.h"

/*!
  \brief Utility functions for HDF5 input/output
*/
namespace hdf5utils {
  void writeAttString(hid_t, const std::string&, const std::string&);
  void writeAttStrings(hid_t, const std::string&, unsigned int,
		       const std::string* const);
  void writeAttStrings(hid_t, const std::string&, 
		       const std::vector<std::string>&);

  void writeAttUnsignedInts(hid_t, const std::string&, unsigned int,
			    const unsigned int* const);
  void writeAttUnsignedInts(hid_t, const std::string&, 
			    const std::vector<unsigned int>&);

  void writeAttBool(hid_t, const std::string&, bool);
  void writeAttBools(hid_t, const std::string&, unsigned int,
		     const bool* const);
  void writeAttBools(hid_t, const std::string&, 
		     const std::vector<bool>&);

  void writeAttFloats(hid_t, const std::string&, unsigned int,
		      const float* const);
  void writeAttFloats(hid_t, const std::string&, 
		      const std::vector<float>&);

  void writeAttDoubles(hid_t, const std::string&, unsigned int,
		       const double* const);
  void writeAttDoubles(hid_t, const std::string&, 
		       const std::vector<double>&);

}

#endif
