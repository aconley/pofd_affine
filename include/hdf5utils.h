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
  /*! \brief Output file types */
  enum outfiletype { UNKNOWN=0, TXT=1, FITS=2, HDF5=3 };

  /*! \brief Determine file type from extension */
  outfiletype getOutputFileType(const std::string& str);

  // Attribute writers

  /*! \brief Write single string as attribute */
  void writeAttString(hid_t, const std::string&, const std::string&);
  /*! \brief Write array of strings as attribute */
  void writeAttStrings(hid_t, const std::string&, unsigned int,
		       const std::string* const);
  /*! \brief Write vector of strings as attribute */
  void writeAttStrings(hid_t, const std::string&, 
		       const std::vector<std::string>&);

  /*! \brief Write array of unsigned ints as attribute */
  void writeAttUnsignedInts(hid_t, const std::string&, unsigned int,
			    const unsigned int* const);
  /*! \brief Write vector of unsigned ints as attribute */
  void writeAttUnsignedInts(hid_t, const std::string&, 
			    const std::vector<unsigned int>&);

  /*! \brief Write single boolean as attribute */
  void writeAttBool(hid_t, const std::string&, bool);
  /*! \brief Write array of booleans as attribute */
  void writeAttBools(hid_t, const std::string&, unsigned int,
		     const bool* const);
  /*! \brief Write vector of booleans as attribute */
  void writeAttBools(hid_t, const std::string&, 
		     const std::vector<bool>&);

  /*! \brief Write array of floats as attribute */
  void writeAttFloats(hid_t, const std::string&, unsigned int,
		      const float* const);
  /*! \brief Write vector of floats as attribute */
  void writeAttFloats(hid_t, const std::string&, 
		      const std::vector<float>&);

  /*! \brief Write array of doubles as attribute */
  void writeAttDoubles(hid_t, const std::string&, unsigned int,
		       const double* const);
  /*! \brief Write vector of doubles as attribute */
  void writeAttDoubles(hid_t, const std::string&, 
		       const std::vector<double>&);

  // Data writers
  // 1D
  /*! \brief Write 1D array of unsigned ints as data */
  void writeDataUnsignedInts(hid_t, const std::string&, unsigned int,
			     const unsigned int* const);
  /*! \brief Write 1D array of floats as data */
  void writeDataFloats(hid_t, const std::string&, unsigned int,
			const float* const);
  /*! \brief Write 1D array of doubles as data */
  void writeDataDoubles(hid_t, const std::string&, unsigned int,
			const double* const);
  
  // 2D
  /*! \brief Write 2D array of floats as data */
  void writeData2DFloats(hid_t, const std::string&, unsigned int,
			 unsigned int, const float* const);
  /*! \brief Write 2D array of doubles as data */
  void writeData2DDoubles(hid_t, const std::string&, unsigned int,
			  unsigned int, const double* const);
}

#endif
