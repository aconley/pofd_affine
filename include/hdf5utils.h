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

  /*! \brief Read a single string as an attribute */
  std::string readAttString(hid_t, const std::string&);
  /*! \brief Read array of strings as an attribute */
  std::vector<std::string> readAttStrings(hid_t, const std::string&);

  /*! \brief Write a single unsigned int as an attribute */
  void writeAttUnsignedInt(hid_t, const std::string&, unsigned int);
  /*! \brief Write array of unsigned ints as attribute */
  void writeAttUnsignedInts(hid_t, const std::string&, unsigned int,
			    const unsigned int* const);
  /*! \brief Write vector of unsigned ints as attribute */
  void writeAttUnsignedInts(hid_t, const std::string&, 
			    const std::vector<unsigned int>&);
  /*! \brief Read unsigned int attribute */
  unsigned int readAttUnsignedInt(hid_t objid, const std::string& name);
  /*! \brief Read vector of unsigned ints attribute */
  std::vector<unsigned int> readAttUnsignedInts(hid_t objid, 
                                                const std::string& name);

  /*! \brief Write a single int as an attribute */
  void writeAttInt(hid_t, const std::string&, int);
  /*! \brief Write array of ints as attribute */
  void writeAttInts(hid_t, const std::string&, unsigned int,
		    const int* const);
  /*! \brief Write vector of ints as attribute */
  void writeAttInts(hid_t, const std::string&, const std::vector<int>&);

  /*! \brief Write single boolean as attribute */
  void writeAttBool(hid_t, const std::string&, bool);
  /*! \brief Write array of booleans as attribute */
  void writeAttBools(hid_t, const std::string&, unsigned int,
		     const bool* const);
  /*! \brief Write vector of booleans as attribute */
  void writeAttBools(hid_t, const std::string&, 
		     const std::vector<bool>&);
  /*! \brief Read bool attribute */
  bool readAttBool(hid_t objid, const std::string& name);
  /*! \brief Read vector of booleans attribute */
  std::vector<bool> readAttBools(hid_t objid, const std::string& name);

  /*! \brief Write single float as attribute */
  void writeAttFloat(hid_t, const std::string&, float);
  /*! \brief Write array of floats as attribute */
  void writeAttFloats(hid_t, const std::string&, unsigned int,
		      const float* const);
  /*! \brief Write vector of floats as attribute */
  void writeAttFloats(hid_t, const std::string&, 
		      const std::vector<float>&);
  /*! \brief Read float attribute */
  float readAttFloat(hid_t objid, const std::string& name);
  /*! \brief Read vector of floats attribute */
  std::vector<float> readAttFloats(hid_t objid, const std::string& name);

  /*! \brief Write single double as attribute */
  void writeAttDouble(hid_t, const std::string&, double);
  /*! \brief Write array of doubles as attribute */
  void writeAttDoubles(hid_t, const std::string&, unsigned int,
		       const double* const);
  /*! \brief Write vector of doubles as attribute */
  void writeAttDoubles(hid_t, const std::string&, 
		       const std::vector<double>&);
  /*! \brief Read double attribute */
  double readAttDouble(hid_t objid, const std::string& name);
  /*! \brief Read vector of doubles attribute */
  std::vector<double> readAttDoubles(hid_t objid, const std::string& name);

  // Data writers and readers
  // 1D
  /*! \brief Write string as data */
  void writeDataString(hid_t, const std::string&, const std::string&);
  /*! \brief Write 1D array of strings as data */
  void writeDataStrings(hid_t, const std::string&, unsigned int,
			const std::string* const);
  /*! \brief Write vector of strings as data */
  void writeDataStrings(hid_t, const std::string&,
			const std::vector<std::string>&);
  /*! \brief Read a single string as data */
  std::string readDataString(hid_t, const std::string&);
  /*! \brief Read a vector of strings string as data */
  std::vector<std::string> readDataStrings(hid_t, const std::string&);

  /*! \brief Write a single boolean as data */
  void writeDataBool(hid_t, const std::string&, bool);
  /*! \brief Write 1D array of bools as data */
  void writeDataBools(hid_t, const std::string&, unsigned int,
		                  const bool* const);
  /*! \brief Write vector of bools as data */
  void writeDataBools(hid_t, const std::string&, const std::vector<bool>&);
  /*! \brief Read a single bool as data */
  bool readDataBool(hid_t, const std::string&);
  /*! \brief Read a vector of bools as data */
  std::vector<bool> readDataBools(hid_t, const std::string&);

  /*! \brief Write a single unsigned int as data */
  void writeDataUnsignedInt(hid_t, const std::string&, unsigned int);
  /*! \brief Write 1D array of unsigned ints as data */
  void writeDataUnsignedInts(hid_t, const std::string&, unsigned int,
			     const unsigned int* const);
  /*! \brief Write vector of unsigned ints as data */
  void writeDataUnsignedInts(hid_t, const std::string&,
			     const std::vector<unsigned int>&);
  /*! \brief Read a single unsigned int as data */
  unsigned int readDataUnsignedInt(hid_t, const std::string&);
  /*! \brief Read a vector of unsigned ints as data */
  std::vector<unsigned int> readDataUnsignedInts(hid_t, const std::string&);

  /*! \brief Write a single int as data */
  void writeDataInt(hid_t, const std::string&, int);
  /*! \brief Write 1D array of ints as data */
  void writeDataInts(hid_t, const std::string&, unsigned int,
		     const int* const);
  /*! \brief Write vector of ints as data */
  void writeDataInts(hid_t, const std::string&, const std::vector<int>&);
  /*! \brief Read a single int as data */
  int readDataInt(hid_t, const std::string&);
  /*! \brief Read a vector of ints as data */
  std::vector<int> readDataInts(hid_t, const std::string&);

  /*! \brief Write a single float as data */
  void writeDataFloat(hid_t, const std::string&, float);
  /*! \brief Write 1D array of floats as data */
  void writeDataFloats(hid_t, const std::string&, unsigned int,
			const float* const);
  /*! \brief Write vector of floats as data */
  void writeDataFloats(hid_t, const std::string&,
		       const std::vector<float>&);
  /*! \brief Read a single float as data */
  float readDataFloat(hid_t, const std::string&);
  /*! \brief Read 1D array of floats as data */
  void readDataFloats(hid_t, const std::string&, unsigned int, float* const);
  
  /*! \brief Write a single double as data */
  void writeDataDouble(hid_t, const std::string&, double);
  /*! \brief Write 1D array of doubles as data */
  void writeDataDoubles(hid_t, const std::string&, unsigned int,
			const double* const);
  /*! \brief Write vector of doubles as data */
  void writeDataDoubles(hid_t, const std::string&,
			const std::vector<double>&);
  /*! \brief Read a single double as data */
  double readDataDouble(hid_t, const std::string&);
  /*! \brief Read 1D array of doubles as data */
  void readDataDoubles(hid_t, const std::string&, unsigned int,
		       double* const);
  
  // 2D
  /*! \brief Write 2D array of floats as data */
  void writeData2DFloats(hid_t, const std::string&, unsigned int,
			 unsigned int, const float* const);
  /*! \brief Write 2D array of doubles as data */
  void writeData2DDoubles(hid_t, const std::string&, unsigned int,
			  unsigned int, const double* const);
}

#endif
