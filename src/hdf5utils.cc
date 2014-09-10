// Utilities for writing things to HDF5 handles
#include<algorithm>

#include "../include/hdf5utils.h"
#include "../include/affineExcept.h"


/*!
  \param[in] str File name
  \returns File type as determined by extension.

  This should only be expected to work for ASCII text, possibly
  not on all compilers.  Ah, the benefits of coding for yourself,
  and only really yourself.
*/
hdf5utils::outfiletype hdf5utils::getOutputFileType(const std::string& str) {
  size_t pos = str.find_last_of('.');
  if (pos == std::string::npos) return UNKNOWN; // Didn't find
  size_t strlen = str.size();
  if (pos == strlen - 1) return UNKNOWN; // . was last character
  std::string extn = str.substr(pos + 1);
  std::transform(extn.begin(), extn.end(), extn.begin(), ::toupper);
  // Switch doesn't work, so big if table
  if (extn == "TXT" || extn == "TEXT")
    return TXT;
  else if (extn == "FITS" || extn == "FIT")
    return FITS;
  else if (extn == "H5" || extn == "HDF5")
    return HDF5;
  else
    return UNKNOWN;
}

// There's a ton of repeated code here, but it seems preferable
// to cut and paste than to add another layer of function calls

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttString(hid_t objid, const std::string& name,
			       const std::string& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttString",
		       "Input handle is not valid when writing " + name);

  // String datatype
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);

  hsize_t adims;
  hid_t mems_id, att_id;

  const char * ctmp;

  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);

  att_id = H5Acreate1(objid, name.c_str(), datatype,
		      mems_id, H5P_DEFAULT);

  ctmp = value.c_str();
  H5Awrite(att_id, datatype, &ctmp);
  H5Aclose(att_id);
  H5Sclose(mems_id);
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttStrings(hid_t objid, const std::string& name,
				unsigned int n, 
				const std::string* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttStrings",
		       "Input handle is not valid");
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttStrings",
		       "Invalid number of elements in value when writing " +
		       name);

  // String datatype
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);

  hsize_t adims;
  hid_t mems_id, att_id;

  const char ** ctmp;
  ctmp = new const char*[n];
  for (unsigned int i = 0; i < n; ++i)
    ctmp[i] = value[i].c_str();

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), datatype,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, datatype, ctmp);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  delete[] ctmp;
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttStrings(hid_t objid, const std::string& name,
				const std::vector<std::string>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttStrings",
		       "Input handle is not valid when writing " + name);
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttStrings",
		       "Invalid number of elements in value when writing " + 
		       name);

  // String datatype
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);

  hsize_t adims;
  hid_t mems_id, att_id;

  const char ** ctmp;
  ctmp = new const char*[n];
  for (unsigned int i = 0; i < n; ++i)
    ctmp[i] = value[i].c_str();

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), datatype,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, datatype, ctmp);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  delete[] ctmp;
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttBool(hid_t objid, const std::string& name,
			     bool value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttBool",
		       "Input handle is not valid when writing " + name);

  hsize_t adims;
  hid_t mems_id, att_id;
  hbool_t bl;

  bl = static_cast<hbool_t>(value);
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &bl);
  H5Aclose(att_id);
  H5Sclose(mems_id);
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttBools(hid_t objid, const std::string& name,
			      unsigned int n, const bool* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttBools",
		       "Input handle is not valid when writing " + name);
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttBools",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, att_id;
  hbool_t *bl;

  bl = new hbool_t[n];
  for (unsigned int i = 0; i < n; ++i)
    bl[i] = static_cast<hbool_t>(value[i]);

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, bl);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  delete[] bl;
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttBools(hid_t objid, const std::string& name,
			      const std::vector<bool>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttBools",
		       "Input handle is not valid when writing " + name);
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttBools",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, att_id;
  hbool_t *bl;

  bl = new hbool_t[n];
  for (unsigned int i = 0; i < n; ++i)
    bl[i] = static_cast<hbool_t>(value[i]);

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, bl);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  delete[] bl;
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttUnsignedInts(hid_t objid, const std::string& name,
				     unsigned int n, 
				     const unsigned int* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttUnsignedInts",
		       "Input handle is not valid when writing " + name);
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttUnsignedInts",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, att_id;

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, value);
  H5Aclose(att_id);
  H5Sclose(mems_id);
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttUnsignedInts(hid_t objid, const std::string& name,
				     const std::vector<unsigned int>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttUnsignedInts",
		       "Input handle is not valid when writing " + name);
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttUnsignedInts",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, att_id;
  unsigned int *v;
  v = new unsigned int[n]; // Don't assume C++ vector.data() exists
  for (unsigned int i = 0; i < n; ++i)
    v[i] = value[i];

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, v);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  delete[] v;
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttFloats(hid_t objid, const std::string& name,
			       unsigned int n, const float* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttFloats",
		       "Input handle is not valid when writing " + name);
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttFloats",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, att_id;

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_FLOAT,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_FLOAT, value);
  H5Aclose(att_id);
  H5Sclose(mems_id);
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttFloats(hid_t objid, const std::string& name,
			       const std::vector<float>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttFloats",
		       "Input handle is not valid when writing " + name);
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttFloats",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, att_id;
  float *v;
  v = new float[n]; // Don't assume C++ vector.data() exists
  for (unsigned int i = 0; i < n; ++i)
    v[i] = value[i];

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_FLOAT,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_FLOAT, v);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  delete[] v;
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttDoubles(hid_t objid, const std::string& name,
				unsigned int n, const double* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttDoubles",
		       "Input handle is not valid when writing " + name);
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttDoubles",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, att_id;

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, value);
  H5Aclose(att_id);
  H5Sclose(mems_id);
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of attribute to write
  \param[in] value Value to write
*/
void hdf5utils::writeAttDoubles(hid_t objid, const std::string& name,
				const std::vector<double>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeAttDoubles",
		       "Input handle is not valid when writing " + name);
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeAttDoubles",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, att_id;
  double *v;
  v = new double[n]; // Don't assume C++ vector.data() exists
  for (unsigned int i = 0; i < n; ++i)
    v[i] = value[i];

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, name.c_str(), H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, v);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  delete[] v;
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of data to write
  \param[in] value Value to write
*/
void hdf5utils::writeDataUnsignedInts(hid_t objid, const std::string& name, 
				      unsigned int n,
				      const unsigned int* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeDataUnsignedInts",
		       "Input handle is not valid when writing " + name);
  if (n == 0)
    throw affineExcept("hdf5utils", "writeDataUnsignedInts",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, dat_id;

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  dat_id = H5Dcreate2(objid, name.c_str(), H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, value);
  H5Dclose(dat_id);
  H5Sclose(mems_id);
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of data to write
  \param[in] value Value to write
*/
void hdf5utils::writeDataFloats(hid_t objid, const std::string& name, 
				unsigned int n, const float* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeDataFloats",
		       "Input handle is not valid when writing " + name);
  if (n == 0)
    throw affineExcept("hdf5utils", "writeDataFloats",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, dat_id;

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  dat_id = H5Dcreate2(objid, name.c_str(), H5T_NATIVE_FLOAT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, value);
  H5Dclose(dat_id);
  H5Sclose(mems_id);
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of data to write
  \param[in] value Value to write
*/
void hdf5utils::writeDataDoubles(hid_t objid, const std::string& name, 
				 unsigned int n, const double* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeDataDoubles",
		       "Input handle is not valid when writing " + name);
  if (n == 0)
    throw affineExcept("hdf5utils", "writeDataDoubles",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims;
  hid_t mems_id, dat_id;

  adims = static_cast<hsize_t>(n);
  mems_id = H5Screate_simple(1, &adims, NULL);
  dat_id = H5Dcreate2(objid, name.c_str(), H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, value);
  H5Dclose(dat_id);
  H5Sclose(mems_id);
}

/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of data to write
  \param[in] value Value to write
*/
void hdf5utils::writeData2DFloats(hid_t objid, const std::string& name, 
				  unsigned int n1, unsigned int n2, 
				  const float* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeDataDoubles",
		       "Input handle is not valid when writing " + name);
  if (n1 * n2 == 0)
    throw affineExcept("hdf5utils", "writeDataDoubles",
		       "Invalid number of elements in value when writing " + 
		       name);

  hsize_t adims[2];
  hid_t mems_id, dat_id;

  adims[0] = static_cast<hsize_t>(n1);
  adims[1] = static_cast<hsize_t>(n2);
  mems_id = H5Screate_simple(2, adims, NULL);
  dat_id = H5Dcreate2(objid, name.c_str(), H5T_NATIVE_FLOAT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, value);
  H5Dclose(dat_id);
  H5Sclose(mems_id);
}
/*!
  \param[in] objid HDF5 handle to write to
  \param[in] name Name of data to write
  \param[in] value Value to write
*/
void hdf5utils::writeData2DDoubles(hid_t objid, const std::string& name, 
				   unsigned int n1, unsigned int n2, 
				   const double* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeDataDoubles",
		       "Input handle is not valid when writing " + name);
  if (n1 * n2 == 0)
    throw affineExcept("hdf5utils", "writeDataDoubles",
		       "Invalid number of elements in value when writing " +
		       name);

  hsize_t adims[2];
  hid_t mems_id, dat_id;

  adims[0] = static_cast<hsize_t>(n1);
  adims[1] = static_cast<hsize_t>(n2);
  mems_id = H5Screate_simple(2, adims, NULL);
  dat_id = H5Dcreate2(objid, name.c_str(), H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, value);
  H5Dclose(dat_id);
  H5Sclose(mems_id);
}
