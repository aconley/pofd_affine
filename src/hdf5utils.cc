// Utilities for writing things to HDF5 handles

#include "../include/hdf5utils.h"
#include "../include/affineExcept.h"

void hdf5utils::writeAttString(hid_t objid, const std::string& name,
			       const std::string& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeString",
		       "Input handle is not valid");

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
  H5Awrite(att_id, datatype, ctmp);
  H5Aclose(att_id);
  H5Sclose(mems_id);
}

void hdf5utils::writeAttStrings(hid_t objid, const std::string& name,
				unsigned int n, 
				const std::string* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeStrings",
		       "Input handle is not valid");
  if (n == 0)
    throw affineExcept("hdf5utils", "writeStrings",
		       "Invalid number of elements in value");

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

void hdf5utils::writeAttStrings(hid_t objid, const std::string& name,
				const std::vector<std::string>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeStrings",
		       "Input handle is not valid");
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeStrings",
		       "Invalid number of elements in value");

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

void hdf5utils::writeAttBool(hid_t objid, const std::string& name,
			     bool value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeBool",
		       "Input handle is not valid");

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

void hdf5utils::writeAttBools(hid_t objid, const std::string& name,
			      unsigned int n, const bool* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeBools",
		       "Input handle is not valid");
  if (n == 0)
    throw affineExcept("hdf5utils", "writeBools",
		       "Invalid number of elements in value");

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

void hdf5utils::writeAttBools(hid_t objid, const std::string& name,
			      const std::vector<bool>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeBools",
		       "Input handle is not valid");
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeBools",
		       "Invalid number of elements in value");

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

void hdf5utils::writeAttUnsignedInts(hid_t objid, const std::string& name,
				     unsigned int n, 
				     const unsigned int* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeUnsignedInts",
		       "Input handle is not valid");
  if (n == 0)
    throw affineExcept("hdf5utils", "writeUnsignedInts",
		       "Invalid number of elements in value");

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

void hdf5utils::writeAttUnsignedInts(hid_t objid, const std::string& name,
				     const std::vector<unsigned int>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeUnsignedInts",
		       "Input handle is not valid");
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeUnsignedInts",
		       "Invalid number of elements in value");

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


void hdf5utils::writeAttFloats(hid_t objid, const std::string& name,
			       unsigned int n, const float* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeFloats",
		       "Input handle is not valid");
  if (n == 0)
    throw affineExcept("hdf5utils", "writeFloats",
		       "Invalid number of elements in value");

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

void hdf5utils::writeAttFloats(hid_t objid, const std::string& name,
			       const std::vector<float>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeFloats",
		       "Input handle is not valid");
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeFloats",
		       "Invalid number of elements in value");

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

void hdf5utils::writeAttDoubles(hid_t objid, const std::string& name,
				unsigned int n, const double* const value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeDoubles",
		       "Input handle is not valid");
  if (n == 0)
    throw affineExcept("hdf5utils", "writeDoubles",
		       "Invalid number of elements in value");

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

void hdf5utils::writeAttDoubles(hid_t objid, const std::string& name,
				const std::vector<double>& value) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("hdf5utils", "writeDoubles",
		       "Input handle is not valid");
  unsigned int n = value.size();
  if (n == 0)
    throw affineExcept("hdf5utils", "writeDoubles",
		       "Invalid number of elements in value");

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
