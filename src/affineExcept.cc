#include <iostream>
#include <sstream> 

#include "../include/affineExcept.h"

using namespace std;

/*!
  \param[in] CLASS Class generating exception, or routine if not in a class
  \param[in] METHOD Method generating exception
  \param[in] ERRSTR Error string 
  \param[in] ERRNUMBER Error number

  Does the work of the constructors, except for setting the flags
*/
void affineExcept::init(const std::string& CLASS,
			const std::string& METHOD,
			const std::string& ERRSTR, int ERRNUMBER) {
  errclass = CLASS;
  errmethod = METHOD;
  errstr = ERRSTR;
  errnum = ERRNUMBER;
}

/*!
  Basic constructor
*/
affineExcept::affineExcept() {
  init("","","",0);
  classset = methodset = strset = errset = false;
}

/*!
  \param[in] errstr Error string
*/
affineExcept::affineExcept(const std::string& errstr) {
  init("", "", errstr, 0);
  strset = true;
}

/*!
  \param[in] errstr Error string
  \paran[in] errnum Error number
*/
affineExcept::affineExcept(const std::string& errstr, int errnum) {
  init("", "", errstr, errnum);
  strset = errset = true;
}

/*!
  \param[in] errclass Class generating error
  \param[in] errmethod Method generating error
  \param[in] errstr Error string
*/
affineExcept::affineExcept(const std::string& errclass,
			   const std::string& errmethod,
			   const std::string& errstr) {
  init(errclass, errmethod, errstr, 0);
  classset = methodset = strset = true;
}

/*!  
  \param[in] errclass Class generating error
  \param[in] errmethod Method generating error
  \param[in] errstr Error string
  \paran[in] errnum Error number

*/
affineExcept::affineExcept(const std::string& errclass,
			   const std::string& errmethod,
			   const std::string& errstr, 
			   int errnum) {
  init(errclass, errmethod, errstr, errnum);
  classset = methodset = strset = errset = true;
}

std::string affineExcept::what() const {
  std::stringstream str;
  bool first; //Last element won't have a linebreak at the end
  first = true;
  if (classset) {
    str << "Error Class/Namespace: " << errclass;
    first = false;
  }
  if (methodset) {
    if (!first) str << std::endl;
    str << "Method: " << errmethod;
    first = false;
  }
  if (strset) {
    if (!first) str << std::endl;
    str << "Error Message: " << errstr;
    first = false;
  }
  if (errset) {
    if (!first) str << std::endl;
    str << "Error Code: " << errnum;
  }
  return str.str();
}

/*
  Provides output capabilities, not outputting stuff not set.
*/
ostream& operator<<(std::ostream& os, const affineExcept& err) {
  os << err.what();
  return os;
}
