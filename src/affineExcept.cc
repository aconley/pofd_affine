#include <iostream>
#include <sstream>

#include "../include/affineExcept.h"

using namespace std;

/*!
  \param[in] CLASS Class generating exception, or routine if not in a class
  \param[in] METHOD Method generating exception
  \param[in] ERRSTR Error string 

  Does the work of the constructors, except for setting the flags
*/
void affineExcept::init(const std::string& CLASS,
			const std::string& METHOD,
			const std::string& ERRSTR) {
  errclass = CLASS;
  errmethod = METHOD;
  errstr = ERRSTR;
}

/*!
  \param[in] errstr Error string
*/
affineExcept::affineExcept(const std::string& errstr) {
  init("", "", errstr);
}

/*!
  \param[in] errclass Class generating error
  \param[in] errmethod Method generating error
  \param[in] errstr Error string
*/
affineExcept::affineExcept(const std::string& errclass,
			   const std::string& errmethod,
			   const std::string& errstr) {
  init(errclass, errmethod, errstr);
}

const char* affineExcept::what() const throw() {
  std::stringstream str;
  bool first; //Last element won't have a linebreak at the end
  first = true;
  if (!errclass.empty()) {
    str << "Error Class/Namespace: " << errclass;
    first = false;
  }
  if (!errmethod.empty()) {
    if (!first) str << std::endl;
    str << "Method: " << errmethod;
    first = false;
  }
  if (!errstr.empty()) {
    if (!first) str << std::endl;
    str << "Error Message: " << errstr;
  }
  whatmsg = str.str();
  return whatmsg.c_str();
}

/*
  Provides output capabilities, not outputting stuff not set.
*/
ostream& operator<<(std::ostream& os, const affineExcept& err) {
  os << err.what();
  return os;
}
