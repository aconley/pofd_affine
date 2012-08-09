#ifndef __affineExcept__
#define __affineExcept__

#include <string>
#include <iostream>

/*! 
  \brief Exception class for mcmc_affine
*/

class affineExcept {

  void init(const std::string& errclass, const std::string& errmethod,
	    const std::string& errstr, int err);  //!< Internal initialization

 public:
  bool classset; //!< Is errclass set
  bool methodset; //!< Is errmethod set
  bool strset; //!< Is errstr set
  bool errset; //!< Is errnum set
  std::string errclass;          //!< Class throwing the exception
  std::string errmethod;         //!< Method throwing the exception
  std::string errstr;            //!< Error string (user consumption)
  int errnum;               //!< Error code (class or method specific)

  // Constructors

  affineExcept(); //!< Basic constructor
  explicit affineExcept(const std::string errstr); //!< Just with errstring
  explicit affineExcept(const std::string errstr,int err); //!< Errstring and number
  explicit affineExcept(const std::string errclass,const std::string errmethod,
			const std::string errstr); //!< Class, method, error string
  explicit affineExcept(const std::string errclass,const std::string errmethod,
			const std::string errstr,int err); //!< Class, method, error string, and number

  std::string what() const; //!< Explain error

};

std::ostream& operator<<(std::ostream& os,const affineExcept& ex);//!< Output operator for affineExcept

#endif
