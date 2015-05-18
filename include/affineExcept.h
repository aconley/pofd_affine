#ifndef __affineExcept__
#define __affineExcept__

#include <string>
#include <ostream>
#include <exception>
/*! 
  \brief Exception class for mcmc_affine

  Adds class/method fields to except
*/

class affineExcept : public std::exception {
 private:
  void init(const std::string&, const std::string&,
            const std::string&);  //!< Initializer

  std::string errclass;          //!< Class throwing the exception
  std::string errmethod;         //!< Method throwing the exception
  std::string errstr;            //!< Error string (user consumption)

  mutable std::string whatmsg;  //!< Holds combined error message
 public:
  // Constructors
  explicit affineExcept(const std::string& errstr); //!< Just with errstring
  /*! \brief Class, method, error string*/
  explicit affineExcept(const std::string& errclass, 
                        const std::string& errmethod,
                        const std::string& errstr); 
  ~affineExcept() throw() {}

  std::string getErrClass() const { return errclass; }
  std::string getErrMethod() const { return errmethod; }
  std::string getErrStr() const { return errstr; }

  const char* what() const throw();
};

std::ostream& operator<<(std::ostream& os, const affineExcept& ex);//!< Output operator for affineExcept

#endif
