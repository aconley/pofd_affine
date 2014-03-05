//1D P(D) case MCMC driver
#include<string>

#include "../include/affineEnsemble.h"
#include "../include/calcLike.h"
#include "../include/specFile.h"

/*!
  \brief Class for doing 1D P(D) fits
*/
class pofdMCMC : public affineEnsemble {
private:
  std::string initfile; //!< File to read initilaization from
  std::string specfile; //!< File to read spec from

  initFileKnots ifile; //!< Stores initial values and limits on parameters
  specFile spec_info; //!< Controls how fit is done

  unsigned int nknots; //!/ Number of knots in model; useful to keep around

  calcLike likeSet; //!< Does likelihood calculation

  bool initChainsMaster(); //!< Initialization routine for master node
  bool initChainsSlave(); //!< Initialization routine for slave node

public:
  pofdMCMC(const std::string& initfile, const std::string& specfile,
	   unsigned int NWALKERS, unsigned int NSAMPLES, 
	   unsigned int INIT_STEPS=50, double INIT_TEMP=2.0,
	   unsigned int MIN_BURN=50, bool FIXED_BURN=false, 
	   float BURN_MULTIPLE=5.0, float SCALEFAC=2.0); //!< Constructor
  ~pofdMCMC() {}; //!< Destructor

  void setVerbosity(unsigned int);

  void initChains(); //!< Initializes data between MPI jobs
  bool areParamsValid(const paramSet&) const; //!< Test against limits
  void generateInitialPosition(const paramSet&); //!< Sets up initial position
  double getLogLike(const paramSet&, bool&); //!< Evaluates log likelihood
  void fillBonusParams(paramSet&, bool rej=false); //!< Add mean flux per area
  void writeToHDF5Handle(hid_t) const; //!< Serialize to HDF5 handle
};
