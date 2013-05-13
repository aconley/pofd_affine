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

  calcLike likeSet; //!< Does likelihood calculation

  bool initChainsMaster(); //!< Initialization routine for master node
  bool initChainsSlave(); //!< Initialization routine for slave node

public:
  pofdMCMC(const std::string& initfile, const std::string& specfile,
	   unsigned int NWALKERS, unsigned int NSAMPLES, 
	   unsigned int INIT_STEPS=50, unsigned int MIN_BURN=50, 
	   bool FIXED_BURN=false, float BURN_MULTIPLE=5.0, 
	   float SCALEFAC=2.0); //!< Constructor
  ~pofdMCMC() {}; //!< Destructor

  void initChains(); //!< Initializes data between MPI jobs
  void generateInitialPosition(const paramSet&); //!< Sets up initial position
  double getLogLike(const paramSet&); //!< Evaluates log likelihood
};
