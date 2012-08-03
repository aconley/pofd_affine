//1D P(D) case MCMC driver
#include<string>

#include<affineEnsemble.h>
#include<calcLike.h>

class pofdMCMC : public affineEnsemble {
private:
  std::string initfile; //!< File to read initilaization from
  std::string specfile; //!< File to read spec from

  calcLike like; //!< Does likelihood calculation
  initFileKnots ifile; //!< Stores initial values and limits on parameters
  
  void initChainsMaster(); //!< Initialization routine for master node
  void initChainsSlave(); //!< Initialization routine for slave node

public:
  pofdMCMC(const std::string& initfile, const std::string& specfile,
	   unsigned int NWALKERS, unsigned int NSAMPLES, 
	   unsigned int INIT_STEPS=50, unsigned int MIN_BURN=50, 
	   double BURN_MULTIPLE=5.0, double SCALEFAC=2.0); //!< Constructor
  ~pofdMCMC() {}; //!< Destructor

  void initChains(); //!< Initializes data between MPI jobs
  double getLogLike(const paramSet&); //!< Evaluates log likelihood
};


