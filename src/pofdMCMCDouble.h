//2D P(D) case MCMC driver
#include<string>

#include<affineEnsemble.h>
#include<calcLikeDouble.h>

/*!
  \brief Class for doing 2D P(D1,D2) MCMC fits
*/
class pofdMCMCDouble : public affineEnsemble {
private:
  std::string initfile; //!< File to read initilaization from
  std::string specfile; //!< File to read spec from

  initFileDoubleLogNormal ifile; //!< Stores initial values and limits on parameters

  calcLikeDouble likeSet; //!< Does likelihood calculation
  
  bool initChainsMaster(); //!< Initialization routine for master node
  bool initChainsSlave(); //!< Initialization routine for slave node

public:
  /*! \brief Constructor */
  pofdMCMCDouble(const std::string& initfile, const std::string& specfile,
		 unsigned int NWALKERS, unsigned int NSAMPLES, 
		 unsigned int INIT_STEPS=50, unsigned int MIN_BURN=50, 
		 double BURN_MULTIPLE=5.0, double SCALEFAC=2.0); 
  ~pofdMCMCDouble() {}; //!< Destructor

  void initChains(); //!< Initializes data between MPI jobs
  double getLogLike(const paramSet&); //!< Evaluates log likelihood
};


