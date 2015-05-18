//2D P(D) case MCMC driver
#include<string>

#include "../include/affineEnsemble.h"
#include "../include/specFileDouble.h"
#include "../include/calcLikeDouble.h"
#include "../include/initFileDoubleLogNormal.h"

/*!
  \brief Class for doing 2D P(D1,D2) MCMC fits
*/
class pofdMCMCDouble : public affineEnsemble {
 private:
  std::string initfile; //!< File to read initilaization from
  std::string specfile; //!< File to read spec from

  initFileDoubleLogNormal ifile; //!< Stores initial values and limits on parameters
  specFileDouble spec_info; //!< Controls how fit is done

  unsigned int ntot; //!< Number of total knots in model; useful to keep around

  calcLikeDouble likeSet; //!< Does likelihood calculation
  
  bool initChainsMaster(); //!< Initialization routine for master node
  bool initChainsSlave(); //!< Initialization routine for slave node

 public:
  /*! \brief Constructor */
  pofdMCMCDouble(const std::string& initfile, const std::string& specfile,
                 unsigned int NWALKERS, unsigned int NSAMPLES, 
                 unsigned int INIT_STEPS=50, double INIT_TEMP=2.0,
                 unsigned int MIN_BURN=50, bool FIXED_BURN=false, 
                 float BURN_MULTIPLE=5.0, float SCALEFAC=2.0); //!< Constructor
  pofdMCMCDouble(const pofdMCMCDouble&)=delete;
  pofdMCMCDouble(pofdMCMCDouble&&)=delete;
  virtual ~pofdMCMCDouble() {}; //!< Destructor

  pofdMCMCDouble& operator=(const pofdMCMCDouble&)=delete;
  pofdMCMCDouble& operator=(pofdMCMCDouble&&)=delete;
  
  void setVerbosity(unsigned int) override;

  void initChains() override; //!< Initializes data between MPI jobs
  bool areParamsValid(const paramSet&) const override; //!< Test against limits
  void generateInitialPosition(const paramSet&) override; //!< Sets up initial position
  double getLogLike(const paramSet&, bool&) override; //!< Evaluates log likelihood
  void fillBonusParams(paramSet& par, bool rej) override;

  void writeToHDF5Handle(hid_t) const override; //!< Serialize to HDF5 handle
};
