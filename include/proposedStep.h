//proposedStep.h

#ifndef __proposedStep__
#define __proposedStep__

#include<mpi.h>

#include "paramSet.h"

/*!
  \brief Convenience class for keeping track of steps
 */
struct proposedStep {
  unsigned int update_idx; //!< Which one to update
  paramSet oldStep; //!< Previous step
  paramSet newStep; //!< New step to try
  double oldLogLike; //!< Log Likelihood of previous step
  double newLogLike; //!< Log Likelihood of new step
  float z; //!< Z value (stretch parameter)

  /*! \brief  Constructor with number of parameters */
  explicit proposedStep(unsigned int);
  proposedStep(const proposedStep&); //!< Copy constructor
  proposedStep(proposedStep&&); //!< Move constructor
  ~proposedStep(); //!< Destructor

  void clear(); //!<Clear contents
  void setNParams(unsigned int); //!< Set number of parameters
  unsigned int getNParams() const { return oldStep.getNParams(); }

  proposedStep& operator=(const proposedStep&); //!< Copy operator
  proposedStep& operator=(proposedStep&&); //!< Move assignment

  /*! MPI copy send operation */
  void sendSelf(MPI_Comm comm, int dest) const;
  /*! MPI copy receive operation */
  void receiveCopy(MPI_Comm comm, int src);
};

std::ostream& operator<<(std::ostream& os, const proposedStep& p);

#endif
