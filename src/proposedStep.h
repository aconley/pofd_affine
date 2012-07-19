//proposedStep.h

#ifndef __proposedStep__
#define __proposedStep__

#include<mpi.h>

#include<paramSet.h>

/*!
  \brief Convenience class for keeping track of steps
 */
struct proposedStep {
  unsigned int update_idx; //!< Which one to update
  paramSet oldStep; //!< Previous step
  paramSet newStep; //!< New step to try
  double   oldLogLike; //!< Log Likelihood of previous step
  double   newLogLike; //!< Log Likelihood of new step
 
  proposedStep(unsigned int); //!< Constructor with number of parameters
  proposedStep(const proposedStep&); //!< Copy constructor
  ~proposedStep(); //!< Destructor

  void clear(); //!<Clear contents

  proposedStep& operator=(const proposedStep&); //!< Copy operator

  /*! MPI copy send operation */
  void sendSelf(MPI::Comm& comm, int dest) const;
  /*! MPI copy recieve operation */
  void recieveCopy(MPI::Comm& comm, int src);
};

std::ostream& operator<<(std::ostream& os, const proposedStep& p);

#endif
