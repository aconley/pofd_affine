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
  paramSet oldStep;
  paramSet newStep;
  double   oldLogLike;
  double   newLogLike;

  proposedStep(unsigned int);
  proposedStep(const proposedStep&);
  ~proposedStep();

  void clear();

  proposedStep& operator=(const proposedStep&);

  void sendSelf(MPI::Comm& comm, int dest) const;
  void recieveCopy(MPI::Comm& comm, int src);
};

std::ostream& operator<<(std::ostream& os, const proposedStep& p);

#endif
