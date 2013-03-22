#include<limits>

#include "../include/global_settings.h"
#include "../include/proposedStep.h"

proposedStep::proposedStep(unsigned int npar) : oldStep(npar), newStep(npar) {
  update_idx = 0;
  oldLogLike = std::numeric_limits<double>::quiet_NaN();
  newLogLike = std::numeric_limits<double>::quiet_NaN();
}

proposedStep::proposedStep(const proposedStep& other) :
  update_idx(other.update_idx),
  oldStep(other.oldStep), newStep(other.newStep),
  oldLogLike(other.oldLogLike), newLogLike(other.newLogLike) { }

proposedStep::~proposedStep() {}

void proposedStep::clear() {
  update_idx = 0;
  oldStep.clear();
  newStep.clear();
  oldLogLike = std::numeric_limits<double>::quiet_NaN();
  newLogLike = std::numeric_limits<double>::quiet_NaN();
}

/*!
  \param[in] n New number of parameters
  
  Generally will not preserve content, but doesn't clear it either
 */
void proposedStep::setNParams(unsigned int n) {
  oldStep.setNParams(n);
  newStep.setNParams(n);
}

proposedStep& proposedStep::operator=(const proposedStep& other) {
  if (this == &other) return *this; //Self copy
  update_idx = other.update_idx;
  oldStep = other.oldStep;
  newStep = other.newStep;
  oldLogLike = other.oldLogLike;
  newLogLike = other.newLogLike;
  return *this;
}

void proposedStep::sendSelf(MPI::Comm& comm, int dest) const {
  comm.Send(&update_idx,1,MPI::UNSIGNED,dest,mcmc_affine::PSTSENDIDX);
  oldStep.sendSelf(comm,dest);
  newStep.sendSelf(comm,dest);
  comm.Send(&oldLogLike,1,MPI::DOUBLE,dest,mcmc_affine::PSTSENDOLIKE);
  comm.Send(&newLogLike,1,MPI::DOUBLE,dest,mcmc_affine::PSTSENDNLIKE);
}

void proposedStep::recieveCopy(MPI::Comm& comm, int src) {
  comm.Recv(&update_idx,1,MPI::UNSIGNED,src,mcmc_affine::PSTSENDIDX);
  oldStep.recieveCopy(comm,src);
  newStep.recieveCopy(comm,src);
  comm.Recv(&oldLogLike,1,MPI::DOUBLE,src,mcmc_affine::PSTSENDOLIKE);
  comm.Recv(&newLogLike,1,MPI::DOUBLE,src,mcmc_affine::PSTSENDNLIKE);
}

std::ostream& operator<<(std::ostream& os, const proposedStep& p) {
  os << "Update idx: " << p.update_idx << std::endl;
  os << "Old step: " << p.oldStep << " logLike: " << p.oldLogLike << std::endl;
  os << "New step: " << p.newStep << " logLike: " << p.newLogLike;
  return os;
}
