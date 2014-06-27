#include<limits>

#include "../include/global_settings.h"
#include "../include/proposedStep.h"

proposedStep::proposedStep(unsigned int npar) : oldStep(npar), newStep(npar) {
  update_idx = 0;
  oldLogLike = std::numeric_limits<double>::quiet_NaN();
  newLogLike = std::numeric_limits<double>::quiet_NaN();
  z = std::numeric_limits<double>::quiet_NaN();
}

proposedStep::proposedStep(const proposedStep& other) :
  update_idx(other.update_idx),
  oldStep(other.oldStep), newStep(other.newStep),
  oldLogLike(other.oldLogLike), newLogLike(other.newLogLike),
  z (other.z) { }

proposedStep::~proposedStep() {}

void proposedStep::clear() {
  update_idx = 0;
  oldStep.clear();
  newStep.clear();
  oldLogLike = std::numeric_limits<double>::quiet_NaN();
  newLogLike = std::numeric_limits<double>::quiet_NaN();
  z = std::numeric_limits<double>::quiet_NaN();
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
  z = other.z;
  return *this;
}

void proposedStep::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&update_idx), 1, MPI_UNSIGNED, 
	   dest, mcmc_affine::PSTSENDIDX, comm);
  oldStep.sendSelf(comm, dest);
  newStep.sendSelf(comm, dest);
  MPI_Send(const_cast<double*>(&oldLogLike), 1, MPI_DOUBLE, dest, 
	   mcmc_affine::PSTSENDOLIKE, comm);
  MPI_Send(const_cast<double*>(&newLogLike), 1, MPI_DOUBLE, dest, 
	   mcmc_affine::PSTSENDNLIKE, comm);
  MPI_Send(const_cast<float*>(&z), 1, MPI_FLOAT, dest, 
	   mcmc_affine::PSTSENDZ, comm);
}

void proposedStep::receiveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  MPI_Recv(&update_idx, 1, MPI_UNSIGNED, src, mcmc_affine::PSTSENDIDX, 
	   comm, &Info);
  oldStep.receiveCopy(comm, src);
  newStep.receiveCopy(comm, src);
  MPI_Recv(&oldLogLike, 1, MPI_DOUBLE, src, mcmc_affine::PSTSENDOLIKE, 
	   comm, &Info);
  MPI_Recv(&newLogLike, 1, MPI_DOUBLE, src, mcmc_affine::PSTSENDNLIKE, 
	   comm, &Info);
  MPI_Recv(&z, 1, MPI_FLOAT, src, mcmc_affine::PSTSENDZ, 
	   comm, &Info);
}

std::ostream& operator<<(std::ostream& os, const proposedStep& p) {
  os << "Update idx: " << p.update_idx << std::endl;
  os << "Old step: " << p.oldStep << " logLike: " << p.oldLogLike << std::endl;
  os << "New step: " << p.newStep << " logLike: " << p.newLogLike << std::endl;
  os << "z value: " << p.z;
  return os;
}
