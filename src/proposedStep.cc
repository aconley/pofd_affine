#include<limits>

#include "../include/global_settings.h"
#include "../include/proposedStep.h"

/*!
  \param[in] npar Number of parameters
*/
proposedStep::proposedStep(unsigned int npar) : oldStep(npar), newStep(npar) {
  update_idx = 0;
  oldLogLike = std::numeric_limits<double>::quiet_NaN();
  newLogLike = std::numeric_limits<double>::quiet_NaN();
  z = std::numeric_limits<double>::quiet_NaN();
}

/*!
  \param[in] other Other proposed step to initialize from
*/
proposedStep::proposedStep(const proposedStep& other) :
  update_idx(other.update_idx),
  oldStep(other.oldStep), newStep(other.newStep),
  oldLogLike(other.oldLogLike), newLogLike(other.newLogLike),
  z (other.z) { }

/*!
  \param[in] other Other proposed step to apply move semantics from
*/
proposedStep::proposedStep(proposedStep&& other) {
  update_idx = other.update_idx;
  oldLogLike = other.oldLogLike;
  newLogLike = other.newLogLike;
  z = other.z;
  oldStep = std::move(other.oldStep);
  newStep = std::move(other.newStep);
  other.update_idx = 0;
  other.oldLogLike = std::numeric_limits<double>::quiet_NaN();
  other.newLogLike = std::numeric_limits<double>::quiet_NaN();
  other.z = std::numeric_limits<double>::quiet_NaN();
}

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

/*!
  \param[in] other proposedStep to copy
*/
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

/*!
  \param[in] other Other proposed step to apply move semantics from
*/
proposedStep& proposedStep::operator=(proposedStep&& other) {
  if (this == &other) return *this; //Self copy
  update_idx = other.update_idx;
  oldLogLike = other.oldLogLike;
  newLogLike = other.newLogLike;
  z = other.z;
  oldStep = std::move(other.oldStep);
  newStep = std::move(other.newStep);
  other.update_idx = 0;
  other.oldLogLike = std::numeric_limits<double>::quiet_NaN();
  other.newLogLike = std::numeric_limits<double>::quiet_NaN();
  other.z = std::numeric_limits<double>::quiet_NaN();
  return *this;
}

/*!
  \param[in] comm Communicator
  \param[in] dest Destination of messages
*/
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

/*!
  \param[in] comm Communicator
  \param[in] src Source of messages
*/
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

/*!
  \param[in] os Output stream to write to
  \param[in] p Step to write to output stream

  \returns A reference to the modifed output stream
*/
std::ostream& operator<<(std::ostream& os, const proposedStep& p) {
  os << "Update idx: " << p.update_idx << std::endl;
  os << "Old step: " << p.oldStep << " logLike: " << p.oldLogLike << std::endl;
  os << "New step: " << p.newStep << " logLike: " << p.newLogLike << std::endl;
  os << "z value: " << p.z;
  return os;
}
