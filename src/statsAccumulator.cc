#include<cmath>
#include<algorithm>

#include "../include/statsAccumulator.h"

/*! \brief 
    \param[in] minval Minimum flux value
    \param[in] maxval Maximum flux value
    \param[in] logspace Space points logarithmically instead of linearly
    \param[in] n Number of points
    \param[out] s On output, loaded with flux values.
 */
void statsAccumulator::setupFluxes(float minval, float maxval, 
                                   bool logspace, unsigned int n,
                                   float* s) {
  if (logspace) {
    float log_minval = log(minval);
    float delta = (log(maxval) - log_minval) /
      static_cast<double>(n - 1);
    for (unsigned int i = 0; i < n; ++i)
      s[i] = exp(log_minval + delta * static_cast<float>(i));
  } else {
    float delta = (maxval - minval) / static_cast<float>(n - 1);
    for (unsigned int i = 0; i < n; ++i)
      s[i] = minval + delta * static_cast<float>(i);
  }
}

/*! 
  \param[in] n1 len of mean, median, 2nd dimension of snake
  \param[in] n2 Second dimension of working
  \param[inout] working Working array, n1 by n2.  Modified on output
  \param[out] mean Mean at each of n1 points here
  \param[out] median Median at each of n1 points
  \param[out] snake Error snake, nprob by n1 by 2, flattened
*/
void statsAccumulator::accumulateStats(unsigned int n1, unsigned int n2,
                                       float *working, float *mean, 
                                       float* median, float *snake) {

  float lowidx[nprob], highidx[nprob]; // Fractional index of interpolation
  for (unsigned int i = 0; i < nprob; ++i)
    lowidx[i] = 0.5 * (1.0 - plevels[i]) * n2;
  for (unsigned int i = 0; i < nprob; ++i)
    highidx[i] = 0.5 * (1 + plevels[i]) * n2;
  
  unsigned int sidx;
  float *wptr;
  for (unsigned int j = 0; j < n1; ++j) {
    wptr = working + j * n2;
    std::sort(wptr, wptr + n2);
    double mnval = wptr[0]; // Accumulate in double
    for (unsigned int k = 1; k < n2; ++k) mnval += wptr[k];
    mean[j] = static_cast<float>(mnval / static_cast<double>(n2));
    median[j] = wptr[n2 / 2];
    for (unsigned int k = 0; k < nprob; ++k) {
      sidx = 2 * (k * n1 + j);
      // Lower value
      unsigned int idxl = static_cast<unsigned int>(lowidx[k]);
      snake[sidx] = wptr[idxl] + 
        (wptr[idxl+1] - wptr[idxl]) * (lowidx[k] - idxl);
      // And upper
      unsigned int idxu = static_cast<unsigned int>(highidx[k]);
      snake[sidx + 1] = wptr[idxu] + 
        (wptr[idxu+1] - wptr[idxu]) * (highidx[k] - idxu);
    }
  }
}
