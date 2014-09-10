//PDDouble.cc
#include<limits>
#include<sstream>

#include<fitsio.h>
#include<fftw3.h>
#include<hdf5.h>

#include "../include/global_settings.h"
#include "../include/PDDouble.h"
#include "../include/affineExcept.h"
#include "../include/hdf5utils.h"

const double PDDouble::lowsigval = 2.5;

/*!
  \param[in] N1 Dimension along band 1
  \param[in] MINFLUX1 Minimum flux density in band 1
  \param[in] DFLUX1 Delta flux density along band 1
  \param[in] N2 Dimension along band 2
  \param[in] MINFLUX2 Minimum flux density in band 2
  \param[in] DFLUX2 Delta flux density along band 2
  \param[in] LOG Assume data is stored as log
*/
PDDouble::PDDouble(unsigned int N1, double MINFLUX1, double DFLUX1,
		   unsigned int N2, double MINFLUX2, double DFLUX2,
		   bool LOG) : n1(N1), n2(N2), capacity(N1*N2),
			       logflat(LOG),
			       minflux1(MINFLUX1), dflux1(DFLUX1),
			       minflux2(MINFLUX2), dflux2(DFLUX2) {
  if (capacity == 0) pd_ = NULL; else
    pd_ = (double *) fftw_malloc(sizeof(double) * capacity);
}

PDDouble::~PDDouble() {
  if (pd_ != NULL) fftw_free(pd_);
}

/*!
  \param[in] N1 New number of elements, band 1
  \param[in] N2 New number of elements, band 2

  Only preserves data if new dimensions are the same as the old.
*/
void PDDouble::resize(unsigned int N1, unsigned int N2) {
  //Doesn't actually resize arrays if it can avoid it
  unsigned int newcap = N1 * N2;
  if (newcap > capacity) {
    if (pd_ != NULL) fftw_free(pd_);
    if (newcap > 0) pd_ = (double *) fftw_malloc(sizeof(double) * newcap);
    else pd_ = NULL;
    capacity = newcap;
  }
  n1 = N1;
  n2 = N2;
}

/*!
  Tries to preserve data
*/
void PDDouble::shrink() {
  unsigned int newcap = n1*n2;
  if (newcap < capacity) {
    if (newcap > 0) {
      double* tmp = (double*) fftw_malloc(sizeof(double) * newcap);
      for (unsigned int i = 0; i < newcap; ++i)
	tmp[i] = pd_[i];
      if (pd_ != NULL) fftw_free(pd_);
      pd_ = tmp;
    } else {
      if (pd_ != NULL) fftw_free(pd_);
      pd_ = NULL;
    }
    capacity = newcap;
  }
}

/*!
  \param[in] N1 New number of elements
  \param[in] N2 New number of elements

  Doesn't preserve data unless new size is the same as old
*/
void PDDouble::strict_resize(unsigned int N1, unsigned int N2) {
  unsigned int newcap = N1*N2;
  if (newcap != capacity) {
    if (pd_ != NULL) fftw_free(pd_);
    if (newcap > 0) pd_ = (double*) fftw_malloc(sizeof(double) * newcap);
    else pd_ = NULL;
    capacity = newcap;
  }
  n1 = N1;
  n2 = N2;
}

/*!
  \returns Total of all elements
*/
double PDDouble::getTotal() const {
  if ((n1 == 0) || (n2 == 0))
    return std::numeric_limits<double>::quiet_NaN();
  double retval;
  unsigned int sz = n1*n2;
  if (logflat) {
    retval = exp2(pd_[0]);
    for (unsigned int i = 1; i < sz; ++i)
      retval += exp2(pd_[i]);
  } else {
    retval = pd_[0];
    for (unsigned int i = 1; i < sz; ++i)
      retval += pd_[i];
  }
  return retval;
}

/*!
  \returns Integral of P(D)

  Uses the trapezoidal rule
*/
double PDDouble::getIntegral() const {
  if ((n1 == 0) || (n2 == 0))
    return std::numeric_limits<double>::quiet_NaN();
  
  double tot, *rowptr;
  if (logflat) {
    tot = 0.5*exp2(pd_[0]);
    for (unsigned int j = 1; j < n2-1; ++j)
      tot += exp2(pd_[j]);
    tot += 0.5*exp2(pd_[n2-1]);
    tot *= 0.5;
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + i*n2;
      tot += 0.5*exp2(rowptr[0]);
      for (unsigned int j = 1; j < n2-1; ++j)
	tot += exp2(rowptr[j]);
      tot += 0.5*exp2(rowptr[n2-1]);
    }
    rowptr = pd_ + (n1-1)*n2;
    tot += 0.25*exp2(rowptr[0]);
    for (unsigned int j = 1; j < n2-1; ++j)
      tot += 0.5*exp2(rowptr[j]);
    tot += 0.25*exp2(rowptr[n2-1]);
  } else {
    tot = 0.5*pd_[0];
    for (unsigned int j = 1; j < n2-1; ++j)
      tot += pd_[j];
    tot += 0.5*pd_[n2-1];
    tot *= 0.5;
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + i*n2;
      tot += 0.5*rowptr[0];
      for (unsigned int j = 1; j < n2-1; ++j)
	tot += rowptr[j];
      tot += 0.5*rowptr[n2-1];
    }
    rowptr = pd_ + (n1-1)*n2;
    tot += 0.25*rowptr[0];
    for (unsigned int j = 1; j < n2-1; ++j)
      tot += 0.5*rowptr[j];
    tot += 0.25*rowptr[n2-1];
  }
  return tot*dflux1*dflux2;
}

/*!
  Normalize the P(D), using the trapezoidal rule
  to integrate
*/
void PDDouble::normalize() {
  if ((n1 == 0) || (n2 == 0))
    throw affineExcept("PDDouble", "normalize",
		       "No information present to normalize");
  //Note, because of the 0.5 edge pieces we don't just use
  // getTotal
  double tot = getIntegral();
  unsigned int sz = n1*n2;
  if (logflat) {
    double lgtot = log2(tot);
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] -= lgtot;
  } else {
    double itot = 1.0/tot;
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] *= itot;
  }
}

/*!
  \param[in] nocheck Don't check whether or not values are positive

  If doing check and value is less than or equal to zero, set to 
  pofd_mcmc::smalllogval
*/
void PDDouble::applyLog(bool nocheck) {
  if (logflat) return;
  unsigned int sz = n1*n2;
  double val;
  if (nocheck)
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] = log2(pd_[i]);
  else 
    for (unsigned int i = 0; i < sz; ++i) {
      val = pd_[i];
      if (val <= 0.0) pd_[i] = pofd_mcmc::smalllogval; 
      else pd_[i] = log2(val);
    }
  logflat = true;
}

void PDDouble::deLog() {
  if (!logflat) return;
  unsigned int sz = n1*n2;
  for (unsigned int i = 0; i < sz; ++i)
    pd_[i] = exp2(pd_[i]);
  logflat = false;
}


/*
  \param[in] donorm Do not assume P(D) is normalized
*/
void PDDouble::edgeFix(bool donorm) {
  //Compute mean and stdev in each row and column.
  if (n1 < 3 || n2 < 3) return; //No point

  if (logflat)
    throw affineExcept("PDDouble", "edgeFix", "Not supported for logged PDs");

  //Get mean and vars
  std::pair<dblpair, dblpair> mnvar;
  double mn1, mn2, var1, var2;
  mnvar = getMeansAndVars(donorm);
  mn1 = mnvar.first.first; var1 = mnvar.first.second;
  mn2 = mnvar.second.first; var2 = mnvar.second.second;
  if (std::isnan(mn1) || std::isinf(mn1) ||
      std::isnan(var1) || std::isinf(var1) ||
      std::isnan(mn2) || std::isinf(mn2) ||
      std::isnan(var2) || std::isinf(var2)) {
    std::stringstream errstr;
    errstr << "Problem with means/vars: " << std::endl;
    if (std::isnan(mn1)) errstr << std::endl << "Mean 1 is NaN";
    if (std::isinf(mn1)) errstr << std::endl<< "Mean 1 is Inf";
    if (std::isnan(mn2)) errstr << std::endl << "Mean 2 is NaN";
    if (std::isinf(mn2)) errstr << std::endl << "Mean 2 is Inf";
    if (std::isnan(var1)) errstr << std::endl << "Var 1 is NaN";
    if (std::isinf(var1)) errstr << std::endl << "Var 1 is Inf";
    if (std::isnan(var2)) errstr << std::endl << "Var 2 is NaN";
    if (std::isinf(var2)) errstr << std::endl << "Var 2 is Inf";
    throw affineExcept("PDDouble", "edgeFix", errstr.str());
  }
  
  double istdev1 = 1.0 / sqrt(var1);
  double istdev2 = 1.0 / sqrt(var2);

  //Figure out what indexes these represent in x and y
  double maxfluxfix1, maxfluxfix2;
  int maxidx1, maxidx2;
  maxfluxfix1 = mn1 - PDDouble::lowsigval * sqrt(var1);
  maxfluxfix2 = mn2 - PDDouble::lowsigval * sqrt(var2);
  maxidx1 = static_cast<int>((maxfluxfix1 - minflux1) / dflux1);
  maxidx2 = static_cast<int>((maxfluxfix2 - minflux2) / dflux2);
  maxfluxfix1 = minflux1 + maxidx1 * dflux1;
  maxfluxfix2 = minflux2 + maxidx2 * dflux2;
  
  //Do edges now
  double pdval, tval, preconst, stepfac, subfac, *rowptr;
  if (maxidx1 > 1) {
    int minidx2 = (maxidx2 > 0) ? maxidx2 : 0;
    for (unsigned int j =  minidx2; j < n2; ++j) {
      pdval = pd_[maxidx1 * n2 + j];
      tval = (maxfluxfix1 - mn1) * istdev1;
      preconst = pdval * exp(0.5 * tval * tval);
      subfac = (minflux1 - mn1) * istdev1;
      stepfac = dflux1 * istdev1;
      for (int i = 0; i < maxidx1; ++i) {
	tval = subfac + i * stepfac;
	pd_[i * n2 + j] = preconst * exp(-0.5 * tval * tval);
      }
    }
  }
  if (maxidx2 > 1) {
    int minidx1 = (maxidx1 > 0) ? maxidx1 : 0;
    for (unsigned int i =  minidx1; i < n1; ++i) {
      rowptr = pd_ + i * n2;
      pdval = rowptr[maxidx2];
      tval = (maxfluxfix2 - mn2) * istdev2;
      preconst = pdval * exp(0.5 * tval * tval);
      subfac = (minflux2 - mn2) * istdev2;
      stepfac = dflux2 * istdev2;
      for (int j = 0; j < maxidx2; ++j) {
	tval = subfac + j * stepfac;
	rowptr[j] = preconst * exp(-0.5 * tval * tval);
      }
    }
  }

  //Corner is tricky.
  // We will extrapolate in from both sides (if available)
  // and take the geometric mean
  if (maxidx1 > 0 && maxidx2 > 0) {
    for (int i = 0; i < maxidx1; ++i) {
      rowptr = pd_ + i*n2;
      pdval = rowptr[maxidx2];
      tval = (maxfluxfix2 - mn2) * istdev2;
      preconst = pdval * exp(0.5 * tval * tval);
      subfac = (minflux2 - mn2) * istdev2;
      stepfac = dflux2 * istdev2;
      for (int j = 0; j < maxidx2; ++j) {
	tval = subfac + j * stepfac;
	rowptr[j] = preconst * exp(-0.5 * tval * tval);
      }
    }
    for (int j =  0; j < maxidx2; ++j) {
      pdval = pd_[maxidx1 * n2 + j];
      tval = (maxfluxfix1 - mn1) * istdev1;
      preconst = pdval * exp(0.5 * tval * tval);
      subfac = (minflux1 - mn1) * istdev1;
      stepfac = dflux1 * istdev1;
      for (int i = 0; i < maxidx1; ++i) {
	tval = subfac + i * stepfac;
	pd_[i*n2+j] = sqrt(pd_[i * n2 + j] * preconst * 
			   exp(-0.5 * tval * tval));
      }
    }
  }
}


/*!
  \param[in] donorm Do not assume that P(D) is normalized.
  \returns Pair of the mean along each axis
*/
dblpair PDDouble::getMeans(bool donorm) const {
  if ((n1 == 0) || (n2 == 0)) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());

  //We use the trapezoidal rule here for the integrals
  // so it isn't quite a simple sum
  double xsum, pval, *rowptr;
  double mean1 = 0.0, mean2 = 0.0;

  if (!logflat) {
    //Integrate over lowest x value.  Note that x = 0
    // here (modulo the minflux, which we add later anyways)
    // so there is no mean1 contribution

    // i = 0 term
    // mean2 = 0.5 * 0 * pd[0];
    for (unsigned int j = 1; j < n2 - 1; ++j)
      mean2 += j * pd_[j]; 
    mean2 += 0.5 * (n2 - 1) * pd_[n2 - 1]; //(y=0 at pd_[0])
    mean2 *= 0.5; //End bit of trap in x, so 1/2 factor

    //Now main body of trap
    for (unsigned int i = 1; i < n1 - 1; ++i) {
      rowptr = pd_ + n2 * i;
      xsum = 0.5 * rowptr[0]; //xsum will be the row sum, mult by x later
      //mean2 += 0.5 * 0 * rowptr[0] obviously not needed
      for (unsigned int j = 1; j < n2 - 1; ++j) {
	pval = rowptr[j];
	xsum += pval;
	mean2 += static_cast<double>(j) * pval;
      }      
      xsum += 0.5 * rowptr[n2 - 1];
      mean1 += xsum * static_cast<double>(i); //Multiply in x value
      mean2 += 0.5 * (n2 - 1) * rowptr[n2 - 1];
    }

    //Endpiece (i=n1-1), all multiplied by 1/2 since last x bit
    rowptr = pd_ + (n1 - 1) * n2;
    xsum = 0.5 * rowptr[0];
    for (unsigned int j = 1; j < n2 - 1; ++j) {
      pval = rowptr[j];
      xsum += pval;
      mean2 += 0.5 * static_cast<double>(j) * pval;
    }
    xsum += 0.5 * rowptr[n2 - 1];
    mean1 += 0.5 * xsum * (n1 - 1);
    mean2 += 0.25 * (n2 - 1) * rowptr[n2 - 1];
  } else {
    for (unsigned int j = 1; j < n2 - 1; ++j)
      mean2 += j * exp2(pd_[j]);
    mean2 += 0.5 * (n2-1) * exp2(pd_[n2 - 1]);
    mean2 *= 0.5; 
    for (unsigned int i = 1; i < n1 - 1; ++i) {
      rowptr = pd_ + n2 * i;
      xsum = 0.5 * exp2(rowptr[0]); 
      for (unsigned int j = 1; j < n2 - 1; ++j) {
	pval = exp2(rowptr[j]);
	xsum += pval;
	mean2 += static_cast<double>(j) * pval;
      }      
      pval = exp2(rowptr[n2-1]);
      xsum += 0.5 * pval;
      mean1 += xsum * static_cast<double>(i);
      mean2 += 0.5 * (n2-1) * pval;
    }
    rowptr = pd_ + (n1 - 1) * n2;
    xsum = 0.5 * exp2(rowptr[0]);
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = exp2(rowptr[j]);
      xsum += pval;
      mean2 += 0.5 * static_cast<double>(j) * pval;
    }
    pval = exp2(rowptr[n2-1]);
    xsum += 0.5 * pval;
    mean1 += 0.5 * xsum * (n1 - 1);
    mean2 += 0.25 * (n2 - 1) * pval;
  }

  //Add on step sizes for each integral,
  // which is both area and step size in x,y
  mean1 *= dflux1 * dflux1 * dflux2;
  mean2 *= dflux1 * dflux2 * dflux2;

  if (donorm) {
    double inorm = 1.0 / getIntegral();
    mean1 *= inorm;
    mean2 *= inorm;
  }
  mean1 += minflux1;
  mean2 += minflux2;
  return std::make_pair(mean1, mean2);
}

/*!
  \param[in] donorm Do not assume that P(D) is normalized.
  \returns A pair of pairs, each component of which is the
     mean and variance along first the 1st then 2nd axis
*/
std::pair<dblpair, dblpair> PDDouble::getMeansAndVars(bool donorm) const {
  if ((n1 == 0) || (n2 == 0)) {
    dblpair v = std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			       std::numeric_limits<double>::quiet_NaN());
    return std::make_pair(v, v);
  }

  double normfac = 1.0;
  if (donorm) normfac = 1.0 / getIntegral();

  //First compute means
  //Why not just call getMeans?  To avoid calling getIntegral
  // twice.  After this, mean1 and mean2 will be the actual
  // means/dflux - minflux.
  double xsum, pval, *rowptr;
  double mean1 = 0.0, mean2 = 0.0;
  if (!logflat) {
    for (unsigned int j = 1; j < n2 - 1; ++j)
      mean2 += j*pd_[j]; 
    mean2 += 0.5 * (n2 - 1) * pd_[n2 - 1]; 
    mean2 *= 0.5; 
    for (unsigned int i = 1; i < n1 - 1; ++i) {
      rowptr = pd_ + n2 * i;
      xsum = 0.5 * rowptr[0];
      for (unsigned int j = 1; j < n2 - 1; ++j) {
	pval = rowptr[j];
	xsum += pval;
	mean2 += static_cast<double>(j) * pval;
      }      
      xsum += 0.5 * rowptr[n2 - 1];
      mean1 += xsum * static_cast<double>(i); 
      mean2 += 0.5 * (n2 - 1) * rowptr[n2 - 1];
    }
    rowptr = pd_ + (n1 - 1) * n2;
    xsum = 0.5 * rowptr[0];
    for (unsigned int j = 1; j < n2-1; ++j) {
      pval = rowptr[j];
      xsum += pval;
      mean2 += 0.5 * static_cast<double>(j) * pval;
    }
    xsum += 0.5 * rowptr[n2 - 1];
    mean1 += 0.5 * xsum * (n1 - 1);
    mean2 += 0.25 * (n2 - 1)*rowptr[n2 - 1];
  } else {
    for (unsigned int j = 1; j < n2 - 1; ++j)
      mean2 += j * exp2(pd_[j]);
    mean2 += 0.5 * (n2 - 1) * exp2(pd_[n2 - 1]);
    mean2 *= 0.5; 
    for (unsigned int i = 1; i < n1 - 1; ++i) {
      rowptr = pd_ + n2 * i;
      xsum = 0.5 * exp2(rowptr[0]); 
      for (unsigned int j = 1; j < n2 - 1; ++j) {
	pval = exp2(rowptr[j]);
	xsum += pval;
	mean2 += static_cast<double>(j) * pval;
      }      
      pval = exp2(rowptr[n2 - 1]);
      xsum += 0.5 * pval;
      mean1 += xsum * static_cast<double>(i);
      mean2 += 0.5 * (n2 - 1) * pval;
    }
    rowptr = pd_ + (n1 - 1) * n2;
    xsum = 0.5 * exp2(rowptr[0]);
    for (unsigned int j = 1; j < n2 - 1; ++j) {
      pval = exp2(rowptr[j]);
      xsum += pval;
      mean2 += 0.5 * static_cast<double>(j) * pval;
    }
    pval = exp2(rowptr[n2 - 1]);
    xsum += 0.5 * pval;
    mean1 += 0.5 * xsum * (n1 - 1);
    mean2 += 0.25 * (n2 - 1) * pval;
  }
  mean1 *= dflux1 * dflux2;
  mean2 *= dflux1 * dflux2;

  if (donorm) {
    mean1 *= normfac;
    mean2 *= normfac;
  }
  
  //Now variances, pretty much the same calculation
  // Recall mean1, mean2 are means/dflux - minflux
  double var1 = 0.0, var2 = 0.0;
  double deltax, deltay;
  if (!logflat) {
    //i=0 bit, 1/2 factor for end
    pval = pd_[0];
    xsum = 0.5 * pval;
    var2 = -0.5 * mean2 * mean2 * pval; // deltay = -mean2
    for (unsigned int j = 1; j < n2 - 1; ++j) {
      pval = pd_[j];
      xsum += pval;
      deltay = j - mean2;
      var2 += pval * deltay * deltay;
    }
    pval = pd_[n2 - 1];
    xsum += 0.5 * pval;
    //1/2 factor for end bit, x-<mean> was -mean1 here
    var1 = 0.5 * mean1 * mean1 * xsum; 
    deltay = (n2 - 1) - mean2;
    var2 += 0.5 * deltay * deltay * pval;
    var2 *= 0.5; //1/2 factor for end
    
    //Now core bit
    for (unsigned int i = 1; i < n1 - 1; ++i) {
      rowptr = pd_ + i * n2;
      pval = rowptr[0];
      deltax = i - mean1;
      xsum = 0.5 * pval;
      var2 += 0.5 * mean2 * mean2 * pval; // deltay = -mean2
      for (unsigned int j = 1; j < n2 - 1; ++j) {
	pval = rowptr[j];
	xsum += pval;
	deltay = j - mean2;
	var2 += deltay * deltay * pval;
      }
      pval = rowptr[n2 - 1];
      xsum += 0.5 * pval;
      var1 += xsum * deltax * deltax;
      deltay = n2 - 1 - mean2;
      var2 += 0.5 * deltay * deltay * pval;
    }

    //And now the top end in x
    rowptr = pd_ + (n1 - 1) * n2;
    pval = rowptr[0];
    deltax = n1 - 1 - mean1;
    xsum = 0.5 * pval;
    var2 += 0.25 * mean2 * mean2 * pval;
    for (unsigned int j = 1; j < n2 - 1; ++j) {
      pval = rowptr[j];
      xsum += pval;
      deltay = j - mean2;
      var2 += 0.5 * deltay * deltay * pval;
    }
    pval = rowptr[n2 - 1];
    xsum += 0.5 * pval;
    var1 += 0.5 * xsum * deltax * deltax; //1/2 for end
    deltay = n2 - 1 - mean2;
    var2 += 0.25 * deltay * deltay * pval;  //1/4 is 1/2*1/2 for double end
  } else {
    pval = exp2(pd_[0]);
    xsum = 0.5 * pval;
    var2 = -0.5 * mean2 * mean2 * pval; // deltay = -mean2
    for (unsigned int j = 1; j < n2 - 1; ++j) {
      pval = exp2(pd_[j]);
      xsum += pval;
      deltay = j - mean2;
      var2 += pval * deltay * deltay;
    }
    pval = exp2(pd_[n2 - 1]);
    xsum += 0.5 * pval;
    var1 = 0.5 * mean1 * mean1 * xsum; 
    deltay = (n2 - 1) - mean2;
    var2 += 0.5 * deltay * deltay * pval;
    var2 *= 0.5;
    for (unsigned int i = 1; i < n1-1; ++i) {
      rowptr = pd_ + i * n2;
      pval = exp2(rowptr[0]);
      deltax = i - mean1;
      xsum = 0.5 * pval;
      var2 += 0.5 * mean2 * mean2 * pval; // deltay = -mean2
      for (unsigned int j = 1; j < n2 - 1; ++j) {
	pval = exp2(rowptr[j]);
	xsum += pval;
	deltay = j - mean2;
	var2 += deltay * deltay * pval;
      }
      pval = exp2(rowptr[n2 - 1]);
      xsum += 0.5 * pval;
      var1 += xsum * deltax * deltax;
      deltay = n2 - 1 - mean2;
      var2 += 0.5 * deltay * deltay * pval;
    }
    rowptr = pd_ + (n1 - 1) * n2;
    pval = exp2(rowptr[0]);
    deltax = n1 - 1 - mean1;
    xsum = 0.5 * pval;
    var2 += 0.25 * mean2 * mean2 * pval;
    for (unsigned int j = 1; j < n2 - 1; ++j) {
      pval = exp2(rowptr[j]);
      xsum += pval;
      deltay = j - mean2;
      var2 += 0.5 * deltay * deltay * pval;
    }
    pval = exp2(rowptr[n2 - 1]);
    xsum += 0.5 * pval;
    var1 += 0.5 * xsum * deltax * deltax; //1/2 for end
    deltay = n2 - 1 - mean2;
    var2 += 0.25 * deltay * deltay * pval;  //1/4 is 1/2*1/2 for double end
  }

  //Integral dx dy
  var1 *= dflux1 * dflux2;
  var2 *= dflux1 * dflux2;

  if (donorm) {
    var1 *= normfac;
    var2 *= normfac;
  }

  //Put in step size in x/y bit
  mean1 *= dflux1;
  mean2 *= dflux2;
  var1  *= dflux1 * dflux1;
  var2  *= dflux2 * dflux2;

  mean1 += minflux1;
  mean2 += minflux2;
  return std::make_pair(std::make_pair(mean1, var1),
			std::make_pair(mean2, var2));
}


/*!
  \param[in] other PDDouble to copy from
*/
PDDouble& PDDouble::operator=(const PDDouble& other) {
  if (this == &other) return *this; //Self-copy
  resize(other.n1, other.n2);
  minflux1 = other.minflux1;
  dflux1 = other.dflux1;
  minflux2 = other.minflux2;
  dflux2 = other.dflux2;
  unsigned int sz = n1 * n2;
  if (sz > 0) {
    if (pd_ == NULL)
      throw affineExcept("PDDouble", "operator=", 
			 "Internal storage not initialized");
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] = other.pd_[i];
  }
  logflat = other.logflat;
  return *this;
}

/*!
  \param[in] N1 New number of elements, band 1
  \param[in] MINFLUX1 New minimum flux density, band 1
  \param[in] DFLUX1 New delta flux density, band 1
  \param[in] N2 New number of elements, band 2
  \param[in] MINFLUX2 New minimum flux density, band 2
  \param[in] DFLUX2 New delta flux density, band 2
  \param[in] PD New values in row-major order
  \param[in] LOG Is the input PD log2?
*/
void PDDouble::fill(unsigned int N1, double MINFLUX1, double DFLUX1,
		    unsigned int N2, double MINFLUX2, double DFLUX2,
		    const double* const PD, bool LOG) {
  logflat = LOG;
  resize(N1, N2);
  minflux1 = MINFLUX1;
  dflux1 = DFLUX1;
  minflux2 = MINFLUX2;
  dflux2 = DFLUX2;
  unsigned int sz = n1*n2;
  if (sz > 0) {
    if (pd_ == NULL)
      throw affineExcept("PDDouble", "operator=", 
			 "Internal storage not initialized");
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] = PD[i];
  }
}

/*!
  \param[in] x Flux density in band 1 to evaluate P(D) at
  \param[in] y Flux density in band 2 to evaluate P(D) at
  \param[in] logval Return log (base e) value
  \returns Interpolated P(D) value.  

  Uses linear interpolation.  Note that this will work differently
  if the P(D) has been converted to log values.
*/
double PDDouble::getPDVal(double x, double y, bool logval) const {
  if (pd_ == NULL) return std::numeric_limits<double>::quiet_NaN();

  //look up the effective indexes
  int idx1 = static_cast<int>((x - minflux1) / dflux1);
  int idx2 = static_cast<int>((y - minflux2) / dflux2);
  int n2idx1 = n2 * idx1;

  double maxflux1 = minflux1 + static_cast<double>(n1 - 1) * dflux1;
  double maxflux2 = minflux2 + static_cast<double>(n2 - 1) * dflux2;

  unsigned int n2n1 = n2 * n1;
  double interp_val;
  //Check to see if we are off the edge
  if (x < minflux1) {
    if (y < minflux2) interp_val = pd_[0];
    else if (y > maxflux2) interp_val = pd_[n2 - 1];
    else interp_val = pd_[idx2];
  } else if (x > maxflux1) {
    if (y < minflux2) interp_val = pd_[n2n1 - n1];
    else if (y > maxflux2) interp_val = pd_[n2n1 - 1];
    else interp_val =  pd_[n2n1 - n1 + idx2];
  } else if (y < minflux2) {
    interp_val =  pd_[n2idx1];
  } else if (y > maxflux2) {
    interp_val = pd_[n2idx1 + n2 - 1];
  } else {
    //Actual interpolation
    double u,t,omu,omt;
    t = (x - minflux1) / dflux1 - static_cast<double>(idx1);
    u = (y - minflux2) / dflux2 - static_cast<double>(idx2);
    omu = 1.0 - u; omt = 1.0-t;
    
    unsigned int baseidx = n2idx1 + idx2;
    interp_val = omt*(omu*pd_[baseidx] + u*pd_[baseidx+1]) +
      t*(omu*pd_[baseidx+n2] + u*pd_[baseidx+n2+1]);
  }
  if (!logval) {
    if (logflat) return exp2(interp_val); else return interp_val;
  } else {
    //Note we return ln, not log2
    if (logflat) return pofd_mcmc::log2toe * interp_val; 
    else return log(interp_val);
  }
}

/*!
  \param[inout] os Stream to output values to
*/
std::ostream& PDDouble::writeToStream(std::ostream& os) const {
  os << n1 << " " << minflux1 << " " << dflux1 << std::endl;
  os << n2 << " " << minflux2 << " " << dflux2 << std::endl;
  double *rowptr;
  if (n1*n2 > 0) {
    for (unsigned int i = 0; i < n1; ++i) {
      rowptr = pd_ + i * n2;
      os << rowptr[0];
      for (unsigned int j = 1; j < n2; ++j)
	os << " " << rowptr[j];
      os << std::endl;
    }
  }
  return os;
}

/*!
  \param[in] outputfile File to write to in FITS format
  \returns 0 on success, an error code (!=0) for anything else
 */
int PDDouble::writeToFits(const std::string& outputfile) const {

  //Make the fits file
  int status = 0;
  fitsfile *fp;

  fits_create_file(&fp, outputfile.c_str(), &status);

  if (status) {
    fits_report_error(stderr,status);
    throw affineExcept("PDDouble", "writeToFits",
		       "Error creating FITS output file");
  }

  long axissize[2];
  axissize[0] = static_cast<long>(n1);
  axissize[1] = static_cast<long>(n2);
  
  fits_create_img(fp, DOUBLE_IMG, 2, axissize, &status);
  
  //Add "WCS" info to hdr
  float crpix = 1;
  double tmpval;
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
		 const_cast<char*>("FLUX1"),
		 const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		 const_cast<char*>("Ref pix of axis 1"), &status);
  tmpval = minflux1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL1"), &tmpval, 
		 const_cast<char*>("val at ref pix"), &status);
  tmpval = dflux1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CDELT1"), &tmpval,
		 const_cast<char*>("delta along axis 1"), &status);

  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE2"),
		 const_cast<char*>("FLUX2"),
		 const_cast<char*>("Type of Data axis 2"),&status);
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX2"), &crpix, 
		 const_cast<char*>("Ref pix of axis 2"), &status);
  tmpval = minflux2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL2"), &tmpval, 
		 const_cast<char*>("val at ref pix"), &status);
  tmpval = dflux2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CDELT2"), &tmpval,
		 const_cast<char*>("delta along axis 2"), &status);
  tmpval = dflux1;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD1_1"), &tmpval,
		 const_cast<char*>("WCS matrix element 1 1"), &status);
  tmpval = 0.0;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD1_2"), &tmpval, 
		 const_cast<char*>("WCS matrix element 1 2"),
		 &status);
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD2_1"), &tmpval, 
		 const_cast<char*>("WCS matrix element 2 1"),
		 &status);
  tmpval = dflux2;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CD2_2"), &tmpval, 
		 const_cast<char*>("WCS matrix element 2 2"),
		 &status);

  int lg = static_cast<int>(logflat);
  fits_write_key(fp, TLOGICAL, const_cast<char*>("LOG"),&lg,
		 const_cast<char*>("Is log P(D) stored?"), &status);

  //Do data writing.  We have to make a transposed copy of the
  // data to do this, which is irritating as hell
  double *tmpdata = new double[ n1 ];
  long fpixel[2] = { 1, 1 };
  for (unsigned int j = 0; j < n2; ++j) {
    for (unsigned int i = 0; i < n1; ++i) tmpdata[i] = pd_[i*n2+j];
    fpixel[1] = static_cast<long>(j+1);
    fits_write_pix(fp, TDOUBLE, fpixel, n1, tmpdata, &status);
  }
  delete[] tmpdata;

  fits_close_file(fp,&status);

  if (status) {
    fits_report_error(stderr,status);
    throw affineExcept("PDDouble", "writeToFits", "Error doing FITS write");
  }
  return status;
}

/*!
  \param[in] outputfile File to write to
*/
void PDDouble::writeToHDF5(const std::string& outputfile) const {
  hid_t file_id;
  file_id = H5Fcreate(outputfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);
  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("PDDouble", "writeToHDF5",
		       "Failed to open HDF5 file to write");
  }

  hsize_t adims;
  hid_t mems_id, att_id;
  
  // Properties
  hdf5utils::writeAttBool(file_id, "isLog", logflat);

  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(file_id, "dflux1", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux1);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "dflux2", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux2);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "minflux1", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &minflux1);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "minflux2", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &minflux2);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  
  // Rfluxes -- by making temporary array
  unsigned int maxn = n1 > n2 ? n1 : n2;
  double *flux = new double[maxn];
  for (unsigned int i = 0; i < n1; ++i) 
    flux[i] = static_cast<double>(i) * dflux1 + minflux1;
  hdf5utils::writeDataDoubles(file_id, "flux1", n1, flux);
  for (unsigned int i = 0; i < n2; ++i) 
    flux[i] = static_cast<double>(i) * dflux2 + minflux2;
  hdf5utils::writeDataDoubles(file_id, "flux2", n2, flux);
  delete[] flux;

  hdf5utils::writeData2DDoubles(file_id, "PD", n1, n2, pd_);

  H5Fclose(file_id);
}

/*!
  \param[in] data Data to evaluate the log Likelihood of
  \returns Log likelihood

  Because this uses interpolation, slightly different values
  will result if the P(D) is stored in log form or not.
*/
double PDDouble::getLogLike(const fitsDataDouble& data) const {
  if (pd_ == NULL) throw affineExcept("PDDouble", "getLogLike",
				      "pd not filled before likelihood calc");
  unsigned int ndata = data.getN();
  if (ndata == 0) throw affineExcept("PDDouble", "getLogLike",
				     "No data present");

  if (data.isBinned()) return getLogLikeBinned(data);
  else return getLogLikeUnbinned(data);
}

/*!
  \param[in] data Data to evaluate the log Likelihood of
  \returns Log likelihood

  Because this uses interpolation, slightly different values
  will result if the P(D) is stored in log form or not.
*/
double PDDouble::getLogLikeUnbinned(const fitsDataDouble& data) const {
  unsigned int ndata = data.getN();

  //Quantities for edge test
  double maxflux1 = minflux1 + static_cast<double>(n1 - 1) * dflux1;
  double maxflux2 = minflux2 + static_cast<double>(n2 - 1) * dflux2;

  int idx1, idx2, n2idx1; //!< Index look up
  unsigned int n2n1, baseidx;
  n2n1 = n2 * n1;

  const double* flux1;
  const double* flux2;
  double cflux1, cflux2, loglike, interp_val, delt1, delt2;
  double u, t, omu, omt;
  double idflux1 = 1.0 / dflux1;
  double idflux2 = 1.0 / dflux2;

  loglike = 0.0;
  flux1 = data.getData1();
  flux2 = data.getData2();

  //Do interpolation.  Note we don't call getPDval, since
  // we do the interpolation here always in log space no matter
  // how it is stored internally, and because it's more efficient
  // to do it in house.
  if (logflat) {
    //Internal information is stored as log2 of P(D)
    for (unsigned int i = 0; i < ndata; ++i) {
      cflux1 = flux1[i]; cflux2 = flux2[i];
      //Get effective indices
      delt1 = (cflux1 - minflux1) * idflux1;
      delt2 = (cflux2 - minflux2) * idflux2;
      idx1 = static_cast<int>(delt1);
      idx2 = static_cast<int>(delt2);
      n2idx1 = n2 * idx1;
      if (cflux1 <= minflux1) {
	if (cflux2 <= minflux2) interp_val = pd_[0];
        else if (cflux2 >= maxflux2) interp_val = pd_[n2 - 1];
        else interp_val = pd_[idx2];
      } else if (cflux1 >= maxflux1) {
        if (cflux2 <= minflux2) interp_val = pd_[n2n1 - n1];
        else if (cflux2 >= maxflux2) interp_val = pd_[n2n1 - 1];
        else interp_val = pd_[n2n1 - n1 + idx2];
      } else if (cflux2 <= minflux2) {
        interp_val = pd_[n2idx1];
      } else if (cflux2 >= maxflux2) {
        interp_val = pd_[n2idx1 + n2 - 1];
      } else {
        //Not off edge
        t = delt1 - static_cast<double>(idx1);
        u = delt2 - static_cast<double>(idx2);
        omu = 1.0 - u; omt = 1.0 - t;
        baseidx = n2idx1 + idx2;
        interp_val = omt * (omu * pd_[baseidx] + u * pd_[baseidx + 1]) +
          t * (omu * pd_[baseidx + n2] + u * pd_[baseidx + n2 + 1]);
      }
      loglike += interp_val;
    }
  } else {
    //Not stored as log2 -- inefficient, but supported
    //Note that it would be really stupid to do this multiplicatively,
    // then take the log.  Also, it's better to interpolate
    // in log space than interpolate, then log
    for (unsigned int i = 0; i < ndata; ++i) {
      cflux1 = flux1[i]; cflux2 = flux2[i];
      delt1 = (cflux1 - minflux1) * idflux1;
      delt2 = (cflux2 - minflux2) * idflux2;
      idx1 = static_cast<int>(delt1);
      idx2 = static_cast<int>(delt2);
      n2idx1 = n2 * idx1;
      
      if (cflux1 < minflux1) {
        if (cflux2 < minflux2) interp_val = log2(pd_[0]);
        else if (cflux2 > maxflux2) interp_val = log2(pd_[n2 - 1]);
        else interp_val = log2(pd_[idx2]);
      } else if (cflux1 > maxflux1) {
        if (cflux2 < minflux2) interp_val = log2(pd_[n2n1 - n1]);
        else if (cflux2 > maxflux2 ) interp_val = log2(pd_[n2n1 - 1]);
        else interp_val = log2(pd_[n2n1 - n1 + idx2]);
      } else if (cflux2 < minflux2) {
        interp_val = log2(pd_[n2idx1]);
      } else if (cflux2 > maxflux2) {
        interp_val = log2(pd_[n2idx1 + n2 - 1]);
      } else {
        //Not off edge
        t = delt1 - static_cast<double>(idx1);
        u = delt2 - static_cast<double>(idx2);
        omu = 1.0 - u; omt = 1.0 - t;
        
        baseidx = n2idx1 + idx2;
        interp_val = 
	  omt * (omu * log2(pd_[baseidx]) + u * log2(pd_[baseidx + 1])) +
          t* (omu * log2(pd_[baseidx + n2]) + u * log2(pd_[baseidx + n2 + 1]));
      }
      loglike += interp_val;
    }
  }
  //This has been base 2 -- convert back to base e
  return pofd_mcmc::log2toe * loglike;
}

/*!
  \param[in] data Data to evaluate the log Likelihood of
  \returns Log likelihood

  Because this uses interpolation, slightly different values
  will result if the P(D) is stored in log form or not.
*/
double PDDouble::getLogLikeBinned(const fitsDataDouble& data) const {
  if (!data.isBinned())
    throw affineExcept("PDDouble", "getLogLikeBinned", "Data is not binned");

  std::pair<unsigned int,unsigned int> nbins = data.getNBins();
  unsigned int nbins1 = nbins.first;
  unsigned int nbins2 = nbins.second;

  //Quantities for edge test
  double maxflux1 = minflux1 + static_cast<double>(n1-1)*dflux1;
  double maxflux2 = minflux2 + static_cast<double>(n2-1)*dflux2;

  int idx1, idx2, n2idx1; //!< Index look up
  unsigned int baseidx, n2n1;
  n2n1 = n2*n1;

  const unsigned int *bins, *binptr;
  unsigned int ninbin;
  double cflux1, cflux2, bincent01, bincent02, bindelta1, bindelta2;
  double loglike, interp_val, idflux1, idflux2;
  double u,t,omu,omt;

  loglike = 0.0;
  std::pair<double,double> pr = data.getBinCent0();
  bincent01 = pr.first;
  bincent02 = pr.second;
  pr = data.getBinDelta();
  bindelta1 = pr.first;
  bindelta2 = pr.second;
  bins = data.getBinnedData();
  idflux1 = 1.0/dflux1;
  idflux2 = 1.0/dflux2;

  //Do interpolation.  Note we don't call getPDval, since
  // we do the interpolation here always in log space no matter
  // how it is stored internally, and because it's more efficient
  // to do it in house.
  if (logflat) {
    //Internal information is stored as log2 of P(D)
    for (unsigned int i = 0; i < nbins1; ++i) {
      cflux1 = bincent01 + static_cast<double>(i)*bindelta1;
      idx1 = static_cast<int>((cflux1 - minflux1) * idflux1);
      n2idx1 = n2 * idx1;
      binptr = bins + i * nbins2;
      for (unsigned int j = 0; j < nbins2; ++j) {
	ninbin = binptr[j];
	if (ninbin == 0) continue;  //skip calculation
	cflux2 = bincent02 + static_cast<double>(j) * bindelta2;
	idx2 = static_cast<int>((cflux2 - minflux2) * idflux2);
	if (cflux1 < minflux1) {
	  if (cflux2 < minflux2) interp_val = pd_[0];
	  else if (cflux2 > maxflux2) interp_val = pd_[n2 - 1];
	  else interp_val = pd_[idx2];
	} else if (cflux1 > maxflux1) {
	  if (cflux2 < minflux2) interp_val = pd_[n2n1 - n1];
	  else if (cflux2 > maxflux2) interp_val = pd_[n2n1 - 1];
	  else interp_val = pd_[n2n1 - n1 + idx2];
	} else if (cflux2 < minflux2) {
	  interp_val = pd_[n2idx1];
	} else if (cflux2 > maxflux2) {
	  interp_val = pd_[n2idx1 + n2 - 1];
	} else {
	  //Not off edge
	  t = (cflux1 - minflux1) * idflux1 - static_cast<double>(idx1);
	  u = (cflux2 - minflux2) * idflux2 - static_cast<double>(idx2);
	  omu = 1.0 - u; omt = 1.0 - t;
	  baseidx = n2idx1 + idx2;
	  interp_val = 
	    omt * (omu * pd_[baseidx] + u * pd_[baseidx+1]) +
	    t * (omu * pd_[baseidx+n2] + u * pd_[baseidx+n2+1]);
	}
	loglike += static_cast<double>(ninbin) * interp_val;
      }
    }
  } else {
    //Not stored as log2 -- inefficient, but supported
    //Note that it would be insane to do this multiplicatively,
    // then take the log.  Also, it's better to interpolate
    // in log space than interpolate, then log
    for (unsigned int i = 0; i < nbins1; ++i) {
      cflux1 = bincent01 + static_cast<double>(i) * bindelta1;
      idx1 = static_cast<int>((cflux1 - minflux1) * idflux1);
      n2idx1 = n2 * idx1;
      binptr = bins + i * nbins2;
      for (unsigned int j = 0; j < nbins2; ++j) {
	ninbin = binptr[j];
	if (ninbin == 0) continue;  //skip calculation
	cflux2 = bincent02 + static_cast<double>(j) * bindelta2;
	idx2 = static_cast<int>((cflux2 - minflux2) * idflux2);
	
	if (cflux1 < minflux1) {
	  if (cflux2 < minflux2) interp_val = log2(pd_[0]);
	  else if (cflux2 > maxflux2) interp_val = log2(pd_[n2 - 1]);
	  else interp_val = log2(pd_[idx2]);
	} else if (cflux1 > maxflux1) {
	  if (cflux2 < minflux2) interp_val = log2(pd_[n2n1 - n1]);
	  else if (cflux2 > maxflux2) interp_val = log2(pd_[n2n1 - 1]);
	  else interp_val = log2(pd_[n2n1 - n1 + idx2]);
	} else if (cflux2 < minflux2) {
	  interp_val = log2(pd_[n2idx1]);
	} else if (cflux2 > maxflux2) {
	  interp_val = log2(pd_[n2idx1 + n2 - 1]);
	} else {
	  //Not off edge
	  t = (cflux1 - minflux1) * idflux1 - idx1;
	  u = (cflux2 - minflux2) * idflux2 - idx2;
	  omu = 1.0 - u; omt = 1.0-t;
	  
	  baseidx = n2idx1 + idx2;
	  interp_val = 
	    omt * (omu * log2(pd_[baseidx]) + u * log2(pd_[baseidx+1])) +
	    t * (omu * log2(pd_[baseidx+n2]) + u * log2(pd_[baseidx+n2+1]));
	}
	loglike += static_cast<double>(ninbin) * interp_val;
      }
    }
  }
  //This has been base 2 -- convert back to base e
  return pofd_mcmc::log2toe * loglike;
}

/*!
  \param[inout] os Stream to write to
  \param[in] b PD to output
*/
std::ostream& operator<<(std::ostream& os, const PDDouble& b) {
  b.writeToStream(os);
  return os;
}
