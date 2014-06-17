//Program for testing the behavior of the ensemble code

//Right now only does affineQueue

#include<iostream>
#include<sstream>
#include<utility>
#include<sstream>
#include<cmath>
#include<limits>

#include "../include/affineEnsemble.h"
#include "../include/affineExcept.h"
#include "../include/ran.h"

//We need a subclass of affineEnsemble to test
class testEnsemble : public affineEnsemble {
public:
  testEnsemble(unsigned int, unsigned int, unsigned int);
  ~testEnsemble() {};

  void initChains();
  void generateInitialPosition(const paramSet&);
  double getLogLike(const paramSet&, bool&);
};

testEnsemble::testEnsemble(unsigned int NWALKERS, unsigned int NPARAMS,
			   unsigned int NSAMPLES) :
  affineEnsemble(NWALKERS,NPARAMS,NSAMPLES) {}

void testEnsemble::initChains() {
  //Just generate random positions from zero to one for
  // initial vectors
  is_init = true;
  if (rank != 0) return;

  unsigned int npar = getNParams();
  paramSet p(npar);
  for (unsigned int i = 0; i < npar; ++i)
    p[i] = 0.5;
  generateInitialPosition(p);
}

void testEnsemble::generateInitialPosition(const paramSet& p) {
  if (rank != 0) return;

  //For master
  chains.clear();
  chains.addChunk(1);

  unsigned int npar = getNParams();
  paramSet p2(npar);
  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    for (unsigned int j = 0; j < npar; ++j)
      p2[j] = rangen.doub() - 0.5 + p[i];
    chains.addNewStep(i, p2, -std::numeric_limits<double>::infinity());
    naccept[i] = 0;
  }
}

double testEnsemble::getLogLike(const paramSet& p, bool& b) {
  b = false;
  return 1.0;
}

/////////////////////////////////

int main(int argc, char** argv) {
  
  unsigned int rank;
  MPI::Init(argc,argv);
  rank = MPI::COMM_WORLD.Get_rank();

  if (rank == 0) {
    try {
      ran ranstruct;
      std::cout << "affineEnsemble tests" << std::endl;
      
      testEnsemble ensemble(250,4,1000);
      
      if (ensemble.getNWalkers() != 250)
	throw affineExcept("ensemble_test", "main",
			   "ensemble should have 250 walkers");
      if (ensemble.getNParams() != 4)
	throw affineExcept("ensemble_test", "main",
			   "ensemble should have 4 params");
      if (ensemble.getNChunks() != 0)
	throw affineExcept("ensemble_test", "main",
			   "ensemble should not have any chunks");
      
      //Initialize
      ensemble.initChains();
      if (ensemble.getNChunks() != 1)
	throw affineExcept("ensemble_test", "main",
			   "ensemble should now have one chunk");
      
      //Test the random number generator
      const unsigned int nsamples = 10000000;
      double *zvals = new double[nsamples];
      for (unsigned int i = 0; i < nsamples; ++i)
	zvals[i] = ensemble.generateZ();
      
      //Test first, second third moment (not central moments!)
      double mom1, mom2, mom3, mom4, mom5, dval;
      dval = zvals[0];
      mom1 = dval; mom2 = dval*dval; 
      mom3 = mom2*dval; mom4 = mom3*dval;
      mom5 = mom4*dval;
      for (unsigned int i=1; i < nsamples; ++i) {
	dval = zvals[i];
	mom1 += dval; mom2 += dval*dval; 
	mom3 += dval*dval*dval; mom4 += dval*dval*dval*dval;
	mom5 += dval*dval*dval*dval*dval;
      }
      delete[] zvals;
      double nfac = 1.0/static_cast<double>(nsamples);
      mom1 *= nfac; mom2 *= nfac; mom3 *= nfac; mom4 *= nfac; mom5 *= nfac;
      
      //Can't expect exact result -- so sometimes these tests
      // will fail even if the code is working!
      double a = ensemble.getScalefac();
      double prefac = 0.5*sqrt(a)/(a-1.0);
      double target = prefac*( std::pow(a,1.5)-std::pow(a,-1.5) )/1.5;
      if ( fabs( (mom1 - target)/mom1 ) > 1e-3 ) {
	std::stringstream errstr;
	errstr << "generateZ failed first moment test -- expected "
	       << target << " got " << mom1 << std::endl;
	throw affineExcept("ensemble_test", "main", errstr.str());
      }
      
      target = prefac*( std::pow(a,2.5)-std::pow(a,-2.5) )/2.5;
      if ( fabs( (mom2 - target)/mom2 ) > 1e-3) {
	std::stringstream errstr;
	errstr << "generateZ failed second moment test -- expected "
	       << target << " got " << mom2 << std::endl;
	throw affineExcept("ensemble_test", "main", errstr.str());
      }
      
      target = prefac*( std::pow(a,3.5)-std::pow(a,-3.5) )/3.5;
      if ( fabs( (mom3 - target)/mom3 ) > 1e-3 ) {
	std::stringstream errstr;
	errstr << "generateZ failed third moment test -- expected "
	       << target << " got " << mom3 << std::endl;
	throw affineExcept("ensemble_test", "main", errstr.str());
      }

      target = prefac*( std::pow(a,4.5)-std::pow(a,-4.5) )/4.5;
      if ( fabs( (mom4 - target)/mom4 ) > 1e-3 ) {
	std::stringstream errstr;
	errstr << "generateZ failed fourth moment test -- expected "
	       << target << " got " << mom4 << std::endl;
	throw affineExcept("ensemble_test", "main", errstr.str());
      }

      target = prefac*( std::pow(a,5.5)-std::pow(a,-5.5) )/5.5;
      if ( fabs( (mom5 - target)/mom5 ) > 1e-3 ) {
	std::stringstream errstr;
	errstr << "generateZ failed fifth moment test -- expected "
	       << target << " got " << mom5 << std::endl;
	throw affineExcept("ensemble_test", "main", errstr.str());
      }
      
      //Resizing tests
      testEnsemble tens(250,2,1000);
      if (tens.getNSteps() != 4) {
	std::stringstream errstr;
	errstr << "Nsteps should be 4, but is: " << tens.getNSteps();
	throw affineExcept("ensemble_test", "main", errstr.str());
      }

      tens.setNWalkers(100);
      if (tens.getNWalkers() != 100) {
	std::stringstream errstr;
	errstr << "Nwalkers should be 100, but is: " << tens.getNWalkers();
	throw affineExcept("ensemble_test", "main", errstr.str());
      }
      if (tens.getNParams() != 2) {
	std::stringstream errstr;
	errstr << "Nparams should be 2, but is: " << tens.getNParams();
	throw affineExcept("ensemble_test", "main", errstr.str());
      }
      if (tens.getNSteps() != 10) {
	std::stringstream errstr;
	errstr << "Nsteps should be 10, but is: " << tens.getNSteps();
	throw affineExcept("ensemble_test", "main", errstr.str());
      }

      //Add some steps, resize again, make sure they go away
      tens.initChains();
      if (tens.getNChunks() != 1)
	throw affineExcept("ensemble_test", "main",
			   "ensemble should now have one chunk");
      tens.setNParams(4);
      if (tens.getNParams() != 4) {
	std::stringstream errstr;
	errstr << "Nparams should be 4, but is: " << tens.getNParams();
	throw affineExcept("ensemble_test", "main", errstr.str());
      }
      if (tens.getNChunks() != 0)
	throw affineExcept("ensemble_test", "main",
			   "After change in nparams, nchunks should be 0");

      std::cout << "All ensemble tests passed" << std::endl;
      
    } catch ( const affineExcept& ex ) {
      std::cerr << "Error encountered" << std::endl;
      std::cerr << ex << std::endl;
      return 1;
    } catch ( const std::bad_alloc& ex ) {
      std::cerr << "Memory allocation problem" << std::endl;
      std::cerr << ex.what() << std::endl;
      return 1;
    }
  }
  MPI::Finalize();
  return 0;
}
