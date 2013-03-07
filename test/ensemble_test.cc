//Program for testing the behavior of the ensemble code

//Right now only does affineQueue

#include<iostream>
#include<sstream>
#include<utility>
#include<sstream>
#include<cmath>

#include "../include/affineEnsemble.h"
#include "../include/affineExcept.h"
#include "../include/ran.h"

//We need a subclass of affineEnsemble to test
class testEnsemble : public affineEnsemble {
public:
  testEnsemble(unsigned int, unsigned int, unsigned int);
  ~testEnsemble() {};

  void initChains();
  double getLogLike(const paramSet&);
};

testEnsemble::testEnsemble(unsigned int NWALKERS, unsigned int NPARAMS,
			   unsigned int NSAMPLES) :
  affineEnsemble(NWALKERS,NPARAMS,NSAMPLES) {}

void testEnsemble::initChains() {
  //Just generate random positions from zero to one for
  // initial vectors

  if (rank != 0) return;

  //For master
  chains.clear();
  chains.addChunk(1);
    
  paramSet p( getNParams() );
  double logLike;
  unsigned int npar = getNParams();
  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    for (unsigned int j = 0; j < npar; ++j)
      p[j] = rangen.doub();
    logLike = getLogLike(p);
    chains.addNewStep( i, p, logLike );
    naccept[i] = 1;
  }
}

double testEnsemble::getLogLike(const paramSet& p) {
  return 1.0;
}

/////////////////////////////////

int main(int argc, char** argv) {
  
  const unsigned int cap = 25;

  unsigned int rank;
  MPI::Init(argc,argv);
  rank = MPI::COMM_WORLD.Get_rank();

  if (rank == 0) {
    try {
      ran ranstruct;

      std::cout << "affineQueue tests" << std::endl;
      affineQueue< unsigned int > queue(cap);
      
      //Insert 25, make sure that's okay
      for (unsigned int i = 0; i < 5; ++i)
	queue.push(i);
      for (unsigned int i = 5; i < cap; ++i)
	queue.push( ranstruct.int32() );
      if (queue.size() != cap)
	throw affineExcept("ensemble_test","main",
			   "Queue should have 25 elems in it",1);
      if (queue.capacity() != cap)
	throw affineExcept("ensemble_test","main",
			   "Queue should have capacity of 25",2);
      
      unsigned int val;
      for (unsigned int i = 0; i < 5; ++i) {
	val = queue.pop();
	if (val != i) {
	  std::stringstream str;
	  str << "Queue pop should have yielded: " << i << " but gave "
	      << val;
	  throw affineExcept("ensemble_test","main",
			     str.str(),3);
	}
      }
      if (queue.size() != 20)
	throw affineExcept("ensemble_test","main",
			   "Queue should have 20 elements",4);
      
      //Pop the rest
      for (unsigned int i = 5; i < cap; ++i)
	val = queue.pop();
      if (queue.size() != 0)
	throw affineExcept("ensemble_test","main",
			   "Queue should have 0 elements",5);
      if (!queue.empty())
	throw affineExcept("ensemble_test","main",
			   "Queue should be empty",6);
      if (queue.capacity() != cap)
	throw affineExcept("ensemble_test","main",
			   "Queue should have capacity of 25",7);
      
      
      //Test clear
      for (unsigned int i = 0; i < cap-3; ++i)
	queue.push( ranstruct.int32() );
      if (queue.size() != cap-3)
	throw affineExcept("ensemble_test","main",
			   "Queue should have 22 elements",8);
      if (queue.capacity() != cap)
	throw affineExcept("ensemble_test","main",
			   "Queue should have capacity of 25",9);
      queue.clear();
      if (queue.size() != 0)
	throw affineExcept("ensemble_test","main",
			   "Queue should have 0 elements",10);
      if (!queue.empty())
	throw affineExcept("ensemble_test","main",
			   "Queue should be empty",11);
      if (queue.capacity() != cap)
	throw affineExcept("ensemble_test","main",
			   "Queue should have capacity of 25",12);
      
      
      //Try a more complicated element
      affineQueue< std::pair<unsigned int, int> > queue2(cap);
      if (queue2.capacity() != cap)
	throw affineExcept("ensemble_test","main",
			   "Queue2 should have capacity of 25",13);
      if (!queue2.empty())
	throw affineExcept("ensemble_test","main",
			   "Queue2 should be empty",14);
      std::pair<unsigned int, int> val2,val3;
      val2.first = 4;
      val2.second = -2;
      queue2.push(val2);
      val2.first = 5;
      val2.second = -4;
      queue2.push(val2);
      for (unsigned int i = 2; i < 10; ++i) {
	val3.first = ranstruct.int32();
	val3.second = ranstruct.int32();
	queue2.push(val3);
      }
      if (queue2.capacity() != cap)
	throw affineExcept("ensemble_test","main",
			   "Queue2 should have capacity of 25",15);
      if (queue2.size() != 10)
	throw affineExcept("ensemble_test","main",
			   "Queue2 should have 10 elements",16);
      val3 = queue2.pop();
      if (queue2.size() != 9)
	throw affineExcept("ensemble_test","main",
			   "Queue2 should have 9 elements",17);
      val3 = queue2.pop();
      if (queue2.size() != 8)
	throw affineExcept("ensemble_test","main",
			   "Queue2 should have 8 elements",18);
      if (val3.first != val2.first)
	throw affineExcept("ensemble_test","main",
			   "First elem of second pop not as expected",19);
      if (val3.second != val2.second)
	throw affineExcept("ensemble_test","main",
			   "Second elem of second pop not as expected",20);
      
      queue2.setCapacity(5);
      if (queue2.capacity() != 5)
	throw affineExcept("ensemble_test","main",
			   "Queue2 should have capacity of 5",21);
      if (!queue2.empty())
	throw affineExcept("ensemble_test","main",
			   "Queue2 should be empty",22);
      

      std::cout << "affineEnsemble tests" << std::endl;
      
      testEnsemble ensemble(250,4,1000);
      
      if (ensemble.getNWalkers() != 250)
	throw affineExcept("ensemble_test","main",
			   "ensemble should have 250 walkers",23);
      if (ensemble.getNParams() != 4)
	throw affineExcept("ensemble_test","main",
			   "ensemble should have 4 params",24);
      if (ensemble.getNChunks() != 0)
	throw affineExcept("ensemble_test","main",
			   "ensemble should not have any chunks",24);
      
      //Initialize
      ensemble.initChains();
      if (ensemble.getNChunks() != 1)
	throw affineExcept("ensemble_test","main",
			   "ensemble should now have one chunk",25);
      
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
	throw affineExcept("ensemble_test","main",errstr.str(),26);
      }
      
      target = prefac*( std::pow(a,2.5)-std::pow(a,-2.5) )/2.5;
      if ( fabs( (mom2 - target)/mom2 ) > 1e-3) {
	std::stringstream errstr;
	errstr << "generateZ failed second moment test -- expected "
	       << target << " got " << mom2 << std::endl;
	throw affineExcept("ensemble_test","main",
			   errstr.str(),27);
      }
      
      target = prefac*( std::pow(a,3.5)-std::pow(a,-3.5) )/3.5;
      if ( fabs( (mom3 - target)/mom3 ) > 1e-3 ) {
	std::stringstream errstr;
	errstr << "generateZ failed third moment test -- expected "
	       << target << " got " << mom3 << std::endl;
	throw affineExcept("ensemble_test","main",
			   errstr.str(),28);
      }

      target = prefac*( std::pow(a,4.5)-std::pow(a,-4.5) )/4.5;
      if ( fabs( (mom4 - target)/mom4 ) > 1e-3 ) {
	std::stringstream errstr;
	errstr << "generateZ failed fourth moment test -- expected "
	       << target << " got " << mom4 << std::endl;
	throw affineExcept("ensemble_test","main",
			   errstr.str(),29);
      }

      target = prefac*( std::pow(a,5.5)-std::pow(a,-5.5) )/5.5;
      if ( fabs( (mom5 - target)/mom5 ) > 1e-3 ) {
	std::stringstream errstr;
	errstr << "generateZ failed fifth moment test -- expected "
	       << target << " got " << mom5 << std::endl;
	throw affineExcept("ensemble_test","main",
			   errstr.str(),30);
      }
      
      //Resizing tests
      testEnsemble tens(250,2,1000);
      if (tens.getNSteps() != 4) {
	std::stringstream errstr;
	errstr << "Nsteps should be 4, but is: " << tens.getNSteps();
	throw affineExcept("ensemble_test", "main", errstr.str(), 31);
      }

      tens.setNWalkers(100);
      if (tens.getNWalkers() != 100) {
	std::stringstream errstr;
	errstr << "Nwalkers should be 100, but is: " << tens.getNWalkers();
	throw affineExcept("ensemble_test", "main", errstr.str(), 32);
      }
      if (tens.getNParams() != 2) {
	std::stringstream errstr;
	errstr << "Nparams should be 2, but is: " << tens.getNParams();
	throw affineExcept("ensemble_test", "main", errstr.str(), 33);
      }
      if (tens.getNSteps() != 10) {
	std::stringstream errstr;
	errstr << "Nsteps should be 10, but is: " << tens.getNSteps();
	throw affineExcept("ensemble_test", "main", errstr.str(), 34);
      }

      //Add some steps, resize again, make sure they go away
      tens.initChains();
      if (tens.getNChunks() != 1)
	throw affineExcept("ensemble_test","main",
			   "ensemble should now have one chunk",35);
      tens.setNParams(4);
      if (tens.getNParams() != 4) {
	std::stringstream errstr;
	errstr << "Nparams should be 4, but is: " << tens.getNParams();
	throw affineExcept("ensemble_test", "main", errstr.str(), 36);
      }
      if (tens.getNChunks() != 0)
	throw affineExcept("ensemble_test", "main",
			   "After change in nparams, nchunks should be 0", 37);

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
