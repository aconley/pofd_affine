//Program for testing the behavior of the chain code
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>

#include<gtest/gtest.h>

#include "../include/paramSet.h"
#include "../include/affineExcept.h"
#include "../include/affineChainSet.h"
#include "../include/ran.h"

//Basic size, resize
TEST(paramSetTest, Sizing) {
  const unsigned int nparams = 3;
  paramSet p(nparams);
  ASSERT_EQ(nparams, p.getNParams()) << "paramSet not right size";

  unsigned int newsize = 5;
  p.setNParams(newsize);
  ASSERT_EQ(newsize, p.getNParams()) << "paramSet resizing failed";
}

//Test paramSet setup and recovery
TEST(paramSetTest, SetRecover) {
  const unsigned int nparams = 3;
  paramSet p(nparams);

  p[0] = 1.0; p[1] = 2.0; p[2] = -4.0;
  //[] check
  ASSERT_FLOAT_EQ(1.0, p[0]) << "First element not recovered";
  ASSERT_FLOAT_EQ(2.0, p[1]) << "Second element not recovered";
  ASSERT_FLOAT_EQ(-4.0, p[2]) << "Third element not recovered";

  //at check
  EXPECT_FLOAT_EQ(1.0, p.at(0)) << "First element not recovered in at";
  EXPECT_FLOAT_EQ(2.0, p.at(1)) << "Second element not recovered in at";
  EXPECT_FLOAT_EQ(-4.0, p.at(2)) << "Third element not recovered";
  EXPECT_THROW(p.at(3), std::range_error) << "Out of range didn't throw";

  //Set with wrong length
  std::vector<float> v1(2);
  EXPECT_THROW(p.setParamValues(v1), affineExcept) 
    << "Set with wrong vector size didn't throw";
  
  //Set with right length, recover
  std::vector<float> v2(3);
  v2[0] = -1.5; v2[1] = -3.0; v2[2] = 10.5;
  ASSERT_NO_THROW(p.setParamValues(v2)) << "Setting param values failed";
  for (unsigned int i = 0; i < 3; ++i)
    ASSERT_FLOAT_EQ(v2[i], p[i]) << "Vector set didn't set right values";

}

//Test operator= and copy constructor
TEST(paramSetTest, Equality) {
  const unsigned int nparams = 3;
  paramSet p(nparams);
  p[0] = 1.0; p[1] = 2.0; p[2] = -4.0;
  ASSERT_EQ(p.getNParams(), nparams) << "Initial paramSet not expected size";

  paramSet p2(p);
  EXPECT_EQ(nparams, p2.getNParams()) << "Copy constructor size not expected";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p2[i], p[i]) << "Copy constructor didn't set right values";

  paramSet p3;
  p3 = p;
  EXPECT_EQ(nparams, p3.getNParams()) << "Operator= size not expected";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p3[i], p[i]) << "Operator= didn't set right values";
 
}

////////////////////////////////

// Basic initialization
TEST(affineChainSetTest, BasicInit) {
  const unsigned int nwalkers = 2;
  const unsigned int nparams  = 3;

  affineChainSet chains(nwalkers,nparams);

  // We should crash out if the basic stuff doesn't work -- not much
  // point in doing more complicated tests
  ASSERT_EQ(nwalkers, chains.getNWalkers()) << "NWalkers doesn't match input";
  ASSERT_EQ(0U, chains.getNIters()) << "NIters should be zero on init";
  ASSERT_EQ(nparams, chains.getNParams()) << "NParams doesn't match input";
  ASSERT_EQ(0U, chains.getNChunks()) << "NChunks should be zero on init";

}

// Ability to insert and get back data
TEST(affineChainSetTest, InsertRecoverClear) {
  const unsigned int nwalkers = 2;
  const unsigned int nparams  = 3;
  const unsigned int niter    = 20;

  affineChainSet chains(nwalkers,nparams);

  // We should crash out if the basic stuff doesn't work -- not much
  // point in doing more complicated tests
  chains.addChunk(niter);
  ASSERT_EQ(1U, chains.getNChunks()) << "NChunks should be one after addChunk";
  
  //Add first step
  paramSet p(nparams);
  p[0] = 1.0; p[1] = 2.0; p[2] = 4.0;
  double logLike = 4.0;
  ASSERT_TRUE(chains.addNewStep(0,p,logLike)) << "Adding step failed";
  ASSERT_EQ(1U, chains.getNIters()) << "NIters should be one after first add";

  //Recovery test
  paramSet p2(nparams);
  double logLike2;
  chains.getLastStep(0, p2, logLike2);
  EXPECT_FLOAT_EQ(logLike, logLike2) <<
    "Recovered step like doesn't match input";
  for (unsigned int i = 0; i < nparams; ++i)
    EXPECT_FLOAT_EQ(p[i], p2[i]) << "Recovered step params don't match input";
  
  //Make sure we can clear
  chains.clear();
  EXPECT_EQ(nwalkers, chains.getNWalkers()) 
    << "Unexpected NWalkers after clear";
  EXPECT_EQ(0U, chains.getNIters()) <<
    "NIters should be zero after clear";
  EXPECT_EQ(nparams, chains.getNParams()) <<
    "NParams doesn't match after clear";
  EXPECT_EQ(0U, chains.getNChunks()) <<
    "NChunks should be zero after clear";

}

// Randomly adding a bunch of data
TEST(affineChainSetTest, RandomInsert) {
  const unsigned int nwalkers = 2;
  const unsigned int nparams  = 3;
  const unsigned int niter    = 20;

  ran ranstruct;
  affineChainSet chains(nwalkers,nparams);
  paramSet p(nparams);
  double logLike, logLike2;

  chains.addChunk(niter);
  for (unsigned int i = 0; i < niter; ++i)
    for (unsigned int j = 0; j < nwalkers; ++j) {
      p[0] = ranstruct.doub();
      p[1] = ranstruct.doub()*2.4;
      p[2] = ranstruct.doub()*4;
      logLike = ranstruct.doub()*100;
      ASSERT_TRUE(chains.addNewStep(j,p,logLike)) 
	<< "Failure to add random step";
    }

  EXPECT_EQ(1U, chains.getNChunks())
    << "NChunks should be one after random insert";
  EXPECT_EQ(niter * nwalkers, chains.getNIters())
    << "niters should be nwalker*niters after random insert";


  //Add a second chunk with fewer elements!
  chains.addChunk(5);
  ASSERT_EQ(2U, chains.getNChunks()) <<
    "NChunks should be two after second addChunk";

  //Max check
  //Add a high log like element to walker 1
  paramSet p2(nparams);
  p2[0] = ranstruct.doub();
  p2[1] = ranstruct.doub()*2.4;
  p2[2] = ranstruct.doub()*4;
  logLike2 = 1000; //Should be larger than anything previous
  ASSERT_TRUE(chains.addNewStep(1,p2,logLike2)) << 
    "Failure to add new step to second chunk, walker 1";
  //And a random one to walker 0
  p[0] = ranstruct.doub();
  p[1] = ranstruct.doub()*2.4;
  p[2] = ranstruct.doub()*4;
  logLike = ranstruct.doub()*100;
  ASSERT_TRUE(chains.addNewStep(0,p,logLike)) << 
    "Failure to add new step to second chunk, walker 0";

  //Fill out 2 more
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < nwalkers; ++j) {
      p[0] = ranstruct.doub();
      p[1] = ranstruct.doub()*2.4;
      p[2] = ranstruct.doub()*4;
      logLike = ranstruct.doub()*100;
      ASSERT_TRUE(chains.addNewStep(j,p,logLike)) 
	<< "Failure to add random step, second chunk";
    }
  EXPECT_EQ(niter*nwalkers + nwalkers*3, chains.getNIters()) <<
    "NIters should be nwalker*niters+nwalkers*3 after last insert";

  //Maximum check
  logLike = -2;
  logLike = chains.getMaxLogLike();
  EXPECT_FLOAT_EQ(logLike, logLike2) << "Didn't recover max logLike";
  //Try the same but also getting the parameters (still in p2)
  logLike = -2;
  chains.getMaxLogLikeParam(logLike,p);
  EXPECT_FLOAT_EQ(logLike, logLike2) <<
    "Didn't recover max logLike w/params";
  for (unsigned int i = 0; i < nparams; ++i)
    EXPECT_FLOAT_EQ(p2[i], p[i]) << 
      "Recovered max logLike params don't match input";
}

// Chain appending
TEST(affineChainSetTest, Append) {
  const unsigned int nwalkers = 3;
  const unsigned int nparams  = 2;
  
  affineChainSet chains(nwalkers, nparams);
  chains.addChunk(4);
  ASSERT_EQ(1U, chains.getNChunks()) << "Chain should have 1 chunk";

  //Add 3 items 
  paramSet p(nparams);
  double logLike, logLike2;
  ran ranstruct;
  for (unsigned int i = 0; i < 3U; ++i)
    for (unsigned int j = 0; j < nwalkers; ++j) {
      p[0] = ranstruct.doub();
      p[1] = ranstruct.doub()*2.4;
      logLike = ranstruct.doub()*100;
      ASSERT_TRUE(chains.addNewStep(j,p,logLike)) << "Failed to add element";
    }
  ASSERT_EQ(chains.getNIters(), 3U * nwalkers) <<
    "Chains had unexpected number of elements";
  //Add one more with a large loglike to last walker
  for (unsigned int j = 0; j < nwalkers-1; ++j) {
    p[0] = ranstruct.doub();
    p[1] = ranstruct.doub()*2.4;
    logLike = ranstruct.doub()*10;
    ASSERT_TRUE(chains.addNewStep(j,p,logLike)) << "Failed to add element";
  }
  p[0] = ranstruct.doub();
  p[1] = ranstruct.doub()*2.4;
  logLike = 2000;
  ASSERT_TRUE(chains.addNewStep(nwalkers-1,p,logLike)) << 
    "Failed to add element";

  ASSERT_EQ(chains.getNIters(), 4U * nwalkers) <<
    "Chains had unexpected number of elements after final add";
  logLike2 = chains.getMaxLogLike();
  //Note that logLike is still set to max value
  ASSERT_FLOAT_EQ(logLike, logLike2) <<
    "Recovered max log like not what was expected";

  // Second chain set
  affineChainSet chains2(nwalkers, nparams);
  chains2.addChunk(3);
  ASSERT_EQ(1U, chains.getNChunks()) << "Chain2 should have 1 chunk";

  // And 2 to the second
  for (unsigned int i = 0; i < 2U; ++i)
    for (unsigned int j = 0; j < nwalkers; ++j) {
      p[0] = ranstruct.doub();
      p[1] = ranstruct.doub()*1.5;
      logLike = ranstruct.doub()*100;
      ASSERT_TRUE(chains2.addNewStep(j,p,logLike)) << "Failed to add element";
    }
   ASSERT_EQ(chains2.getNIters(), 2U * nwalkers) <<
    "Chains2 had unexpected number of elements";

   //Add a second chunk, fill
   chains2.addChunk(1);
   ASSERT_EQ(2U, chains2.getNChunks()) << "Chain2 should have 2 chunks";
   for (unsigned int j = 0; j < nwalkers; ++j) {
     p[0] = ranstruct.doub();
     p[1] = ranstruct.doub()*1.5;
     logLike = ranstruct.doub()*100;
     ASSERT_TRUE(chains2.addNewStep(j,p,logLike)) << "Failed to add element";
   }
   ASSERT_EQ(chains2.getNIters(), 3U * nwalkers) <<
    "Chains2 had unexpected number of elements";

   //Append and test
   chains += chains2;
   EXPECT_EQ(chains.getNChunks(), 3U) <<
     "Appended chain should have three chunks";
   EXPECT_EQ(chains.getNIters(), 7U * nwalkers) <<
     "Appended chain had wrong number of entries";

   //Make sure max log like still as expected
   // The old max is still in logLike2
   logLike = chains.getMaxLogLike();
   EXPECT_FLOAT_EQ(logLike, logLike2) <<
     "Recovered max log like not what was expected after append";

   //Make sure this clears
   chains2.clear();
   EXPECT_EQ(0U, chains2.getNChunks()) <<
     "Chains2 not empty after clear";
   
}

// Averaging
TEST(affineChainSetTest, Averaging) {
  const unsigned int nwalkers = 2;
  const unsigned int nparams  = 3;

  affineChainSet chains(nwalkers,nparams);
  ASSERT_EQ(nwalkers, chains.getNWalkers()) << "NWalkers doesn't match input";
  ASSERT_EQ(0U, chains.getNIters()) << "NIters should be zero on init";
  ASSERT_EQ(nparams, chains.getNParams()) << "NParams doesn't match input";
  ASSERT_EQ(0U, chains.getNChunks()) << "NChunks should be zero on init";

  chains.addChunk(2);
  ASSERT_EQ(1U, chains.getNChunks()) << "NChunks should be 1";

  paramSet p(nparams);
  double logLike = 10.0;
  ASSERT_EQ(nparams, p.getNParams()) << 
    "paramSet.getNParams should be " << nparams;
  p[0] = 1.0; p[1] = 2.0; p[2] = 4.0;
  ASSERT_TRUE(chains.addNewStep(0,p,logLike)) << "Failed to add step";
  p[0] = 0.5; p[1] = 1.0; p[2] = 2.0;
  ASSERT_TRUE(chains.addNewStep(1,p,logLike)) << "Failed to add step";
  p[0] = 3.0; p[1] = 1.0; p[2] = -2.0;
  ASSERT_TRUE(chains.addNewStep(0,p,logLike)) << "Failed to add step";
  p[0] = 1.5; p[1] = 1.7; p[2] = 2.0;
  ASSERT_TRUE(chains.addNewStep(1,p,logLike)) << "Failed to add step";
  float p1av0 = 0.5*(2.0+1.0);  //First step
  float p1av1 = 0.5*(1.0+1.7);  //Second step

  //First, check recovery of specific entry
  // This should get the values of parameter 2 from walker 1.  We
  // have entered two values, so we should get two back
  std::vector<float> parvec;
  chains.getParamVector(1, 2, parvec);
  ASSERT_EQ(parvec.size(), 2U) <<
    "Vector from chains::getParamVector not of expected length (2)";
  //This is the one we just stuck in using p
  EXPECT_FLOAT_EQ(parvec[1], p[2]) <<
    "chains::getParamVector didn't give back right value";

  //Now, check averages
  parvec.resize(3); //Resize to check getAverage does resizing as well
  // Parameter 1
  chains.getAverageParamVector(1,parvec);
  ASSERT_EQ(parvec.size(), 2U) <<
    "chains::getAverageParamVector didn't resize vector correctly";
  EXPECT_FLOAT_EQ(parvec[0], p1av0) <<
    "getAverageParamVector didn't return right elem0";
  EXPECT_FLOAT_EQ(parvec[1], p1av1) <<
    "getAverageParamVector didn't return right elem1";
}

//Getting last point
TEST(affineChainSetTest, GetLastPoint) {
  const unsigned int nwalkers = 2;
  const unsigned int nparams  = 3;
  affineChainSet chains(nwalkers, nparams);
  paramSet p(nparams), p2(nparams);
  double logLike, logLike2;

  logLike = 3012.5;
  chains.addChunk(2);
  p[0] = 0.0; p[1] = 0.0; p[2] = 0.0;
  ASSERT_TRUE(chains.addNewStep(0,p,logLike)) << "Failed to add point";
  p[0] = 0.5; p[1] = 1.0; p[2] = 2.0;
  ASSERT_TRUE(chains.addNewStep(0,p,logLike)) << "Failed to add second point";
  p2[0] = 1.5; p2[1] = 2.0; p2[2] = -2.0;
  ASSERT_TRUE(chains.addNewStep(1,p2,logLike)) << "Failed to add third point";

  ASSERT_EQ(chains.getNIters(), 3U) <<
    "chain should have 3 steps after insert";
  
  //Actual recovery test
  paramSet pold(nparams), pnew(nparams);
  chains.getLastStep(0,pold,logLike2); //Last step on walker 0
  EXPECT_FLOAT_EQ(logLike2, logLike) << 
    "Didn't recover last step from walker 0";
  for (unsigned int i = 0; i < nparams; ++i)
    EXPECT_FLOAT_EQ(pold[i], p[i]) <<
      "Didn't recover last step from walker 0";
  chains.getLastStep(1,pold,logLike2); //Walker 1
  EXPECT_FLOAT_EQ(logLike2, logLike) << 
    "Didn't recover last step from walker 1";
  for (unsigned int i = 0; i < nparams; ++i)
    EXPECT_FLOAT_EQ(pold[i], p2[i]) <<
      "Didn't recover last step from walker 0";
}

//Accessing older chunks after adding new ones
TEST(affineChainSetTest, AccessOldChunks) {
  const unsigned int nwalkers = 1;
  const unsigned int nparams  = 3;

  affineChainSet chains(nwalkers, nparams);
  paramSet pnew(nparams), pold(nparams);
  double logLike, logLike2;
  ran ranstruct;

  logLike = 1000.0;
  pold[0] = 1.0; pold[1] = 2.0; pold[2] = -3.5;
  chains.addChunk(1);
  ASSERT_TRUE(chains.addNewStep(0, pold, logLike)) << 
    "Failed to add first step";
   
  //Now add 10 new chunks and a bunch of steps
  for (unsigned int idx = 0; idx < 10; ++idx) {
    chains.addChunk(1000);
    for (unsigned int i = 0; i < 1000; ++i) {
      pnew[0] = 100*ranstruct.doub();
      pnew[1] = 50*ranstruct.doub();
      pnew[2] = 20*ranstruct.doub();
      ASSERT_TRUE(chains.addNewStep(0, pnew, 2.0)) << "Failed to add new step";
    }
  }
  
  //First -- do we have the right number of chunks!
  ASSERT_EQ(chains.getNChunks(), 11U) << "Wrong number of chunks";
  //And total number of steps
  ASSERT_EQ(chains.getNIters(), 10U * 1000U + 1U) <<
    "Wrong number of steps";

  //And try to access the original one (pold) into pnew
  chains.getStep(0, 0, 0, pnew,logLike2);
  EXPECT_FLOAT_EQ(logLike, logLike2) <<
    "Got back wrong logLike from 0, 0, 0 step";
  for (unsigned int i = 0; i < nparams; ++i)
    EXPECT_FLOAT_EQ(pold[i], pnew[i]) <<
      "Got back wrong 0 0 0 step value for param " << i;
}

//Make sure clearing while preserving the last step works
TEST(affineChainSetTest, ClearPreserveLast) {
    const unsigned int nwalkers = 5;
    const unsigned int nparams  = 3;
    paramSet pnew(nparams), pold(nparams);
    double logLike, logLike2;
    ran ranstruct;

    affineChainSet chains(nwalkers, nparams);
    chains.addChunk(1000);
    for (unsigned int i = 0; i < 1000; ++i) 
      for (unsigned int j = 0; j < nwalkers; ++j) {
	pnew[0] = 100*ranstruct.doub();
	pnew[1] = 50*ranstruct.doub();
	pnew[2] = 20*ranstruct.doub();
        ASSERT_TRUE(chains.addNewStep(j,pnew,2.0)) << 
	  "Failure to add random step";
      }
    //And, do it again, but with a smaller chunk we don't totally fill
    chains.addChunk(100);
    for (unsigned int i = 0; i < 30; ++i) 
      for (unsigned int j = 0; j < nwalkers; ++j) {
	pnew[0] = 100*ranstruct.doub();
	pnew[1] = 50*ranstruct.doub();
	pnew[2] = 20*ranstruct.doub();
        ASSERT_TRUE(chains.addNewStep(j,pnew,2.0)) << 
	  "Failure to add random step";
      }
    ASSERT_EQ(chains.getNIters(), nwalkers * (1000U + 30U) ) <<
      "Got unexpected number of total steps";

    //Now, add one more to all walkers
    pold[0] = 1.0; pold[1] = 2.0; pold[2] = 3.0;
    logLike = -11.7;
    for (unsigned int j = 0; j < nwalkers; ++j) 
      ASSERT_TRUE(chains.addNewStep(j,pold,logLike)) <<
	"Failure to add last step";

    chains.clearPreserveLast();
    ASSERT_EQ(chains.getNChunks(), 1U) <<
      "nchunks should be one after clearPreserveLast";
    ASSERT_EQ(chains.getNIters(), nwalkers) <<
      "NIters should be nwalkers after clearPreserveLast";
    for (unsigned int j = 0; j < nwalkers; ++j) {
      chains.getLastStep(j,pnew,logLike2);
      EXPECT_FLOAT_EQ(logLike2, logLike) <<
	"didn't recover logLike after clearPreserveLast in walker: " << j;
      for (unsigned int i = 0; i < nparams; ++i)
	EXPECT_FLOAT_EQ(pnew[i], pold[i]) <<
	  "didn't recover param " << i << 
	  "after clearPreserveLast in walker: " << j;
    }
}

//Test Acor using test vector
TEST(affineChainSetTest, Acor) {
  const unsigned int ntest = 8815;
  std::ifstream ifs("testdata/testvec.txt");
  ASSERT_TRUE(ifs) << "Can't open testvec.txt";
  affineChainSet chains3(1,1);
  chains3.addChunk(ntest);
  paramSet p(1);
  float val;
  double logLike = 1.0;
  for (unsigned int i = 0; i < ntest; ++i) {
    ifs >> val;
    p[0] = val;
    ASSERT_TRUE(chains3.addNewStep(0, p, logLike)) << 
      "Failed to add step " << i;
  }
  ifs.close();
  double tau, mean, sigma;
  bool succ;
  tau = chains3.getAcor(0,mean,sigma,succ);
  ASSERT_TRUE(succ) << "Couldn't compute acor";

  //Values computed using original acor code
  EXPECT_NEAR(tau, 130.052, 0.002) << "Incorrect tau: " << tau;
  EXPECT_NEAR(mean, 5.82089, 0.0002) << "Incorrect mean: " << mean;
  EXPECT_NEAR(sigma, 0.0134451, 0.0002) << "Incorrect sigma: " << sigma;
}

////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running chain_tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

