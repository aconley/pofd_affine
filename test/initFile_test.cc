#include<iostream>
#include<cmath>

#include<gtest/gtest.h>

#include "../include/numberCountsKnots.h"
#include "../include/numberCountsDoubleLogNormal.h"
#include "../include/affineExcept.h"
#include "../include/paramSet.h"

//This is a pretty simple test.  There are two files in the test
// directory.  We read them and check that the resulting object
// is what we expect
TEST(initFile1DTest, ReadCheck1) {
  //File 1
  //Simple file, no limits, etc.
  initFileKnots test1("testdata/fiducial_1Dmodel.txt",false,false);
  
  unsigned int nknots = 10;
  ASSERT_EQ(test1.getNKnots(), nknots) <<
    "test file1 had wrong number of knots";

  //Test knot positions
  const double knotpos[] = { 0.0001, 0.002, 0.005, 0.010, 0.02, 0.045,
			     0.100, 0.2, 0.45, 1.0 };
  std::vector<double> kp;
  test1.getKnotPos(kp);
  for (unsigned int i = 0; i < nknots; ++i) {
    EXPECT_FLOAT_EQ(kp[i], knotpos[i]) <<
      "Knot position problem in test file; expected " << 
      knotpos[i] << " got " << kp[i] << " on knot " << i;
  }

  //Test knot values
  const double knotval[] = { 10.75, 7.0, 6.3, 5.9, 5.1, 4.0, 2.6, 1.4, 
			     0.6, -0.4 };
  std::vector<double> kv;
  test1.getKnotVals(kv);
  for (unsigned int i = 0; i < nknots; ++i) {
    EXPECT_FLOAT_EQ(kp[i], knotpos[i]) <<
      "Knot value problem in test file; expected " << 
      knotval[i] << " got " << kv[i] << " on knot " << i;
  }


  //No knots should be fixed, since we didn't read sigma
  for (unsigned int i = 0; i < nknots; ++i) 
    EXPECT_FALSE(test1.isKnotFixed(i)) <<
      "No knots should be fixed in test file, but " <<
      "knot " << i << " is";

  //Test getting parameters
  paramSet p(10);
  test1.getParams(p);
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_FLOAT_EQ(p[i], knotval[i]) <<
      "Parameter get in test file; expected " << knotval[i] << " got " << 
      p[i] << " on knot " << i;
}

TEST(initFile1DTest, ReadCheck2) {
  //Test file 2 -- a more complicated case with sigmas and limits
  const double knotpos[] = { 0.002, 0.005, 0.045, 0.100, 0.45, 1.0 };
  const double knotval[] = { 7.0, 6.3, 4.0, 2.6, 0.6, -0.4 };
  const bool knotfixed[] = { false, false, false, false, false, true };
  const double knotsigma[] = { 0.02, 0.03, 0.02, 0.01, 0.05, 0.0 };
  const bool haslowlim[] = { true, false, true, false, true, false };
  const double lowlim[] = { 3.0, 0, 1.0, 0, -2.0, 0 };
  const bool hasuplim[] = { true, false, false, false, true, false };
  const double uplim[] = { 8.0, 0.0, 0.0, 0.0, 3.0, 0.0 };
  
  initFileKnots test2("testdata/initfile_test1D.txt",true,true);
    
  //First make sure it has the right number of entries
  unsigned int nknots = test2.getNKnots();
  ASSERT_EQ(nknots, 6U) << "test file had wrong number of knots;"
			<< " expected " << 6 << " got " << nknots;

  //Knot position test
  for (unsigned int i = 0; i < nknots; ++i) {
    double kp = test2.getKnotPos(i);
    EXPECT_FLOAT_EQ(kp, knotpos[i]) <<
      "Knot position problem in test file; expected " << 
      knotpos[i] << " got " << kp << " on knot " << i;

    //Knot value test
    double kv = test2.getKnotValue(i);
    EXPECT_FLOAT_EQ(kv, knotval[i]) <<
      "Knot value problem in test file; expected " <<
      knotval[i] << " got " << kv << " on knot " << i;

    //Knot fixed
    EXPECT_EQ(test2.isKnotFixed(i), knotfixed[i]) <<
      "Knot fixed problem; expected " << knotfixed[i] <<
      " got " << test2.isKnotFixed(i) << " on knot " << i;

    //Knot sigma
    double ks = test2.getKnotSigma(i);
    EXPECT_FLOAT_EQ(ks, knotsigma[i]) <<
      "Knot sigma problem in test file; expected " <<
      knotsigma[i] << " got " << ks  << " on knot " << i;

    //Lower limit
    bool chaslowlim = test2.knotHasLowerLimit(i);
    EXPECT_EQ(chaslowlim, haslowlim[i]) <<
      "Knot has lower limit problem; expected " << haslowlim[i] <<
      " got " << chaslowlim;

    if (chaslowlim) {
      double kl = test2.getLowerLimit(i);
      EXPECT_FLOAT_EQ(kl, lowlim[i]) <<
	"Knot lower limit problem in test file; expected " <<
	lowlim[i] << " got " << kl << " on knot " << i;
    }

    //Upper limit
    bool chasuplim = test2.knotHasUpperLimit(i);
    EXPECT_EQ(chasuplim, hasuplim[i]) <<
      "Knot has upper limit problem; expected " << hasuplim[i] <<
      " got " << chasuplim;

    if (chasuplim) {
      double ku = test2.getUpperLimit(i);
      EXPECT_FLOAT_EQ(ku, uplim[i]) <<
	"Knot upper limit problem in test file; expected " <<
	uplim[i] << " got " << ku << " on knot " << i;
    }
  }

  //Test getting parameters
  paramSet p(10);
  test2.getParams(p);
  for (unsigned int i = 0; i < nknots; ++i) 
    EXPECT_FLOAT_EQ(p[i], knotval[i]) <<
      "Parameter get in test file; expected " <<
      knotval[i] << " got " << p[i] << " on knot " << i;
  
}

TEST(initFile2DTest, ReadCheck) {
  initFileDoubleLogNormal test1("testdata/initfile_test2D.txt",true,true);

  unsigned int nknots = test1.getNKnots();
  ASSERT_EQ(nknots, 6U) << "test file had wrong number of knots";

  unsigned int nsigmas = test1.getNSigmas();
  ASSERT_EQ(nsigmas, 3U) << "test file had wrong number of sigmas";
    
  unsigned int noffsets = test1.getNOffsets();
  ASSERT_EQ(noffsets, 2U) << "test file had wrong number of offsets";

  unsigned int ntot = test1.getNTot();
  ASSERT_EQ(ntot, 11U) << "test file had wrong number of total params";
      
  //Knots tests
  const double knotpos[] = { 0.002, 0.005, 0.045, 0.100, 0.45, 1.0 };
  const double knotval[] = { 7.0, 6.3, 4.0, 2.6, 0.6, -0.4 };
  const bool knotfixed[] = { false, false, false, false, false, true };
  const double knotsigma[] = { 0.02, 0.03, 0.02, 0.01, 0.05, 0.0 };
  const bool haslowlim[] = { true, false, true, false, true, false };
  const double lowlim[] = { 3.0, 0, 1.0, 0, -2.0, 0 };
  const bool hasuplim[] = { true, false, false, false, true, false };
  const double uplim[] = { 8.0, 0.0, 0.0, 0.0, 3.0, 0.0 };
  
  std::vector<double> kp, kv;
  test1.getKnotPos(kp);
  test1.getKnotVals(kv);
  for (unsigned int i = 0; i < nknots; ++i) {
    EXPECT_FLOAT_EQ(knotpos[i], kp[i]) << 
      "Knot position problem in test file on knot " << i;

    EXPECT_FLOAT_EQ(knotval[i], kv[i]) << 
      "Knot value problem in test file on knot " << i;

    //Knot fixed
    EXPECT_EQ(test1.isKnotFixed(i), knotfixed[i]) << 
      "Knot fixed problem on knot " << i;

    //Knot sigma
    EXPECT_FLOAT_EQ(test1.getKnotSigma(i), knotsigma[i]) <<
      "Knot sigma problem in test file on knot " << i;
    
    //Lower limit
    EXPECT_EQ(test1.knotHasLowerLimit(i), haslowlim[i]) <<
      "Knot has lower limit problem on knot " << i;

    if (haslowlim[i]) {
      EXPECT_FLOAT_EQ(test1.getLowerLimit(i), lowlim[i]) <<
	"Knot lower limit problem in test file on knot " << i;
    }

    //Upper limit
    EXPECT_EQ(test1.knotHasUpperLimit(i), hasuplim[i]) <<
      "Knot has upper limit problem on knot " << i;

    if (hasuplim[i]) 
      EXPECT_FLOAT_EQ(test1.getUpperLimit(i), uplim[i]) <<
	"Knot upper limit problem in test file on knot " << i;
  
  }

  //Same thing, but sigmas
  const double sigmapos[] = { 0.004, 0.03, 0.10 };
  const double sigmaval[] = { 0.3, 0.35, 0.45 };
  const bool sigmafixed[] = { false, false, true };
  const double sigmasigma[] = { 0.03, 0.03, 0.0 };
  const bool haslowlimsigma[] = { false, true, false };
  const double lowlimsigma[] = { 0.0, 0.1, 0.0 };
  const bool hasuplimsigma[] = { false, true, false };
  const double uplimsigma[] = { 0.0, 0.6, 0.0 };
  test1.getSigmaPos(kp);
  test1.getSigmaVals(kv);
  for (unsigned int i = 0; i < nsigmas; ++i) {
    //Knot position test
    EXPECT_FLOAT_EQ(kp[i], sigmapos[i]) << 
      "Sigma knot position problem in test file; expected " <<
      sigmapos[i] << " got " << kp[i] << " on knot " << i;
    
    //Knot value test
    EXPECT_FLOAT_EQ(kv[i], sigmaval[i]) << 
      "Sigma knot value problem in test file; expected " <<
      sigmaval[i] << " got " << kv[i] << " on knot " << i;
    
    //Sigma of sigma knots
    EXPECT_EQ(test1.isKnotFixed(i+nknots), sigmafixed[i]) <<
      "Sigma knot fixed problem; expected " << sigmafixed[i] <<
      " got " << test1.isKnotFixed(i+nknots) << std::endl <<
      " on sigma knot " << i;
    
    EXPECT_FLOAT_EQ(test1.getKnotSigma(i+nknots), sigmasigma[i]) <<
      "Sigma knot sigma problem in test file; expected " << 
      sigmasigma[i] << " got " << test1.getKnotSigma(i+nknots) << 
      std::endl << " on sigma knot " << i;
    
    //Lower limit
    EXPECT_EQ(test1.knotHasLowerLimit(i+nknots), haslowlimsigma[i]) <<
      "Sigma knot has lower limit problem; expected " << haslowlimsigma[i] <<
      " got " << test1.knotHasLowerLimit(i+nknots) << std::endl << 
      " on sigma knot " << i;
    if (haslowlimsigma[i]) 
      EXPECT_FLOAT_EQ(test1.getLowerLimit(i+nknots), lowlimsigma[i]) <<
	"Sigma knot lower limit problem in test file; expected " <<
	lowlimsigma[i] << " got " << test1.getLowerLimit(i+nknots) << 
	" on sigma knot " << i;
    
    //Upper limit
    EXPECT_EQ(test1.knotHasUpperLimit(i+nknots), hasuplimsigma[i]) <<
      "Sigma Knot has upper limit problem; expected " << hasuplimsigma[i] <<
      " got " << test1.knotHasUpperLimit(i+nknots) << " on sigma knot " << i;
    if (hasuplimsigma[i]) 
      EXPECT_FLOAT_EQ(test1.getUpperLimit(i+nknots), uplimsigma[i]) <<
	"Sigma knot upper limit problem in test file; expected " << 
	uplimsigma[i] << " got " << test1.getUpperLimit(i+nknots) << 
	" on sigma knot " << i;
  }
}

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running initFile tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
