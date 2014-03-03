//Testing models
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<utility>

#include<gtest/gtest.h>

#include "../include/global_settings.h"
#include "../include/paramSet.h"
#include "../include/affineExcept.h"
#include "../include/beam.h"
#include "../include/doublebeam.h"
#include "../include/numberCountsKnotsSpline.h"
#include "../include/numberCountsDoubleLogNormal.h"



//////////////////////////////////////////////////
// numberCountsKnotsSpline

//Basic Instantiation
TEST(model1DTest, Init) {
  
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  std::vector<double> kposvec(nknots);
  for (unsigned int i = 0; i < nknots; ++i) kposvec[i] = knotpos[i];

  numberCountsKnotsSpline model1(nknots);
  EXPECT_EQ(nknots, model1.getNKnots()) << 
    "Unexpected number of knots in numberCountsKnotsSpline(int) constructor";
  EXPECT_FALSE(model1.isValid()) <<
    "Model should not be valid after numberCountsKnotsSpline(int) constructor";

  numberCountsKnotsSpline model2(kposvec);
  EXPECT_EQ(nknots, model2.getNKnots()) << 
    "Unexpected number of knots in numberCountsKnotsSpline(vec) constructor";
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_FLOAT_EQ(kposvec[i], model2.getKnotPos(i)) <<
      "Didn't get expected knot position after numberCountsKnotsSpline(vec)" <<
      " at position: " << i;
  EXPECT_FLOAT_EQ(kposvec[0], model2.getMinFlux()) <<
    "Didn't get expected maximum model knot after vector constructor";
  EXPECT_FLOAT_EQ(kposvec[nknots-1], model2.getMaxFlux()) <<
    "Didn't get expected maximum model knot after vector constructor";
  EXPECT_FALSE(model2.isValid()) <<
    "Model should not be valid after numberCountsKnotsSpline(vec) constructor";
  
  numberCountsKnotsSpline model3(nknots, knotpos);
  EXPECT_EQ(nknots, model3.getNKnots()) << 
    "Unexpected number of knots in numberCountsKnotsSpline(int, double*)" <<
    " constructor";
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_FLOAT_EQ(kposvec[i], model3.getKnotPos(i)) <<
      "Didn't get expected knot position after " <<
      "numberCountsKnotsSpline(int, double*)" << " at position: " << i;
  EXPECT_FLOAT_EQ(kposvec[0], model3.getMinFlux()) <<
    "Didn't get expected maximum model knot after c array constructor";
  EXPECT_FLOAT_EQ(kposvec[nknots-1], model3.getMaxFlux()) <<
    "Didn't get expected maximum model knot after c array constructor";
  EXPECT_FALSE(model3.isValid()) <<
    "Model should not be valid after c array constructor";

  //Copy constructor
  numberCountsKnotsSpline model4(model2);
   EXPECT_EQ(nknots, model4.getNKnots()) << 
    "Unexpected number of knots after copy constructor";
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_FLOAT_EQ(kposvec[i], model4.getKnotPos(i)) <<
      "Didn't get expected knot position after " <<
      "copy constructor" << " at position: " << i;
  EXPECT_FLOAT_EQ(kposvec[0], model4.getMinFlux()) <<
    "Didn't get expected maximum model knot after copy constructor";
  EXPECT_FLOAT_EQ(kposvec[nknots-1], model4.getMaxFlux()) <<
    "Didn't get expected maximum model knot after copy constructor";
  EXPECT_FALSE(model4.isValid()) <<
    "Model should not be valid after copy constructor";

}

//Test set positions
TEST(model1DTest, SetPositions) {
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };

  //Make it with fewer knots -- set should resize
  numberCountsKnotsSpline model1(nknots-1);
  ASSERT_EQ(nknots-1, model1.getNKnots()) << 
    "Unexpected number of knots in numberCountsKnotsSpline(int) constructor";
  model1.setKnotPositions(nknots, knotpos);
  ASSERT_EQ(nknots, model1.getNKnots()) <<
    "setKnotPositions failed to resize";
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_FLOAT_EQ(knotpos[i], model1.getKnotPos(i)) <<
      "After setKnotPositions, didn't find expected knot position at " << i;

  
  //Reject negative knot positions
  const double knotpos2[nknots] = { -0.002, 0.005, 0.010, 0.020, 0.040 };
  EXPECT_THROW(model1.setKnotPositions(nknots, knotpos2), affineExcept) <<
    "Trying to set a negative knot position should throw an exception";
}

//Test set params
TEST(model1DTest, SetParams) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const float knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  ASSERT_EQ(nknots, model.getNKnots()) << "Model has wrong number of knots";
  EXPECT_FALSE(model.isValid()) << "Model should not be valid";
  paramSet p(nknots, knotval);
  ASSERT_EQ(nknots, p.getNParams()) << "paramSet has wrong number of values";

  model.setParams(p);
  ASSERT_EQ(nknots, model.getNKnots()) <<
    "Model has wrong number of knots after setParams";
  EXPECT_TRUE(model.isValid()) << "Model should be valid after setParams";
  std::pair<double, double> pr;
  for (unsigned int i = 0; i < nknots; ++i) {
    pr = model.getLogKnot(i);
    EXPECT_FLOAT_EQ(knotpos[i], pr.first) << 
      "Recovered wrong knot position after setParams";
    EXPECT_FLOAT_EQ(pofd_mcmc::logfac * knotval[i], pr.second) << 
      "Recovered wrong knot log2 value after setParams";
  }

  //Try setting the wrong length
  paramSet p2(nknots-1);
  ASSERT_EQ(nknots-1, p2.getNParams()) <<
    "paramSet has wrong number of knots";
  EXPECT_THROW(model.setParams(p2), affineExcept) <<
    "Should throw exception if number of elements in paramSet is wrong";
}

//Test operator=
TEST(model1DTest, Equals) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const float knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model1(nknots, knotpos);
  ASSERT_EQ(nknots, model1.getNKnots()) << "Model has wrong number of knots";
  paramSet p(nknots, knotval);
  ASSERT_EQ(nknots, p.getNParams()) << "paramSet has wrong number of values";
  model1.setParams(p);
  ASSERT_EQ(nknots, model1.getNKnots()) <<
    "Model has wrong number of knots after setParams";

  numberCountsKnotsSpline model2;
  model2 = model1;
  ASSERT_EQ(nknots, model2.getNKnots()) << 
    "Model2 has wrong number of knots after =";
  std::pair<double, double> pr;
  for (unsigned int i = 0; i < nknots; ++i) {
    pr = model2.getLogKnot(i);
    EXPECT_FLOAT_EQ(knotpos[i], pr.first) << 
      "Recovered wrong knot position after operator=";
    EXPECT_FLOAT_EQ(pofd_mcmc::logfac * knotval[i], pr.second) << 
      "Recovered wrong knot log2 value after setParams";
  }
}

//Test getting number counts; this depends on the details of
// the spline model, so will change if GSL changes anything
// Hence, it's a little hard to say exactly what the results should be.
TEST(model1DTest, NumberCounts) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const float knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  paramSet p(nknots, knotval);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should consider itself valid";

  //Evaluate at knot positions -- these should be very close
  for (unsigned int i = 1; i < nknots-1; ++i)
    EXPECT_NEAR(knotval[i], log10(model.getNumberCounts(knotpos[i])), 1e-4) << 
      "Didn't get expected number counts at " << knotpos[i];
  
  const unsigned int ntest = 3;
  const double testpos[ntest] = {0.003, 0.007, 0.015};
  const double testval[ntest] = {8.69422, 5.88217, 4.53303};
  //Don't require as tight tolerances for this test
  for (unsigned int i = 0; i < ntest; ++i)
    EXPECT_NEAR(testval[i], log10(model.getNumberCounts(testpos[i])), 1e-3) << 
      "Didn't get expected number counts at " << testpos[i];

}

//Model integration -- also hard to ensure the values are right
TEST(model1DTest, Integration) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const float knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  paramSet p(nknots, knotval);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  EXPECT_NEAR(3.12022e6, model.getNS(), 1e4) <<
    "Unexpected total number of sources for test case";
  EXPECT_NEAR(7386.67, model.getFluxPerArea(), 1.0) <<
    "Unexpected flux per area";
  EXPECT_NEAR(18.132, model.getFluxSqPerArea(), 0.01) <<
    "Unexpected flux^2 per area";
}


//Test getting R range
TEST(model1DTest, getRRange) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const float knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  paramSet p(nknots, knotval);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should consider itself valid";
  beam bm("testdata/band1_beam.fits");
  
  ASSERT_NEAR(knotpos[0], model.getMinFlux(), 1e-5)
    << "Unexpected model minimum flux";
  ASSERT_NEAR(knotpos[nknots-1], model.getMaxFlux(), 1e-5)
    << "Unexpected model minimum flux";

  dblpair rng = model.getRRange(bm);
  EXPECT_FLOAT_EQ(0, rng.first) << "Model lower R range should be zero";
  EXPECT_NEAR(0.0488993, rng.second, 1e-4) << "Unexpected model upper R range";
}


//Test getting R; this depends on the details of
// the spline model, so will change if GSL changes anything
// Hence, it's a little hard to say exactly what the results should be.
TEST(model1DTest, getR) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const float knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  paramSet p(nknots, knotval);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should consider itself valid";

  beam bm("testdata/band1_beam.fits");
  const unsigned int ntest = 3;
  const double fluxdens[ntest] = {0.003, 0.011, 0.03};
  const double rexpect[ntest] = {1891.52, 0.817009, 0.000313631};
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens[i], bm);
    reldiff = fabs(rval - rexpect[i]) / rexpect[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R value not as expected -- wanted: " << rexpect[i] <<
      " got " << rval << " relative difference: " << reldiff;
  }

  EXPECT_FLOAT_EQ(0.0, model.getR(-0.001, bm))
    << "R should be zero for negative flux density";
  EXPECT_FLOAT_EQ(0.0, model.getR(0.1, bm))
    << "R should be zero above model";

  //Vectorized version
  double rarr[ntest];
  model.getR(ntest, fluxdens, bm, rarr);
  for (unsigned int i = 0; i < ntest; ++i) {
    reldiff = fabs(rarr[i] - rexpect[i]) / rexpect[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R value not as expected -- wanted: " << rexpect[i] <<
      " got " << rarr[i] << " relative difference: " << reldiff;
  }
}

//Test getting R with histogrammed beam
TEST(model1DTest, getRHist) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const float knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  paramSet p(nknots, knotval);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should consider itself valid";

  beam bm("testdata/band1_beam.fits", true, 120, 1e-5);
  ASSERT_TRUE(bm.isPosHist()) << "Beam should have weights after binning";

  const unsigned int ntest = 3;
  const double fluxdens[ntest] = {0.003, 0.011, 0.03};
  const double rexpect[ntest] = {1891.52, 0.817009, 0.000313631};
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens[i], bm);
    reldiff = fabs(rval - rexpect[i]) / rexpect[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R value not as expected -- wanted: " << rexpect[i] <<
      " got " << rval << " relative difference: " << reldiff;
  }

  //Vectorized version
  double rarr[ntest];
  model.getR(ntest, fluxdens, bm, rarr);
  for (unsigned int i = 0; i < ntest; ++i) {
    reldiff = fabs(rarr[i] - rexpect[i]) / rexpect[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R value not as expected -- wanted: " << rexpect[i] <<
      " got " << rarr[i] << " relative difference: " << reldiff;
  }
}


//Same tests, but beam has negative components
TEST(model1DTest, getRNeg) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const float knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  paramSet p(nknots, knotval);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should consider itself valid";

  // First, non-histogrammed beam
  beam bm("testdata/band12_beam_f100.fits[0]", false, 120, 1e-5);
  ASSERT_FALSE(bm.isPosHist()) << "Pos beam should not have weights";
  ASSERT_FALSE(bm.isNegHist()) << "Neg beam should not have weights";
  ASSERT_FLOAT_EQ(8.0, bm.getPixSize()) << "Wrong pixel size";

  const unsigned int ntest = 5;
  const double fluxdens[ntest] = {-0.03, -0.003, -1e-4, 0.003, 0.03};
  const double rexpect[ntest] = {0.0, 0.019075718, 11916638.0, 969.44687,
				 7.2612855e-5};
				 
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens[i], bm);
    if (rexpect[i] > 0)
      reldiff = fabs(rval - rexpect[i]) / rexpect[i];
    else
      reldiff = fabs(rval - rexpect[i]);
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R value not as expected -- wanted: " << rexpect[i] <<
      " got " << rval << " relative difference: " << reldiff <<
      " for flux: " << fluxdens[i];
  }

  //Vectorized version
  double rarr[ntest];
  model.getR(ntest, fluxdens, bm, rarr);
  for (unsigned int i = 0; i < ntest; ++i) {
    if (rexpect[i] > 0)
      reldiff = fabs(rarr[i] - rexpect[i]) / rexpect[i];
    else
      reldiff = fabs(rarr[i] - rexpect[i]);
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R value not as expected -- wanted: " << rexpect[i] <<
      " got " << rarr[i] << " relative difference: " << reldiff <<
      " for flux density: " << fluxdens[i];
  }

  // Repeat with histogramming
  bm.makeHistogram(200);
  ASSERT_TRUE(bm.isPosHist()) << "Pos beam should have weights";
  ASSERT_TRUE(bm.isNegHist()) << "Neg beam should have weights";
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens[i], bm);
    if (rexpect[i] > 0)
      reldiff = fabs(rval - rexpect[i]) / rexpect[i];
    else
      reldiff = fabs(rval - rexpect[i]);
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R value not as expected -- wanted: " << rexpect[i] <<
      " got " << rval << " relative difference: " << reldiff <<
      " for flux: " << fluxdens[i];
  }
  model.getR(ntest, fluxdens, bm, rarr);
  for (unsigned int i = 0; i < ntest; ++i) {
    if (rexpect[i] > 0)
      reldiff = fabs(rarr[i] - rexpect[i]) / rexpect[i];
    else
      reldiff = fabs(rarr[i] - rexpect[i]);
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R value not as expected -- wanted: " << rexpect[i] <<
      " got " << rarr[i] << " relative difference: " << reldiff <<
      " for flux density: " << fluxdens[i];
  }
}


//////////////////////////////////////////////////
// numberCountsDoubleLogNormal

//Basic Instantiation
TEST(model2DTest, Init) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const float sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 2;
  const float offsetpos[noffsets] = { 0.03, 0.05 };
  
  //Just set size version, no positions set
  numberCountsDoubleLogNormal model1(nknots, nsigmas, noffsets);
  EXPECT_EQ(nknots, model1.getNKnots()) << 
    "Unexpected number of knots in int, int, int constructor";
  EXPECT_EQ(nsigmas, model1.getNSigmas()) << 
    "Unexpected number of sigmas in int, int, int constructor";
  EXPECT_EQ(noffsets, model1.getNOffsets()) << 
    "Unexpected number of Offsets in int, int, int constructor";
  EXPECT_EQ(nknots + nsigmas + noffsets, model1.getNParams()) <<
    "Unexpected number of total parameters after int, int, int constructor";
  EXPECT_FALSE(model1.isValid()) <<
    "Model should not be valid after int, int, int constructor";
  
  //C array version
  numberCountsDoubleLogNormal model2(nknots, knotpos, nsigmas, sigmapos,
				     noffsets, offsetpos);
  EXPECT_EQ(nknots, model2.getNKnots()) << 
    "Unexpected number of knots in c-array constructor";
  EXPECT_EQ(nsigmas, model2.getNSigmas()) << 
    "Unexpected number of sigmas in c-array constructor";
  EXPECT_EQ(noffsets, model2.getNOffsets()) << 
    "Unexpected number of Offsets in c-array constructor";
  EXPECT_FALSE(model2.isValid()) <<
    "Model should not be valid after c-array constructor";
  std::vector<double> kp, sp, op;
  model2.getPositions(kp, sp, op);
  for (unsigned int i = 0; i < nknots; ++i) 
    EXPECT_FLOAT_EQ(knotpos[i], kp[i]) <<
      "Knot position mismatch after c-array constructor at index: " << i;
  for (unsigned int i = 0; i < nsigmas; ++i) 
    EXPECT_FLOAT_EQ(sigmapos[i], sp[i]) <<
      "Sigma position mismatch after c-array constructor at index: " << i;
  for (unsigned int i = 0; i < noffsets; ++i) 
    EXPECT_FLOAT_EQ(offsetpos[i], op[i]) <<
      "Offset position mismatch after c-array constructor at index: " << i;

  //Vector version -- note that kp, sp, op are now set to the kpositions
  numberCountsDoubleLogNormal model3(kp, sp, op);
  EXPECT_EQ(nknots, model3.getNKnots()) << 
    "Unexpected number of knots in vector constructor";
  EXPECT_EQ(nsigmas, model3.getNSigmas()) << 
    "Unexpected number of sigmas in vector constructor";
  EXPECT_EQ(noffsets, model3.getNOffsets()) << 
    "Unexpected number of Offsets in vector constructor";
  EXPECT_FALSE(model3.isValid()) <<
    "Model should not be valid after vector constructor";
  std::vector<double> kp2, sp2, op2;
  model3.getPositions(kp2, sp2, op2);
  for (unsigned int i = 0; i < nknots; ++i) 
    EXPECT_FLOAT_EQ(kp[i], kp2[i]) <<
      "Knot position mismatch after vector constructor at index: " << i;
  for (unsigned int i = 0; i < nsigmas; ++i) 
    EXPECT_FLOAT_EQ(sp[i], sp2[i]) <<
      "Sigma position mismatch after vector constructor at index: " << i;
  for (unsigned int i = 0; i < noffsets; ++i) 
    EXPECT_FLOAT_EQ(op[i], op2[i]) <<
      "Offset position mismatch after vector constructor at index: " << i;

  //And, finally, copy constructor
  numberCountsDoubleLogNormal model4(model3);
  EXPECT_EQ(nknots, model4.getNKnots()) << 
    "Unexpected number of knots in copy constructor";
  EXPECT_EQ(nsigmas, model4.getNSigmas()) << 
    "Unexpected number of sigmas in copy constructor";
  EXPECT_EQ(noffsets, model4.getNOffsets()) << 
    "Unexpected number of Offsets in copy constructor";
  EXPECT_FALSE(model4.isValid()) <<
    "Model should not be valid after copy constructor";
  model4.getPositions(kp2, sp2, op2);
  for (unsigned int i = 0; i < nknots; ++i) 
    EXPECT_FLOAT_EQ(kp[i], kp2[i]) <<
      "Knot position mismatch after copy constructor at index: " << i;
  for (unsigned int i = 0; i < nsigmas; ++i) 
    EXPECT_FLOAT_EQ(sp[i], sp2[i]) <<
      "Sigma position mismatch after copy constructor at index: " << i;
  for (unsigned int i = 0; i < noffsets; ++i) 
    EXPECT_FLOAT_EQ(op[i], op2[i]) <<
      "Offset position mismatch after copy constructor at index: " << i;

}

//Test set positions
TEST(model2DTest, SetPositions) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const float sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 2;
  const float offsetpos[noffsets] = { 0.03, 0.05 };

  //Make it with fewer knots -- set should resize
  numberCountsDoubleLogNormal model1(nknots-1, nsigmas, noffsets+1);
  ASSERT_EQ(nknots-1, model1.getNKnots()) << 
    "Unexpected number of knots after int, int, int constructor";
  ASSERT_EQ(nsigmas, model1.getNSigmas()) << 
    "Unexpected number of sigma knots after int, int, int constructor";
  ASSERT_EQ(noffsets+1, model1.getNOffsets()) << 
    "Unexpected number of offset knots after int, int, int constructor";
  model1.setKnotPositions(nknots, knotpos);
  model1.setSigmaPositions(nsigmas, sigmapos);
  model1.setOffsetPositions(noffsets, offsetpos);
  ASSERT_EQ(nknots, model1.getNKnots()) << "setKnotPositions failed to resize";
  ASSERT_EQ(nsigmas, model1.getNSigmas()) << 
    "setSigmaPositions failed to resize";
  ASSERT_EQ(noffsets, model1.getNOffsets()) << 
    "setOffsetPositions failed to resize";
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_FLOAT_EQ(knotpos[i], model1.getKnotPosition(i)) <<
      "After setKnotPositions, didn't find expected knot position at " << i;
  for (unsigned int i = 0; i < nsigmas; ++i)
    EXPECT_FLOAT_EQ(sigmapos[i], model1.getSigmaPosition(i)) <<
      "After setSigmaPositions, didn't find expected knot position at " << i;
  for (unsigned int i = 0; i < noffsets; ++i)
    EXPECT_FLOAT_EQ(offsetpos[i], model1.getOffsetPosition(i)) <<
      "After setOffsetPositions, didn't find expected offset position at " << i;
  
  //Reject negative knot positions
  const float knotpos2[nknots] = { -0.002, 0.005, 0.010, 0.020, 0.040 };
  EXPECT_THROW(model1.setKnotPositions(nknots, knotpos2), affineExcept) <<
    "Trying to set a negative knot position should throw an exception";
  EXPECT_THROW(model1.setSigmaPositions(nknots, knotpos2), affineExcept) <<
    "Trying to set a negative sigma position should throw an exception";
  EXPECT_THROW(model1.setOffsetPositions(nknots, knotpos2), affineExcept) <<
    "Trying to set a negative offset position should throw an exception";
}

//Test set params
TEST(model2DTest, SetGetParams) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const float sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const float offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const float pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  ASSERT_EQ(nknots, model.getNKnots()) << "Model has wrong number of knots";
  ASSERT_EQ(nsigmas, model.getNSigmas()) << "Model has wrong number of sigmas";
  ASSERT_EQ(noffsets, model.getNOffsets()) << 
    "Model has wrong number of offsets";
  EXPECT_FALSE(model.isValid()) << "Model should not be valid";

  paramSet p(ntot, pvals);
  ASSERT_EQ(ntot, p.getNParams()) << "paramSet has wrong number of values";

  model.setParams(p);
  ASSERT_EQ(nknots, model.getNKnots()) <<
    "Model has wrong number of knots after setParams";
  EXPECT_TRUE(model.isValid()) << "Model should be valid after setParams";

  paramSet p2(1);
  EXPECT_THROW(model.getParams(p2), affineExcept) <<
    "Should throw exception if wrong number of params in getParams paramSet";
  p2.setNParams(ntot);
  model.getParams(p2);
  ASSERT_EQ(ntot, p2.getNParams()) << 
    "Output params had wrong size after getParams";
  for (unsigned int i = 0; i < ntot; ++i)
    EXPECT_FLOAT_EQ(p[i], p2[i]) << 
      "Recovered wrong parameter value after set/get params at index: " << i;
}


//Test model min/max
TEST(model2DTest, getModelRange) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const float sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const float offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const float pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  dblpair minmax;
  minmax = model.getMinFlux();
  EXPECT_NEAR(knotpos[0], minmax.first, 1e-5)
    << "Unexpected band 1 minimum flux density";
  minmax = model.getMaxFlux();
  EXPECT_NEAR(0.040, minmax.first, 1e-5) 
    << "Band 1 model max not as expected";
  EXPECT_NEAR(0.0636028, minmax.second, 1e-3) 
    << "Band 2 model max not as expected";
}

//Test getting R range
TEST(model2DTest, getRRange) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const float sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const float offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const float pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  doublebeam bm("testdata/band1_beam.fits",
		"testdata/band2_beam.fits");
  
  ASSERT_NEAR(bm.getMinMax1(0).second, 0.97798554, 1e-5)
    << "Unexpected beam band 1 maximum";
  
  dblpair bmrange = bm.getMinMax1(0);
  ASSERT_NEAR(0.97798554, bmrange.second, 1e-5)
    << "Unexpected beam band 1 maximum";
  bmrange = bm.getMinMax2(0);
  ASSERT_NEAR(0.98850347, bmrange.second, 1e-5)
    << "Unexpected beam band 2 maximum";

  std::pair<dblpair, dblpair> rng = model.getRRange(bm);
  EXPECT_FLOAT_EQ(0, rng.first.first) << 
    "Model band 1 lower R range should be zero";
  EXPECT_NEAR(0.039119421, rng.first.second, 1e-4) 
    << "Unexpected model upper R range band 1";
  EXPECT_FLOAT_EQ(0, rng.second.first) << 
    "Model band 2 lower R range should be zero";
  EXPECT_NEAR(0.062871593, rng.second.second, 1e-4) 
    << "Unexpected model upper R range band 2";
}


//For the following tests (NumberCounts, Integration, getR)
// it's hard to say what the 'correct' value is because the models
// are complicated.  But we can still record values from the time
// they were written as a check against things changing unexpectedly later.
TEST(model2DTest, NumberCounts) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const float sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const float offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const float pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  const unsigned int ntest = 4;
  const double fluxdens1[ntest] = {0.003, 0.011, 0.011, 0.03};
  const double fluxdens2[ntest] = {0.003, 0.010, 0.015, 0.028};
  const double numexp[ntest] = { 2.09837e+10, 1.67999e+06, 94125.8,
				 29739.2 };

  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getNumberCounts(fluxdens1[i], fluxdens2[i]);
    reldiff = fabs(rval - numexp[i])/numexp[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "Number counts not as expected -- wanted: " << numexp[i] <<
      " got: " << rval << " for " << fluxdens1[i] << " " << fluxdens2[i];
  }
}

TEST(model2DTest, Integration) {
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const float sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const float offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const float pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  EXPECT_NEAR(57266.8, model.getNS(), 10.0) <<
    "Total number of sources not as expected";
  EXPECT_NEAR(152.262, model.getFluxPerArea(0), 0.1) <<
    "Flux per area in band 1 not as expected";
  EXPECT_NEAR(148.585, model.getFluxPerArea(1), 0.1) <<
    "Flux per area in band 2 not as expected";
}

TEST(model2DTest, getR) {

  // Simple model that can be tested against IDL code
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 1;
  const float sigmapos[nsigmas] = { 0.015 };
  const unsigned int noffsets = 1;
  const float offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const float pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.15, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  doublebeam bm("testdata/band1_beam.fits",
		"testdata/band2_beam.fits");
  
  // Comparison is direct computation from IDL without histogramming
  const unsigned int ntest = 4;
  const double fluxdens1[ntest] = {0.003, 0.011, 0.011, 0.03};
  const double fluxdens2[ntest] = {0.003, 0.010, 0.015, 0.028};
  const double rexp[ntest] = { 76926.027, 9.6701326, 11.520839, 0.12455423 };
  
  //Scalar version
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens1[i], fluxdens2[i], bm);
    reldiff = fabs(rval - rexp[i]) / rexp[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "Rnot as expected -- wanted: " << rexp[i] <<
      " got: " << rval << " for " << fluxdens1[i] << " " << fluxdens2[i];
  }

  // Vector version
  const unsigned int ntest_1 = 3;
  const unsigned int ntest_2 = 2;
  const double fluxdens1_2d[ntest_1] = {0.003, 0.011, 0.02};
  const double fluxdens2_2d[ntest_2] = {0.003, 0.015};
  
  // Based on IDL direct evaluation
  const double rexp_2d[ntest_1 * ntest_2] = 
    {76926.027, 0.024892799, 1.2063052e-15, 11.520839,
     1.9838477e-34, 0.21958590};
  double rval_2d[ntest_1 * ntest_2];
  model.getR(ntest_1, fluxdens1_2d, ntest_2, fluxdens2_2d,
	     bm, rval_2d);
  unsigned int idx;
  for (unsigned int i = 0; i < ntest_1; ++i)
    for (unsigned int j = 0; j < ntest_2; ++j) {
      idx = i * ntest_2 + j;
      reldiff = fabs(rval_2d[idx] - rexp_2d[idx]) / rexp_2d[idx];
      EXPECT_NEAR(0.0, reldiff, 1e-3) <<
	"Rnot as expected -- wanted: " << rexp_2d[idx] <<
	" got: " << rval_2d[idx] << " for " << fluxdens1[i] << " " << 
	fluxdens2[j];
    }
}

//Same test, but with histogrammed beam
TEST(model2DTest, getRHist) {

  // Simple model that can be tested against IDL code
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 1;
  const float sigmapos[nsigmas] = { 0.020 };
  const unsigned int noffsets = 1;
  const float offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const float pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.15, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  doublebeam bm("testdata/band1_beam.fits", "testdata/band2_beam.fits",
		true, 150);
  ASSERT_TRUE(bm.isHistogrammed(0)) << "Beam should be histogrammed";

  // The comparison values are from a non-histogrammed beam
  // with 0 minflux, so shouldn't be precisely the same.
  const unsigned int ntest = 4;
  const double fluxdens1[ntest] = {0.003, 0.011, 0.011, 0.03};
  const double fluxdens2[ntest] = {0.003, 0.010, 0.015, 0.028};
  const double rexp[ntest] = { 76926.027, 9.6701326, 11.520839, 0.12455423 };
  
  //Scalar version
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens1[i], fluxdens2[i], bm);
    reldiff = fabs(rval - rexp[i])/rexp[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "Rnot as expected -- wanted: " << rexp[i] <<
      " got: " << rval << " for " << fluxdens1[i] << " " << fluxdens2[i];
  }

  // Array version
  const unsigned int ntest_1 = 3;
  const unsigned int ntest_2 = 2;
  const double fluxdens1_2d[ntest_1] = {0.003, 0.011, 0.02};
  const double fluxdens2_2d[ntest_2] = {0.003, 0.015};
  
  // Based on IDL direct evaluation
  const double rexp_2d[ntest_1 * ntest_2] = 
    {76926.027, 0.024892799, 1.2063052e-15, 11.520839,
     1.9838477e-34, 0.21958590};
  double rval_2d[ntest_1 * ntest_2];
  model.getR(ntest_1, fluxdens1_2d, ntest_2, fluxdens2_2d,
	     bm, rval_2d);
  unsigned int idx;
  for (unsigned int i = 0; i < ntest_1; ++i)
    for (unsigned int j = 0; j < ntest_2; ++j) {
      idx = i * ntest_2 + j;
      reldiff = fabs(rval_2d[idx] - rexp_2d[idx]) / rexp_2d[idx];
      EXPECT_NEAR(0.0, reldiff, 1e-3) <<
	"Rnot as expected -- wanted: " << rexp_2d[idx] <<
	" got: " << rval_2d[idx] << " for " << fluxdens1[i] << " " << 
	fluxdens2[j];
    }
}

// Check with beam that has negative bits
TEST(model2DTest, getRNeg) {

  // Simple model that can be tested against IDL code
  const unsigned int nknots = 5;
  const float knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 1;
  const float sigmapos[nsigmas] = { 0.015 };
  const unsigned int noffsets = 1;
  const float offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const float pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.15, -0.3 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  doublebeam bm("testdata/band12_beam_f100.fits[0]",
		"testdata/band12_beam_f100.fits[1]", false, 150, 1e-7);
  
  // Comparison is direct computation from IDL without histogramming
  const unsigned int ntest = 4;
  const double fluxdens1[ntest] = {-2e-4, -2e-4, 0.003, 0.011};
  const double fluxdens2[ntest] = {-3e-4, 1e-4, 0.002, 0.010};
  const double rexp[ntest] = { 1.596908e8, 12727124, 75840.624, 20.409118 };
  
  //Scalar version
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens1[i], fluxdens2[i], bm);
    reldiff = fabs(rval - rexp[i])/rexp[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R not as expected -- wanted: " << rexp[i] <<
      " got: " << rval << " for " << fluxdens1[i] << " " << fluxdens2[i];
  }

  EXPECT_FLOAT_EQ(0, model.getR(0.050, 0.030, bm))
    << "R should be zero outside model coverage";
  EXPECT_FLOAT_EQ(0, model.getR(0.050, -0.030, bm))
    << "R should be zero outside model coverage";
  EXPECT_FLOAT_EQ(0, model.getR(-0.050, 0.030, bm))
    << "R should be zero outside model coverage";
  EXPECT_FLOAT_EQ(0, model.getR(-0.050, -0.030, bm))
    << "R should be zero outside model coverage";

  // Vector version
  const unsigned int ntest_1 = 3;
  const unsigned int ntest_2 = 2;
  const double fluxdens1_2d[ntest_1] = {-3e-4, -1e-4, 1e-3};
  const double fluxdens2_2d[ntest_2] = {-2e-4, 0.5e-3};
  
  // Based on IDL direct evaluation
  const double rexp_2d[ntest_1 * ntest_2] = 
    {37281355, 771.54407, 45388812, 49.613048, 0.0, 7762.3524 };
  double rval_2d[ntest_1 * ntest_2];
  model.getR(ntest_1, fluxdens1_2d, ntest_2, fluxdens2_2d,
	     bm, rval_2d);
  unsigned int idx;
  for (unsigned int i = 0; i < ntest_1; ++i)
    for (unsigned int j = 0; j < ntest_2; ++j) {
      idx = i * ntest_2 + j;
      if (rexp_2d[idx] > 0)
	reldiff = fabs(rval_2d[idx] - rexp_2d[idx]) / rexp_2d[idx];
      else
	reldiff = fabs(rval_2d[idx] - rexp_2d[idx]);
      EXPECT_NEAR(0.0, reldiff, 1e-3) <<
	"R not as expected -- wanted: " << rexp_2d[idx] <<
	" got: " << rval_2d[idx] << " for " << fluxdens1[i] << " " << 
	fluxdens2[j];
    }

  // Histogram and try again; pos-neg shouldn't be histed (too few elems)
  bm.makeHistogram(200);
  ASSERT_TRUE(bm.isHistogrammed(0)) << "Pos-pos should be histogrammed";
  ASSERT_FALSE(bm.isHistogrammed(1)) << "Pos-neg should not be histogrammed";
  ASSERT_TRUE(bm.isHistogrammed(2)) << "Neg-pos should be histogrammed";
  ASSERT_TRUE(bm.isHistogrammed(3)) << "Neg-neg should be histogrammed";

  // Scalar
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens1[i], fluxdens2[i], bm);
    reldiff = fabs(rval - rexp[i])/rexp[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "R not as expected -- wanted: " << rexp[i] <<
      " got: " << rval << " for " << fluxdens1[i] << " " << fluxdens2[i];
  }
  // Vector
  model.getR(ntest_1, fluxdens1_2d, ntest_2, fluxdens2_2d,
	     bm, rval_2d);
  for (unsigned int i = 0; i < ntest_1; ++i)
    for (unsigned int j = 0; j < ntest_2; ++j) {
      idx = i * ntest_2 + j;
      if (rexp_2d[idx] > 0)
	reldiff = fabs(rval_2d[idx] - rexp_2d[idx]) / rexp_2d[idx];
      else
	reldiff = fabs(rval_2d[idx] - rexp_2d[idx]);
      EXPECT_NEAR(0.0, reldiff, 1e-3) <<
	"R not as expected -- wanted: " << rexp_2d[idx] <<
	" got: " << rval_2d[idx] << " for " << fluxdens1[i] << " " << 
	fluxdens2[j];
    }

}


////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running model tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
