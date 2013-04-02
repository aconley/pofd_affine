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
//beam
TEST(beam1DTest, Read) {
  beam bm("testdata/band1_beam.fits");
  
  EXPECT_TRUE(bm.hasPos()) << "Beam should have positive pixels";
  EXPECT_EQ(625U, bm.getNPos()) << "Beam had wrong number of positive pixels";
  EXPECT_FALSE(bm.hasPosWeights()) << "Beam should not have positive weights";
  EXPECT_FALSE(bm.hasNeg()) << "Beam should have no negative pixels";
  EXPECT_EQ(0U, bm.getNNeg()) << "Beam had wrong number of negative pixels";
  EXPECT_FALSE(bm.hasNegWeights()) << "Beam should not have negative weights";
  EXPECT_NEAR(4.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";
  EXPECT_NEAR(2.8327253e-5, bm.getEffectiveArea(), 1e-8) <<
    "Beam had wrong effective area";
  EXPECT_NEAR(2.8327253e-5, bm.getEffectiveAreaPos(), 1e-8) <<
    "Beam had wrong positive effective area";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaNeg()) <<
    "Beam should have zero negative area";
  EXPECT_NEAR(22.945073, bm.getEffectiveAreaPix(), 1e-4) << 
    "Didn't get expected beam area in pixels";
  EXPECT_NEAR(1.3852891e-5, bm.getEffectiveAreaSq(), 1e-8) <<
    "Beam^2 had wrong effective area";
  EXPECT_NEAR(0.97798554, bm.getMaxPos(), 1e-5) <<
    "Beam had unexpected maximum pixel value";
}

TEST(beam1DTest, CopyConstructor) {
  beam bm("testdata/band1_beam.fits");
  beam bm2(bm);

  EXPECT_TRUE(bm2.hasPos()) << "Beam should have positive pixels";
  EXPECT_EQ(bm.getNPos(), bm2.getNPos()) << 
    "Beam had wrong number of positive pixels";
  EXPECT_FALSE(bm2.hasNeg()) << "Beam should have no negative pixels";
  EXPECT_EQ(bm.getNNeg(), bm2.getNNeg()) << 
    "Beam had wrong number of negative pixels";
  EXPECT_FLOAT_EQ(bm.getPixSize(), bm2.getPixSize()) << "Wrong pixel size";
  EXPECT_FLOAT_EQ(bm.getEffectiveArea(), bm2.getEffectiveArea()) <<
    "Beam had wrong effective area";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaPos(), bm2.getEffectiveAreaPos()) <<
    "Beam had wrong positive effective area";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaNeg(), bm2.getEffectiveAreaNeg()) <<
    "Beam should have zero negative area";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaPix(), bm2.getEffectiveAreaPix()) <<
    "Didn't get expected beam area in pixels";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSq(), bm2.getEffectiveAreaSq()) <<
    "Beam^2 had wrong effective area";
  EXPECT_FLOAT_EQ(bm.getMaxPos(), bm2.getMaxPos()) <<
    "Beam had unexpected maximum pixel value";
}

TEST(beam1DTest, Equals) {
  beam bm("testdata/band1_beam.fits");
  beam bm2("testdata/band2_beam.fits");
  bm2 = bm;

  EXPECT_TRUE(bm2.hasPos()) << "Beam should have positive pixels";
  EXPECT_EQ(bm.getNPos(), bm2.getNPos()) << 
    "Beam had wrong number of positive pixels";
  EXPECT_FALSE(bm2.hasNeg()) << "Beam should have no negative pixels";
  EXPECT_EQ(bm.getNNeg(), bm2.getNNeg()) << 
    "Beam had wrong number of negative pixels";
  EXPECT_FLOAT_EQ(bm.getPixSize(), bm2.getPixSize()) << "Wrong pixel size";
  EXPECT_FLOAT_EQ(bm.getEffectiveArea(), bm2.getEffectiveArea()) <<
    "Beam had wrong effective area";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaPos(), bm2.getEffectiveAreaPos()) <<
    "Beam had wrong positive effective area";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaNeg(), bm2.getEffectiveAreaNeg()) <<
    "Beam should have zero negative area";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaPix(), bm2.getEffectiveAreaPix()) <<
    "Didn't get expected beam area in pixels";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSq(), bm2.getEffectiveAreaSq()) <<
    "Beam^2 had wrong effective area";
  EXPECT_FLOAT_EQ(bm.getMaxPos(), bm2.getMaxPos()) <<
    "Beam had unexpected maximum pixel value";
}

TEST(beam1DTest, Histogram) {
  beam bm("testdata/band1_beam.fits", true, 0.2);
  
  ASSERT_TRUE(bm.hasPosWeights()) << "Beam should have positive weights";
  EXPECT_TRUE(bm.hasPos()) << "Beam should have positive pixels";
  EXPECT_EQ(76U, bm.getNPos()) << "Beam had wrong number of positive pixels";
  EXPECT_FALSE(bm.hasNeg()) << "Beam should have no negative pixels";
  EXPECT_EQ(0U, bm.getNNeg()) << "Beam had wrong number of negative pixels";
  EXPECT_NEAR(4.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";
  EXPECT_NEAR(2.8327253e-5, bm.getEffectiveArea(), 1e-8) <<
    "Beam had wrong effective area";
  EXPECT_NEAR(2.8327253e-5, bm.getEffectiveAreaPos(), 1e-8) <<
    "Beam had wrong positive effective area";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaNeg()) <<
    "Beam should have zero negative area";
  EXPECT_NEAR(22.945073, bm.getEffectiveAreaPix(), 1e-4) << 
    "Didn't get expected beam area in pixels";
  EXPECT_NEAR(1.3852891e-5, bm.getEffectiveAreaSq(), 1e-8) <<
    "Beam^2 had wrong effective area";
  EXPECT_NEAR(0.97798554, bm.getMaxPos(), 1e-5) <<
    "Beam had unexpected maximum pixel value";

  //Check histogram
  const double *wts = bm.getPosWeights();
  double exp_posweights[5] = {4.0, 8.0, 8.0, 4.0, 8.0};
  for (unsigned int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(exp_posweights[i], wts[i]) << "Didn't get expected weight";
  const double* iparr = bm.getPosInvPixArr();
  const double exp_inv[5] = {1.99712, 1.74686, 1.33647, 1.169, 1.02251};
  for (unsigned int i = 71; i < 76; ++i)
    EXPECT_NEAR(exp_inv[i-71], iparr[i], 0.001) << 
      "Got unexpected binned inverse pixel value in index: " << i;

}

//////////////////////////////////////////////////
// doublebeam
TEST(beam2DTest, Read) {
  doublebeam bm("testdata/band1_beam.fits",
		"testdata/band2_beam.fits");
  
  EXPECT_TRUE(bm.hasSign(0)) << "Beam should have pos-pos pixels";
  EXPECT_FALSE(bm.hasSign(1)) << "Beam should not have pos-neg pixels";
  EXPECT_FALSE(bm.hasSign(2)) << "Beam should not have neg-pos pixels";
  EXPECT_FALSE(bm.hasSign(3)) << "Beam should not have neg-neg";
  EXPECT_EQ(625U, bm.getTotalNPix()) << "Wrong number of total pixels in beam";
  EXPECT_EQ(625U, bm.getNPix(0)) << "Beam had wrong number of pos-pos pixels";
  EXPECT_EQ(0U, bm.getNPix(1)) << "Beam had wrong number of pos-neg pixels";
  EXPECT_EQ(0U, bm.getNPix(2)) << "Beam had wrong number of neg-pos pixels";
  EXPECT_EQ(0U, bm.getNPix(3)) << "Beam had wrong number of neg-neg pixels";
  EXPECT_NEAR(4.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";

  //Beam1 tests
  EXPECT_NEAR(2.8327253e-5, bm.getEffectiveArea1(), 1e-8) <<
    "Beam had wrong effective area in band1";
  EXPECT_NEAR(2.8327253e-5, bm.getEffectiveAreaSign1(0), 1e-8) <<
    "Beam had wrong pos-pos effective area in band1";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign1(1)) <<
    "Beam had wrong pos-neg effective area in band1";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign1(2)) <<
    "Beam had wrong neg-pos effective area in band1";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign1(3)) <<
    "Beam had wrong neg-neg effective area in band1";
  EXPECT_NEAR(0.97798554, bm.getMax1(0), 1e-5) <<
    "Beam had unexpected maximum pixel value in band1";

  //Beam2 tests
  EXPECT_NEAR(5.4643351e-5, bm.getEffectiveArea2(), 1e-8) <<
    "Beam had wrong effective area in band2";
  EXPECT_NEAR(5.4643351e-5, bm.getEffectiveAreaSign2(0), 1e-8) <<
    "Beam had wrong pos-pos effective area in band2";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign2(1)) <<
    "Beam had wrong pos-neg effective area in band2";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign2(2)) <<
    "Beam had wrong neg-pos effective area in band2";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign2(3)) <<
    "Beam had wrong neg-neg effective area in band2";
  EXPECT_NEAR(0.98850347, bm.getMax2(0), 1e-5) <<
    "Beam had unexpected maximum pixel value in band 2";
}

TEST(beam2DTest, CopyConstructor) {
  doublebeam bm("testdata/band1_beam.fits",
		"testdata/band2_beam.fits");
  doublebeam bm2(bm);
  
  EXPECT_TRUE(bm2.hasSign(0)) << "Beam should have pos-pos pixels";
  EXPECT_FALSE(bm2.hasSign(1)) << "Beam should not have pos-neg pixels";
  EXPECT_FALSE(bm2.hasSign(2)) << "Beam should not have neg-pos pixels";
  EXPECT_FALSE(bm2.hasSign(3)) << "Beam should not have neg-neg";
  EXPECT_EQ(bm.getTotalNPix(), bm2.getTotalNPix()) << 
    "Wrong number of total pixels in beam";
  EXPECT_EQ(bm.getNPix(0), bm2.getNPix(0)) << 
    "Beam had wrong number of pos-pos pixels";
  EXPECT_EQ(bm.getNPix(1), bm2.getNPix(1)) << 
    "Beam had wrong number of pos-neg pixels";
  EXPECT_EQ(bm.getNPix(2), bm2.getNPix(2)) << 
    "Beam had wrong number of neg-pos pixels";
  EXPECT_EQ(bm.getNPix(3), bm2.getNPix(3)) << 
    "Beam had wrong number of neg-neg pixels";
  EXPECT_FLOAT_EQ(bm.getPixSize(), bm2.getPixSize()) << "Wrong pixel size";

  //Beam1 tests
  EXPECT_FLOAT_EQ(bm.getEffectiveArea1(), bm2.getEffectiveArea1()) <<
    "Beam had wrong effective area in band1";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSign1(0), bm2.getEffectiveAreaSign1(0)) <<
    "Beam had wrong pos-pos effective area in band1";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSign1(1), bm2.getEffectiveAreaSign1(1)) <<
    "Beam had wrong pos-neg effective area in band1";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSign1(2), bm2.getEffectiveAreaSign1(2)) <<
    "Beam had wrong neg-pos effective area in band1";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSign1(3), bm2.getEffectiveAreaSign1(3)) <<
    "Beam had wrong neg-neg effective area in band1";
  EXPECT_FLOAT_EQ(bm.getMax1(0), bm2.getMax1(0)) <<
    "Beam had unexpected maximum pixel value in band1";

  //Beam2 tests
  EXPECT_FLOAT_EQ(bm.getEffectiveArea2(), bm2.getEffectiveArea2()) <<
    "Beam had wrong effective area in band2";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSign2(0), bm2.getEffectiveAreaSign2(0)) <<
    "Beam had wrong pos-pos effective area in band2";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSign2(1), bm2.getEffectiveAreaSign2(1)) <<
    "Beam had wrong pos-neg effective area in band2";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSign2(2), bm2.getEffectiveAreaSign2(2)) <<
    "Beam had wrong neg-pos effective area in band2";
  EXPECT_FLOAT_EQ(bm.getEffectiveAreaSign2(3), bm2.getEffectiveAreaSign2(3)) <<
    "Beam had wrong neg-neg effective area in band2";
  EXPECT_FLOAT_EQ(bm.getMax2(0), bm2.getMax2(0)) <<
    "Beam had unexpected maximum pixel value in band 2";
}

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
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

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
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

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
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  paramSet p(nknots, knotval);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should consider itself valid";

  //Evaluate at knot positions -- these should be very close
  for (unsigned int i = 1; i < nknots-1; ++i)
    EXPECT_NEAR(knotval[i], log10(model.getNumberCounts(knotpos[i])),0.0001) << 
      "Didn't get expected number counts at " << knotpos[i];
  
  const unsigned int ntest = 3;
  const double testpos[ntest] = {0.003, 0.007, 0.015};
  const double testval[ntest] = {8.69422, 5.88217, 4.53303};
  //Don't require as tight tolerances for this test
  for (unsigned int i = 0; i < ntest; ++i)
    EXPECT_NEAR(testval[i], log10(model.getNumberCounts(testpos[i])),0.001) << 
      "Didn't get expected number counts at " << testpos[i];

}

//Model integration -- also hard to ensure the values are right
TEST(model1DTest, Integration) {
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

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

//Test getting R; this depends on the details of
// the spline model, so will change if GSL changes anything
// Hence, it's a little hard to say exactly what the results should be.
TEST(model1DTest, getR) {
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

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
    rval = model.getR(fluxdens[i], bm, numberCounts::BEAMBOTH);
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

//Test getting R with histogrammed beam
TEST(model1DTest, getRHist) {
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.050 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCountsKnotsSpline model(nknots, knotpos);
  paramSet p(nknots, knotval);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should consider itself valid";

  beam bm("testdata/band1_beam.fits", true, 0.2);
  ASSERT_TRUE(bm.hasPosWeights()) << "Beam should have weights after binning";

  const unsigned int ntest = 3;
  const double fluxdens[ntest] = {0.003, 0.011, 0.03};
  const double rexpect[ntest] = {1891.52, 0.817009, 0.000313631};
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens[i], bm, numberCounts::BEAMBOTH);
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


//////////////////////////////////////////////////
// numberCountsDoubleLogNormal

//Basic Instantiation
TEST(model2DTest, Init) {
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const double sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 2;
  const double offsetpos[noffsets] = { 0.03, 0.05 };
  
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
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const double sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 2;
  const double offsetpos[noffsets] = { 0.03, 0.05 };

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
  const double knotpos2[nknots] = { -0.002, 0.005, 0.010, 0.020, 0.040 };
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
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const double sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const double offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const double pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

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

//For the following tests (NumberCounts, Integration, getR)
// it's hard to say what the 'correct' value is because the models
// are complicated.  But we can still record values from the time
// they were written as a check against things changing unexpectedly later.
TEST(model2DTest, NumberCounts) {
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const double sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const double offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const double pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

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
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const double sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const double offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const double pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

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
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const double sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const double offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const double pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  doublebeam bm("testdata/band1_beam.fits",
		"testdata/band2_beam.fits");

  const unsigned int ntest = 4;
  const double fluxdens1[ntest] = {0.003, 0.011, 0.011, 0.03};
  const double fluxdens2[ntest] = {0.003, 0.010, 0.015, 0.028};
  const double rexp[ntest] = { 97571.4, 10.9065, 9.95045, 0.102312 };
  
  //Scalar version
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens1[i], fluxdens2[i], bm);
    reldiff = fabs(rval - rexp[i])/rexp[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "Rnot as expected -- wanted: " << rexp[i] <<
      " got: " << rval << " for " << fluxdens1[i] << " " << fluxdens2[i];
  }

  const unsigned int ntest_1 = 3;
  const unsigned int ntest_2 = 2;
  const double fluxdens1_2d[ntest_1] = {0.003, 0.011, 0.02};
  const double fluxdens2_2d[ntest_2] = {0.003, 0.015};
  
  const double rexp_2d[ntest_1 * ntest_2] = 
    {97571.4, 0.277285, 8.66458e-13,
     9.95045, 1.30238e-19, 0.46487};
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
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const unsigned int nsigmas = 2;
  const double sigmapos[nsigmas] = { 0.002, 0.020 };
  const unsigned int noffsets = 1;
  const double offsetpos[noffsets] = { 0.03 };
  const unsigned int ntot = nknots + nsigmas + noffsets;
  const double pvals[ntot] = { 8.0, 6.0, 4.0, 3.3, 2.0, 0.1, 0.2, -0.03 };

  numberCountsDoubleLogNormal model(nknots, knotpos, nsigmas, sigmapos,
				    noffsets, offsetpos);
  paramSet p(ntot, pvals);
  model.setParams(p);
  ASSERT_TRUE(model.isValid()) << "Model should be valid";

  doublebeam bm("testdata/band1_beam.fits", "testdata/band2_beam.fits",
		true, 0.2);
  ASSERT_TRUE(bm.hasWeights(0)) << "Beam should have weights after binning";

  const unsigned int ntest = 4;
  const double fluxdens1[ntest] = {0.003, 0.011, 0.011, 0.03};
  const double fluxdens2[ntest] = {0.003, 0.010, 0.015, 0.028};
  const double rexp[ntest] = { 97571.4, 10.9065, 9.95045, 0.102312 };
  
  //Scalar version
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(fluxdens1[i], fluxdens2[i], bm);
    reldiff = fabs(rval - rexp[i])/rexp[i];
    EXPECT_NEAR(0.0, reldiff, 1e-3) <<
      "Rnot as expected -- wanted: " << rexp[i] <<
      " got: " << rval << " for " << fluxdens1[i] << " " << fluxdens2[i];
  }

  const unsigned int ntest_1 = 3;
  const unsigned int ntest_2 = 2;
  const double fluxdens1_2d[ntest_1] = {0.003, 0.011, 0.02};
  const double fluxdens2_2d[ntest_2] = {0.003, 0.015};
  
  const double rexp_2d[ntest_1 * ntest_2] = 
    {97571.4, 0.277285, 8.66458e-13,
     9.95045, 1.30238e-19, 0.46487};
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


////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running model tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


