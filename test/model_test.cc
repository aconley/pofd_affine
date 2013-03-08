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

//////////////////////////////////////////////////
//beam
TEST(beam1DTest, Read) {
  beam bm("testdata/band1_beam.fits");
  
  EXPECT_TRUE(bm.hasPos()) << "Beam should have positive pixels";
  EXPECT_EQ(625U, bm.getNPos()) << "Beam had wrong number of positive pixels";
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

////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running model tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


