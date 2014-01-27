//Testing models
#include<iostream>
#include<fstream>
#include<string>

#include<gtest/gtest.h>

#include "../include/global_settings.h"
#include "../include/affineExcept.h"
#include "../include/beam.h"
#include "../include/doublebeam.h"

//////////////////////////////////////////////////
//beam
TEST(beam1DTest, Read) {
  beam bm;
  ASSERT_FALSE(bm.hasData()) << "Beam should not have data before reading";

  bm.readFile("testdata/band1_beam.fits");
  ASSERT_TRUE(bm.hasData()) << "Beam should have data";

  EXPECT_TRUE(bm.hasPos()) << "Beam should have positive pixels";
  EXPECT_EQ(277U, bm.getNPos()) << "Beam had wrong number of positive pixels";
  EXPECT_FALSE(bm.isPosHist()) << "Beam should not have positive weights";
  EXPECT_FALSE(bm.hasNeg()) << "Beam should have no negative pixels";
  EXPECT_EQ(0U, bm.getNNeg()) << "Beam had wrong number of negative pixels";
  EXPECT_FALSE(bm.isNegHist()) << "Beam should not have negative weights";
  EXPECT_NEAR(4.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";
  EXPECT_NEAR(2.8327039e-5, bm.getEffectiveArea(), 1e-8) <<
    "Beam had wrong effective area";
  EXPECT_NEAR(2.8327039e-5, bm.getEffectiveAreaPos(), 1e-8) <<
    "Beam had wrong positive effective area";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaNeg()) <<
    "Beam should have zero negative area";
  EXPECT_NEAR(22.944899, bm.getEffectiveAreaPix(), 1e-4) << 
    "Didn't get expected beam area in pixels";
  EXPECT_NEAR(1.3852891e-5, bm.getEffectiveAreaSq(), 1e-8) <<
    "Beam^2 had wrong effective area";
  EXPECT_NEAR(0.97798554, bm.getMinMaxPos().second, 1e-5) <<
    "Beam had unexpected maximum pixel value";
}

TEST(beam1DTest, ReadConstructor) {
  beam bm("testdata/band1_beam.fits");
  ASSERT_TRUE(bm.hasData()) << "Beam should have data";

  EXPECT_TRUE(bm.hasPos()) << "Beam should have positive pixels";
  EXPECT_EQ(277U, bm.getNPos()) << "Beam had wrong number of positive pixels";
  EXPECT_FALSE(bm.isPosHist()) << "Beam should not have positive weights";
  EXPECT_FALSE(bm.hasNeg()) << "Beam should have no negative pixels";
  EXPECT_EQ(0U, bm.getNNeg()) << "Beam had wrong number of negative pixels";
  EXPECT_FALSE(bm.isNegHist()) << "Beam should not have negative weights";
  EXPECT_NEAR(4.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";
}

TEST(beam1DTest, CopyConstructor) {
  beam bm("testdata/band1_beam.fits");
  beam bm2(bm);

  ASSERT_TRUE(bm2.hasData()) << "Beam should have data";
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
  EXPECT_FLOAT_EQ(bm.getMinMaxPos().first, bm2.getMinMaxPos().first) <<
    "Beam had unexpected minimum pixel value";
  EXPECT_FLOAT_EQ(bm.getMinMaxPos().second, bm2.getMinMaxPos().second) <<
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
  EXPECT_FLOAT_EQ(bm.getMinMaxPos().first, bm2.getMinMaxPos().first) <<
    "Beam had unexpected minimum pixel value";
  EXPECT_FLOAT_EQ(bm.getMinMaxPos().second, bm2.getMinMaxPos().second) <<
    "Beam had unexpected maximum pixel value";
}

TEST(beam1DTest, Histogram) {
  beam bm("testdata/band1_beam.fits", true, 120, 1e-5);
  
  ASSERT_TRUE(bm.hasPos()) << "Beam should have positive pixels";
  ASSERT_TRUE(bm.isPosHist()) << "Beam should have positive weights";
  EXPECT_EQ(39U, bm.getNHistPos()) 
    << "Beam had wrong number of positive pixels";
  EXPECT_FALSE(bm.hasNeg()) << "Beam should have no negative pixels";
  EXPECT_EQ(0U, bm.getNNeg()) << "Beam had wrong number of negative pixels";
  EXPECT_NEAR(4.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";
  EXPECT_NEAR(2.8327039e-5, bm.getEffectiveArea(), 1e-8) <<
    "Beam had wrong effective area";
  EXPECT_NEAR(2.8327039e-5, bm.getEffectiveAreaPos(), 1e-8) <<
    "Beam had wrong positive effective area";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaNeg()) <<
    "Beam should have zero negative area";
  EXPECT_NEAR(22.944899, bm.getEffectiveAreaPix(), 1e-6) << 
    "Didn't get expected beam area in pixels";
  EXPECT_NEAR(1.3852891e-5, bm.getEffectiveAreaSq(), 1e-8) <<
    "Beam^2 had wrong effective area";
  EXPECT_NEAR(0.97798554, bm.getMinMaxPos().second, 1e-5) <<
    "Beam had unexpected maximum pixel value";
  EXPECT_FLOAT_EQ(1e-5, bm.getMinval()) << "Unexpected minimum acceptance val";

  //Check histogram
  const double *wts = bm.getPosHistWeights();
  const double* iparr = bm.getPosHist();
  double exp_posweights[5] = {8.0, 4.0, 4.0, 4.0, 1.0};
  for (unsigned int i = 34; i < 39; ++i)
    EXPECT_FLOAT_EQ(exp_posweights[i-34], wts[i]) << 
      "Didn't get expected weight";
  const double exp_inv[5] = {1.99712, 1.74686, 1.33647, 1.169, 1.02251};
  for (unsigned int i = 34; i < 39; ++i)
    EXPECT_NEAR(exp_inv[i-34], iparr[i], 0.0001) << 
      "Got unexpected binned inverse pixel value in index: " << i;

}

// Now do some similar tests on beams with negative components
TEST(beam1DTest, NegBeam) {
  beam bm("testdata/band12_beam_f100.fits[0]", false, 120, 1e-7);
  ASSERT_TRUE(bm.hasData()) << "Beam should have data";
  EXPECT_FLOAT_EQ(1e-7, bm.getMinval()) << "Unexpected minimum acceptance val";
  ASSERT_TRUE(bm.hasPos()) << "Beam should have positive pixels";
  EXPECT_EQ(477U, bm.getNPos()) << "Beam had wrong number of positive pixels";
  ASSERT_TRUE(bm.hasNeg()) << "Beam should have negative pixels";
  EXPECT_EQ(364U, bm.getNNeg()) << "Beam had wrong number of positive pixels";
  EXPECT_FALSE(bm.isPosHist()) << "Beam should not have positive weights";
  EXPECT_FALSE(bm.isNegHist()) << "Beam should not have negative weights";
  EXPECT_NEAR(8.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";
  EXPECT_NEAR(1.0765914e-4, bm.getEffectiveArea(), 1e-8) <<
    "Beam had wrong effective area";
  EXPECT_NEAR(5.3829571e-5, bm.getEffectiveAreaPos(), 1e-8) <<
    "Beam had wrong positive effective area";
  EXPECT_FLOAT_EQ(5.3829571e-5, bm.getEffectiveAreaNeg()) <<
    "Beam had wrong negative effective area";
  EXPECT_NEAR(0.82171527, bm.getMinMaxPos().second, 1e-5) <<
    "Beam had unexpected positive maximum pixel value";
  EXPECT_NEAR(0.082822058, bm.getMinMaxNeg().second, 1e-6) <<
    "Beam had unexpected negative maximum pixel value";
}


//////////////////////////////////////////////////
// doublebeam
TEST(beam2DTest, Read) {
  doublebeam bm;
  ASSERT_FALSE(bm.hasData()) << "Beam should not have data before read";

  bm.readFiles("testdata/band1_beam.fits", "testdata/band2_beam.fits", 0.0);
  ASSERT_TRUE(bm.hasData()) << "Beam should have data after read";

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
  EXPECT_FLOAT_EQ(0.0, bm.getMinval()) << "Unexpected minimum value";

  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_FALSE(bm.isHistogrammed(i)) 
      << "Beam should not be histogrammed for component: " << i;

  // Re read with standard minval
  bm.readFiles("testdata/band1_beam.fits", "testdata/band2_beam.fits", 1e-6);
  EXPECT_TRUE(bm.hasSign(0)) << "Beam should have pos-pos pixels";
  EXPECT_FALSE(bm.hasSign(1)) << "Beam should not have pos-neg pixels";
  EXPECT_FALSE(bm.hasSign(2)) << "Beam should not have neg-pos pixels";
  EXPECT_FALSE(bm.hasSign(3)) << "Beam should not have neg-neg";
  EXPECT_EQ(325U, bm.getTotalNPix()) << "Wrong number of total pixels in beam";
  EXPECT_EQ(325U, bm.getNPix(0)) << "Beam had wrong number of pos-pos pixels";
  EXPECT_EQ(0U, bm.getNPix(1)) << "Beam had wrong number of pos-neg pixels";
  EXPECT_EQ(0U, bm.getNPix(2)) << "Beam had wrong number of neg-pos pixels";
  EXPECT_EQ(0U, bm.getNPix(3)) << "Beam had wrong number of neg-neg pixels";
  EXPECT_NEAR(4.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";
  EXPECT_FLOAT_EQ(1e-6, bm.getMinval()) << "Unexpected minimum value";

  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_FALSE(bm.isHistogrammed(i)) 
      << "Beam should not be histogrammed for component: " << i;

  //Beam1 tests
  EXPECT_NEAR(2.8327225e-5, bm.getEffectiveArea1(), 1e-8) <<
    "Beam had wrong effective area in band1";
  EXPECT_NEAR(2.8327225e-5, bm.getEffectiveAreaSign1(0), 1e-8) <<
    "Beam had wrong pos-pos effective area in band1";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign1(1)) <<
    "Beam had wrong pos-neg effective area in band1";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign1(2)) <<
    "Beam had wrong neg-pos effective area in band1";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign1(3)) <<
    "Beam had wrong neg-neg effective area in band1";
  EXPECT_NEAR(0.97798554, bm.getMinMax1(0).second, 1e-5) <<
    "Beam had unexpected maximum pixel value in band1";

  //Beam2 tests
  EXPECT_NEAR(5.4604548e-5, bm.getEffectiveArea2(), 1e-8) <<
    "Beam had wrong effective area in band2";
  EXPECT_NEAR(5.4604548e-5, bm.getEffectiveAreaSign2(0), 1e-8) <<
    "Beam had wrong pos-pos effective area in band2";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign2(1)) <<
    "Beam had wrong pos-neg effective area in band2";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign2(2)) <<
    "Beam had wrong neg-pos effective area in band2";
  EXPECT_FLOAT_EQ(0.0, bm.getEffectiveAreaSign2(3)) <<
    "Beam had wrong neg-neg effective area in band2";
  EXPECT_NEAR(0.98850347, bm.getMinMax2(0).second, 1e-5) <<
    "Beam had unexpected maximum pixel value in band 2";
}

TEST(beam2DTest, Histogram) {
  doublebeam bm("testdata/band1_beam.fits", "testdata/band2_beam.fits",
		true, 150, 1e-6);
  ASSERT_TRUE(bm.hasData()) << "Beam should have data after read";
  EXPECT_FLOAT_EQ(1e-6, bm.getMinval()) << "Unexpected minimum accepted val";
  ASSERT_TRUE(bm.isHistogrammed(0)) << "pp beam should be histogrammed";
  for (unsigned int i = 1; i < 4; ++i)
    EXPECT_FALSE(bm.isHistogrammed(i)) 
      << "Beam should not be histogrammed for component: " << i;
  EXPECT_EQ(150, bm.getNBins()) << "Unexpected number of histogram bins";
  EXPECT_EQ(45, bm.getNHist(0)) << "Unexpected number of hist bins in 0 (pp)"
				<< " component";
  for (unsigned int i = 1; i < 4; ++i)
    EXPECT_EQ(0, bm.getNHist(i)) << "Unexpected number of hist bins in "
				 << "component " << i;
}


// Now do some similar tests on beams with negative components
TEST(beam2DTest, NegBeam) {
  doublebeam bm("testdata/band12_beam_f100.fits[0]",
		"testdata/band12_beam_f100.fits[1]", false, 150, 1e-7);
  ASSERT_TRUE(bm.hasData()) << "Beam should have data";
  EXPECT_FLOAT_EQ(1e-7, bm.getMinval()) << "Unexpected minimum acceptance val";
  ASSERT_TRUE(bm.hasSign(0)) << "Beam should have pos-pos pixels";
  ASSERT_TRUE(bm.hasSign(1)) << "Beam should have pos-neg pixels";
  ASSERT_TRUE(bm.hasSign(2)) << "Beam should have neg-pos pixels";
  ASSERT_TRUE(bm.hasSign(3)) << "Beam should have neg-neg pixels";
  
  EXPECT_EQ(841U, bm.getTotalNPix()) << "Wrong number of total pixels";
  EXPECT_EQ(469U, bm.getNPix(0)) << "Wrong number of pos-pos pixels";
  EXPECT_EQ(8U, bm.getNPix(1)) << "Wrong number of pos-neg pixels";
  EXPECT_EQ(16U, bm.getNPix(2)) << "Wrong number of neg-pos pixels";
  EXPECT_EQ(348U, bm.getNPix(3)) << "Wrong number of neg-neg pixels";

  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_FALSE(bm.isHistogrammed(i)) << "Beam should not be histogrammed" <<
      " in component: " << i;

  EXPECT_NEAR(8.0, bm.getPixSize(), 0.0001) << "Wrong pixel size";

  EXPECT_NEAR(5.3796224e-5, bm.getEffectiveAreaSign1(0), 1e-8) <<
    "pos-pos beam had wrong effective area";
  EXPECT_NEAR(3.3357658e-8, bm.getEffectiveAreaSign1(1), 1e-8) <<
    "pos-neg beam had wrong effective area";
  EXPECT_NEAR(3.3918851e-6, bm.getEffectiveAreaSign1(2), 1e-8) <<
    "neg-pos beam had wrong effective area";
  EXPECT_NEAR(5.0437686e-5, bm.getEffectiveAreaSign1(3), 1e-8) <<
    "pos-pos beam had wrong effective area";

  dblpair minmax;
  minmax = bm.getMinMax1(0);
  EXPECT_NEAR(0.0010022319, minmax.first, 1e-6) <<
    "Unexpected pos-pos minimum value, band 1";
  EXPECT_NEAR(0.82171527, minmax.second, 1e-6) <<
    "Unexpected pos-pos maximum value, band 1";
  minmax = bm.getMinMax2(0);
  EXPECT_NEAR(0.0022131729, minmax.first, 1e-6) <<
    "Unexpected pos-pos minimum value, band 2";
  EXPECT_NEAR(0.74992756, minmax.second, 1e-6) <<
    "Unexpected pos-pos maximum value, band 2";
  minmax = bm.getMinMax1(2);
  EXPECT_NEAR(0.019993591, minmax.first, 1e-6) <<
    "Unexpected neg-pos minimum value, band 1";
  EXPECT_NEAR(0.055385745, minmax.second, 1e-6) <<
    "Unexpected neg-pos maximum value, band 1";
  minmax = bm.getMinMax2(2);
  EXPECT_NEAR(0.0015057219, minmax.first, 1e-6) <<
    "Unexpected neg-pos minimum value, band 2";
  EXPECT_NEAR(0.065584647, minmax.second, 1e-6) <<
    "Unexpected neg-pos maximum value, band 2";
  
  // Histogram; pos-neg shouldn't be because it has too few elements
  bm.makeHistogram(200);
  EXPECT_TRUE(bm.isHistogrammed(0)) << "Pos-pos should be histogrammed";
  EXPECT_FALSE(bm.isHistogrammed(1)) << "Pos-neg should not be histogrammed";
  EXPECT_TRUE(bm.isHistogrammed(2)) << "Neg-pos should be histogrammed";
  EXPECT_TRUE(bm.isHistogrammed(3)) << "Neg-neg should be histogrammed";
}


////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running model tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
