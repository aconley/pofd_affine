//Testing PD and PDFactory classes
#include<cmath>

#include<gtest/gtest.h>

#include "../include/global_settings.h"
#include "../include/affineExcept.h"
#include "../include/PD.h"
#include "../include/PDDouble.h"

////////////////////////////////////////////
// PD tests
TEST(pd1DTest, Init) {
  PD pd(1024, -0.3, 1e-3);
  EXPECT_EQ(1024U, pd.getDim()) << "Unexpected dimension of PD";
  EXPECT_FLOAT_EQ(-0.3, pd.minflux) << "Unexpected minimum flux";
  EXPECT_FLOAT_EQ(1e-3, pd.dflux) << "Unexpected minimum flux";
}

TEST(pd1DTest, Data) {
  const unsigned int n = 1024;
  PD pd;
  
  double dat[n];
  for (unsigned int i = 0; i < n; ++i) 
    dat[i] = exp(-0.5 * (i-200.) * (i-200.) / 2500.0);
  pd.fill(n, -0.3, 1e-3, dat, false);
  EXPECT_EQ(n, pd.getDim()) << "Unexpected size";
  EXPECT_FLOAT_EQ(-0.3, pd.minflux) << "Unexpected minimum flux";
  EXPECT_FLOAT_EQ(1e-3, pd.dflux) << "Unexpected minimum flux";
  EXPECT_FALSE(pd.isLog()) << "Shouldn't be log";
  
  unsigned int tstidx[4] = {0, 20, 133, 1023};
  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_NEAR(-0.3 + 1e-3 * tstidx[i], pd.getFluxVal(tstidx[i]), 1e-6)
      << "Unexpected x axis value for test index: " << tstidx[i];
  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(dat[tstidx[i]], pd[tstidx[i]])
      << "Unexpected pd value at index: " << tstidx[i];
  // This tests interpolation at the index
  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(dat[tstidx[i]], pd.getPDVal(pd.getFluxVal(tstidx[i])))
      << "Unexpected pd interpolant value at index: " << tstidx[i];
  // This tests at the midpoint; can't use the 4th point as it's off the end
  for (unsigned int i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(0.5*(dat[tstidx[i]] + dat[tstidx[i]+1]),
		    pd.getPDVal(pd.getFluxVal(tstidx[i]) + 0.5*pd.dflux))
      << "Unexpected pd interpolant value between indices: "
      << tstidx[i] << " and " << tstidx[i]+1;
  }

  // Log and re-revaluate.  Use a smaller range of test indices
  //  because some of the values will be effectively 0 in the pre-log PD
  pd.applyLog();
  ASSERT_TRUE(pd.isLog()) << "Should be logged now";
  unsigned int tstidxlog[4] = { 160, 200, 230, 255 };
  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(log2(dat[tstidxlog[i]]), pd[tstidxlog[i]])
      << "Unexpected pd value at index: " << tstidxlog[i] << " after log2ing";
}

////////////////////////////////////////////
// PD2D tests
TEST(pd2DTest, Init) {
  PDDouble pd(1024, -0.3, 1e-3, 512, -0.2, 2e-3, false);
  EXPECT_EQ(1024U, pd.getDim1()) << "Unexpected dimension 1 of PD";
  EXPECT_FLOAT_EQ(-0.3, pd.minflux1) << "Unexpected minimum flux band 1";
  EXPECT_FLOAT_EQ(1e-3, pd.dflux1) << "Unexpected minimum flux band 1";
  EXPECT_EQ(512U, pd.getDim2()) << "Unexpected dimension 2 of PD";
  EXPECT_FLOAT_EQ(-0.2, pd.minflux2) << "Unexpected minimum flux band 2";
  EXPECT_FLOAT_EQ(2e-3, pd.dflux2) << "Unexpected minimum flux band 2";
}

TEST(pd2DTest, Data) {
  const unsigned int n = 256;
  PDDouble pd;
  
  double distsq;
  double dat[n * n];
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j) {
      distsq = (i-50.0) * (i - 50.0) + (j - 25) * (j - 25);
      dat[i * n + j] = exp(- distsq / 800.0);
    }
  pd.fill(n, -0.3, 1e-3, n, -0.1, 3e-3, dat, false);
  EXPECT_EQ(n, pd.getDim1()) << "Unexpected size dim 1";
  EXPECT_EQ(n, pd.getDim2()) << "Unexpected size dim 2";
  EXPECT_FLOAT_EQ(-0.3, pd.minflux1) << "Unexpected minimum flux band 1";
  EXPECT_FLOAT_EQ(1e-3, pd.dflux1) << "Unexpected minimum flux band 1";
  EXPECT_FLOAT_EQ(-0.1, pd.minflux2) << "Unexpected minimum flux band 2";
  EXPECT_FLOAT_EQ(3e-3, pd.dflux2) << "Unexpected minimum flux band 2";
  EXPECT_FALSE(pd.isLog()) << "Shouldn't be log";
  
  unsigned int tstidx[4] = {10, 20, 30};
  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_NEAR(-0.3 + 1e-3 * tstidx[i], pd.getFluxVal1(tstidx[i]), 1e-6)
      << "Unexpected x axis value for test index: " << tstidx[i] 
      << " band 1";
  for (unsigned int i = 0; i < 4; ++i)
    EXPECT_NEAR(-0.1 + 3e-3 * tstidx[i], pd.getFluxVal2(tstidx[i]), 1e-6)
      << "Unexpected x axis value for test index: " << tstidx[i] 
      << " band 2";
  for (unsigned int i = 0; i < 4; ++i)
    for (unsigned int j = 0; j < 4; ++j)
      EXPECT_FLOAT_EQ(dat[tstidx[i]*n + tstidx[j]], pd[tstidx[i]][tstidx[j]])
	<< "Unexpected pd value at index: " << tstidx[i] << ", "
	<< tstidx[j];

  // Log and re-revaluate.  Use a smaller range of test indices
  //  because some of the values will be effectively 0 in the pre-log PD
  pd.applyLog();
  ASSERT_TRUE(pd.isLog()) << "Should be logged now";
  unsigned int tstidxlog1[4] = { 45, 50, 55, 60 };
  unsigned int tstidxlog2[4] = { 18, 23, 27, 29 };
  for (unsigned int i = 0; i < 4; ++i)
    for (unsigned int j = 0; j < 4; ++j) {
      unsigned int idx = tstidxlog1[i] * n + tstidxlog2[j];
      EXPECT_FLOAT_EQ(log2(dat[idx]), pd[tstidxlog1[i]][tstidxlog2[j]])
		<< "Unexpected pd value at index: " << tstidxlog1[i]
		<< ", " << tstidxlog2[j] << " after log2ing";
    }
}


////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running PD tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
