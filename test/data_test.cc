//Program for testing fitsData
#include<iostream>
#include<utility>
#include<gtest/gtest.h>

#include "../include/fitsData.h"
#include "../include/fitsDataDouble.h"
#include "../include/affineExcept.h"

TEST(fitsDataTest, Read) {
  const std::string testfile="testdata/testmodel2D_band1.fits";
  fitsData dat;

  EXPECT_FALSE(dat.hasData()) << "Not expected to have data";
  EXPECT_FALSE(dat.isBinned()) << "Data should not be binned";

  dat.readData(testfile);
  EXPECT_TRUE(dat.hasData()) << "After read, should have data";
  EXPECT_FALSE(dat.isBinned()) << "Data should not be binned after read";
  EXPECT_EQ(38431U, dat.getN()) << "Unexpected number of data points";

  EXPECT_NEAR(-0.011066, dat.getMin(), 1e-5) << "Didn't get expected minimum";
  EXPECT_NEAR(0.0659472, dat.getMax(), 1e-5) << "Didn't get expected maximum";

  dblpair minmax = dat.getMinMax();
  EXPECT_NEAR(-0.011066, minmax.first, 1e-5) << "Didn't get expected minimum";
  EXPECT_NEAR(0.0659472, minmax.second, 1e-5) << "Didn't get expected maximum";

  //Same, but without masking
  dat.readData(testfile, true);
  EXPECT_TRUE(dat.hasData()) << "After re read, should have data";
  EXPECT_FALSE(dat.isBinned()) << "Data should not be binned after re read";
  EXPECT_EQ(40401U, dat.getN()) << "Unexpected number of data points";

  EXPECT_NEAR(-0.011066, dat.getMin(), 1e-5) << "Didn't get expected minimum";
  EXPECT_NEAR(0.0659472, dat.getMax(), 1e-5) << "Didn't get expected maximum";

  //Read constructor
  fitsData dat2(testfile);
  EXPECT_TRUE(dat2.hasData()) << "After constructor read, should have data";
  EXPECT_FALSE(dat2.isBinned()) << "Data should not be binned after constructor read";
  EXPECT_EQ(38431U, dat2.getN()) << "Unexpected number of data points";
  minmax = dat2.getMinMax();
  EXPECT_NEAR(-0.011066, minmax.first, 1e-5)  << "Didn't get expected minimum";
  EXPECT_NEAR(0.0659472, minmax.second, 1e-5) << "Didn't get expected maximum";

}

TEST(fitsDataDoubleTest, Read) {
  const std::string testfile1="testdata/testmodel2D_band1.fits";
  const std::string testfile2="testdata/testmodel2D_band2.fits";
  fitsDataDouble dat;

  EXPECT_FALSE(dat.hasData()) << "Not expected to have data";
  EXPECT_FALSE(dat.isBinned()) << "Data should not be binned";

  dat.readData(testfile1, testfile2);
  EXPECT_TRUE(dat.hasData()) << "After read, should have data";
  EXPECT_FALSE(dat.isBinned()) << "Data should not be binned after read";
  EXPECT_EQ(38431U, dat.getN()) << "Unexpected number of data points";

  std::pair<double, double> pr;
  pr = dat.getMin();
  EXPECT_NEAR(-0.011066, pr.first, 1e-5) << 
    "Didn't get expected minimum, band 1";
  EXPECT_NEAR(-0.0176356, pr.second, 1e-5) << 
    "Didn't get expected minimum, band 2";
  pr = dat.getMax();
  EXPECT_NEAR(0.0659472, pr.first, 1e-5) << 
    "Didn't get expected maximum, band 1";
  EXPECT_NEAR(0.0776496, pr.second, 1e-5) << 
    "Didn't get expected maximum, band 2";

  //Same, but without masking
  dat.readData(testfile1, testfile2, true);
  EXPECT_TRUE(dat.hasData()) << "After re read, should have data";
  EXPECT_FALSE(dat.isBinned()) << "Data should not be binned after re read";
  EXPECT_EQ(40401U, dat.getN()) << "Unexpected number of data points";
  pr = dat.getMin();
  EXPECT_NEAR(-0.011066, pr.first, 1e-5) << 
    "Didn't get expected minimum, band 1";
  EXPECT_NEAR(-0.0176356, pr.second, 1e-5) << 
    "Didn't get expected minimum, band 2";
  pr = dat.getMax();
  EXPECT_NEAR(0.0659472, pr.first, 1e-5) << 
    "Didn't get expected maximum, band 1";
  EXPECT_NEAR(0.0776496, pr.second, 1e-5) << 
    "Didn't get expected maximum, band 2";

  // Read on construction
  fitsDataDouble dat2(testfile1, testfile2);
  EXPECT_TRUE(dat2.hasData()) << "After read, should have data";
  EXPECT_FALSE(dat2.isBinned()) << "Data should not be binned after read";
  EXPECT_EQ(38431U, dat2.getN()) << "Unexpected number of data points";

  pr = dat2.getMin();
  EXPECT_NEAR(-0.011066, pr.first, 1e-5) << 
    "Didn't get expected minimum, band 1";
  EXPECT_NEAR(-0.0176356, pr.second, 1e-5) << 
    "Didn't get expected minimum, band 2";
  pr = dat2.getMax();
  EXPECT_NEAR(0.0659472, pr.first, 1e-5) << 
    "Didn't get expected maximum, band 1";
  EXPECT_NEAR(0.0776496, pr.second, 1e-5) << 
    "Didn't get expected maximum, band 2";
}

      
////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running data tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

