#include<iostream>
#include<string>
#include<cmath>
#include<gtest/gtest.h>

#include "../include/specFile.h"
#include "../include/specFileDouble.h"
#include "../include/affineExcept.h"

//In both cases the test is to read in test files and make sure they
// look like what is expected
TEST(specFile1DTest, ReadCheck) {

  try {
    specFile spec_info("testdata/specfile_test1D.txt");
  
    ASSERT_EQ(spec_info.datafiles.size(), 2U) << 
      "Unexpected number of datasets";
    EXPECT_EQ(spec_info.datafiles[0], std::string("file1.fits")) <<
      "Expected first data file to be file1.fits, but it is " <<
      spec_info.datafiles[0];
    EXPECT_EQ(spec_info.datafiles[1], std::string("file2.fits")) <<
      "Expected second data file to be file2.fits, but it is " <<
      spec_info.datafiles[1];
    
    ASSERT_EQ(spec_info.sigmas.size(), 2U) <<
      "Unexpected number of sigmas";
    EXPECT_NEAR(spec_info.sigmas[0], 0.0014, 0.0001) <<
      "Expected first instrument sigma to be 0.0014";
    EXPECT_NEAR(spec_info.sigmas[1], 0.002, 0.0001) <<
      "Expected second instrument sigma to be 0.002";
    
    ASSERT_EQ(spec_info.psffiles.size(), 2U) <<
      "Unexpected number of psffiles";
    EXPECT_EQ(spec_info.psffiles[0], std::string("beam1.fits")) <<
      "Expected first beam file to be beam1.fits, but it is " << 
      spec_info.psffiles[0];
    EXPECT_EQ(spec_info.psffiles[1], std::string("beam3.fits")) <<
      "Expected second beam file to be beam3.fits, but it is " << 
      spec_info.psffiles[1];
    
    ASSERT_EQ(spec_info.like_norm.size(), 2U) << 
      "Unexpected number of like_norms";
    EXPECT_NEAR(spec_info.like_norm[0], 1.0, 0.0001) <<
      "Expected first like norm to be 1, but it is " << spec_info.like_norm[0];
    EXPECT_NEAR(spec_info.like_norm[1], 1.1, 0.0001) <<
      "Expected second like norm to be 1.1, but it is " << spec_info.like_norm[1];
    
    EXPECT_TRUE(spec_info.bin_data) << "Expected data binning to be on";
    EXPECT_EQ(spec_info.nbins, 500U) << "Expected nbins to be 500 but it is " << 
      spec_info.nbins;
    
    EXPECT_TRUE(spec_info.mean_sub) << "Expected mean_sub binning to be on";
    
    EXPECT_FALSE(spec_info.ignore_mask) << "Expected ignore_mask to be off";
    
    EXPECT_EQ(spec_info.fftsize, 2048U) <<
      "Expected fftsize to be 2048 but it is " << spec_info.fftsize;
    
    EXPECT_EQ(spec_info.ninterp, 1200U) <<
      "Expected ninterp to be 1200 but it is " << spec_info.ninterp;
    
    EXPECT_TRUE(spec_info.beam_histogram) << 
      "Expected beam histogramming to be on";
    
    EXPECT_TRUE(spec_info.fit_sigma) << "Expected fit_sigma to be on";
    
    ASSERT_TRUE(spec_info.has_sigprior) << "Expected to have sigma prior";
    EXPECT_NEAR(spec_info.sigprior_stdev, 0.15, 0.0001) <<
      "Unexpected sigma prior stdev value";
    
    EXPECT_EQ(spec_info.exp_conf, 0.0) <<
      "Expected zero expected confusion noise";
    
    EXPECT_FLOAT_EQ(spec_info.minbeamval, 1e-6) <<
      "Unexpected minimum beam value";
    
    EXPECT_EQ(spec_info.nbeamhist, 120) <<
      "Unexpected number of beam histogram bins";
    
    ASSERT_TRUE(spec_info.has_cfirbprior) << "Expected to have CFIRB prior";
    EXPECT_NEAR(spec_info.cfirbprior_mean, 1.05, 0.0001) <<
      "Expected cfirb prior mean to be 1.05, but it is " <<
      spec_info.cfirbprior_mean;
    EXPECT_NEAR(spec_info.cfirbprior_stdev, 0.1, 0.0001) <<
      "Expected cfirb prior stdev to be 0.1, but it is " <<
      spec_info.cfirbprior_stdev;
    
    ASSERT_TRUE(spec_info.has_wisdom_file) << "Expected to have wisdom file";
    EXPECT_EQ(spec_info.wisdom_file, std::string("a_wisdom_file.txt")) <<
      "Expected wisdom file name 'a_wisdom_file.txt', got " << 
      spec_info.wisdom_file;
    
    EXPECT_EQ(1, spec_info.verbosity) << "Expected level 1 verbosity";
  } catch (const affineExcept& ex) {
    std::cout << "Error processing specFile1D" << std::endl;
    std::cout << ex << std::endl;
    throw ex;
  }
}


TEST(specFile2DTest, ReadCheck) {
  
  specFileDouble spec_info("testdata/specfile_test2D.txt");
    
  ASSERT_EQ(spec_info.datafiles1.size(), 2U) <<
    "Unexpected number of datasets1";
  ASSERT_EQ(spec_info.datafiles2.size(), 2U) <<
    "Unexpected number of datasets2";
  EXPECT_EQ(spec_info.datafiles1[0], std::string("file1_1.fits")) <<
    "Expected first data file band 1 to be file1_1.fits, but it is " <<
    spec_info.datafiles1[0];
  EXPECT_EQ(spec_info.datafiles2[0], std::string("file1_2.fits")) <<
    "Expected first data file band 2 to be file1_2.fits, but it is " <<
    spec_info.datafiles2[0];

  EXPECT_EQ(spec_info.datafiles1[1], std::string("file2.fits")) <<
      "Expected second data file band 1 to be file2.fits, but it is " << 
    spec_info.datafiles1[1];
  EXPECT_EQ(spec_info.datafiles2[1], std::string("file4.fits")) <<
      "Expected second data file band 2 to be file4.fits, but it is " << 
    spec_info.datafiles2[1];

  ASSERT_EQ(spec_info.sigmas1.size(), 2U) << "Unexpected number of sigmas1";
  EXPECT_NEAR(spec_info.sigmas1[0], 0.0014, 0.00001) <<
    "Expected first instrument sigma band 1 to be 0.0014, " <<
    "but it is " << spec_info.sigmas1[0];
  EXPECT_NEAR(spec_info.sigmas2[0], 0.003, 0.000001) <<
    "Expected first instrument sigma (band 2) to be 0.003, " << 
    "but it is " << spec_info.sigmas2[0];
  EXPECT_NEAR(spec_info.sigmas1[1], 0.002, 0.000001) << 
    "Expected second instrument sigma (band 1) to be 0.002";
  EXPECT_NEAR(spec_info.sigmas2[1], 0.04, 0.000001) << 
    "Expected second instrument sigma (band 2) to be 0.04";

  ASSERT_EQ(spec_info.psffiles1.size(), 2U) << "Unexpected number of psffiles1";
  ASSERT_EQ(spec_info.psffiles2.size(), 2U) << "Unexpected number of psffiles2";
  EXPECT_EQ(spec_info.psffiles1[0], std::string("beam1.fits")) <<
    "Expected first beam file (band 1) to be beam1.fits, but it is " <<
    spec_info.psffiles1[0];
  EXPECT_EQ(spec_info.psffiles2[0], std::string("beam5.fits")) <<
    "Expected first beam file (band 2) to be beam5.fits, but it is " << 
    spec_info.psffiles2[0];
  EXPECT_EQ(spec_info.psffiles1[1], std::string("beam3.fits")) <<
    "Expected second beam file (band 1) to be beam3.fits, " << 
    "but it is " << spec_info.psffiles1[1];
  EXPECT_EQ(spec_info.psffiles2[1], std::string("funny_beam.fits")) <<
    "Expected second beam file (band 2) to be funny_beam.fits, " << 
    "but it is " << spec_info.psffiles2[1];

  ASSERT_EQ(spec_info.like_norm.size(), 2U) << 
    "Unexpected number of like_norms";
  EXPECT_NEAR(spec_info.like_norm[0], 1.0, 0.000001) <<
    "Expected first like norm to be 1";
  EXPECT_NEAR(spec_info.like_norm[1], 1.1, 0.0000001) <<
    "Expected second like norm to be 1.1";
  
  ASSERT_TRUE(spec_info.bin_data) << "Expected data binning to be on";
  EXPECT_EQ(spec_info.nbins, 300U) <<
    "Expected nbins to be 300 but it is " << spec_info.nbins;

  EXPECT_TRUE(spec_info.mean_sub) << "Expected mean_sub binning to be on";

  EXPECT_TRUE(spec_info.ignore_mask) << "Expected ignore_mask to be on";

  EXPECT_EQ(spec_info.fftsize, 2048U) <<
      "Expected fftsize to be 2048 but it is " << spec_info.fftsize;

  ASSERT_TRUE(spec_info.edge_set) << "Expected edge_set to be on";
  EXPECT_EQ(spec_info.nedge, 100U) <<
    "Expected nedge to be 100 but it is " << spec_info.nedge;

  EXPECT_TRUE(spec_info.beam_histogram) << 
    "Expected beam histogramming to be on";

  EXPECT_TRUE(spec_info.fit_sigma1) << "Expected fit_sigma1 to be on";
  ASSERT_TRUE(spec_info.has_sigprior1) <<
    "Expected to have sigma prior in band 1";
  EXPECT_NEAR(spec_info.sigprior_stdev1, 0.15, 0.000001) <<
    "Expected sigma prior1 stdev to be 0.15";
  EXPECT_FALSE(spec_info.fit_sigma2) << "Expected fit_sigma2 to be off";
  EXPECT_FALSE(spec_info.has_sigprior2) << "Expected no band 2 sigma prior";

  EXPECT_EQ(spec_info.exp_conf1, 0.0) <<
    "Expected zero expected confusion noise, band 1";
  EXPECT_EQ(spec_info.exp_conf2, 0.006) <<
    "Expected non-zero expected confusion noise, band 2";

  EXPECT_EQ(spec_info.nbeamhist, 150) <<
    "Unexpected number of beam histogram bins";

  EXPECT_FALSE(spec_info.has_cfirbprior1) << "Expected no band 1 cfirb prior";
  ASSERT_TRUE(spec_info.has_cfirbprior2) << 
    "Expected to have CFIRB prior in band 2";
  EXPECT_NEAR(spec_info.cfirbprior_mean2, 1.05, 0.000001) <<
    "Expected cfirb prior2 mean to be 1.05";
  EXPECT_NEAR(spec_info.cfirbprior_stdev2, 0.1, 0.0000001) << 
    "Expected cfirb prior2 stdev to be 0.1";

  ASSERT_TRUE(spec_info.has_wisdom_file) << "Expected to have wisdom file";
  EXPECT_EQ(spec_info.wisdom_file, std::string("a_wisdom_file.txt")) <<
    "Expected wisdom file name 'a_wisdom_file.txt', got " << 
    spec_info.wisdom_file;

  EXPECT_EQ(1, spec_info.verbosity) << "Expected level 1 verbosity";
}

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running specFile tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

