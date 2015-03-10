#include<iostream>

#include<gtest/gtest.h>

#include "../include/paramSet.h"

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

  // Set not using []
  p.setParamValue(2, 7);
  EXPECT_FLOAT_EQ(p[2], 7.0) << "setParamValue failed";

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

  // Set array
  float v3[3] = {-2, 11.5, 21.7};
  ASSERT_NO_THROW(p.setParamValues(3, v3)) << "Setting param values failed";
  for (unsigned int i = 0; i < 3; ++i)
    ASSERT_FLOAT_EQ(v3[i], p[i]) << "Array set didn't set right values";

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

//Test distance
TEST(paramSetTest, Distance) {
  const unsigned int nparams = 2;
  paramSet p1(nparams), p2(nparams);
  p1[0] = 2.0; p1[1] = 2.0;
  p2[0] = -1.0; p2[1] = 6.0;
  ASSERT_FLOAT_EQ(5.0, p1.getDist(p2)) << "Distance not correct";
}


////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running chain_tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
