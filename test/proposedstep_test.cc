#include<iostream>

#include<gtest/gtest.h>

#include "../include/paramSet.h"
#include "../include/proposedStep.h"

//Basic size, resize
TEST(proposedStepTest, Sizing) {
  const unsigned int nparams = 3;
  proposedStep p(nparams);
  ASSERT_EQ(nparams, p.getNParams()) << "proposedStep not right size";

  unsigned int newsize = 5;
  p.setNParams(newsize);
  ASSERT_EQ(newsize, p.getNParams()) << "proposedStep resizing failed";
}


//Test operator= and copy constructor
TEST(proposedStepTest, Equality) {
  const unsigned int nparams = 3;
  proposedStep p(nparams);
  const float pvalsold[3] = {1.0, 2.0, -4.0};
  const float pvalsnew[3] = {0.5, -1.1, 3.0};
  for (unsigned int i = 0; i < nparams; ++i)
    p.oldStep[i] = pvalsold[i];
  for (unsigned int i = 0; i < nparams; ++i)
    p.newStep[i] = pvalsnew[i];
  p.oldLogLike = 4.0;
  p.newLogLike = 10.0;
  p.z = 3;
  p.update_idx = 20;
  ASSERT_EQ(p.getNParams(), nparams) << "Initial proposedStep not expected size";

  proposedStep p2(p);
  EXPECT_EQ(nparams, p2.getNParams()) << "Copy constructor size not expected";
  ASSERT_FLOAT_EQ(3.0, p2.z) << "Copy constructor didn't copy z";
  ASSERT_FLOAT_EQ(4.0, p2.oldLogLike) << "Copy constructor didn't copy oldLogLike";
  ASSERT_FLOAT_EQ(10.0, p2.newLogLike) << "Copy constructor didn't copy newLogLike";
  ASSERT_FLOAT_EQ(3.0, p2.z) << "Copy constructor didn't copy z";
  ASSERT_EQ(20, p2.update_idx) << "Copy constructor didn't preserve update_idx";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p2.oldStep[i], pvalsold[i]) 
      << "Copy constructor didn't set right values in oldStep";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p2.newStep[i], pvalsnew[i]) 
      << "Copy constructor didn't set right values in oldStep";

  proposedStep p3(1);
  p3 = p;
  EXPECT_EQ(nparams, p3.getNParams()) << "Assignment size not expected";
  ASSERT_FLOAT_EQ(3.0, p3.z) << "Assignment didn't copy z";
  ASSERT_FLOAT_EQ(4.0, p3.oldLogLike) << "Assignment didn't copy oldLogLike";
  ASSERT_FLOAT_EQ(10.0, p3.newLogLike) << "Assignment didn't copy newLogLike";
  ASSERT_FLOAT_EQ(3.0, p3.z) << "Assignment didn't copy z";
  ASSERT_EQ(20, p3.update_idx) << "Assignment didn't preserve update_idx";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p3.oldStep[i], pvalsold[i]) 
      << "Assignment didn't set right values in oldStep";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p3.newStep[i], pvalsnew[i]) 
      << "Assignment didn't set right values in oldStep";
}

//Test move semantics
TEST(proposedStep, Move) {
  const unsigned int nparams = 3;
  proposedStep p(nparams);
  const float pvalsold[3] = {1.0, 2.0, -4.0};
  const float pvalsnew[3] = {0.5, -1.1, 3.0};
  for (unsigned int i = 0; i < nparams; ++i)
    p.oldStep[i] = pvalsold[i];
  for (unsigned int i = 0; i < nparams; ++i)
    p.newStep[i] = pvalsnew[i];
  p.oldLogLike = 4.0;
  p.newLogLike = 10.0;
  p.z = 3;
  p.update_idx = 20;
  

  proposedStep p2(std::move(p));
  EXPECT_EQ(nparams, p2.getNParams()) << "Move constructor size not expected";
  ASSERT_FLOAT_EQ(3.0, p2.z) << "Move constructor didn't copy z";
  ASSERT_FLOAT_EQ(4.0, p2.oldLogLike) << "Move constructor didn't copy oldLogLike";
  ASSERT_FLOAT_EQ(10.0, p2.newLogLike) << "Move constructor didn't copy newLogLike";
  ASSERT_FLOAT_EQ(3.0, p2.z) << "Move constructor didn't copy z";
  ASSERT_EQ(20, p2.update_idx) << "Move constructor didn't preserve update_idx";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p2.oldStep[i], pvalsold[i]) 
      << "Move constructor didn't set right values in oldStep";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p2.newStep[i], pvalsnew[i]) 
      << "Move constructor didn't set right values in oldStep";
  ASSERT_EQ(0, p.getNParams()) << "Move constructor didn't empty other";

  // Move assign back
  p = std::move(p2);
  EXPECT_EQ(nparams, p.getNParams()) << "Move assignemnt size not expected";
  ASSERT_FLOAT_EQ(3.0, p.z) << "Move assignment didn't copy z";
  ASSERT_FLOAT_EQ(4.0, p.oldLogLike) << "Move assignment didn't copy oldLogLike";
  ASSERT_FLOAT_EQ(10.0, p.newLogLike) << "Move assignment didn't copy newLogLike";
  ASSERT_FLOAT_EQ(3.0, p.z) << "Move assignment didn't copy z";
  ASSERT_EQ(20, p.update_idx) << "Move assignment didn't preserve update_idx";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p.oldStep[i], pvalsold[i]) 
      << "Move assignment didn't set right values in oldStep";
  for (unsigned int i = 0; i < nparams; ++i)
    ASSERT_FLOAT_EQ(p.newStep[i], pvalsnew[i]) 
      << "Move assignment didn't set right values in oldStep";
  ASSERT_EQ(0, p2.getNParams()) << "Move assignment didn't empty other";
}

////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running proposedStep tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
