//Program for testing affineQueue
#include<iostream>
#include<utility>
#include<gtest/gtest.h>

#include "../include/affineQueue.h"
#include "../include/affineExcept.h"
#include "../include/ran.h"

TEST(queueTest, Initialization) {
  const unsigned int cap = 25;

  affineQueue<unsigned int> queue(cap);
  EXPECT_EQ(cap, queue.capacity()) << "Unexpected capacity";
  EXPECT_TRUE(queue.empty()) << "Queue should be empty";
  EXPECT_EQ(0U, queue.size()) << "Queue should have 0 elements";
}

TEST(queueTest, InsertRecover) {
  const unsigned int cap = 25;
  ran ranstruct;

  affineQueue<unsigned int> queue(cap);

  //Insert 25, make sure that's okay
  for (unsigned int i = 0; i < 5; ++i) queue.push(i);
  for (unsigned int i = 5; i < cap; ++i) queue.push(ranstruct.int32());
  EXPECT_EQ(cap, queue.capacity()) << "Unexpected capacity";
  EXPECT_FALSE(queue.empty()) << "Queue should not be empty";
  EXPECT_EQ(25U, queue.size()) << "Queue should have 25 elements";
  
  //Shouldn't be able to insert more
  EXPECT_THROW(queue.push(0), affineExcept) <<
    "Trying to add another element to a full queue should throw an exception";

  //Pop first few
  for (unsigned int i = 0; i < 5; ++i)
    EXPECT_EQ(i, queue.pop()) << "Didn't get expected pop value on " <<
      i << "th pop";

  EXPECT_EQ(cap, queue.capacity()) << "Unexpected capacity";
  EXPECT_EQ(20U, queue.size()) << "Queue should have 20 elements after 5 pops";

  //Pop the rest
  for (unsigned int i = 5; i < cap; ++i) queue.pop();
  EXPECT_EQ(0U, queue.size()) << "Queue should have 0 elements now";
  EXPECT_TRUE(queue.empty()) << "Queue should now be empty";
  EXPECT_EQ(cap, queue.capacity()) << "Unexpected capacity after full pop";

  // Shouldn't be able to pop more
  EXPECT_THROW(queue.pop(), affineExcept) << 
    "Popping an empty queue should throw exception";


  //Try a more complicated element
  affineQueue< std::pair<unsigned int, int> > queue2(cap);
  EXPECT_EQ(cap, queue.capacity()) << "Unexpected capacity in pair queue";
  EXPECT_TRUE(queue.empty()) << "Pair queue should be empty";
  
  std::pair<unsigned int, int> val2, val3;
  val2.first = 4;
  val2.second = -2;
  queue2.push(val2);
  val2.first = 5;
  val2.second = -4;
  queue2.push(val2);
  for (unsigned int i = 2; i < 10; ++i) {
    val3.first = ranstruct.int32();
    val3.second = ranstruct.int32();
    queue2.push(val3);
  }
  EXPECT_EQ(cap, queue2.capacity()) << 
    "Unexpected capacity in pair after 10 push";
  EXPECT_EQ(10U, queue2.size()) << "Unexpected size for pair queue";
  val3 = queue2.pop();
  EXPECT_EQ(9U, queue2.size()) << "Unexpected size for pair queue after pop";
  val3 = queue2.pop();
  EXPECT_EQ(8U, queue2.size()) << 
    "Unexpected size for pair queue after second pop";
  //val2 should be val3 after this series of operations
  EXPECT_EQ(val2.second, val3.second) <<
    "Second elem of second pair pop not as expected";
  EXPECT_EQ(val2.second, val3.second) <<
    "Second elem of second pair pop not as expected";
}

TEST(queueTest, Clear) {
  const unsigned int cap = 25;
  ran ranstruct;

  affineQueue<unsigned int> queue(cap);
  for (unsigned int i = 0; i < cap-3; ++i) queue.push(ranstruct.int32());
  EXPECT_EQ(cap, queue.capacity()) << "Unexpected capacity after partial fill";
  EXPECT_EQ(cap-3, queue.size()) << "Unexpected number of queue elements";
  
  queue.clear();
  EXPECT_EQ(cap, queue.capacity()) << "Unexpected capacity after clear";
  EXPECT_TRUE(queue.empty()) << "Queue should now be empty";
  EXPECT_EQ(0U, queue.size()) << "Queue should have 0 elements now";
}

TEST(queueTest, setCapacity) {
  const unsigned int cap = 25;
  affineQueue<int> queue(cap);
  
  for (unsigned int i = 0; i < 8; ++i)
    queue.push(i);

  EXPECT_EQ(cap, queue.capacity()) << "Unexpected capacity";
  EXPECT_EQ(8U, queue.size()) << "Unexpected size";
  EXPECT_FALSE(queue.empty()) << "Queue should not be empty after pushes";

  queue.setCapacity(5);
  EXPECT_EQ(5U, queue.capacity()) << "Unexpected capacity after setCapacity";
  EXPECT_EQ(0U, queue.size()) << "Queue should have 0 elements after resize";
  EXPECT_TRUE(queue.empty()) << "Queue should be empty after resize";

  for (unsigned int i = 0; i < 3; ++i)
    queue.push(i);
  EXPECT_EQ(3U, queue.size()) << "Unexpected size";
  EXPECT_FALSE(queue.empty()) << "Queue should not be empty after pushes";

  //Try set capacity without changing size
  queue.setCapacity(queue.size());
  EXPECT_EQ(3U, queue.capacity()) << "Unexpected capacity after setCapacity";
  EXPECT_EQ(0U, queue.size()) << "Queue should have 0 elements after resize";
  EXPECT_TRUE(queue.empty()) << "Queue should be empty after resize";
}
      
      
////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running queue tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

