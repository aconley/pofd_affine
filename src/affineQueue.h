//affineQueue.h

#ifndef __affineQueue__
#define __affineQueue__

#include<affineExcept.h>

//So, why roll my own queue class?  Why not just use std::queue?
// First, we have a fixed maximum size, so can do this in a more
// efficient, faster fashion.  But, much more importantly, I am
// experiencing memory allocation bugs with std::queue, at least
// the gcc version!

/*!
  \brief Queue (FIFO stack) of fixed capacity.  Implemented as a circular list.
 */
template <class Item> class affineQueue {
 private :
  Item* data; //!< Holds actual data
  unsigned int cap; //!< Capacity of this queue
  unsigned int nelem; //!< Number currently in queue
  unsigned int head; //!< Points at head of queue (place to remove items)
  unsigned int tail; //!< Points at where new item will be inserted

 public :

  affineQueue(); //!< Constructor
  affineQueue(unsigned int); //!< Constructor with specified capacity
  ~affineQueue(); //!< Destructor

  unsigned int capacity() const { return cap; } //!< Get capactity
  unsigned int size() const { return nelem; } //!< Get number of elements in queue
  bool empty() const { return nelem == 0; } //!< Is queue empty?

  void clear(); //!< Clear out contents (but preserve capacity)
  void setCapacity( unsigned int ); //!< Change capacity -- does not preserve contents

  void push(const Item&) throw (affineExcept);  //!< Push something onto queue
  Item pop() throw (affineExcept); //!< Pull item off of queue
};

template< class Item > affineQueue<Item>::affineQueue() {
  cap      = 0;
  nelem    = 0;
  head     = 0;
  tail     = 0;
  data = NULL;
}

/*!
  \param[in] CAP capacity
*/
template< class Item > affineQueue<Item>::affineQueue(unsigned int CAP) {
  cap      = CAP;
  nelem    = 0;
  head     = 0;
  tail     = 0;
  if (cap > 0)
    data = new Item[cap];
  else 
    data = NULL;
}

template< class Item > affineQueue<Item>::~affineQueue() {
  if (data != NULL) delete[] data;
}

/*!
  \param[in] CAP capacity
*/
template< class Item > void affineQueue<Item>::setCapacity(unsigned int CAP) {
  if (cap != CAP) {
    if (data != NULL) delete[] data;
    if (CAP > 0) data = new Item[CAP]; else data = NULL;
    cap = CAP;
  }
  nelem = 0;
  head = 0;
  tail = 0;
  return;
}

template< class Item > void affineQueue<Item>::clear() {
  nelem = 0;
  head  = 0;
  tail  = 0;
  return;
}

/*!
  \param[in] item Item to push onto queue
 */
template< class Item > void affineQueue<Item>::push(const Item& item)
  throw (affineExcept) {
  
  //Note this also checks for zero cap
  if (nelem >= cap) 
    throw affineExcept("affineQueue","push",
		       "No room to add another element",1);  
  data[tail] = item;
  tail += 1;
  if (tail == cap) tail = 0;
  nelem += 1;
}

/*!
  \returns Item pulled off of queue
 */
template< class Item > Item affineQueue<Item>::pop()
  throw (affineExcept) {
  if (nelem == 0) 
    throw affineExcept("affineQueue","pop",
		       "Queue is empty",1);
  unsigned int oldhead = head;
  head += 1;
  if (head == cap) head = 0; //Wrap
  nelem -= 1;
  return data[oldhead];
}

#endif
