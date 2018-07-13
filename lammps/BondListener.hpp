/* BondListener.hpp
 *
 * An interface for handling bond events
 */

#ifndef BONDLISTENER_HPP
#define BONDLISTENER_HPP

#include <memory>

using std::shared_ptr;

class Bond; 

class BondListener {

public:
  virtual void bondCreated(const shared_ptr<Bond>& bond);
  virtual void bondRemoved(const shared_ptr<Bond>& bond);
  virtual void bondTypeChanged(const shared_ptr<Bond>& bond, 
			       int oldType, int newType) = 0;

};

#endif
