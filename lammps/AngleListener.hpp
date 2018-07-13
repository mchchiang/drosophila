/* AngleListener.hpp
 *
 * An interface for handling bond events
 */

#ifndef ANGLELISTENER_HPP
#define ANGLELISTENER_HPP

#include <memory>

using std::shared_ptr;

class Angle; 

class AngleListener {

public:
  virtual void angleCreated(const shared_ptr<Angle>& angle);
  virtual void angleRemoved(const shared_ptr<Angle>& angle);
  virtual void angleTypeChanged(const shared_ptr<Angle>& angle, 
			       int oldType, int newType) = 0;

};

#endif
