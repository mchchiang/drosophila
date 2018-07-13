/* Angle.hpp
 *
 * A class for storing the info of three beads that are angle-bonded
 */

#ifndef ANGLE_HPP
#define ANGLE_HPP

#include <memory>

using std::weak_ptr;
using std::shared_ptr;

class Bead; // Forward declaration to break cyclic references

class Angle : public std::enable_shared_from_this<Angle> {
	
private:
  const int numOfBeads {3};
  int type;
  weak_ptr<Bead> beads[3];
  
public:
  // Constructor
  Angle(int type, shared_ptr<Bead> bead1, 
	shared_ptr<Bead> bead2, shared_ptr<Bead> bead3);
  
  // Accessor methods
  shared_ptr<Bead> getBead(int id);
  void setType(int type);
  int getType();
};

#endif
