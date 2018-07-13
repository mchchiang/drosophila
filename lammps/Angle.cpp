// Angle.cpp

#include <memory>
#include "Bead.hpp"
#include "Angle.hpp"

using std::weak_ptr;
using std::shared_ptr;

// Constructor
Angle::Angle(int t, shared_ptr<Bead> b1, 
	     shared_ptr<Bead> b2, shared_ptr<Bead> b3) :
  type {t}, beads {weak_ptr<Bead>(b1), weak_ptr<Bead>(b2), 
		weak_ptr<Bead>(b3)} {}

// Accessor methods
shared_ptr<Bead> Angle::getBead(int id){
  return beads[id].lock();
}

void Angle::setType(int t){
  type = t;
}

int Angle::getType(){
  return type;
}

