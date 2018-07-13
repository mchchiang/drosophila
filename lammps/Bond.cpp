// Bond.cpp

#include <memory>
#include "Bead.hpp"
#include "Bond.hpp"

using std::weak_ptr;
using std::shared_ptr;

// Constructor
Bond::Bond(int t, shared_ptr<Bead> b1, shared_ptr<Bead> b2) :
  type {t}, beads {weak_ptr<Bead>(b1), weak_ptr<Bead>(b2)} {}

// Accessor methods
shared_ptr<Bead> Bond::getBead(int id){
  return beads[id].lock();
}

void Bond::setType(int t){
  type = t;
}

int Bond::getType(){
  return type;
}

