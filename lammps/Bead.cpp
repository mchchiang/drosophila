// Bead.cpp

#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include "Bead.hpp"
#include "Bond.hpp"
#include "Angle.hpp"

using std::cout;
using std::endl;
using std::pair;
using std::weak_ptr;
using std::shared_ptr;
using std::make_shared;
using BondIterator = vector<shared_ptr<Bond> >::iterator;
using AngleIterator = vector<shared_ptr<Angle> >::iterator;

// Constructors
Bead::Bead(double x, double y, double z,
		   double vx, double vy, double vz,
		   int nx, int ny, int nz, int t, int l) :
  position {x, y, z}, velocity {vx, vy, vz},
  boundaryCount {nx, ny, nz}, type {t}, label {l} {
	bondList.reserve(0);
	angleList.reserve(0);
}

Bead::Bead(double x, double y, double z) :
  Bead {x, y, z, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0} {}

Bead::Bead() :
  Bead {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0} {}


// Accessor methods
void Bead::setPosition(int dim, double value){
  position[dim] = value;
}

double Bead::getPosition(int dim){
  return position[dim];
}

void Bead::setVelocity(int dim, double value){
  velocity[dim] = value;
}

double Bead::getVelocity(int dim){
  return velocity[dim];
}

void Bead::setBoundaryCount(int dim, int count){
  boundaryCount[dim] = count;
}

int Bead::getBoundaryCount(int dim){
  return boundaryCount[dim];
}

void Bead::setType(int t){
  type = t;
}

int Bead::getType(){
  return type;
}

void Bead::setLabel(int l){
  label = l;
}

int Bead::getLabel(){
  return label;
}

int Bead::getNumOfBonds(){
  return bondList.size();
}

int Bead::getNumOfAngles(){
  return angleList.size();
}

shared_ptr<Bond> Bead::getBondWith(shared_ptr<Bead> bead){
  for (BondIterator it = bondList.begin(); it != bondList.end(); it++){
    shared_ptr<Bond> bond = *it;
    if (bond->getBead(0) == bead || bond->getBead(1) == bead){
      return bond;
    }
  }
  return nullptr;
}

shared_ptr<Angle> Bead::getAngleWith(shared_ptr<Bead> bead1,
				     shared_ptr<Bead> bead2){
  for (AngleIterator it = angleList.begin(); it != angleList.end(); it++){
    shared_ptr<Angle> angle = *it;
    shared_ptr<Bead> b1 = angle->getBead(0);
    shared_ptr<Bead> b2 = angle->getBead(1);
    shared_ptr<Bead> b3 = angle->getBead(2);
    if ((bead1 == b1 && (bead2 == b2 || bead2 == b3)) ||
	(bead1 == b2 && (bead2 == b1 || bead2 == b3)) ||
	(bead1 == b3 && (bead2 == b1 || bead2 == b2))){
      return angle;
    }
  }
  return nullptr;
}

vector<shared_ptr<Bond> >& Bead::getBonds(){
  return bondList;
}

vector<shared_ptr<Angle> >& Bead::getAngles(){
  return angleList;
}

void Bead::addBond(shared_ptr<Bond> bond){
  bondList.push_back(bond);
}

void Bead::addBondWith(int t, shared_ptr<Bead> bead){
  shared_ptr<Bond> bond {make_shared<Bond>(t, shared_from_this(),bead)};
  bead->addBond(bond);
  addBond(bond);
}

void Bead::removeBond(shared_ptr<Bond> bond){
  BondIterator it = std::find(bondList.begin(), bondList.end(), bond);
  if (it != bondList.end())
    bondList.erase(it);
}

void Bead::removeBondWith(shared_ptr<Bead> bead){
  for (BondIterator it = bondList.begin(); it != bondList.end(); it++){
    shared_ptr<Bond> bond = *it;
    if (bond->getBead(0) == bead || bond->getBead(1) == bead){
      bead->removeBond(bond);
      bondList.erase(it);
      break;
    }
  }
}

void Bead::addAngle(shared_ptr<Angle> angle){
  angleList.push_back(angle);
}

void Bead::addAngleWith(int t, 
			shared_ptr<Bead> bead1, shared_ptr<Bead> bead2){
  shared_ptr<Angle> angle {
    make_shared<Angle>(t, shared_from_this(), bead1, bead2)};
  bead1->addAngle(angle);
  bead2->addAngle(angle);
  addAngle(angle);
}

void Bead::removeAngle(shared_ptr<Angle> angle){
  AngleIterator it = std::find(angleList.begin(), angleList.end(), angle);
  if (it != angleList.end())
    angleList.erase(it);
}


void Bead::removeAngleWith(shared_ptr<Bead> bead1, shared_ptr<Bead> bead2){
  for (AngleIterator it = angleList.begin(); it != angleList.end(); it++){
    shared_ptr<Angle> angle = *it;
    shared_ptr<Bead> b1 = angle->getBead(0);
    shared_ptr<Bead> b2 = angle->getBead(1);
    shared_ptr<Bead> b3 = angle->getBead(2);
    if ((bead1 == b1 && (bead2 == b2 || bead2 == b3)) ||
	(bead1 == b2 && (bead2 == b1 || bead2 == b3)) ||
	(bead1 == b3 && (bead2 == b1 || bead2 == b2))){
      bead1->removeAngle(angle);
      bead2->removeAngle(angle);
      angleList.erase(it);
      break;
    }
  }
}

void Bead::removeAllBonds(){
  for (auto const& bond : bondList){
    shared_ptr<Bead> bead1 = bond->getBead(0);
    shared_ptr<Bead> bead2 = bond->getBead(1);
    
    if (bead1 != shared_from_this())
      bead1->removeBond(bond);
    else
      bead2->removeBond(bond);
  }
  bondList.clear();
}

void Bead::removeAllAngles(){
  for (auto const& angle : angleList){
    shared_ptr<Bead> bead1 = angle->getBead(0);
    shared_ptr<Bead> bead2 = angle->getBead(1);
    shared_ptr<Bead> bead3 = angle->getBead(2);
    
    if (bead1 == shared_from_this()){
      bead2->removeAngle(angle);
      bead3->removeAngle(angle);
    } else if (bead2 == shared_from_this()){
      bead1->removeAngle(angle);
      bead3->removeAngle(angle);
    } else {
      bead1->removeAngle(angle);
      bead2->removeAngle(angle);
    }
  }
  angleList.clear();
}
