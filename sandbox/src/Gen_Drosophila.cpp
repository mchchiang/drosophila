/* Gen_Drosophila.cpp
 * This is a code that generaetes a chromosome with heterchromatin
 * (determined by the signal of H3K9me) and laminar contact (LaminB signal)
 * info encoded
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <random>
#include <cstdlib>
#include <cmath>
#include <map>
#include "Bead.hpp"
#include "Polymer.hpp"
#include "LAMMPS.hpp"
#include "DataManager.hpp"

using std::cout;
using std::cin;
using std::endl;
using std::map;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;

double getBoxLength(int numOfBeads, double density);
int getBeadType(vector<double>& fracOfStates);
void shapePolymer(shared_ptr<Polymer> polymer, double maxRadius, 
		  double x0, double y0, double z0);
double distance(double x, double y, double z);

int main(int argc, char * argv[]){
  if (argc != 8){
    cout << "Usage: Gen_Drosophila [chromoFile] [hmmFile] [ploidy] "
	 << "[fraction] [buffer] [outfile] [outMapFile]" << endl;
    return 1;
  }

  int argi {};
  string chromoFile (argv[++argi]);
  string hmmFile (argv[++argi]);
  int ploidy {stoi(string(argv[++argi]), nullptr, 10)};
  double fraction {stod(string(argv[++argi]), nullptr)};
  double buffer {stod(string(argv[++argi]), nullptr)};
  string outFile (argv[++argi]);
  string outMapFile (argv[++argi]);

  bool isHaploid {true};
  if (ploidy > 1){
    isHaploid = false;
  }

  const int haploidNum {6}; // For drosophila (2L/2R and 3L/3R)
  const int numOfFibres {haploidNum*ploidy};

  vector<string> chrNumToName = {"chr2L","chr2R","chr3L",
				 "chr3R","chr4","chrX"};
  map<string,int> chrNameToNum;
  for (size_t i {}; i < chrNumToName.size(); i++){
    chrNameToNum[chrNumToName[i]] = i;
  }

  // Determine the length of the fibre representing each chromosome
  int fibreLength [numOfFibres];
  const int bpPerBead {3000}; // Coarse-graining (resolution)

  // Read the number of bp in each chromosome
  string header {};
  ifstream chromoReader;
  chromoReader.open(chromoFile);
  
  if (!chromoReader){
    cout << "Unable to read chromo file \""
         << chromoFile << "\"\n"
         << "Aborting reading process" << endl;
    return 1;
  }

  // Skip the header of the file
  getline(chromoReader, header);
  getline(chromoReader, header);

  // Read the number of bp in each chromosome
  string index {};
  double length {};
  for (int i {}; i < numOfFibres; i++){
    chromoReader >> index >> length;
    cout << index << " " << length << endl;
    fibreLength[chrNameToNum[index]] = 
      static_cast<int>(ceil(length / bpPerBead));
  }
  chromoReader.close();

  // Read the chromatin state file
  cout << "Reading hmm file ..." << endl;
  ifstream hmmReader;
  hmmReader.open(hmmFile);

  // Skip the header of the file
  getline(hmmReader, header);

  int chrNum {}, state {};
  long start {}, end {};
  double frac {};

  vector<double> zero (9, 0.0);
  vector<vector<vector<double> > >* fracOfStates
  {new vector<vector<vector<double> > >()};
  for (int i {}; i < numOfFibres; i++){
    fracOfStates->push_back(vector<vector<double> >());
    for (int j {}; j < fibreLength[i]; j++){
      (*fracOfStates)[i].push_back(zero);
    }
  }
  while (!hmmReader.eof()){
    hmmReader >> index >> start >> end >> state;
    int startBead = start / bpPerBead;
    int endBead = end / bpPerBead;
    state--; // Make sure state runs from 0 to 8 (not 1 to 9)
    
    // Skip unrelated chromosomal data
    try {
      chrNum = chrNameToNum.at(index);
    } catch (const std::out_of_range& oor){
      continue;
    }
    
    // For content contained within the same bead
    if (startBead == endBead){
      frac = static_cast<double>(end-start)/bpPerBead;
      (*fracOfStates)[chrNum][startBead][state] += frac;
      if (!isHaploid){
	(*fracOfStates)[chrNum+haploidNum][startBead][state] += frac;
      }
    } else { // For content spread over multiple beads
      // Start bead content
      frac = 1.0-(static_cast<double>(start)/bpPerBead-
		  static_cast<double>(startBead));
      (*fracOfStates)[chrNum][startBead][state] += frac;
      if (!isHaploid){
	(*fracOfStates)[chrNum+haploidNum][startBead][state] += frac;
      }
      // End bead content
      frac = (static_cast<double>(end)/bpPerBead-
	      static_cast<double>(endBead));
      (*fracOfStates)[chrNum][endBead][state] += frac;
      if (!isHaploid){
	(*fracOfStates)[chrNum+haploidNum][endBead][state] += frac;
      }
      // Other beads in between are completely coded by the content
      for (int i {startBead+1}; i < endBead; i++){
	(*fracOfStates)[chrNum][i][state] = 1.0;
	if (!isHaploid){
	  (*fracOfStates)[chrNum+haploidNum][i][state] = 1.0;
	}
      }
    }
  }
  hmmReader.close();
  cout << "Finish reading hmm file." << endl;

  // Determine the simulation box size from the total number of beads
  int numOfBeads {};
  for (int i {}; i < numOfFibres; i++){
    numOfBeads += fibreLength[i];
  }
  numOfBeads++; // Add in the tethering bead
  const double L {getBoxLength(numOfBeads, fraction)};
  const double lx {L-buffer*2.0}, ly {L-buffer*2.0}, lz {L-buffer*2.0};

  // Generate the chromosome polymers
  cout << "Generating chromosome polymers ..." << endl;
  shared_ptr<LAMMPS> lammps = make_shared<LAMMPS>(L, L, L);
  shared_ptr<Polymer> polymer {};
  shared_ptr<Bead> bead {};

  
  lammps->setTypesOfBeads(3);
  lammps->setTypesOfBonds(2);
  lammps->setTypesOfAngles(1);
  
  double radius {0.15*lx};
  vector<double> xshift {-0.2,-0.3,0.3,0.2,0.0,0.0};
  vector<double> yshift {-0.3,0.0,0.0,-0.3,0.0,0.3};
  vector<double> zshift {0.0,0.0,0.0,0.0,-0.3,0.0};

  for (int i {}; i < numOfFibres; i++){
    cout << "Generating " << chrNumToName[i] << endl;
    polymer = lammps->createPolymer(i,fibreLength[i]);
    shapePolymer(polymer,radius,xshift[i]*lx,
		 yshift[i]*ly,zshift[i]*lz);
    for (int j {}; j < fibreLength[i]; j++){
      bead = polymer->getBead(j);
      bead->setLabel(i);
      bead->setType(getBeadType((*fracOfStates)[i][j]));
    }
    cout << "Done generating " << chrNumToName[i] << endl;
  }
  
  
  // Add a tethering bead to create a Rabl-like configuration
  shared_ptr<Bead> tether = lammps->createBead(numOfFibres);
  
  // Set the tethering bead to locate at the bottom of the box
  tether->setPosition(2,-0.45*lz);
  tether->setLabel(numOfFibres+1);
  tether->setType(3);

  // Find the ending bead of each polymer and bond it with the tethering bead
  // End of chr2L
  polymer = lammps->getPolymer(chrNameToNum["chr2L"]);
  bead = polymer->getBead(polymer->getNumOfBeads()-1);
  tether->addBondWith(2,bead);
  // Start of chr2R
  polymer = lammps->getPolymer(chrNameToNum["chr2R"]);
  bead = polymer->getBead(0);
  tether->addBondWith(2,bead);
  // End of chr3L
  polymer = lammps->getPolymer(chrNameToNum["chr3L"]);
  bead = polymer->getBead(polymer->getNumOfBeads()-1);
  tether->addBondWith(2,bead);
  // Start of chr3R
  polymer = lammps->getPolymer(chrNameToNum["chr3R"]);
  bead = polymer->getBead(0);
  tether->addBondWith(2,bead);
  // End of chr4
  polymer = lammps->getPolymer(chrNameToNum["chr4"]);
  bead = polymer->getBead(polymer->getNumOfBeads()-1);
  tether->addBondWith(2,bead);
  // End of chrX
  polymer = lammps->getPolymer(chrNameToNum["chrX"]);
  bead = polymer->getBead(polymer->getNumOfBeads()-1);
  tether->addBondWith(2,bead);
  
  // Write the input file
  lammps->exportData(outFile, outMapFile);
  
  // Delete resources
  delete fracOfStates;
}

double getBoxLength(int numOfBeads, double density){
  return pow(numOfBeads*M_PI/6.0/density,1.0/3.0);
}


int getBeadType(vector<double>& fracOfStates){
  double activeFrac {}; // State 1-5
  double inactiveFrac {}; // State 6-9
  for (int i {}; i < 5; i++){
    activeFrac += fracOfStates[i];
  }
  for (int i {5}; i < 9; i++){
    inactiveFrac += fracOfStates[i];
  }
  if (inactiveFrac >= activeFrac){
    return 2;
  } else {
    return 1;
  }
}


void shapePolymer(shared_ptr<Polymer> polymer, double maxRadius, 
		     double x0, double y0, double z0){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> randDouble(0.0,1.0);
  
  bool ok {true};
  const double pi {M_PI};
  double x, y, z;
  double r, costheta, sintheta, phi;
  shared_ptr<Bead> bead, previous, b;
  double minSep {0.5}; // Minimum separation between any beads
  double sep {0.8}; // Separation between consecutive beads

  int nBeads {polymer->getNumOfBeads()};

  // The first bead is at the origin
  previous = polymer->getBead(0);
  previous->setPosition(0, 0.0);
  previous->setPosition(1, 0.0);
  previous->setPosition(2, 0.0);

  for (int i {1}; i < nBeads; i++){
    bead = polymer->getBead(i);
    do {
      ok = true;
      r = randDouble(mt);
      costheta = 1.0-2.0*r;
      sintheta = sqrt(1.0-costheta*costheta);
      r = randDouble(mt);
      phi = 2.0*pi*r;
      x = previous->getPosition(0) + sep * sintheta * cos(phi);
      y = previous->getPosition(1) + sep *sintheta * sin(phi);
      z = previous->getPosition(2) + sep * costheta;
      
      // Check if the bead's position exceeds the boundary
      if (distance(x,y,z) >= maxRadius){
	ok = false;
	continue;
      }

      // Check if the bead overlaps with previous beads
      for (int j {}; j < i; j++){
	b = polymer->getBead(j);
	if (distance(x-b->getPosition(0),y-b->getPosition(1),
		     z-b->getPosition(2)) < minSep){
	  ok = false;
	  continue;
	}
      }
    } while (!ok);
    
    bead->setPosition(0,x);
    bead->setPosition(1,y);
    bead->setPosition(2,z);
    
    previous = bead;
  }

  // Get polymer's centre of mass 
  vector<double> cm {polymer->getCentreOfMass(0,0,0)};

  // Shift the beads such that the centre of mass of
  // the polymer is at (x0,y0,z0)
  for (int i {}; i < nBeads; i++){
    bead = polymer->getBead(i);
    bead->setPosition(0, bead->getPosition(0)-cm[0]+x0);
    bead->setPosition(1, bead->getPosition(1)-cm[1]+y0);
    bead->setPosition(2, bead->getPosition(2)-cm[2]+z0);
  }
}

double distance(double x, double y, double z){
  return sqrt(x*x+y*y+z*z);
}
