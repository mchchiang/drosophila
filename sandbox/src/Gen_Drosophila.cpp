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
  
  // Make each chromosome arm in a smaller box first
  double lxmin {lx/2.0};
  double lymin {ly/2.0};
  double lzmin {lz/2.0};
  vector<double> cm; // Store the centre of mass of a polymer

  /*
   * Strategy to place the chromosome arms:
   * Divide simulation box into eight compartments (each of same size)
   * Randomly place the chromosome arms in these compartments 
   */
  const int numOfBoxes {8};
  const double xs {lx/4.0}, ys {ly/4.0}, zs {lz/4.0};
  vector<double> xshift {-xs,xs,xs,-xs,-xs,xs,xs,-xs};
  vector<double> yshift {-ys,-ys,ys,ys,-ys,-ys,ys,ys};
  vector<double> zshift {-zs,-zs,-zs,-zs,zs,zs,zs,zs};

  std::random_device rd;
  std::mt19937 mt(rd());

  // Create a list of the indices of chromosomes to be picked (i.e. 0,1,...,N)
  vector<int> chromoToPick (numOfFibres);
  std::iota(chromoToPick.begin(), chromoToPick.end(), 0);
  
  // Create a list of the indices of boxes to be picked (i.e 0,1,...,8)
  vector<int> boxToPick (numOfBoxes);
  std::iota(boxToPick.begin(), boxToPick.end(), 0);

  int chr, chrIndex, box, boxIndex;
  while (chromoToPick.size() > 0){
    std::uniform_int_distribution<int> pickChromo(0,chromoToPick.size()-1);
    chrIndex = pickChromo(mt);
    chr = chromoToPick[chrIndex];
    cout << "Generating chromosome " << chrNumToName[chr] << endl;
    polymer = lammps->createRandomWalkPolymer(chr, fibreLength[chr], 0, 
					      0.0, 0.0, 0.0, 
					      lxmin, lymin, lzmin);

    // Pick compartment
    std::uniform_int_distribution<int> pickBox(0,boxToPick.size()-1);
    boxIndex = pickBox(mt);
    box = boxToPick[boxIndex];
    
    // Shift the polymer to a particular place
    cm = polymer->getCentreOfMass(lxmin,lymin,lzmin);
    vector<double> shift {xshift[box]-cm[0], yshift[box]-cm[1], 
	zshift[box]-cm[2]};

    for (int i {}; i < fibreLength[chr]; i++){
      bead = polymer->getBead(i);
      for (int j {}; j < 3; j++){
	bead->setPosition(j, bead->getPosition(j)+shift[j]);
      }
      bead->setLabel(chr);
      bead->setType(getBeadType((*fracOfStates)[chr][i]));
    }
    cout << "Done generating " << chrNumToName[chr] << endl;
    chromoToPick.erase(chromoToPick.begin()+chrIndex);
    boxToPick.erase(boxToPick.begin()+boxIndex);
  }
  
  // Add a tethering bead to create a Rabl-like configuration
  shared_ptr<Bead> tether = lammps->createBead(numOfFibres+1);
  
  // Set the tethering bead to locate at the bottom of the box
  tether->setPosition(2,-zs);
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
