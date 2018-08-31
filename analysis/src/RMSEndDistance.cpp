/* RMSEndDistance.cpp
 * A program that reads the position file and computes
 * the root mean square end-to-end distance of the polymer
 * as a function of lienar separation
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "PositionReader.hpp"
#include "LAMMPSMapReader.hpp"

using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::vector;
using std::map;

double distance2(const double (& v1)[3], const double (& v2)[3]);

int main(int argc, char* argv[]){
  
  if (argc != 11){
    cout << "Usage: [numOfBeads] [lx] [ly] [lz] [startTime] [endTime] "
	 << "[timeInc] [posFile] [mapFile] [outFile]" << endl;
    return 1;
  }

  int argi {};
  const int numOfBeads {stoi(string(argv[++argi]), nullptr, 10)};
  const double lx {stod(string(argv[++argi]), nullptr)};
  const double ly {stod(string(argv[++argi]), nullptr)};
  const double lz {stod(string(argv[++argi]), nullptr)};
  const long startTime {stol(string(argv[++argi]), nullptr, 10)};
  const long endTime {stol(string(argv[++argi]), nullptr, 10)};
  const long timeInc {stol(string(argv[++argi]), nullptr, 10)};
  const string posFile (argv[++argi]);
  const string mapFile (argv[++argi]); // For mapping different chromosomes
  const string outFile (argv[++argi]);

  PositionReader reader;
  reader.open(posFile, numOfBeads, lx, ly, lz, timeInc);
  if (!reader.isOpen()){
    cout << "Problem with reading the position file!" << endl;
    return 1;
  }
  
  LAMMPSMapReader mapReader;
  mapReader.open(mapFile);
  if (!reader.isOpen()){
    cout << "Problem with reading the map file!" << endl;
    return 1;
  }
  map<int,vector<int> > polymerMap {mapReader.getPolymerMap()};
  size_t numOfPolymers {polymerMap.size()};

  // For storing the distance data (for each chromosome)
  vector<vector<double> > distances;
  for (size_t i {}; i < numOfPolymers; i++){
    vector<double> zeros (polymerMap[i][1]-polymerMap[i][0]+1,0.0);
    distances.push_back(zeros);
  }
  
  long time;
  double r1 [3], r2 [3];
  int startBead, endBead, s, d;
  while (reader.nextFrame()){
    time = reader.getTime();
    cout << "Reading frame: " << time << endl;
    if (time >= startTime && time <= endTime){
      for (size_t i {}; i < numOfPolymers; i++){
	cout << "Doing polymer " << i << endl;
	startBead = polymerMap[i][0];
	endBead = polymerMap[i][1];
	for (int j {startBead+1}; j <= endBead; j++){
	  for (int k {startBead}; k < j; k++){
	    s = abs(k-j);
	    for (int l {}; l < 3; l++){
	      r1[l] = reader.getUnwrappedPosition(j,l);
	      r2[l] = reader.getUnwrappedPosition(k,l);
	    }
	    d = sqrt(distance2(r1,r2));
	    distances[i][s] += d;
	  }
	}
      }
    } 
  }

  // Normalise
  long numOfFrames {(endTime-startTime)/timeInc};
  for (size_t i {}; i < numOfPolymers; i++){
    size_t polymerSize {distances[i].size()};
    for (size_t j {}; j < polymerSize; j++){
      distances[i][j] /= static_cast<double>((polymerSize-j)*numOfFrames);
    }
  }

  // Output the results
  ofstream writer;
  for (size_t i {}; i < numOfPolymers; i++){
    ostringstream oss;
    oss << outFile << "." << i;
    writer.open(oss.str());
    if (!writer){
      cout << "Problem with opening output file!" << endl;
      return 1;
    }
    size_t polymerSize {distances[i].size()};
    for (size_t j {}; j < polymerSize; j++){
      writer << j << " " << distances[i][j] << endl;
    }
    writer.close();
  }
}

double distance2(const double (& v1)[3], const double (& v2)[3]){
  double dx, dy, dz;
  dx = v1[0] - v2[0];
  dy = v1[1] - v2[1];
  dz = v1[2] - v2[2];
  return dx*dx+dy*dy+dz*dz;
}
