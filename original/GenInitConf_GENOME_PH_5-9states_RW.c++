//This program "COLORS" all the Chromosomes of Drosophila in such a way that
//there are strong/weak binding sites for PolyComb (PH) (see Ajaz)
//AND it takes into account the 9-states models (see Kharchenko, Nature 480 2011) (R2)
// It colors Euchromatin Strong (TSS), weak (transcribed), heterochrom/blue chromatin(HP1) and deserts/black chromatin (H1)
//
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<sstream>
#include<vector>
#include<ctime>
using namespace std;


int InitialConfg=0;
//0 == RANDOM WALK - RW; 1 == Asymmetric RW (stretched) - S; 2 == CRUMPLED inside SPHERE - C; 3 == Fractal globule - FG; 5 == mitotic - M

int Nmax= 30000;
int Nfibermax= 1;
int Ntypemax= 4;
int Nbondmax= Nmax;
int Nanglemax= Nmax;
int nfiber= 1;
int bond1= 2;
int ang1= 3;
int nbeadsmax= 30000;
int linesRed= 226;
int linesBlack= 6183;

const int BasePairsCh2L= 22985362;
const int BasePairsCh2R= 21146687;
const int BasePairsCh3L= 24543600;
const int BasePairsCh3R= 27905053;
const int BasePairsCh4=  1351940;
const int BasePairsChX=  22422294;
const int Genome=BasePairsCh2R+BasePairsCh2L+BasePairsCh3R+BasePairsCh3L+BasePairsChX+BasePairsCh4;

const int LengthChromosome[6]={22985362,21146687,24543600,27905053,1351940,22422294};
const int LChr[6]={7662,7049,8181,9302,451,7475};
const int GenomeCG=7662+7049+8181+9302+451+7475;
// total = 40119

string Chromosome;
const char *ChrNAME[6]={"chr2L","chr2R","chr3L","chr3R","chr4","chrX"};
//
const int MAXbeads=10000;
int ColorChr[6][MAXbeads];
//
int PCprotein=0;
int HeterChrProtein=0;
int EuChrProtein=0;
int H1Protein=0;

int nproteins=0;
//
//STATES
//
double Het[6][MAXbeads];
double PC[6][MAXbeads];
double Eu[6][MAXbeads];
double Desert[6][MAXbeads];
double AbundancePC[6][MAXbeads];
double AbundanceHet[6][MAXbeads];
double AbundanceEu[6][MAXbeads];
double AbundanceDesert[6][MAXbeads];
//
double center[3]={0.0,0.0,0.0};
double distance(double v1[3], double v2[3]);
void InitialiseColors();
void Paint(int i);
//
double position[6][MAXbeads][3];


int main(int argc, char* argv[]){

  long int nbead;
  int npoly,ntype,nbond,nangle,type;
  int nbeads,ntypepol;

  cout << "Input argv: [1] red (polym-PH) sites; [2] black (non-polym-PH) sites ; [3] 9-states model; [4] Type (RW,S,C,FG,M); [5] Replicas; [6] NPHp [7] NHP1p; [8] NEUp; [9] Config (1=RW; 2=STR; 3=CRUMPL; 4=FG; 5=M)"<<endl;
  ntypepol=19;

  int coarsegrain;
  coarsegrain=3000;

  int datatype=0;
  datatype=1;


  PCprotein=atoi(argv[6]);
  cout << "PHp " << PCprotein<<endl;

  HeterChrProtein=atoi(argv[7]);
  cout << "HP1p " << HeterChrProtein<<endl;

  EuChrProtein=atoi(argv[8]);
  cout << "EUp " << EuChrProtein<<endl;

  H1Protein=HeterChrProtein/2;
  cout << "H1p " << H1Protein<<endl;

  nproteins=PCprotein+HeterChrProtein+EuChrProtein+H1Protein;
  cout << GenomeCG<< " + " << nproteins<< endl;


  double density=0.18;   // Monomer density
  double L=pow((GenomeCG+nproteins)*1.0/density,1.0/3.0);
  double Lx=L;
  double Ly=L;
  double Lz=L;

  string Type=argv[4];
  cout << "Type " << Type<<endl;
  int RepsMax=atoi(argv[5]);
  cout <<  "Replicas " <<RepsMax<<endl;
  InitialConfg=atoi(argv[9]);
  cout <<  "Init confg " << InitialConfg<<endl;

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  for(int reps=1;reps<=RepsMax;reps++){
    srand(time(NULL));
    InitialiseColors();
    //
    ofstream writeW;
    stringstream writeFileW;
    writeFileW << "LammpsInput.GENOME.Wild"<<".N"<<GenomeCG<<Type<<reps<<".PHp"<<PCprotein<<"HP1p"<<HeterChrProtein<<"EUp"<<EuChrProtein;
    writeW.open(writeFileW.str().c_str());

    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    //start reading and go through chromosomes
    for(int index_ch=0; index_ch<6; index_ch++){   // Loop over chromosomes

      for(int cc=0; cc<LChr[index_ch]; cc++)ColorChr[index_ch][cc]=1;

      ifstream readRed;
      readRed.open(argv[1]);

      ifstream readBlack;
      readBlack.open(argv[2]);

      ifstream read9;
      read9.open(argv[3]);

      /////////////////////////////////////////////////////////////
      ////////		 POLYCOMB            //////////////////////
      /////////////////////////////////////////////////////////////
      string dummy;

      getline(readRed,dummy);
      //cout << dummy <<endl;

      double thrR=2.0;
      double thrB=1.0;
      long int start,end;
      double Awild;
      int nb=0;
      while(!readRed.eof()){
        readRed >> Chromosome >> start >> end;
        if(Chromosome.compare(ChrNAME[index_ch])==0){
          int nstart = floor(start*1.0/coarsegrain);
          int nend = ceil(end*1.0/coarsegrain);
	  for(int nnb=nstart; nnb<=nend; nnb++){
	    PC[index_ch][nnb]=1;
	    AbundancePC[index_ch][nnb]=1;
	    if(nnb==nstart) AbundancePC[index_ch][nnb]= fabs(1.0-(start*1.0/coarsegrain-nstart));
	    if(nnb==nend) AbundancePC[index_ch][nnb]= fabs(1.0-(nend-end*1.0/coarsegrain));
	  }
        }
      }
      readRed.close();

      getline(readBlack,dummy);
      while(!readBlack.eof()){
        readBlack>> Chromosome >> start >> end;
	if(Chromosome.compare(ChrNAME[index_ch])==0){
          int nstart = floor(start*1.0/coarsegrain);
          int nend = ceil(end*1.0/coarsegrain);
          for(int nnb=nstart; nnb<=nend; nnb++){ 
            PC[index_ch][nnb]=2;
            AbundancePC[index_ch][nnb]=1;
            if(nnb==nstart) AbundancePC[index_ch][nnb]= fabs(1.0-(start*1.0/coarsegrain-nstart));
            if(nnb==nend) AbundancePC[index_ch][nnb]= fabs(1.0-(nend-end*1.0/coarsegrain));
          }
        }
      }
      readBlack.close();

      ///////////////////////////////////////////////////////////////////
      ////////      9-STATES MODEL  --  EUCHROMATIN        //////////////
      ///////////////////////////////////////////////////////////////////
      for(int n=0;n<LChr[index_ch];n++){
        Het[index_ch][n]=Eu[index_ch][n]=Desert[index_ch][n]=0;
        AbundanceHet[index_ch][n]=0;
        AbundanceEu[index_ch][n]=0;
        AbundanceDesert[index_ch][n]=0;
      }

      int nread=0;
      int typeChr;
      int Euchr1a=1;
      int Euchr1b=3;
      int Euchr2a=2;
      int Euchr2b=4;
      while(true){
        read9 >> Chromosome >> start >> end >> typeChr;
        cout << Chromosome <<" " <<   start<< " " << end << " " <<typeChr<<endl;
        if(read9.eof())break;
    
        if(Chromosome.compare(ChrNAME[index_ch])==0){
          if(typeChr==Euchr1a or typeChr==Euchr1b){
            int nstart = floor(start*1.0/coarsegrain);
            int nend = ceil(end*1.0/coarsegrain);
            for(int nnb=nstart; nnb<nend; nnb++){
              Eu[index_ch][nnb]=1;
              AbundanceEu[index_ch][nnb]=1;
              if(nnb==nstart) AbundanceEu[index_ch][nnb]= fabs(1.0-(start*1.0/coarsegrain-nstart));
              if(nnb==nend-1) AbundanceEu[index_ch][nnb]= fabs(1.0-(nend-end*1.0/coarsegrain));
            }
          }
          if(typeChr==Euchr2a or typeChr==Euchr2b){
            int nstart = floor(start*1.0/coarsegrain);
            int nend = ceil(end*1.0/coarsegrain);
            for(int nnb=nstart; nnb<nend; nnb++){
              Eu[index_ch][nnb]=2;
              AbundanceEu[index_ch][nnb]=1;
              if(nnb==nstart) AbundanceEu[index_ch][nnb]= fabs(1.0-(start*1.0/coarsegrain-nstart));
              if(nnb==nend-1) AbundanceEu[index_ch][nnb]= fabs(1.0-(nend-end*1.0/coarsegrain));
            }
          }	
          if(typeChr==7 || typeChr==8){ //7-8 heterochrom
            int nstart = floor(start*1.0/coarsegrain);
            int nend = ceil(end*1.0/coarsegrain);
            for(int nnb=nstart; nnb<nend; nnb++){
              Het[index_ch][nnb]=1;
              AbundanceHet[index_ch][nnb]=1;
              if(nnb==nstart) AbundanceHet[index_ch][nnb]= fabs(1.0-(start*1.0/coarsegrain-nstart));
              if(nnb==nend-1) AbundanceHet[index_ch][nnb]= fabs(1.0-(nend-end*1.0/coarsegrain));
            }
          }
          if(typeChr==9){ //9 low transcription
            int nstart = floor(start*1.0/coarsegrain);
            int nend = ceil(end*1.0/coarsegrain);
            for(int nnb=nstart; nnb<nend; nnb++){
              Desert[index_ch][nnb]=1;
              AbundanceDesert[index_ch][nnb]=1;
              if(nnb==nstart) AbundanceDesert[index_ch][nnb]= fabs(1.0-(start*1.0/coarsegrain-nstart));
              if(nnb==nend-1) AbundanceDesert[index_ch][nnb]= fabs(1.0-(nend-end*1.0/coarsegrain));
            }
          }
        }

        nread++;
      }
      cout << "I'VE FINISHED PAINTING..." <<endl;

      //////////////////////////////////////////////////////////////
      ////////          PICK COLOR            //////////////////////
      //////////////////////////////////////////////////////////////
      Paint(index_ch);

    } //close loop over chromosomes


    ////////////////////////////////////////////////////////////////
    ////////           MAKE POLYMER         ////////////////////////
    ////////////////////////////////////////////////////////////////
    cout << "MAKE THE POLYMERS ..." <<endl;

    int tether=1;

    writeW<< "LAMMPS data file from restart file: timestep = 0,\t procs = 1"<<endl;
    writeW << GenomeCG+nproteins+tether  << " atoms "<<endl;
    writeW << GenomeCG-6+6<< " bonds "<<endl; // 6 tethers
    writeW << GenomeCG-12 << " angles "<<endl;
    writeW << "\n";
    writeW << ntypepol+1 << " atom types "<<endl;
    writeW << 2 << " bond types "<<endl;
    writeW << 1 << " angle types "<<endl;
    writeW << "\n";
    writeW << -Lx/2.0 << " " << (Lx-Lx/2.0) << " xlo xhi"<<endl;
    writeW << -Ly/2.0 << " " << (Ly-Ly/2.0) << " ylo yhi"<<endl;
    writeW << -Lz/2.0 << " " << (Lz-Lz/2.0) << " zlo zhi"<<endl;
    writeW << "\nMasses \n"<<endl;

    for(int j=1; j<=ntypepol+1;j++) writeW << j << " " << 1 << endl;
 
    writeW << "\nAtoms \n"<<endl;
    int atomcounter=0;

    for(int index_ch=0;index_ch<6;index_ch++){
      ///////////////////////////////////////////////////////////
      //make polymer
      double theta[2];
      double phi[4];
      double r = 0.8;
      ///////////////////////////////////////////////////////////
      //initial shift of chromosomes
      double xsh=0;
      double ysh=0;
      if(index_ch==0){xsh=-0.7;ysh=-0.75;}
      if(index_ch==1){xsh=+0.7;ysh=-0.75;}
      if(index_ch==2){xsh=-0.7;ysh=0;}
      if(index_ch==3){xsh=+0.7;ysh=0;}
      if(index_ch==4){xsh=-0.7;ysh=+0.75;}
      if(index_ch==5){xsh=+0.7;ysh=+0.75;}
      ///////////////////////////////////////////////////////////
      double com[3]={0.0,0.0,0.0};
      ///////////////////////////////////////////////////////////
      if(InitialConfg==1){
        double noise=0.9;
        double R1=20;
        for(int m=0;m<LChr[index_ch];m++){
	  if(m==0){
	    theta[0] = ((double)(rand()))/((double)(RAND_MAX))*2.0*M_PI;
	    phi[0] = ((double)(rand()))/((double)(RAND_MAX))*M_PI;
	    position[index_ch][0][0]=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(10./2.0);
	    position[index_ch][0][1]=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(10./2.0);
	    position[index_ch][0][2]=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(10./2.0);
	  }
          else{
	    int check=0;
	    while(distance(position[index_ch][m],center)>R1 or check==0){
		
	      theta[1]=theta[0]+noise*(((double)(rand()))/((double)(RAND_MAX))-0.5)*2.0*M_PI;
	      phi[1]=phi[0]+noise*(((double)(rand()))/((double)(RAND_MAX))-0.5)*M_PI;	
	
	      position[index_ch][m][0]=position[index_ch][m-1][0]+r*cos(theta[1])*sin(phi[1]);
	      position[index_ch][m][1]=position[index_ch][m-1][1]+r*sin(theta[1])*sin(phi[1]);
	      position[index_ch][m][2]=position[index_ch][m-1][2]+r*cos(phi[1]);
	      check=1;

	      if(m%100==0&&m>0)for(int mm=0;mm<m-1;mm++)if(distance(position[index_ch][mm],position[index_ch][m])<0.5)check=0;
	      if(position[index_ch][m][0]>Lx/2. or position[index_ch][m][0]<-Lx/2.) check=0;
	      if(position[index_ch][m][1]>Ly/2. or position[index_ch][m][1]<-Ly/2.) check=0;
	      if(position[index_ch][m][2]>Lz/2. or position[index_ch][m][2]<-Lz/2.) check=0;

	    }
	    theta[0]=theta[1];
	    phi[0]=phi[1];
	    com[0]+=position[index_ch][m][0]*1.0/LChr[index_ch];
            com[1]+=position[index_ch][m][1]*1.0/LChr[index_ch];
            com[2]+=position[index_ch][m][2]*1.0/LChr[index_ch];
          }
        }

        for(int m=0;m<LChr[index_ch];m++){
          writeW << atomcounter+1 << " " << 1 << " " << ColorChr[index_ch][m] << " " << position[index_ch][m][0]-com[0]+xsh*L/4.<<" " << position[index_ch][m][1]-com[1]+ysh*L/4. << " " << position[index_ch][m][2]-com[2] << " " << 0 << " " << 0 << " " << 0 << endl;
          atomcounter++;
        }

        cout << "DONE RW "<< index_ch << endl;

      }//close if over type initial confg
    }//loop over chromosomes

    /*
    if(InitialConfg==5){ //Mitotic
      double R=L/2.;
      double Lx=L;
      double Ly=L;
      double Lz=L;

      double com[3]={0.0,0.0,0.0};

      double s=1.0;
      double rchr=21*s;
      double p=s;
      int kk=6;
      double xc=0.3;
      //cout << " x " << x <<endl;
      //cin.get();
      for(int m=0;m<nbeads;m++){
        //in a single rosette there are 200 beads
        double phi=m*1.0/300.0*2.0*M_PI;

        position[index_ch][m][0]=rchr*(xc+(1.0-xc)*pow(cos(kk*phi),2.0))*cos(phi);
        position[index_ch][m][1]=rchr*(xc+(1.0-xc)*pow(cos(kk*phi),2.0))*sin(phi);
        position[index_ch][m][2]=p*phi/2.0/M_PI;
    
        com[0]+=position[index_ch][m][0];
        com[1]+=position[index_ch][m][1];
        com[2]+=position[index_ch][m][2];
    
        cout << m << " "  << position[index_ch][m][0]<< " " << position[index_ch][m][1] <<" " <<  position[index_ch][m][2]<< " "<< distance(position[index_ch][m],position[index_ch][m-1])<<endl;
      }

      for(int m=0;m<nbeads;m++){
        writeW << m+1 << " " << 1 << " " << colorWild[m] << " " << position[index_ch][m][0]-com[0]*1.0/nbeads<<" " << position[index_ch][m][1]-com[1]*1.0/nbeads << " " << position[index_ch][m][2]-com[2]*1.0/nbeads << " " << 0 << " " << 0 << " " << 0 << endl;
        if(m==1) cout << "distance " <<distance(position[m],position[m-1]) <<endl;
      }

      cout << "DONE MITOTIC " <<reps <<endl;
    }
    */

    ///////////////////////////////////////////////////////////
    //   FINISHED POLYMER
    ///////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////
    //   MAke PROTEINS
    ///////////////////////////////////////////////////////////
    if(InitialConfg!=2){
      for(int m=0;m<PCprotein;m++){
        double rx=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Lx-5);
        double ry=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Ly-5);
        double rz=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Lz-5);
        writeW<< GenomeCG+m+1 << " " << 2 << " " << ntypepol-3 << " " << rx<<" " << ry << " " << rz << " " << 0 << " " << 0 << " " << 0 << endl;
      }
      for(int m=0;m<H1Protein;m++){
        double rx=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Lx-5);
        double ry=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Ly-5);
        double rz=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Lz-5);
        writeW<< GenomeCG+PCprotein+m+1 << " " << 2 << " " << ntypepol-2 << " " << rx<<" " << ry << " " << rz << " " << 0 << " " << 0 << " " << 0 << endl;
      }
      for(int m=0;m<EuChrProtein;m++){
        double rx=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Lx-5);
        double ry=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Ly-5);
        double rz=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Lz-5);
        writeW<< GenomeCG+PCprotein+H1Protein+m+1 << " " << 2 << " " << ntypepol-1 << " " << rx<<" " << ry << " " << rz << " " << 0 << " " << 0 << " " << 0 << endl;
      }
      for(int m=0;m<HeterChrProtein;m++){
        double rx=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Lx-5);
        double ry=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Ly-5);
        double rz=(((double)(rand()))/((double)(RAND_MAX))-0.5)*(Lz-5);
        writeW<< GenomeCG+PCprotein+H1Protein+EuChrProtein+m+1 << " " << 2 << " " << ntypepol << " " << rx<<" " << ry << " " << rz << " " << 0 << " " << 0 << " " << 0 << endl;
      }
    }
    for(int t=0;t<tether;t++)writeW <<GenomeCG+nproteins+t+1<< " "<< 1 << " " << ntypepol+1 << " 0 0 " << -Lx/2.+2.0 << " 0 0 0 "<<endl;

    writeW << endl;
    writeW << "\n Velocities \n" <<endl;
    for(int j=0;j<GenomeCG+nproteins+tether; j++) writeW<<j+1 << " "<<0 << " "<<0<<" "<<0 <<endl;
    writeW << endl;

    cout << "BONDS" <<endl;
    //
    int bond[2];
    int angle[3];
    int beadcounter=0;
    int nbond=0;
    writeW << "\n Bonds \n"<<endl;
    //
    for(int index_ch=0;index_ch<6;index_ch++){
      for(int i=0;i<LChr[index_ch];i++){
        //LAST BOND IF RING
        /*if(i==nbeads-1){		
	  bond[nbond][0] = chainbead[i][j];
	  bond[nbond][1] = chainbead[0][j];
	}*/
	if(i<LChr[index_ch]-1){
	  bond[0] = i+beadcounter;
	  bond[1] = i+1+beadcounter;
          writeW << i+1+beadcounter <<" "<< 1 <<" "<< bond[0]+1<<" "<<bond[1]+1 << endl;
	}
        if(i==LChr[index_ch]-1)beadcounter+=LChr[index_ch];
      }
    }
    /////  ADD TETHERS FROM CHROMSOMES  /////////////////////////
    //end of chr2L
    writeW << beadcounter << " "<< 2 << " " << GenomeCG+nproteins+tether << " "<<LChr[0]<<endl;
    //start for chr2R
    writeW << beadcounter+1 << " "<< 2 << " " << GenomeCG+nproteins+tether << " "<<LChr[0]+1<<endl;
    //end of chr 3L
    writeW << beadcounter+2 << " "<< 2 << " " << GenomeCG+nproteins+tether << " "<<LChr[0]+LChr[1]+LChr[2]<<endl;
    //start of chr 3R
    writeW << beadcounter+3 << " "<< 2 << " " << GenomeCG+nproteins+tether << " "<<LChr[0]+LChr[1]+LChr[2]+1<<endl;
    //end of chrom 4
    writeW << beadcounter+4 << " "<< 2 << " " << GenomeCG+nproteins+tether << " "<<LChr[0]+LChr[1]+LChr[2]+LChr[3]+LChr[4]<<endl;
    //end of chrom X
    writeW << beadcounter+5 << " "<< 2 << " " << GenomeCG+nproteins+tether << " "<<LChr[0]+LChr[1]+LChr[2]+LChr[3]+LChr[4]+LChr[5]<<endl;
    //////////////////////////////////////////////////////////

    cout << "ANGLES" <<endl;
    writeW << "\n Angles \n"<<endl;
    nangle=0;
    beadcounter=0;
    for(int index_ch=0;index_ch<6;index_ch++){
      for(int i=0;i<LChr[index_ch];i++){
        /*if(i==nbeads-2){
	  angle[nangle][0] = chainbead[i][j];
	  angle[nangle][1] = chainbead[i+1][j];
	  angle[nangle][2] = chainbead[0][j];
        }
        if(i==nbeads-1){
	  angle[nangle][0] = chainbead[i][j];
	  angle[nangle][1] = chainbead[0][j];
	  angle[nangle][2] = chainbead[1][j];
        }*/
        if(i<LChr[index_ch]-2){
	  angle[0] = i+beadcounter;
	  angle[1] = i+1+beadcounter;
	  angle[2] = i+2+beadcounter;
          writeW << i+1 <<" "<< 1 <<" "<< angle[0]+1<<" "<<angle[1]+1 <<" "<< angle[2]+1<< endl;
	}
        if(i==LChr[index_ch]-1)beadcounter+=LChr[index_ch];
      }
    }

    writeW.close();
  }//close loop over replicas
  return 0;
}



double distance(double v1[3], double v2[3]){
  double d=0;
  for(int n=0;n<3;n++) d+=(v1[n]-v2[n])*(v1[n]-v2[n]);
  return sqrt(d);
}

void InitialiseColors(){
  for(int c=0;c<6;c++){
    for(int n=0;n<MAXbeads;n++){
      Het[c][n]=PC[c][n]=Eu[c][n]=Desert[c][n]=0;
      ColorChr[c][n]=0;
      AbundancePC[c][n]=0;
      AbundanceHet[c][n]=0;
      AbundanceEu[c][n]=0;
      AbundanceDesert[c][n]=0;
    }
  }
}

void Paint(int index_ch){
  for(int n=0;n<LChr[index_ch];n++){
    cout << "Chr" << index_ch << " Bead " << n << " has Het=" << Het[index_ch][n] << " Eu=" << Eu[index_ch][n] << " Desert="<<Desert[index_ch][n] << "and PC="<< PC[index_ch][n]<<endl;
    cout << "in respective abundances " << AbundanceHet[index_ch][n] << " " << AbundanceEu[index_ch][n] << " "<<AbundanceDesert[index_ch][n]<<" "<<AbundancePC[index_ch][n] <<endl;
    
    
    //RESOLVE
    //CONFLICTS -- Het VS Eu
    if(Het[index_ch][n]==1 && Eu[index_ch][n]==1){
      cout << "CONFLICT1 " << n << endl;
      cout << "Abundances " << AbundanceHet[index_ch][n]<< " " <<AbundanceEu[index_ch][n]<<endl;
      if(AbundanceHet[index_ch][n]>=AbundanceEu[index_ch][n])Eu[index_ch][n]=0;
      else Het[index_ch][n]=0;
      cout << "SOLVED"<<endl;
      //cin.get();
    }
    if(Het[index_ch][n]==1 && Eu[index_ch][n]==2){
      cout << "CONFLICT2 " << n << endl;
      cout << "Abundances " << AbundanceHet[index_ch][n]<< " " <<AbundanceEu[index_ch][n]<<endl;
      if(AbundanceHet[index_ch][n]>=AbundanceEu[index_ch][n])Eu[index_ch][n]=0;
      else Het[index_ch][n]=0;
      cout << "SOLVED"<<endl;
      //cin.get();
    }
    if(Het[index_ch][n]==1 && Desert[index_ch][n]==1){
      cout << "CONFLICT3 " << n << endl;
      cout << "Abundances " << AbundanceHet[index_ch][n]<< " " <<AbundanceDesert[index_ch][n]<<endl;
      if(AbundanceHet[index_ch][n]>=AbundanceDesert[index_ch][n])Desert[index_ch][n]=0;
      else Het[index_ch][n]=0;
      cout << "SOLVED"<<endl;
      //cin.get();
    }
    if(Desert[index_ch][n]==1 && Eu[index_ch][n]==1){
      cout << "CONFLICT4 " << n << endl;
      cout << "Abundances " << AbundanceDesert[index_ch][n]<< " " <<AbundanceEu[index_ch][n]<<endl;
      if(AbundanceDesert[index_ch][n]>=AbundanceEu[index_ch][n])Eu[index_ch][n]=0;
      else Desert[index_ch][n]=0;
      cout << "SOLVED"<<endl;
      //cin.get();
    }
    if(Desert[index_ch][n]==1 && Eu[index_ch][n]==2){
      cout << "CONFLICT4 " << n << endl;
      cout << "Abundances " << AbundanceDesert[index_ch][n]<< " " <<AbundanceEu[index_ch][n]<<endl;
      if(AbundanceDesert[index_ch][n]>=AbundanceEu[index_ch][n])Eu[index_ch][n]=0;
      else Desert[index_ch][n]=0;
      cout << "SOLVED"<<endl;
      //cin.get();
    }

    //Now I should have that beads can have state either
    // Blue or Desert/Black or Eu1 or Eu2 AND PC: 0 or 1 or 2
    // 
    //Easy stuff: PC==0
    if(PC[index_ch][n]==0 && Het[index_ch][n]==0 && Eu[index_ch][n]==0 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=1; //normal
    if(PC[index_ch][n]==0 && Het[index_ch][n]==1 && Eu[index_ch][n]==0 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=4; //only Blue
    if(PC[index_ch][n]==0 && Het[index_ch][n]==0 && Eu[index_ch][n]==1 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=7; //only Eu-high
    if(PC[index_ch][n]==0 && Het[index_ch][n]==0 && Eu[index_ch][n]==2 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=10; //only Eu-Low
    if(PC[index_ch][n]==0 && Het[index_ch][n]==0 && Eu[index_ch][n]==0 && Desert[index_ch][n]==1) ColorChr[index_ch][n]=13; //only Desert/Black
    // PC==1
    if(PC[index_ch][n]==1 && Het[index_ch][n]==0 && Eu[index_ch][n]==0 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=2; //PC-high
    if(PC[index_ch][n]==1 && Het[index_ch][n]==1 && Eu[index_ch][n]==0 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=5; //Blue-PC1
    if(PC[index_ch][n]==1 && Het[index_ch][n]==0 && Eu[index_ch][n]==1 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=8; //Eu1-PC1
    if(PC[index_ch][n]==1 && Het[index_ch][n]==0 && Eu[index_ch][n]==2 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=11; //Eu2-PC1
    if(PC[index_ch][n]==1 && Het[index_ch][n]==0 && Eu[index_ch][n]==0 && Desert[index_ch][n]==1) ColorChr[index_ch][n]=14; //Desert-PC1
    // PC==2
    if(PC[index_ch][n]==2 && Het[index_ch][n]==0 && Eu[index_ch][n]==0 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=3; //PC-low
    if(PC[index_ch][n]==2 && Het[index_ch][n]==1 && Eu[index_ch][n]==0 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=6; //Blue-PC2
    if(PC[index_ch][n]==2 && Het[index_ch][n]==0 && Eu[index_ch][n]==1 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=9; //Eu1-PC2
    if(PC[index_ch][n]==2 && Het[index_ch][n]==0 && Eu[index_ch][n]==2 && Desert[index_ch][n]==0) ColorChr[index_ch][n]=12; //Eu2-PC2
    if(PC[index_ch][n]==2 && Het[index_ch][n]==0 && Eu[index_ch][n]==0 && Desert[index_ch][n]==1) ColorChr[index_ch][n]=15; //Desert-PC2
    
    cout<< "fin color "<<n <<" "<<ColorChr[index_ch][n]<<endl;
  }
}

