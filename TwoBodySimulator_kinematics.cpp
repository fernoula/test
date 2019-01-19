/////// TwoBodySimulator.cpp.
/////// Author: Maria Anastasiou.
/////// Date:November15_2015.


// ----------bourdes--------------

#include <iostream>
#include <fstream>
#include <string.h>
#ifndef __Two_Body_Simulator__cpp
#define __Two_Body_Simulator__cpp



#include <string>
#include <math.h>   
#include <TGraph.h> 
#include <TGraph2D.h>
#include <TRandom3.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TCanvas.h>
#include <TMath.h> 
#include <TROOT.h> 
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <vector> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TTree.h>


//#include "test.cpp"
#include "cross_section.cpp"
#include "straggling.cpp"
#include "thetaphi.cpp"
#include "qqq_hit.cpp"
#include "TwoBodySimulator.hpp"


using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Energy Loss Class
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


EnergyLoss::EnergyLoss(string Eloss_file, double IonMass /*MeV/c^2*/) 
{

  double IonEnergy;
  double dEdx_e,dEdx_n; 
 
  string aux;

  ifstream Read(Eloss_file.c_str()); // reads the SRIM files

  last_point = 0;

  //cout << " Opening " << Eloss_file <<endl;
  if(!Read.is_open()) {
    cout << "*** EnergyLoss Error: File " << Eloss_file << " was not found." << endl;
    GoodELossFile = 0;
  } 
  
  else {
    GoodELossFile = 1;        
    Read >> aux >> aux >> aux;     // The first line has three strings (columns' description).
    points = 0;                       // Cout the number of points.

    do{
      Read >> IonEnergy >> dEdx_e >> dEdx_n;
      points++;
      // cout << IonEnergy << " " << dEdx_e << " " << dEdx_n << " " << points << endl ; 
    }
    while(!Read.eof());
    // while(!Read.eof() && (points<200));
    Read.close();

    // cout << points << endl ; 
    
    // Create the arrays depending on the number rows in the file.
    this->IonEnergy = new double[points];
    this->dEdx_e = new double[points];
    this->dEdx_n = new double[points]; 
    
    // Go to the begining of the file and read it again to now save the info in the newly created arrays.
    
    Read.open(Eloss_file.c_str());
    Read >> aux >> aux >> aux;

    for(int p=0; p<points; p++){
      Read >> IonEnergy >> dEdx_e >> dEdx_n;
      
      this->IonEnergy[p] = IonEnergy;
      this->dEdx_e[p] = dEdx_e;
      this->dEdx_n[p] = dEdx_n;
    }     

    Energy_in_range = 1;
    this->IonMass = IonMass;        
    c = 29.9792458;                  
    EvD = new TGraph();
    
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t EnergyLoss::GetDistance_new(Double_t InitialE, Double_t FinalE, Double_t StepSize)
{
  
  Double_t dist = 0;
  Double_t E = 0, Elast=0;
  //Double_t  Tolerance = 0.01;

  E = InitialE;
  
  
  while(E>FinalE){
    dist += StepSize;
    Elast=E;
    E = E - GetEnergyLoss(E,StepSize);
  }

  return ((dist-StepSize)-(StepSize*(Elast-FinalE)/(E-Elast)));
}


double EnergyLoss::GetDistance(double InitialE, double FinalE, double StepSize, int MaxSteps)
{

  double dist = 0;
  double E = 0;
  double  Tolerance = 0.01;

  for (int i=0; i<MaxSteps; i++) {

    dist += StepSize;

    E = GetFinalEnergy(InitialE, dist, StepSize);
    // cout << " E = " << E << endl;   

    if (((fabs(E-FinalE)/FinalE)<Tolerance) || (E<FinalE)) 
      break;
  }

 return dist;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

// Get the energy loss of the gas for an initial ion energy and a certain distance through the target.

double EnergyLoss::GetEnergyLoss(double energy /*MeV*/, double distance /*cm*/)
{

  Float_t a11=0.0, a12=0.0, a21=0.0, a22=0.0, a23=0.0, a32=0.0, a33=0.0;
  Float_t b11=0.0, b22=0.0, b33=0.0;
  Float_t a1=0.0, a2=0.0, b1=0.0, b2=0.0;
  Float_t K0=0.0, K1=0.0, K2=0.0; 
  Float_t N1=0.0, N2=0.0, N3=0.0;
  Float_t T1=0.0,T2=0.0,q1=0.0,q2=0.0;

  int i = -1;  
  if(energy < 0.01)
    return(0);
  
  // uses p points that were defined in EnergyLoss above when reading the SRIM files

  for(int p=0; p<points-1; p++){
   
    if(energy>=IonEnergy[p]  && energy<IonEnergy[p+1]){

      i = p+1;
      last_point = p;

      break;
    }
  } 

  if(i==-1){
    cout << "*** EnergyLoss Error: energy not within range: " << energy << endl;
    Energy_in_range = 0;
    return 0;
  }
  
  // If the initial energy is within the range of the function, get the stopping power for the initial energy. 
  
  //Ion Energy
  double x0=IonEnergy[i-1];
  double x1=IonEnergy[i];
  double x2=IonEnergy[i+1];

  //Total Energy Loss (electric + nuclear) for one step
  double y0=dEdx_e[i-1]+dEdx_n[i-1];
  double y1=dEdx_e[i]+dEdx_n[i];
  double y2=dEdx_e[i+1]+dEdx_n[i+1];

  a11=2/(x1-x0);
  a12=1/(x1-x0);
  a21=1/(x1-x0);
  a22=2*((1/(x1-x0))+(1/(x2-x1)));
  a23=1/(x2-x1);
  a32=1/(x2-x1);
  a33=2/(x2-x1);
  
  b11=3*((y1-y0)/((x1-x0)*(x1-x0)));
  b22=3*(((y1-y0)/((x1-x0)*(x1-x0)))+((y2-y1)/((x2-x1)*(x2-x1))));
  b33=3*((y2-y1)/((x2-x1)*(x2-x1)));  

  
  N1=(a21*a33*a12-a11*(a22*a33-a23*a32))/(a33*a12);
  N2=(b22*a33-a23*b33)/a33;
  N3=b11*(a22*a33-a23*a32)/(a33*a12);
  //cout<<"N1="<<N1<<" N2="<<N2<<" N3="<<N3<<endl;
  

  K0=(N2-N3)/N1;
  K1=(b11-a11*K0)/a12;
  K2=(b33-a32*K1)/a33;
  // cout<<"K0="<<K0<<" K1="<<K1<<" K2="<<K2<<endl;
  
  a1=K0*(x1-x0)-(y1-y0);  
  b1=-K1*(x1-x0)+(y1-y0);  
  a2=K1*(x2-x1)-(y2-y1);  
  b2=-K2*(x2-x1)+(y2-y1);  
  //cout<<"a1="<<a1<<" b1="<<b1<<" a2="<<a2<<" b2="<<b2<<endl;
  
  T1=(energy-x0)/(x1-x0);
  T2=(energy-x1)/(x2-x1);
  
 
  q1=(1-T1)*y0+T1*y1+T1*(1-T1)*(a1*(1-T1)+b1*T1);
  q2=(1-T2)*y1+T2*y2+T2*(1-T2)*(a2*(1-T2)+b2*T2);

  //cout<<"q1="<<q1<<" q2="<<q2<<endl; 
  
  return (q1*10*distance);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double EnergyLoss::GetInitialEnergy(double FinalEnergy /*MeV*/, double PathLength /*cm*//*dist*/, double StepSize/*cm*/)
{

  double Energy = FinalEnergy;
  int Steps = (int)floor(PathLength/StepSize);
  last_point = 0;

  // The function starts by assuming FinalEnergy is within the energy range, but this could be changed in the GetEnergyLoss() function.

  Energy_in_range = 1;

  for (int s=0; s<Steps; s++) {

    Energy = Energy + GetEnergyLoss(Energy,PathLength/Steps);
    // Energy = Energy + GetSpline(Energy,PathLength/Steps);

    if (!Energy_in_range)

      break;
  } 
  Energy = Energy + GetEnergyLoss(Energy,PathLength-Steps*StepSize);
  //Energy = Energy + GetSpline(Energy,PathLength-Steps*StepSize);

  if (!Energy_in_range)
    Energy = -1000; // Return an unrealistic value.
  
  // cout << "d: K_lf=" << FinalEnergy << "  K_lr=" << Energy << "  l=" << PathLength <<endl;

  return Energy;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double EnergyLoss::GetFinalEnergy(double InitialEnergy /*MeV*/, double PathLength /*cm*/, double StepSize/*cm*/)
{

  double Energy = InitialEnergy;
  int Steps = (int)floor(PathLength/StepSize);
  // The function starts by assuming InitialEnergy is within the energy range, but
  // this could be changes in the GetEnergyLoss() function.
  Energy_in_range = 1;

  for (int s=0; s<Steps; s++) {
    Energy = Energy - GetEnergyLoss(Energy,PathLength/Steps);  

    if (!Energy_in_range)
      break;
  }  
  Energy = Energy - GetEnergyLoss(Energy,PathLength-Steps*StepSize);
 
  if (!Energy_in_range) 
    Energy = -1000;
  //cout << "O: K_bw=" << InitialEnergy << "  K_br=" << Energy << "  l=" << PathLength <<endl;
  return Energy;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calulates the ion's time of flight in ns.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double EnergyLoss::GetTimeOfFlight(double InitialEnergy, double PathLength, double StepSize)
{

  double TOF = 0;
  double Kn = InitialEnergy;
  int Steps = (int)(PathLength/StepSize);  
  if (IonMass==0) {
    cout << "Error: Time of flight cannot be calculated because mass is zero." << endl;
  }
  else {
    for (int n=0; n<Steps; n++) {
      TOF += sqrt(IonMass/(2*Kn))*StepSize/c;               // DeltaT going from point n to n+1.
      Kn -= GetEnergyLoss(Kn, StepSize); //After the TOF is added the K.E. at point n+1 is calc.} 
    }
  }
    return TOF;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EnergyLoss::SetIonMass(double IonMass)
{

  this->IonMass = IonMass;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////
////////////////////------add LookUp Funtion----------------/////////////////////////////
void EnergyLoss::InitializeLookupTables(Double_t MaximumEnergy, Double_t MaximumDistance, 
				    Double_t DeltaE, Double_t DeltaD){

  this->MaximumEnergy = MaximumEnergy;
  this->MaximumDistance = MaximumDistance;
  this->DeltaD = DeltaD;
  this->DeltaE = DeltaE; 

  int noE = (int)ceil(MaximumEnergy / DeltaE );
  int noD = (int)ceil(MaximumDistance / DeltaD );

  
  
  if (!(EtoDtab = new Double_t[noE])||
      !(DtoEtab = new Double_t[noD])){
    cerr << "Could not allocate memory for " << noE << " " << noD << endl;
  }

  //Double_t D;
  int i;
  //-----------------------------------------------------------
  DtoEtab[0] = MaximumEnergy;
  cout << " Number of distance entries " << noD << endl;
  
  for (i=1; i<noD; i++){
    DtoEtab[i] = GetFinalEnergy(DtoEtab[i-1],DeltaD,0.05*DeltaD);
    
    //cout<<"2: "<<DtoEtab[i]<< " D="<<D<<endl;
    //if (i%1000 == 0 ){
    // cout << " Passed  " << i << endl;
      //}
  }

  int j;
  
  cout << " Number of Energy entries " << noE << endl;
  
  EtoDtab[0] = 0.;
  for (j=1;j<noE;j++){
    EtoDtab[j] = EtoDtab[j-1] + GetDistance_new((MaximumEnergy-(j-1)*DeltaE),
					    (MaximumEnergy-(j)*DeltaE),
					    (0.05*DeltaD));     
    //if (j%1000 == 0 ){
    //cout << " Passed  " << j << endl;
      //}
    
  }
  //cout << " Passed  " << j << endl;
  //-----------------------------------------------------------

}
//=======================================================
void EnergyLoss::PrintLookupTables(){
  int noE = (int)ceil(MaximumEnergy / DeltaE );
  int noD = (int)ceil(MaximumDistance / DeltaD );
  int i;
  cout << "Maximum Energy = " << MaximumEnergy << " " << DeltaE << noE << endl;
  for (i=0;i<noE;i++){
    cout << "E,D= "<< MaximumEnergy - DeltaE*i << "," <<EtoDtab[i] <<endl; 
  }
  cout << "Maximum Distance = "<< MaximumDistance << " " << DeltaD << noD << endl;
  for (i=0;i<noD;i++){
    cout << "D,E= "<< DeltaD*i << "," <<DtoEtab[i]<<endl; 
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t EnergyLoss::GetLookupEnergy(Double_t InitialEnergy, Double_t distance){

  Double_t D1,D2,D;
  Double_t E1,E2,E;
  int index;
  int imhere;

  if (InitialEnergy<0 || InitialEnergy>MaximumEnergy){
    return(-1.);
  }

  // Find the distance for which the initial energy is matched, interpolating 
  index = (int)floor((MaximumEnergy-InitialEnergy) / DeltaE) ;

  D1=EtoDtab[index];D2=EtoDtab[index+1];

  E1=MaximumEnergy-index*DeltaE; E2=MaximumEnergy-(index+1)*DeltaE; 

  D = ( InitialEnergy - E1 ) / ( E2 - E1 ) * (D2-D1) + D1 ;
  //cout << " 1:  "<< index << ' ' <<E1 << ' ' << E2 << ' ' << D1 << ' ' << D2 << ' ' << D << endl;


  //Still in the table ?
  if(((D+distance <=0)) || ((D+distance)> MaximumDistance)){
    cout<<"i m here"<<endl;
    //imhere++;
    return(0.);
  }

  // Lookup what energy is reached for (D + distance) and interpolate
  index = (int)floor((D + distance) / DeltaD );
  //cout << endl;
  //cout << "D: " << index << " " << distance << " " << DeltaD << " " << D << endl;

  E1=DtoEtab[index];E2=DtoEtab[index+1];
  //cout << "E1: " << E1 << " E2: " << E2 << endl;
  D1=index*DeltaD;D2=(index+1)*DeltaD;

  E = ( (D+distance) - D1 ) / ( D2 - D1 ) * (E2-E1) + E1 ;
  //cout <<" 2:   "<< E1 << ' ' << E2 << ' ' << D1 << ' ' << D2 << ' ' << E << endl;

  //cout<<"E == "<<E<<endl;
  //cout<<imhere<<endl;
  return(E);
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//--------------------Simulator Class--------------------------------/////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

Simulator::Simulator(string FileELoss1,string FileELoss2,string FileELoss3,string FileELoss4,
		     string FileELoss5,
		     string WorldCoordinatesFilename, string PCWirePCRFile) 
{

  rand = new TRandom3();
  rand->SetSeed();
  ConvFromUtoMeV = 931.494061;

  //initializing mass of atoms in amu
  mB=0.0, mH=0.0, mL=0.0, mT=0.0, mH2=0.0;// mass of the Ion of the beam, of the heavy recoil, of the light recoil,of the target respectively

  //mass of atoms in MeV
  MB=mB;;
  MT=mT;
  ML=mL;
  MH=mH;
  MH2=mH2;

  //pointers
  IonInGas = new EnergyLoss(FileELoss1,MB); // loads the EnergyLoss class for the SRIM file of the Ion of the beam in the gas. 
  // cout << "loaded ion in gas"<<endl ;
 
  proton = new EnergyLoss(FileELoss2,ML); 
  // cout << "loaded proton in gas"<<endl ; 

  Heavy_ion = new EnergyLoss(FileELoss3,MH);
  //cout << "loaded Heavy_ion in gas"<<endl ;

  Heavy_ion2 = new EnergyLoss(FileELoss4,MH2);
  //cout << "loaded Heavy_ion in gas"<<endl ;
 
  proton_Si = new EnergyLoss(FileELoss5,ML); 
  //cout << "loaded Proton_Si in Silicon"<<endl ; 

   //Distances:
  dist = 0.0, D=0.0;
 
 
 //Energies:
  FinalE = 0.0,InitialEnergy=0.0;
  Ep_in=0.0,Eh_in=0;
  Ep_final=0.0,Eh_final=0,Ebeam_final=0;
  ResidualEn = 0.0;
  Ex = 0.0;

  //................................
  Hi_e=0.0;//Heavy ion total energy.
  Hi_p=0.0;//Heavy ion momentum.

  
  ee_L_lab=0.0,ee_H_lab=0.0,ke_L_lab=0.0,ke_H_lab=0.0;
  /////
  EH_cm=0.0, Px_H_cm=0.0, Py_H_cm=0.0, Pz_H_cm=0.0;
  Px_H_lab=0.0, Py_H_lab=0.0, Pz_H_lab=0.0;
  Px_P_lab=0.0, Py_P_lab=0.0, Pz_P_lab=0.0;
  /////
  ee_1_lab=0.0,ee_2_lab=0.0,ke_1_lab=0.0,ke_2_lab=0.0;
  
  //................................
  //Angles:
  theta_L=0.0, theta_H=0.0, phi_L=0.0, phi_H=0.0;
  ////
  theta_h1=0, theta_h2=0, phi_h1 = 0.0, phi_h2 = 0.0; 

  
  //Time of Flights:
  TOF_rxn=0.0,TOF_light=0.0,TOF_heavy=0.0,TOF_p=0.0,TOF_beam=0.0;

  //for the Geometry of Detector plates:
  x11=0,x22=0,x33=0,x44=0,y11=0,y22=0,y33=0,y44=0;
  xr = 0.0,yr = 0.0, zr=0.0;
 
  //Radius of the ANASEN Rings of SX3 Detectors 
  Ra=8.9;//cm
  La=55.0545;//cm

  if (WorldCoordinatesFilename!="")
    InitWorldCoordinates(WorldCoordinatesFilename);

  if (PCWirePCRFile!="")
    InitPCWirePCR(PCWirePCRFile);

}


Simulator::~Simulator()  
{
  //delete rand;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

bool Simulator::SetMasses(Float_t mB,Float_t mT,Float_t mL,Float_t mH, Float_t mH2)
{

  this->mB = mB;
  this->mT = mT;
  this->mL = mL;
  this->mH = mH;
  this->mH2 = mH2;
  
  MB = mB;
  MT = mT;
  ML = mL;
  MH = mH;
  MH2 = mH2;


  IonInGas->SetIonMass(MB); //uses SetIonMass function from EnergyLoss class
  proton->SetIonMass(ML);
  Heavy_ion->SetIonMass(MH);
  Heavy_ion2->SetIonMass(MH2);
  proton_Si->SetIonMass(ML);
 
  return 1;
}



// This event generator is based on a flat energy distribution from a MinE to a MaxE.
// Generates a random final energy at the reaction point after the Ion of the beam has lost energy in the gas. 

bool Simulator::GenerateEvent(CS *cross, Strag *straggle)
{

  FinalE = rand->Uniform(MinBeamEnergy, MaxBeamEnergy);
  //cout << "FinalE = " << FinalE  << " " << MinBeamEnergy << " " <<MaxBeamEnergy << endl;

  Float_t crossSection = rand->Uniform(0,1);

  if (crossSection > cross->CrossSection<Simulator>(FinalE,this))
    return false;
  
  
  // excited states of 21Na:
  Float_t Ex_Heavy[10]={0.0,0.3319,1.7161,2.4238,3.544,4.169,4.984,5.457,6.070,6.468};
  Int_t RandomIndex = fmod(rand->Uniform(0,10),10);

  Ex = Ex_Heavy[RandomIndex];

  //cout << " Index: " << RandomIndex << " Ex: " << Ex << endl;
  
  dist = IonInGas->GetDistance(MaxBeamEnergy, FinalE, 0.1, 550); // uses GetDistance function from EnergyLoss class for the Ion of the beam 
 
  Double_t y_str =0.0, z_str =0.0;  
  straggle->Range_Straggling(dist,y_str,z_str);
  y_str = y_str/10.0; // convert them to cm from mm
  z_str = z_str/10.0;
 
  dist = rand->Gaus(dist,z_str);
 
  if (dist > La){
    return 0;
  }
  zr =(La-dist); 

  // add xr, yr to be created with a random distribution +/- 0.5cm for our beam spot
  xr = rand->Uniform(-0.5,+0.5);
  yr = rand->Uniform(-0.5,+0.5);

  xr = rand->Gaus(xr,y_str);
  yr = rand->Gaus(yr,y_str);
  
  
  return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

bool Simulator::SetBeamEnergyRange(Float_t MinE, Float_t MaxE)
{

  MinBeamEnergy = MinE;
  MaxBeamEnergy = MaxE;
  // cout << "Beam energy from " << MinE << " to " << MaxE << endl;

  return 1;
}


//////////////////////////////////////////////////////////////////////////////////////////////////

bool Simulator::InitWorldCoordinates(string WorldCoordinatesFilename)
{
  ifstream WorldCoordFile;
  Int_t line_counter, NumDets;
  string line; 

  WorldCoordFile.open(WorldCoordinatesFilename.c_str());
  
  if (!WorldCoordFile.is_open()) {
    cout << ">\tTERROR: The file \"" << WorldCoordinatesFilename << "\" couldn't be opened." << endl;
    WorldCoordinatesLoaded = 0;
  }
  else {
    cout << ">\tThe file \"" << WorldCoordinatesFilename << "\" was opened successfully." << endl;
    
    // Count the number of lines with data in the text file.
    line_counter = 0;
    do{
      getline(WorldCoordFile, line);
      // Only count non-empty lines.
      if (!line.empty())          
	line_counter++; 
    } while (!WorldCoordFile.eof());
    WorldCoordFile.close();
    
    NumDets = line_counter-1;
    DetectorID = new Int_t[NumDets];
    ZOffset = new Float_t[NumDets];
    XAt0 = new Float_t[NumDets];
    XAt4 = new Float_t[NumDets];
    YAt0 = new Float_t[NumDets];
    YAt4 = new Float_t[NumDets];
    
    WorldCoordFile.open(WorldCoordinatesFilename.c_str());
    getline(WorldCoordFile, line);    // The 1st line is for column description.
    for (Int_t i=0; i<NumDets; i++) {
      WorldCoordFile >> DetectorID[i] >> ZOffset[i] >> XAt0[i] >> XAt4[i] >> YAt0[i] >> YAt4[i];
      getline(WorldCoordFile, line);  // The last column is for comments
    }
    WorldCoordFile.close();
    WorldCoordinatesLoaded = 1;   // Ready to use the world coordinates.
    cout << ">\tWorld coordinates loaded for " << NumDets << " SX3 detectors." << endl;
  }
  return 1;
}

//------------------- read PC Radius for each wire --------------------------------////////////////////

bool Simulator::InitPCWirePCR(string PCWirePCRFile)
{
  ifstream pcr_file;
  string line;
  Double_t pcr_dum;
  Int_t WireNumber=0;
  Int_t WireNum = 24;
  
  Wire_ID = new Int_t[WireNum];
  PC_Radius = new Float_t[WireNum];
  

  for (Int_t i=0; i<WireNum; i++) {
    Wire_ID[i]=0;
    PC_Radius[i]=0.0;
  }
  
  pcr_file.open(PCWirePCRFile.c_str());
    
  if (pcr_file.is_open()) {
    cout << "File for different PC Radius ";
    cout << PCWirePCRFile << " opened successfully. " << endl;
    getline (pcr_file,line);//Skips the first line in PCWireCalFilename
    //cout<<"line = "<<line<<endl;

    while (!pcr_file.eof()) {
      pcr_file >> WireNumber >> pcr_dum;
      Wire_ID[WireNumber] = WireNumber;
      PC_Radius[WireNumber] = pcr_dum;
  
      //cout << Wire_ID[WireNumber] << " " << PC_Radius[WireNumber] << endl;
    }  
  }
  else cout << "DIDN'T OPEN FILE!" << endl;
  
  
  return 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////--------------SimTheta----------------//////////////////////////////////////

//-------It doesn't use MinE and MaxE (only as function input), just FinalE of Beam calculated in SetBeamEnergyRange & GenerateEvent-----////////////////


bool Simulator::SimThetaPhiE_in(Float_t MinE, Float_t MaxE) 
{

  MinBeamEnergy = MinE;
  MaxBeamEnergy = MaxE;
  rand->SetSeed();
  
  //for light particle: for first reaction:
  Float_t theta_L_cm_deg = rand->Uniform(0.0, 180.0);

  Float_t theta_L_cm = theta_L_cm_deg*TMath::Pi()/180;

  Float_t phi_L_cm_deg = rand->Uniform(0.0,360.0);
  Float_t phi_L_cm = phi_L_cm_deg*TMath::Pi()/180;
  //==================================================================================================


  //Energy & Momentum of combined parents before rxn in lab frame:

  Float_t Parent_M = MB+MT;

  Float_t Beam_E = FinalE + MB; // beam energy  where the FinalE is taken from GenerateEvent function where it is created for the projectile
  Float_t Beam_P = sqrt(Beam_E*Beam_E-MB*MB); // beam momentum
 

  //if(!(Beam_P==0 || Beam_E==0) return 0;
  Float_t Target_E = MT; // target energy just the rest mass
 
  //initializing the four-momentum vectors for the Parent, Light & Heavy particles in Lab
  TLorentzVector Beam_LV(0.,0.,0.,MB);
  TLorentzVector Target_LV(0.,0.,0.,MT);
  TLorentzVector Parent_LV(0.,0.,0.,Parent_M);  // initialization so use only rest masses
  TLorentzVector L_LV(0.,0.,0.,ML);
  TLorentzVector H_LV(0.,0.,0.,MH);

  //---used for testing----------
  // same as H_LV but used to test the heavyE at the CM before boosting back to lab frame
  TLorentzVector H1_LV(0.,0.,0.,MH);  

  //////////////////////////////////
  
  Float_t bx = 0.0;    // beam momentum components for 4-vector
  Float_t by = 0.0;
  Float_t bz = Beam_P;

  Float_t tx = 0.0;   // target momentum components for 4-vector
  Float_t ty = 0.0;
  Float_t tz = 0.0;

  Beam_LV.SetPxPyPzE(bx,by,bz,Beam_E);
  Target_LV.SetPxPyPzE(tx,ty,tz,Target_E);

  Parent_LV = Beam_LV + Target_LV;

  //-- if you want to extract later a boostVector at the compound system
  ParentE_lab = Parent_LV.E();
  Px_P_lab = Parent_LV.Px();
  Py_P_lab = Parent_LV.Py();
  Pz_P_lab = Parent_LV.Pz();
  //------------------------------------
  
  boostv = Parent_LV.BoostVector();

  //boosting back 4-momentum vectors for parent from Lab to CM frame.
  Parent_LV.Boost(-boostv);

  // components of 4-momentum vectors in CM frame for Parents. 
  Float_t Qe_cm = Parent_LV.E();
  //cout<<"Qe_cm="<<Qe_cm<<endl;
  if((Parent_LV.M()-ML-MH-Ex)<=0){
    //cout << "maria watch out for Ex " << Ex << " BeamE: " << FinalE << endl;
    return 0;
    }
  
  //calculating  energy & magnitude of momentum for one of the daughter particle in the cm-frame.                         
  Float_t EL=(ML*ML-(MH+Ex)*(MH+Ex)+Qe_cm*Qe_cm)/(2*Qe_cm);

   //----testing different formulas for EL_at_CM see if it matches the one of EL right above-----
  Double_t Qv = (MB+MT-ML-MH-Ex);
  //iliadis appendix formula: Elight at cm
  Float_t EL_test = (MH/(ML+MH))*(Qv + FinalE*(1-(MB/(ML+MH)))) + ML;
  Float_t EL_test_2 = ((MH+Ex)*(Qe_cm-MH-ML-Ex))/(ML+MH+Ex) + ML;

  //cout << " EL: " << EL << " EL_test: " << EL_test << " EL_test_2: " << EL_test_2 << endl;
  //////////////////////////////
  
  Float_t PL= sqrt(EL*EL-ML*ML);
  if(PL <= 0.0)
  cout<<"PL="<<PL<<endl;

  if (!((PL >0)&(EL>0)))
    return 0;

  if (PL == 0) return 0;

  //components of 4-momentum vectors in CM frame for lighter particle.
  Float_t q1x=PL*sin(theta_L_cm)*cos(phi_L_cm);
  Float_t q1y=PL*sin(theta_L_cm)*sin(phi_L_cm);
  Float_t q1z=PL*cos(theta_L_cm);

  L_LV.SetPxPyPzE(q1x,q1y,q1z,EL); 

  //--testing Eheavy_at_cm to compare Lorentz vector result with formula from Iliadis-----------
  H1_LV = Parent_LV-L_LV;
  EH_cm = H1_LV.E();  
  //iliadis appendix formula: Eheavy at cm
  Double_t EH_test = (ML/(ML+MH))*(Qv + FinalE*(1-(MB/(ML+MH)))) + MH;

  //cout << " EH_test_Lor: " << EH_cm << " EH_Lor - Ex: " << EH_cm - Ex << " EH_test " << EH_test << " Ex: " <<  Ex << endl;
  
  ////////////////////////////////////
  
  //components of 4-momentum vectors in lab frame for light particle.
  L_LV.Boost(boostv);
 
  ee_L_lab = L_LV.E();
  // cout<<"ee_L_lab="<<ee_L_lab<<endl;

  //boosting back parents from CM to Lab frame.
  Parent_LV.Boost(boostv);

  //calculating Heavy particle four momentum in Lab frame.
  H_LV = Parent_LV - L_LV;

 //components of 4-momentum vectors in lab frame for heavy particle.
  ee_H_lab = H_LV.E();
  //cout << " ee_H_lab= " << ee_H_lab << endl;  
   
  // Need the momenta as well to create the initial Heavy TLorentzVector for
  // the SecondProton Function
  Px_H_lab = H_LV.Px();
  Py_H_lab = H_LV.Py();
  Pz_H_lab = H_LV.Pz();

   //////////
 
  //calculating the angles for Lighter & Heavier particles in Lab frame.:   
   phi_L = L_LV.Phi(); 


   while (phi_L < 0) 
    {
      // cout << phi_L << endl;
      phi_L += 2* M_PI;
    }
      
   while (phi_L > 2*M_PI) 
    { 
      // cout << phi_L << endl;
      phi_L-= 2*M_PI;
    } 

   phi_H = H_LV.Phi(); 

   while (phi_H < 0) 
    {
      // cout << phi_H << endl;
      phi_H += 2* M_PI;
    }
      
   while (phi_H > 2*M_PI) 
    { 
      // cout << phi_H << endl;
      phi_H -= 2*M_PI;
    } 

   theta_L = L_LV.Theta();
   theta_H = H_LV.Theta();
   //cout<<"phi_L="<<phi_L<<" phi_H="<<phi_H<<" theta_L="<<theta_L<<" theta_H="<<theta_H<<endl;

   ke_L_lab=ee_L_lab-ML;
   ke_H_lab=ee_H_lab-MH;
   //cout<<"ke_L_lab="<<ke_L_lab<<"ke_H_lab="<<ke_H_lab<<endl;

   
   //---testing the formula at the lab frame from Iliadis to compare with the TLorentz Vector for ke_L_lab---------------
   Double_t a = (1+ML/MH);
   Double_t b = -2*sqrt(MB*ML*FinalE)*cos(theta_L)/MH;
   Double_t c = -Qv - FinalE*(1-MB/MH);
   Double_t s1 = (-b + sqrt(b*b-4*a*c))/2./a;  // sqrt(Ea) = sqrt(Beam_KE)
   Double_t s2 = (-b - sqrt(b*b-4*a*c))/2./a;

   if(s2>s1)
     cout << " s1 " << s1 << " s2 " << s2 <<  endl;
   if(s2>s1)
     std::swap(s1,s2);

   Double_t ke_L_lab_test = s1*s1;

   //compare TLorentz vectors results with formulas from iliadis result
   //cout << " ke_L_lab= " << ke_L_lab << "  " << ke_L_lab_test << endl;
   
   //////////////////
   //from energy conservation:
   Double_t ke_H_lab_test = Qv + FinalE - ke_L_lab_test;

   //compare TLorentz vectors results with formulas from iliadis result
   //cout << " ke_H_lab: " << ke_H_lab << " ke_H_lab - Ex: " << ke_H_lab - Ex << " ke_H_lab_test: " << ke_H_lab_test << " ke_H_lab_test_add_MH_and_Ex: " << ke_H_lab_test + MH + Ex << endl;

   ///////////////
  
   return 1;
} 
///////////////////////////////////////////////////////////////////
//bool Simulator::SecondProton(const Float_t ee_H_lab, const Float_t Ex)
bool Simulator::SecondProton(const Float_t ee_H_lab, const Double_t Px_H_lab,
			     const Double_t Py_H_lab, const Double_t Pz_H_lab,
			     const Float_t Ex) 
{

  //for second reaction:
  Float_t theta_h1_cm_deg = rand->Uniform(0.0, 180.0);
  Float_t theta_h1_cm = theta_h1_cm_deg*TMath::Pi()/180;

  Float_t phi_h1_cm_deg = rand->Uniform(0.0,360.0);
  Float_t phi_h1_cm = phi_h1_cm_deg*TMath::Pi()/180;

  //initializing the four-momentum vectors for the Heavy particle & two decay particles.
  
  TLorentzVector Hi_LV(0.,0.,0.,MH);
  TLorentzVector H1_LV(0.,0.,0.,ML);
  TLorentzVector H2_LV(0.,0.,0.,MH2);
  
  //--used for tests same way as H1_LV at the previous function
  TLorentzVector H2_2_LV(0.,0.,0.,MH2);
  //////////////////////////////

  Hi_LV.SetPxPyPzE(Px_H_lab,Py_H_lab,Pz_H_lab,ee_H_lab);

  boostv = Hi_LV.BoostVector();

  //boosting back 4-momentum vectors for heavy particle to CM frame
  Hi_LV.Boost(-boostv);

  // components of 4-momentum vectors in CM frame for the heavy particle.
 
  Float_t Ee_cm = Hi_LV.E();
  //cout << " Ee_cm= " << Ee_cm << endl;
  if((MH+Ex-ML-MH2)<=0){
    //cout << "maria watch out " << " Ex " << Ex << endl;
    return 0;
    }
  
  //calculating  energy & magnitude of momentum for one of the daughter particle in the cm-frame.                         
  Float_t E1=(ML*ML-MH2*MH2+(Ee_cm)*(Ee_cm))/(2*(Ee_cm));

  //--test--------------
  //----testing different formulas for E1_at_CM see if it matches the one of E1 right above-----
  Double_t Qv2 = (MH+Ex-ML-MH2);
  //iliadis appendix formula: Elight at cm
  Float_t E1_test = (MH2/(ML+MH2))*(Qv2 + ke_H_lab*(1-(MH/(ML+MH2)))) + ML;
  
  //cout << " E1: " << E1 << " E1_test: " << E1_test <<  endl;
  /////////////////////////////////
  
  Float_t P1= sqrt(E1*E1-ML*ML);                       
  //cout<<"P1="<<P1<<endl;
  if (P1 == 0)
    return 0;

  //components of 4-momentum vectors in CM frame for first daughter particle.
  Float_t p1x=P1*sin(theta_h1_cm)*cos(phi_h1_cm);
  Float_t p1y=P1*sin(theta_h1_cm)*sin(phi_h1_cm);
  Float_t p1z=P1*cos(theta_h1_cm);

  H1_LV.SetPxPyPzE(p1x,p1y,p1z,E1); 

  //--test-----------
  //--testing Eheavy2_at_cm to compare Lorentz vector result with formula from Iliadis-----------
  H2_2_LV = Hi_LV-H1_LV;
  Double_t EH2_test_CM_Lor = H2_2_LV.E();
  //iliadis appendix formula: Eheavy at cm
  Double_t EH2_CM_test = (ML/(ML+MH2))*(Qv2 + ke_H_lab*(1-(MH/(ML+MH2)))) + MH2;

  //cout << " EH2_test_CM_Lor: " << EH2_test_CM_Lor <<  " EH2_CM_test " << EH2_CM_test << " Ex: " <<  Ex << endl;

  ///////////////////////////////////////
  
  //components of 4-momentum vectors in lab frame for first daughter particle.
  H1_LV.Boost(boostv);
  
  ee_1_lab = H1_LV.E();
  //cout<<"ee_1_lab="<<ee_1_lab<<endl;
 
  //boosting back heavy ion 4-momentum in lab frame.
  Hi_LV.Boost(boostv);

  //calculating 4-momentum vector for second daughter particle.
  H2_LV = Hi_LV - H1_LV;

 //components of 4-momentum vectors in lab frame for second daughter particle.
 
  ee_2_lab = H2_LV.E();
  //cout << " ee_2_lab: " << ee_2_lab << endl; 

  //calculating the angles of Daughter nucleus: 
  phi_h1 = H1_LV.Phi(); 
  phi_h2 = H2_LV.Phi();

  while (phi_h1 < 0) 
    {
      phi_h1 += 2* M_PI;
    }
      
   while (phi_h1 > 2*M_PI) 
    { 
      phi_h1-= 2*M_PI;
    } 

   while (phi_h2 < 0) 
    {
      phi_h2 += 2* M_PI;
    }
      
   while (phi_h2 > 2*M_PI) 
    { 
      phi_h2 -= 2*M_PI;
    } 
  
  theta_h1 = H1_LV.Theta();
  theta_h2 = H2_LV.Theta(); 
  // cout<<"phi_h1="<<phi_h1<<" phi_h2="<<phi_h2<<" theta_h1="<<theta_h1<<" theta_h2="<<theta_h2<<endl;

  ke_1_lab=ee_1_lab-ML;
  ke_2_lab=ee_2_lab-MH2;
  //cout << " ke_1_lab= " << ke_1_lab <<" ke_2_lab= " << ke_2_lab << " Ex = " << Ex << endl;


  //---test---------------
   Double_t a = (1+ML/MH2);
   Double_t b = -2*sqrt(MH*ML*ke_H_lab)*cos(theta_h1)/MH2;
   Double_t c = -Qv2 - ke_H_lab*(1-MH/MH2);
   Double_t s1 = (-b + sqrt(b*b-4*a*c))/2./a;  // sqrt(Ea) = sqrt(Beam_KE)
   Double_t s2 = (-b - sqrt(b*b-4*a*c))/2./a;

   if(s2>s1)
     cout << " s1 " << s1 << " s2 " << s2 <<  endl;
   if(s2>s1)
     std::swap(s1,s2);

   Double_t ke_1_lab_test = s1*s1;
   
   //cout << " ke_1_lab: " << ke_1_lab << " ke_1_lab_test: " << ke_1_lab_test << " Ex: " << Ex << endl;

   //////////////////
   //from energy conservation:
   Double_t ke_2_lab_test = Qv2 + ke_H_lab - ke_1_lab_test;

   //compare TLorentz vectors results with formulas from iliadis result
   //cout << " ke_2_lab: " << ke_2_lab << " ke_2_lab_test: " << ke_2_lab_test << " ee_2_lab_add_MH2: " << ke_2_lab_test + MH2 << endl;

   ///////////////
   
  return 1;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////------------------------Needle Simulation-------------------////////////////////////////////////////
Float_t Simulator::NeedleEnergy(Float_t heavyEn, Float_t heavyTheta, Float_t IntPoint)
{
  heavyEn = ke_H_lab;
  heavyTheta = theta_H;
  IntPoint = zr;
 
  Float_t Path_H_Full = 0.0;
  Float_t Path_H_IntPtoNeedle = 0.0;
  Float_t E_FullPath = 0.0;
  Float_t BeamAtNeedle = 0.0;
  Float_t deltaBeam = 0.0;
  Float_t y_dist = 0.0;

  Float_t E_atNeedleStart = 0.0;
  Float_t Efinal = 0.0;

  //cout << heavyEn << " "  << heavyTheta*180/3.14 << " " << zr << endl; 

  if(heavyTheta != 0.0 || heavyTheta != TMath::Pi()/2 || heavyTheta != TMath::Pi() || 
     heavyTheta != - TMath::Pi()/2 || heavyTheta != -TMath::Pi())
    {
      Path_H_Full = 2.2 / sin(heavyTheta);  // 2.2 cm is the radius of the IC, inner radius of the PC
      if(IntPoint > 17.8)
	Path_H_IntPtoNeedle = (IntPoint - 17.8)/cos(heavyTheta);  // 17.8 cm is the position the needle extends upstream of QQQ
    }


  if(IntPoint <= 17.8)
    {
      E_FullPath = Heavy_ion->GetLookupEnergy(heavyEn, Path_H_Full);
      BeamAtNeedle = IonInGas->GetLookupEnergy(MaxBeamEnergy, La-17.8);
      deltaBeam = BeamAtNeedle - FinalE;

      //y_dist = sin(heavyTheta)*(Heavy_ion->GetDistance_new(heavyEn, 0.01, 0.1));

      if(E_FullPath <= 0.0)
	ResidualEn = heavyEn + deltaBeam;
      else
	ResidualEn = (heavyEn - E_FullPath) + deltaBeam;

      /*if(E_FullPath > 0.0)
	{
	  cout << " y_dist: " << y_dist << " E_FullPath: " << E_FullPath << " ResidualEn: " << ResidualEn << endl;
	  cout << " ResidualEn: " << ResidualEn << " BeamE: " << FinalE << " BeamAtNeedle: " << BeamAtNeedle << " IntP: " << zr 
	       << " Theta: " << heavyTheta << " EnergyHeavy: " << heavyEn <<  " deltaBeam: " << deltaBeam << endl;
	       }*/
      
     }
  else
    {
      E_atNeedleStart = Heavy_ion->GetLookupEnergy(heavyEn,Path_H_IntPtoNeedle);
        
      if(E_atNeedleStart <= 0.0)
	ResidualEn = -10.0;
      else if(Path_H_Full < Path_H_IntPtoNeedle)
	ResidualEn = -5.0;
      else
	{
	  Efinal = Heavy_ion->GetLookupEnergy(E_atNeedleStart,(Path_H_Full-Path_H_IntPtoNeedle));
	  
	  if(Efinal <= 0.0)
	    ResidualEn = E_atNeedleStart;
	  else
	    ResidualEn = E_atNeedleStart - Efinal;
	}
    }


  /* if(E_FullPath < 0.0) // tested for E_FullPath < 0.0 (or ==-1000.0) and E_FullPath == 0.0 - seems fine
    {
      cout << " y_dist: " << y_dist << " E_FullPath: " << E_FullPath << " ResidualEn: " << ResidualEn << endl;
      cout << " ResidualEn: " << ResidualEn << " BeamE: " << FinalE << " BeamAtNeedle: " << BeamAtNeedle << " IntP: " << zr 
	 << " Theta: " << heavyTheta << " EnergyHeavy: " << heavyEn << " deltaBeam: " << deltaBeam << endl;
      cout << " E_atNeedleStart: " << E_atNeedleStart << " Efinal: " << Efinal 
	   <<  " Path_Full: " << Path_H_Full << " Path_IntPtoNeedle: " << Path_H_IntPtoNeedle << endl;
      cout << " /////////////////////////////////////////////////////// " << endl;
      }*/

  return ResidualEn;

}



///////////////////////////////------SX3 Coord------///////////////////////////////////////////////////////////////////////////////////////////////////////////

bool Simulator::GetSX3Coord(Float_t& X, Float_t& Y, Float_t& Z, Float_t& Xpc, Float_t& Ypc, Float_t& Zpc, 
 const Float_t Theta, const Float_t PHI, int ring)
{
 
  Float_t x1,x2,x3,y1,y2,y3,z1,z2,z3,t,p,q,r,L=100.0; 
  Float_t rho=0.0; //PC Radius 
  Float_t Xw=0.0, Yw=0.0, Zw=0.0;    
  bool GoodHit = 0;
   
 

  for (Int_t Det_ID=6; Det_ID<30; Det_ID++) {
    if (Det_ID == 17 || Det_ID == 25 || Det_ID == 26 || Det_ID == 27 || Det_ID == 28 || Det_ID == 29) continue;
    // mask the detectors that are not working : DetID= R1_2F, R2_2B,2C,2D,2E,2F (or DetID = 10,18,19,20,21,22 on Main_new.cpp)
    
    if(PHI < 0)
      cout << PHI*180/M_PI << " detector " << Det_ID << " " << endl;

    // To find the point (X,Y,Z) where the particle hit on the detector we need to find where 
    // the line defined by the particle's starting point (xr,yr,zr) and 
    // a point P(p,q,r) somewhere in space (L=100 cm) at the same direction as the direction of the particle
    // INTERSECT with the plane of the detector. 

       p= xr + (L*sin(Theta)*cos(PHI));
       q= yr + (L*sin(Theta)*sin(PHI));
       r= zr - (L*cos(Theta));     // assuming wherever you are cos(Theta) = (zr-r)/L // for Theta<90 cos>0 for Theta>90 cos<0
    
       // cout <<"p="<<p<<" q= "<<q<<" r= "<<r<<endl;

       // P1(x1,y1,z1), P2(x2,y2,z2), P3(x3,y3,z3) points to define the detector's plane

       x1= XAt0[Det_ID];
       y1= YAt0[Det_ID];
       z1= ZOffset[Det_ID];
       // cout <<"x1="<<x1<<" y1= "<<y1<<" z1= "<<z1<<endl;

       x2= XAt4[Det_ID]; 
       y2= YAt4[Det_ID];
       z2= ZOffset[Det_ID];
       // cout <<"x2="<<x2<<" y2= "<<y2<<" z2= "<<z2<<endl;

       x3=XAt0[Det_ID];
       y3=YAt0[Det_ID];
       z3=ZOffset[Det_ID]+7.5;
       // cout <<"x3="<<x3<<" y3= "<<y3<<" z3= "<<z3<<endl;

       // vectors for the plane of the detector P1P2, P1P3:
       TVector3 P1P2(0.0,0.0,0.0);
       TVector3 P1P3(0.0,0.0,0.0);
       TVector3 Pd(0.0,0.0,0.0); //point on the detector
       TVector3 P(0.0,0.0,0.0);  //point somewhere in space 
       TVector3 Pr(0.0,0.0,0.0); //reaction point of the particle
       TVector3 n(0.0,0.0,0.0);  //normal vector to the plane

       P1P2.SetXYZ(x2-x1,y2-y1,z2-z1); 
       P1P3.SetXYZ(x3-x1,y3-y1,z3-z1);
       n = P1P2.Cross(P1P3);
       Pd.SetXYZ(x1,y1,z1);
       P.SetXYZ(p,q,r);
       Pr.SetXYZ(xr,yr,zr);
       
       // Equation of line: P(s) = Pr + t*(P-Pr) 
       // @Intersection with plane the dot product of the normal vector and the vector between
       // the point Pd on the detector and the point of intersection is zero: n (dotproduct) (P(t_inter)-Pd) = 0. 

       t =  (n.Dot(Pd-Pr))/(n.Dot(P-Pr));
    
       //cout << " t: " << t << endl;
           
       //In order to Intersect 1ST: t needs to be greater than zero
       
       if(t>0.0)
	 {
	   // We substitute t in the equation of line with P and Pr points and find
	   // the point of the particle hit on the detector X,Y,Z:
	   X= xr+t*(p-xr);
	   Y= yr+t*(q-yr);
	   Z= zr+t*(r-zr);
	   
	   TVector3 P_hit(0.0,0.0,0.0);
	   P_hit.SetXYZ(X,Y,Z);
	   
	   //2ND: The vector on the plane of the detector (P_hit - Pd) has to be within the coordinate limits of the detector
	   // so 0<t1<1 AND 0<t2<1: with P1P2, P1P3 the vectors that defined before for the limits of the detector

	   Double_t t1 = ((P_hit - Pd).Dot(P1P2))/(P1P2.Mag2());
	   Double_t t2 = ((P_hit - Pd).Dot(P1P3))/(P1P3.Mag2());

	   if((t1>0 && t1<1) && (t2>0 && t2<1))
	     { 
	       GoodHit = 1;
	       /*cout << " t1: " << t1 << " t2: " << t2 << " t: " << t << " X: " << X << " Y: " << Y << " Z: " <<  Z
		    << " zr: " << zr << " Z_detector: " << (zr - 8.9/tan(Theta))
		    << " Det_ID: " << Det_ID << " Theta: " << Theta*180/TMath::Pi() << " Phi: " << PHI*180/TMath::Pi() << endl;*/
	     }
	 }
	 
	    
	
       //=============== Find-Assign a hit on a PCwire ================

      
       if(GoodHit == 1)
       	 {
	   Float_t phi_wire = 0.0;      
	   Float_t rho_offset = 0.0;     
	   Float_t phi_woffset = 0.0;

	   Float_t phi_diff = 0.0, phi_diff_mod = 0.0, min_diff = 380.0, min_rho=0.0 ;
	   Int_t min_index = -1;

    	   for(Int_t Wire_ID=0; Wire_ID<24; Wire_ID++)
	     {
	       // the phi angle of each wire made so that Wire0 is at 90 deg matching the analysis config
	       rho = PC_Radius[Wire_ID];	      
	       phi_wire = (24-Wire_ID)*2*TMath::Pi()/24 + TMath::Pi()/2;  
	       if(phi_wire >= 2*TMath::Pi())
		 phi_wire = phi_wire - 2*TMath::Pi();

	       Xw = rho*cos(phi_wire);   // using the physical radius rho of the detector for each wire's coordinates
	       Yw = rho*sin(phi_wire);

	       // need to find the "rho_offset" and "phi_woffset" if the particle is created at an offset (xr,yr) point 
	       // and not at the origin of the coordinate system
	       Float_t dXw = Xw-xr;
	       Float_t dYw = Yw-yr;
	       rho_offset = sqrt((dXw)*(dXw) + (dYw)*(dYw));  
	       phi_woffset = phi_from_cart(dXw, dYw);

	       // In order to assign a hit on a wire we check for the PHI of the particle that is closest to a Wire_Phi
	       // doing so by defining the following phi_difference and looking for the minimum value: 

	       phi_diff = abs(PHI-phi_woffset);
	       if(phi_diff > TMath::Pi() && phi_diff <= 2*TMath::Pi())
		 phi_diff = 2*TMath::Pi() - phi_diff;

	       phi_diff_mod = (fmodf((fabs(PHI-phi_woffset) + 2*TMath::Pi()), 2*TMath::Pi()));

	       if(phi_diff < min_diff)
		 {
		   min_diff = phi_diff;
		   min_index = Wire_ID;
		   min_rho = rho;
		   Xpc = Xw;
		   Ypc = Yw;
		   if(Theta == TMath::Pi()/2)
		     Zpc = zr;
		   else if(Theta == 0 || Theta == TMath::Pi())
		     Zpc = -10.0;
		   else
		     Zpc = zr - (rho_offset/tan(Theta));
		 }

	       /*cout << " wire: " << Wire_ID << " rho: " << rho << " phi_wire: " << phi_wire*180/TMath::Pi() 
		    << " xr: " << xr << " yr: " << yr
		    << " Xw: " << Xw << " Yw: " << Yw 
		    << " phi_woffset: " << phi_woffset*180/TMath::Pi() 
		    << " PHI: " << PHI*180/TMath::Pi() << " phi_diff: " << phi_diff*180/TMath::Pi() 
		    << " phi_diff_mod: " << phi_diff_mod*180/TMath::Pi() << endl;*/
	     	
	     } // end of 24-Wire loop 
    
	   /*
	   cout << " min_diff: " << min_diff*180/TMath::Pi() << " Wire_min: " << min_index << " PCR: " << min_rho
		  << " Xpc: " << Xpc << " Ypc: " << Ypc << " Zpc: " << Zpc
		  << " Det_ID: " << Det_ID 
		  << " X: " << X << " Y: " << Y << " Z: " <<  Z
		  << " zr: " << zr << " Z_detector: " << (zr - 8.9/tan(Theta))
		  << " Theta: " << Theta*180/TMath::Pi() << " Phi: " << PHI*180/TMath::Pi() << endl;*/

	 } // end of: if GoodHit==1 find the PC-hit
	    
	
       if(GoodHit)
	 return GoodHit;
  }// detector loop

  return GoodHit;
}

//////////////////////////////-------QQQ Coord-------/////////////////////////////////////////////////////////////////////////////////

bool Simulator::GetQQQ3FwdCoord(Float_t& X, Float_t& Y, Float_t& Z, Float_t& Xpc, Float_t& Ypc, Float_t& Zpc,Int_t& detid,
				const Float_t Theta,const Float_t PHI)
{

  Float_t L,rF,RF,g;
  Float_t x1,y1,z1,x2,y2,p,q,r;
  Float_t rho =0; //PC Radius
  Float_t Xw=0.0, Yw=0.0, Zw=0.0;
  Float_t phi1_particle = 0.0;
  Float_t phi2_particle = 0.0;

  RF=9.90,rF=5.01,g=.501; //outer radius, inner radius, gap between each QQQ            

  bool HitInDet=0;

  for (Int_t Det_ID=0; Det_ID<4; Det_ID++) {
   
   

    if(PHI < 0)
      cout << PHI*180/M_PI << " detector " << Det_ID << endl;

  
    if(Theta >= TMath::Pi()/2. || zr <= ZOffset[Det_ID]) // since we know that anything that can hit the QQQs is forward
      continue;
    else{
      L = (zr - ZOffset[Det_ID])/cos(Theta);   
      p= xr + (L*sin(Theta)*cos(PHI));
      q= yr + (L*sin(Theta)*sin(PHI));
      r= ZOffset[Det_ID];  
    }
    
    if(Theta >= TMath::Pi()/2. || zr <= ZOffset[Det_ID])
      cout <<"Theta= " << Theta*180/TMath::Pi() << " X= "<< p <<" Y= "<< q <<" Z= "<< r << " zr: " << zr << endl;

    // look at the qqq_hit.cpp file for the functions below:
    // qqq_hit(x,y) takes care of the hits that end up in the 4 gaps between the QQQs 
    // and the hits that are outside the r_min = rF = 5.01 cm and r_max = RF = 9.9 cm
    if(!qqq_hit(p,q))
      continue;
    
    x1= XAt0[Det_ID];
    y1= YAt0[Det_ID];
   
    x2= XAt4[Det_ID];
    y2= YAt4[Det_ID];

    Float_t dX1= x1-xr;
    Float_t dY1= y1-yr;

    Float_t dX2= x2-xr;
    Float_t dY2= y2-yr;

    phi1_particle = phi_from_cart(dX1,dY1);
    phi2_particle = phi_from_cart(dX2,dY2);

    // this function checks if which of the 4 QQQs the particle hit based on: 
    // the coordinates given for each quadrant and the particle's phi
    if(0>(detid=qqq_phi_hit(x1,y1,PHI,phi1_particle,phi2_particle))) // no detector found
      continue;
    //if (detid != Det_ID){
    //  cout <<" Inconsistent: x1="<<x1<<" y1= "<<y1<<" z1= "<<z1<<endl;
    //   cout <<" Inconsistent: x1="<<p<<" y1= "<<q<<" z1= "<<r<<endl;
    // }


    HitInDet = 1;
    
    X = p;
    Y = q;
    Z = r;
    
 
    //=============== Find-Assign a hit on a PCwire ================
    if(HitInDet==1)
      {
    	 
	Float_t phi_wire = 0.0;      
	Float_t rho_offset = 0.0;     
	Float_t phi_woffset = 0.0;

	Float_t phi_diff = 0.0, min_diff = 380.0, min_rho = 0.0;
	Int_t min_index = -1;
    
	for(Int_t Wire_ID=0; Wire_ID<24; Wire_ID++)
	  {
	    
	    // the phi angle of each wire made so that Wire0 is at 90 deg matching the analysis config
	    rho = PC_Radius[Wire_ID];
	    phi_wire = (24-Wire_ID)*2*TMath::Pi()/24 + TMath::Pi()/2;  
	    if(phi_wire > 2*TMath::Pi())
	      phi_wire = phi_wire - 2*TMath::Pi();

	    Xw = rho*cos(phi_wire);   // using the physical radius rho of the detector for each wire's coordinates
	    Yw = rho*sin(phi_wire);

	    // need to find the "rho_offset" and "phi_woffset" if the particle is created at an offset (xr,yr) point 
	    // and not at the origin of the coordinate system
	    Float_t dXw = Xw-xr;
	    Float_t dYw = Yw-yr;
	    rho_offset = sqrt((dXw)*(dXw) + (dYw)*(dYw));  
	    phi_woffset = phi_from_cart(dXw, dYw);

	    // In order to assign a hit on a wire we check for the PHI of the particle that is closest to a Wire_Phi
	    // doing so by defining the following phi_difference and looking for the minimum value: 
	    phi_diff = abs(PHI-phi_woffset); 
	    if(phi_diff > TMath::Pi() && phi_diff <= 2*TMath::Pi())
	      phi_diff = 2*TMath::Pi() - phi_diff;

	    if(phi_diff < min_diff)
	      {
		min_diff = phi_diff;
		min_index = Wire_ID;
		min_rho = rho;
		Xpc = Xw;
		Ypc = Yw;
		if(Theta == TMath::Pi()/2)
		  Zpc = zr;
		else if(Theta == 0 || Theta == TMath::Pi())
		  Zpc = -10.0;
		else
		  Zpc = zr - (rho_offset/tan(Theta));
	      }
	    /*cout << " wire: " << Wire_ID << " rho: " << rho << " phi_wire: " << phi_wire*180/TMath::Pi() 
		<< " xr: " << xr << " yr: " << yr
		<< " Xw: " << Xw << " Yw: " << Yw 
		<< " phi_woffset: " << phi_woffset*180/TMath::Pi() 
		<< " PHI: " << PHI*180/TMath::Pi() << " phi_diff: " << phi_diff*180/TMath::Pi() << endl;*/
	     	
	  } // end of 24-Wire (i) loop 

	/*if(Det_ID != detid)
	//if( PHI*180/TMath::Pi()>180.0  && PHI*180/TMath::Pi()<270.0 && !(X<0.0 && Y<0.0))
	  {  
	    cout << " min_diff: " << min_diff*180/TMath::Pi() << " Wire_min: " << min_index << " PCR: " << min_rho
		 << " Xpc: " << Xpc << " Ypc: " << Ypc << " Zpc: " << Zpc
		 << " Det_ID: " << Det_ID << " assigned DetID: " << detid 
		 << " X: " << X << " Y: " << Y << " Z: " << Z
		 << " zr: " << zr << " xr: " << xr 
		 << " y_app: " << (zr-ZOffset[Det_ID])*tan(Theta)<< " yr: " << yr 
		 << " Theta: " << Theta*180/TMath::Pi() << " Phi: " << PHI*180/TMath::Pi() <<  " phi1: " << phi1_particle*180/TMath::Pi() << " phi2: " << phi2_particle*180/TMath::Pi() << endl; 
	      }*/

      } // end if(HitInDet==1) for PCWire assignment 
    break;

  }// end Det_ID loop

  return HitInDet;  

}

/////////////////////////////--------TOF-----------///////////////////////////////////////////////////////////////////////
// use GetTimeOfFlight(InitialEnergy,Pathlength,StepSize) function from EnergyLoss class

bool Simulator::TimeOfFlight()
{

  TOF_rxn = IonInGas->GetTimeOfFlight(MaxBeamEnergy,dist,.1); // reaction point n so needs to start at "dist" 
  TOF_beam = IonInGas->GetTimeOfFlight(MaxBeamEnergy,La,.1); // beam at window so needs to start at "La"
  TOF_light = proton->GetTimeOfFlight(Ep_in,D,.1);
  TOF_heavy = Heavy_ion->GetTimeOfFlight(Eh_in,D,.1);
  TOF_p = TOF_rxn+TOF_light;
  TOF_h = TOF_rxn+TOF_heavy;

  return 1;
}

/////////////////////////////////-------Draw Lines--------///////////////////////////////////////////////////////////////////

bool Simulator::DrawLines(Int_t Det_ID, Int_t m, Float_t& Wx, Float_t& Wy, Float_t& Wz)
{  
  
  if (m==0){
    Wx = (XAt0[Det_ID]);
    Wy = (YAt0[Det_ID]);
    Wz = (ZOffset[Det_ID]);
  }

  else if (m==1){
    Wx = (XAt0[Det_ID]); 
    Wy = (YAt0[Det_ID]);
    Wz = (ZOffset[Det_ID]+7.5);
  }

  else if(m==2){
    Wx = (XAt4[Det_ID]);
    Wy = (YAt4[Det_ID]);
    Wz = (ZOffset[Det_ID]+7.5);
  }
    
  else if (m==3){
    Wx = (XAt4[Det_ID]);
    Wy = (YAt4[Det_ID]);
    Wz = (ZOffset[Det_ID]);
  }
    
  else if (m==4){
    Wx = (XAt0[Det_ID]);
    Wy = (YAt0[Det_ID]);
    Wz = (ZOffset[Det_ID]);
  }
  //cout << "Det_ID=" << Det_ID <<" m= "<<m<<" Wx="<<Wx<<" Wy="<<Wy<<" Wz="<<Wz<<endl;  
 
  return 1;
}


////////////////////////////--------Reconstruction Point-----------//////////////////////////////////////////////////////////////////////////////

bool Simulator::Simulation_Parameters(Float_t& z_SIM, Float_t& Path_sim, Float_t& EL_final,
				      const Float_t& Ep_kin, const Float_t& X, const Float_t& Y, const Float_t& Z)
{
    z_SIM = zr; // Simulation generated IntPoint

    //Need to calculate the actual path using the Si-coord (X,Y,Z) 
    //and the simulation IntP (xr,yr,zr)
    //since our Silicon Detectors can measure the Energy of the particle PRETTY ACCURATELY.
    Path_sim = sqrt(pow(X-xr,2) + pow(Y-yr,2) + pow(Z-zr,2));

    //The energy measured at the detector using the initial energy from kinematics generator
    //and by calculating the "real" path of the particle based on 
    //the geometrically generated parameters.
    //P.S. we have to assume that we know from SRIM the Eloss of the protons in our gas
    EL_final = proton->GetLookupEnergy(Ep_kin, Path_sim);
    
    return 1;
  
}

				      
//get the Reconstructed IntPoint and the Reconstructed parameters based on this IntPoint 
//NOT the one generated from the simulation
bool Simulator::Reac_point_Reconstruction(Float_t& x_br, Float_t& y_br, Float_t& z_br,
					  Float_t& Theta_Rec, Float_t& Path_Rec, Float_t& Phi_Rec,
					  const Float_t& X, const Float_t& Y, const Float_t& Z, 
					  const Float_t& Xa, const Float_t& Ya, const Float_t& Za)
{

  //(X,Y,Z) are the hits on the Outer Silicon detector.
  //(Xa,Ya,Za) can be EITHER hits on the Inner Silicon detector OR hits in the Proportional Counter.
  //depending on your choice
  //Your choice should be (if the Energy of the particle is High, Inner Silicon Detector hits
  //& if the Energy of thr particle is Low, Proportional Counter hits)
  
  Float_t t1=0.0,t2=0.0;
  Float_t rho=0.0, rho_a=0.0;
  
  if((X-Xa) != 0 || (Y-Ya) != 0){

    rho=sqrt(X*X + Y*Y);
    rho_a=sqrt(Xa*Xa + Ya*Ya);
    
    t1=-Xa/(X-Xa);
    t2=-Ya/(Y-Ya);
    
    x_br=Xa+t1*(X-Xa);
    y_br=Ya+t2*(Y-Ya);
    z_br=Za-(rho_a/(rho-rho_a))*(Z-Za);   // reconstructed IntPoint
    
  
    if ((z_br - Z) > 0){	 
      Theta_Rec = atan(rho/(z_br - Z));
      Path_Rec = rho/sin(Theta_Rec);

      if(Theta_Rec*180/TMath::Pi() > 90.0)
	cout << " Theta_Rec is MORE than 90: " << Theta_Rec*180/TMath::Pi() << endl;
    }
    else if((z_br - Z) < 0){	 
      Theta_Rec = TMath::Pi() + atan(rho/(z_br - Z));
      Path_Rec = rho/sin(Theta_Rec);

      if(Theta_Rec*180/TMath::Pi() < 90.0  || Theta_Rec*180/TMath::Pi() > 180.0)
	cout << " Theta_Rec is LESS than 90: " << Theta_Rec*180/TMath::Pi() << endl;
    }
    else{
      Theta_Rec = TMath::Pi()/2;
      Path_Rec = rho;
    }

    // get Phi from the Si-detector coord X,Y as in the experiment
    Phi_Rec = phi_from_cart(X,Y);
    
  }  
  
  /*if(Z > 2.0 && abs(z_br - z_SIM)>3.5)
    cout<<"  z_br= "<<z_br<<"  z_SIM= " << z_SIM << " Z= " << Z << " Zpc: " << Za << " Theta_zbr: " 
	<< atan(8.9/(z_br-Z))*180/TMath::Pi() << " Theta_zsim: " << atan(8.9/(z_SIM-Z))*180/TMath::Pi() 
	<< endl;*/
  
  return 1;  

}


bool Simulator::Kinematics_Rec(Float_t& Ex_Heavy_Rec, const Float_t& Ep_Rec,
			       const Float_t& Theta_Rec, const Float_t& Path_Rec, 
			       const Float_t& Phi_Rec, const Float_t& IntPoint)
{
  
  TLorentzVector Light_LV(0.,0.,0.,0.);  // proton 
  TLorentzVector Beam_LV(0.,0.,0.,0.);   // 18Ne 
  TLorentzVector Target_LV(0.,0.,0.,0.); // 4He
  TLorentzVector Parent_LV(0.,0.,0.,0.); // 22Mg 
  TLorentzVector Heavy_LV(0.,0.,0.,0.);  // 21Na 

  Float_t E_light_si = 0.0;  // regarding the light particle in general
  Float_t d_light = 0.0;
  Float_t E_light_rxn = 0.0;
  Float_t theta_light = 0.0;
  Float_t phi_light = 0.0;
  Float_t BeamEnergy = 0.0;
  
  Float_t M_Light = ML;  
  Float_t M_Beam = MB;
  Float_t M_Target = MT;
  Float_t M_Heavy = MH;

  //total energy of proton.
  Float_t E_light_tot =0.0;

  //Total momentum of proton.
  Float_t P_light =0.0;

  //Components of 4-vectors for protons.
  Float_t light_x =0.0;
  Float_t light_y =0.0;
  Float_t light_z =0.0;

  //Kinetic energy and Excitation energy of the 18Ne.

  Float_t Beam_E_tot =0.0;
  Float_t Beam_Pz =0.0;
  Float_t Target_E =0.0;	 

  // FOR PROTON OR IN GENERAL THE LIGHT OUTGOING PARTICLE
  //d_light, theta and phi based on the Reconstructed IntPoint not the one generated from simulation
  d_light = Path_Rec;  
  theta_light = Theta_Rec;
  phi_light = Phi_Rec;

  // E_Light_final as detected on the Si-det using the geometrical path created in the simulation
  E_light_si = Ep_Rec; 

  // proton energy at the reaction point  
  // Watch out!!! : the Initial Energy of the light particle at this point
  // IS NOT the one from the kinematics generator
  // NEEDS to be calculated based on:
  // a) the Final Energy that is detected on the Silicon (using the "real geometrical" path) 
  // b) the reconstructed Path of the particle from the Rec IntPoint (NOT the simulation generated one)
  E_light_rxn = proton->GetLookupEnergy(E_light_si,(-d_light));

  // same logic for BeamEnergy:
  // we have to use the Recons IntP to find the BeamEnergy
  // NOT use the FinalE from the simulation
  // because if the Rec IntP differs from the Sim IntP so will the BeamEnergy
  // keep in mind BeamEnergy is based on the SRIM Eloss 
  BeamEnergy = IonInGas->GetLookupEnergy(MaxBeamEnergy,(La-IntPoint));

  //total energy of protons.
  // the initial energy (@interaction point) and the rest mass
  E_light_tot = E_light_rxn + M_Light;

  //Total momentum of protons.
  P_light = sqrt(E_light_tot*E_light_tot - M_Light*M_Light);
  //Components of 4-vectors for proton.
  light_x = P_light*sin(theta_light)*cos(phi_light);
  light_y = P_light*sin(theta_light)*sin(phi_light);
  light_z = P_light*cos(theta_light);

  //Four vectors for proton.
  Light_LV.SetPxPyPzE(light_x,light_y,light_z,E_light_tot); 
  
  //Proton reconstruction
  Beam_E_tot = BeamEnergy + M_Beam;
  Beam_Pz = sqrt(Beam_E_tot*Beam_E_tot -M_Beam*M_Beam);
  Target_E = M_Target;

  //four vectors for beam & Target
  Beam_LV.SetPxPyPzE(0,0,Beam_Pz,Beam_E_tot);
  Target_LV.SetPxPyPzE(0,0,0,Target_E);
  
  //four vectors for parents 
  Parent_LV = Beam_LV + Target_LV;  // Parent(compound) = Beam + target (e.g. 22Mg = 18Ne + 4He)

  //four vectors for Sodium 21Na
  Heavy_LV = Parent_LV-Light_LV;        // Heavy outgoing = Parent(compound)-light (e.g. 21Na = 22Mg - p)

  Float_t Theta_Heavy_Rec = 0.0, Phi_Heavy_Rec=0.0;

  Theta_Heavy_Rec = Heavy_LV.Theta()*180/TMath::Pi();
  Phi_Heavy_Rec = Heavy_LV.Phi()*180/TMath::Pi();
  
  //Kinetic energy and Excitation energy of the Heavy Outgoing.
  Float_t KE_Heavy_Rec =-10.0;
  
  KE_Heavy_Rec = Heavy_LV.E()- Heavy_LV.M(); // .M()=sqrt(E^2-p^2) .E()=4th component of 4-vector, Energy
  Ex_Heavy_Rec = Heavy_LV.M() - M_Heavy;

  if((Parent_LV.M()-M_Light-M_Heavy-Ex_Heavy_Rec + BeamEnergy)<=0)
    cout << " Ex problem: " << Ex_Heavy_Rec << " Beam: " << BeamEnergy << " p_energy: " << ke_L_lab << endl;

  /*if(abs(Ex_Heavy_Rec - Ex) > 2.0)
    cout << " Ex_Heavy_Rec: " << Ex_Heavy_Rec << " BeamEnergy: " << BeamEnergy 
	 << " Theta: " << theta_light*180/TMath::Pi() << " p energy at det: " << E_light_si 
	 << " p at IntP: " << E_light_rxn << " IntP: " << IntPoint << endl;*/

   return 1;
  
  
}

bool Simulator::RecBeam_SecProton(Float_t& BeamRec_20Ne, const Float_t& Ep1_Rec, const Float_t& Ep2_Rec,
				  const Float_t& Theta1_Rec, const Float_t& Theta2_Rec,
				  const Float_t& Path1_Rec, const Float_t& Path2_Rec, 
				  const Float_t& Phi1_Rec, const Float_t& Phi2_Rec){

   Double_t E_p1_si = 0.0;
   Double_t E_p2_si = 0.0;  
   Double_t d_p1 = 0.0;
   Double_t d_p2 = 0.0;
   Double_t E_p1_rxn = 0.0;
   Double_t E_p2_rxn = 0.0;
   Double_t theta_p1 = 0.0;
   Double_t phi_p1 = 0.0;
   Double_t theta_p2 = 0.0;
   Double_t phi_p2 = 0.0;


  //total energy of proton.
  Double_t E_p1_tot =0.0;
  Double_t E_p2_tot =0.0;

  //Total momentum of proton.
  Double_t P_p1 =0.0;
  Double_t P_p2 =0.0;

  //Components of 4-vectors for protons.
  Double_t p1_x =0.0;
  Double_t p1_y =0.0;
  Double_t p1_z =0.0;
  Double_t p2_x =0.0;
  Double_t p2_y =0.0;
  Double_t p2_z =0.0;

   
  TVector3 p1_V(0.,0.,0.);  
  TVector3 p2_V(0.,0.,0.);  
  TVector3 SUM_p1p2_V(0.,0.,0.);  

  ///////////////////////////////////////////////////////
  E_p1_si = Ep1_Rec;
  E_p2_si = Ep2_Rec;

  d_p1 = Path1_Rec;
  d_p2 = Path2_Rec;

  theta_p1 = Theta1_Rec;
  theta_p2 = Theta2_Rec;

  phi_p1 = Phi1_Rec;
  phi_p2 = Phi2_Rec; 

  E_p1_rxn = proton->GetLookupEnergy(E_p1_si,(-d_p1));
  E_p2_rxn = proton->GetLookupEnergy(E_p2_si,(-d_p2));
   
  E_p1_tot = E_p1_rxn + ML;
  E_p2_tot = E_p2_rxn + ML;

 
  //Total momentum of protons.
  P_p1 = sqrt(E_p1_tot*E_p1_tot - ML*ML);
  //Components of 4-vectors for proton.
  p1_x = P_p1*sin(theta_p1)*cos(phi_p1);
  p1_y = P_p1*sin(theta_p1)*sin(phi_p1);
  p1_z = P_p1*cos(theta_p1);

  P_p2 = sqrt(E_p2_tot*E_p2_tot - ML*ML);
  //Components of 4-vectors for proton.
  p2_x = P_p2*sin(theta_p2)*cos(phi_p2);
  p2_y = P_p2*sin(theta_p2)*sin(phi_p2);
  p2_z = P_p2*cos(theta_p2);

  
  //Four vectors for proton.
  p1_V.SetXYZ(p1_x,p1_y,p1_z); 
  p2_V.SetXYZ(p2_x,p2_y,p2_z); 
  SUM_p1p2_V = p1_V + p2_V;


  //Double_t SUM_p1p2_Angle = p1_V.Angle(p2_V);  // angle between p1 and p2 vector
  Double_t SUM_p1p2_Theta = SUM_p1p2_V.Theta(); // angle between the sum of p1 and p2 with respect to the beam axis z 

  Double_t SUM_p1p2_Mag = SUM_p1p2_V.Mag();

  
  //cout << "SUM_p1p2_Theta " << SUM_p1p2_Theta*180/TMath::Pi() << " " << "SUM_p1p2_Mag " << SUM_p1p2_Mag << endl;
 
  ///////////////////////////////////////////////////
  
  Double_t ma = MB;
  Double_t mb = ML;
  Double_t mA = MT;
  Double_t mB = MH2;
  Double_t Q = MB+MT-(2*ML)-MH2;

  //cout << " Qvalue: " << Q << endl;
  
  Double_t cost = cos(SUM_p1p2_Theta);
  
  Double_t a = (ma/mB-1);
  Double_t b = -(sqrt(2*ma)/(mB))*SUM_p1p2_Mag*cost;
  Double_t c = E_p1_rxn + E_p2_rxn + (SUM_p1p2_Mag*SUM_p1p2_Mag)/(2*mB) - Q;
  
  Double_t s1 = (-b + sqrt(b*b-4*a*c))/2./a;  // sqrt(Ea) = sqrt(Beam_KE)
  Double_t s2 = (-b - sqrt(b*b-4*a*c))/2./a;

  //cout << " s1: " << s1 << " s2: " << s2 << endl;
  if(s1>s2)
    std::swap(s1,s2);

  BeamRec_20Ne = s2*s2;
 
   return 1;
 }


bool Simulator::Reconstruct_20Ne(Float_t& Ex_20Ne, const Float_t& Ep1_Rec, const Float_t& Ep2_Rec,
				  const Float_t& Theta1_Rec, const Float_t& Theta2_Rec,
				  const Float_t& Path1_Rec, const Float_t& Path2_Rec, 
				 const Float_t& Phi1_Rec, const Float_t& Phi2_Rec, const Float_t IntP){

   Double_t E_p1_si = 0.0;
   Double_t E_p2_si = 0.0;  
   Double_t d_p1 = 0.0;
   Double_t d_p2 = 0.0;
   Double_t E_p1_rxn = 0.0;
   Double_t E_p2_rxn = 0.0;
   Double_t theta_p1 = 0.0;
   Double_t phi_p1 = 0.0;
   Double_t theta_p2 = 0.0;
   Double_t phi_p2 = 0.0;


  //total energy of proton.
  Double_t E_p1_tot =0.0;
  Double_t E_p2_tot =0.0;

  //Total momentum of proton.
  Double_t P_p1 =0.0;
  Double_t P_p2 =0.0;

  //Components of 4-vectors for protons.
  Double_t p1_x =0.0;
  Double_t p1_y =0.0;
  Double_t p1_z =0.0;
  Double_t p2_x =0.0;
  Double_t p2_y =0.0;
  Double_t p2_z =0.0;
  
  Double_t Beam_E_tot =0.0;
  Double_t Beam_Pz =0.0;
  Double_t Target_E =0.0;	 


  TLorentzVector Beam_LV(0.,0.,0.,0.);  
  TLorentzVector Target_LV(0.,0.,0.,0.); 
  TLorentzVector Parent_LV(0.,0.,0.,0.); 

  TLorentzVector Heavy_LV(0.,0.,0.,0.); 
   
  TLorentzVector p1_LV(0.,0.,0.,0.);  
  TLorentzVector p2_LV(0.,0.,0.,0.);  
  TLorentzVector SUM_p1p2_LV(0.,0.,0.,0.);  

  Double_t Theta_20Ne = 0.0, Phi_20Ne = 0.0, KE_20Ne = 0.0;
  ///////////////////////////////////////////////////////
 
  E_p1_si = Ep1_Rec;
  E_p2_si = Ep2_Rec;

  d_p1 = Path1_Rec;
  d_p2 = Path2_Rec;

  theta_p1 = Theta1_Rec;
  theta_p2 = Theta2_Rec;

  phi_p1 = Phi1_Rec;
  phi_p2 = Phi2_Rec;
  
  E_p1_rxn = proton->GetLookupEnergy(E_p1_si,(-d_p1));
  E_p2_rxn = proton->GetLookupEnergy(E_p2_si,(-d_p2));
   
  E_p1_tot = E_p1_rxn + ML;
  E_p2_tot = E_p2_rxn + ML;

  //Total momentum of protons.
  P_p1 = sqrt(E_p1_tot*E_p1_tot - ML*ML);
  //Components of 4-vectors for proton.
  p1_x = P_p1*sin(theta_p1)*cos(phi_p1);
  p1_y = P_p1*sin(theta_p1)*sin(phi_p1);
  p1_z = P_p1*cos(theta_p1);

  P_p2 = sqrt(E_p2_tot*E_p2_tot - ML*ML);
  //Components of 4-vectors for proton.
  p2_x = P_p2*sin(theta_p2)*cos(phi_p2);
  p2_y = P_p2*sin(theta_p2)*sin(phi_p2);
  p2_z = P_p2*cos(theta_p2);

  
  //Four vectors for proton.
  p1_LV.SetPxPyPzE(p1_x,p1_y,p1_z,E_p1_tot); 
  p2_LV.SetPxPyPzE(p2_x,p2_y,p2_z,E_p2_tot); 
  SUM_p1p2_LV = p1_LV + p2_LV;


  Double_t BeamEnergy = IonInGas->GetLookupEnergy(MaxBeamEnergy,(La-IntP));
  Beam_E_tot = BeamEnergy + MB;
  Beam_Pz = sqrt(Beam_E_tot*Beam_E_tot - MB*MB);
  Target_E = MT;

  //four vectors for beam & Target
  Beam_LV.SetPxPyPzE(0,0,Beam_Pz,Beam_E_tot);
  Target_LV.SetPxPyPzE(0,0,0,Target_E);
  
  //four vectors for parents 
  Parent_LV = Beam_LV + Target_LV;  // Parent(compound) = Beam + target (e.g. 22Mg = 18Ne + 4He)

  //four vectors for Sodium 21Na
  Heavy_LV = Parent_LV-p1_LV-p2_LV;        // Heavy outgoing = Parent(compound)-light (e.g. 20Ne = 22Mg - 2p)

 
  Theta_20Ne = Heavy_LV.Theta()*180/TMath::Pi();
  Phi_20Ne = Heavy_LV.Phi()*180/TMath::Pi();
  
  //Kinetic energy and Excitation energy of the Heavy Outgoing.

  KE_20Ne = Heavy_LV.E()- Heavy_LV.M(); // .M()=sqrt(E^2-p^2) .E()=4th component of 4-vector, Energy
  Ex_20Ne = Heavy_LV.M() - MH2;

 
  ///////////////////////////////////////////////////

  
   return 1;
 }

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
