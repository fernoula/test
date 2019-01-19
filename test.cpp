#ifndef __cross__section__cpp__
#define __cross__section__cpp__

///--------not use anymore------------

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>


using namespace std;

class CS
{
public:

  Int_t val=0;
  vector <Double_t> Energy;
  vector <Double_t> Cross;

  Double_t maxCrossSection;
  Double_t maxCSEnergy;

  
  CS()
  {
    ;
  }


  Int_t ReadFile(const char* CrossSectionFile)
  {
    ifstream cross_file;
    string line;
    Int_t line_counter;
   

    cross_file.open(CrossSectionFile);
    if(!cross_file.is_open())
      {
	cout << " cross section file: " << CrossSectionFile << " is NOT open " << endl;
      }
    else
      {
	cout << " cross section file: " << CrossSectionFile << " is being read " << endl;
	line_counter = 0;
	do
	  {
	    getline(cross_file, line);
	    if(!line.empty())
	      line_counter++;
	  } while(!cross_file.eof());
	cross_file.close();

	val = line_counter-1;
	Energy.resize(val,0);
	Cross.resize(val,0);

	cout << val << endl;

	cross_file.open(CrossSectionFile);
	getline(cross_file,line);
	cout << line << endl;

	for(Int_t i=0; i<val; i++)
	  {
	    cross_file >> Energy[i] >> Cross[i];
	    //cout << Energy [i] << " " << Cross[i] << endl;
	  }

	// std::vector<Double_t>::iterator it
	auto it = std::max_element(Cross.begin(),Cross.end());
	maxCrossSection = *it; 

	Int_t index = std::distance(Cross.begin(),it);
	maxCSEnergy = Energy[index];

	//cout << " index: " << index << " maxCrossSection: " << maxCrossSection << " maxCSEnergy: " << maxCSEnergy << endl;
	
      }

    return 1;
  }


  
  template<class T> 
  Double_t Delta_x(Double_t En, T* S)
  {
    En = En*(4./22);

    Double_t dx = 0.0;
    Double_t FinalE = 0.0, InitialE = 0.0;

    auto it = std::lower_bound(Energy.begin(),Energy.end(),En);

    if(En > Energy[Energy.size()-1])
      return 0;
    
    Int_t index = std::distance(Energy.begin(),it);
    //Int_t low_index = index - 1;

    Int_t high_index = index + 41;
    Int_t low_index = index - 41;

    if(index > 3859) high_index = 3901;
    if(index < 42) low_index = 1;

    /*cout << " En: " << En << " Energy_high: " << Energy[high_index] << " " << Cross[high_index] << " index: " << index
	 << " Energy_low: " << Energy[low_index] << " " << Cross[low_index]
	 << " low_index: " << low_index << " high_index: " << high_index << endl;*/
    
    FinalE = Energy[low_index]*22./4;
    InitialE = Energy[high_index]*22./4;
    
    dx = S->IonInGas->GetDistance_new(InitialE,FinalE,0.001);

    //cout << " FinalE: " << FinalE << " InitialE: " << InitialE << " dx: " << dx
    //	 << " S->FinalE: " << S->FinalE << endl;
	
    
    return dx;   
  }


  template<class T> 
  Double_t CrossSection(Double_t En, T* S)
  {    

    Double_t density = (1.1617e-4)*(6.022e23)*0.96/4.0;
    Double_t Nin = 9.18e8;

    Double_t Nout = 0.0;
    Double_t cross = 0.0;
    
    Double_t dx = 0.0;
    Double_t dx_max = 0.0;
    dx_max = Delta_x<T>(maxCSEnergy*22./4, S);

    //cout << " dx_max: " << dx_max << endl;
    
    dx = Delta_x<T>(En, S);
    
    En = En*(4./22);
    
    // returns an iterator pointing to the first element in the range
    // [first,last) which doesn't compare less than val
    auto it = std::lower_bound(Energy.begin(),Energy.end(),En);

    if(En > Energy[Energy.size()-1])
      return 0;
    
    Int_t index = std::distance(Energy.begin(),it);

    //cout << " S->FinalE: " << S->FinalE << " En: " << En << " Energy: " << Energy[index] << " "
    //	 << Cross[index] << " dx: " << dx << " dx_max: " << dx_max << endl;
    
    //cross = Cross[index]/maxCrossSection;
    //return cross;
    
    
    cross = Cross[index];
    Nout = (dx*cross)/(dx_max*maxCrossSection);

    /*if(Nout == 1.0)
      cout << " S->FinalE: " << S->FinalE << " En: " << En << " Energy: " << Energy[index] << " "
	   << Cross[index] << " dx: " << dx
    	   << " maxCSEnergy: " << maxCSEnergy << " maxCS: " << maxCrossSection << " dx_max: " << dx_max << endl;*/
    
    return Nout;
  }

  
};

#endif
