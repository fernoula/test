#ifndef __cross__section__cpp__
#define __cross__section__cpp__


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

  
  //function Thickness is made specifically in 
  template<class T>
  Double_t Thickness(T* S)
  {
    ofstream out("cross_product_thick.txt");
    
    Int_t length = Energy.size();
    //cout << length << endl;
    
    vector <Double_t> thickness(length,0.0);
    vector <Double_t> product(length,0.0);

    vector <Double_t> Enlab(length,0.0);
    
    Double_t maxProduct = 0.0;
    
    for(Int_t i=0; i<length; i++)
      {
	Enlab[i] = Energy[i]*22./4;
      }

    
    for(Int_t i=0; i<length-1; i++)
      {
	thickness[i] = S->IonInGas->GetDistance_new(Enlab[i+1],Enlab[i],0.001);
	product[i] = Cross[i]*thickness[i];
	

	//cout << " Energy: " << Enlab[i] << " Energy[i+1]: " << Enlab[i+1]
	//     << " dx: " << thickness[i] << " cross: " << Cross[i] << " product: " << product[i] << endl;

	out << " Energy: " << Energy[i] 
	     << " dx: " << thickness[i] << " cross: " << Cross[i] << " product: " << product[i] << endl;
	
      }

    auto it = std::max_element(product.begin(),product.end());
    maxProduct = *it;
    Int_t index = std::distance(product.begin(),it);
    
    //cout << " maxProduct: " << maxProduct << " index: " << index << endl;
    
    return maxProduct;
   
    }
 
  
  template<class T> 
  Double_t Delta_x(Double_t En, T* S)
  {
    En = En*(4./22);

    Double_t dx = 0.0;
    Double_t FinalE = 0.0, InitialE = 0.0;

    // Res: corresponds to resolution at CM 0.083MeV (120,0,10), so since in the file we have
    // Energy gaps of 1keV then to create the 0.083MeV resolution we need steps of
    // (left and right of given En, as in analysis) 0.083/2 = 0.0415 MeV = 41.5 keV or 42 steps
    //Int_t Res = 41;

    Int_t Res = 1;  // corresponds to 1keV gap between the energies at the txt file
    
    auto it = std::lower_bound(Energy.begin(),Energy.end(),En);

    if(En > Energy[Energy.size()-1])
      return 0;
    
    Int_t index = std::distance(Energy.begin(),it);
    //by choosing the high_index and low_index the following way we have the same energy gap
    //for the thickness calculation, as the function Thickness above, so the products cross*dx
    //for the final cross section calculation can be comparable.
    
    Int_t high_index = index;      
    Int_t low_index = index - Res;

    if(index > 13048) high_index = 13050;
    if(index < 2) low_index = 1;
      
    FinalE = Energy[low_index]*22./4;
    InitialE = Energy[high_index]*22./4;
    
    dx = S->IonInGas->GetDistance_new(InitialE,FinalE,0.001);

    /*cout << " En: " << En << " Energy_high: " << Energy[high_index] << " " << Cross[high_index]
	 << " index: " << index
	 << " Energy_low: " << Energy[low_index] << " " << Cross[low_index]
	 << " low_index: " << low_index << " high_index: " << high_index
	 << " FinalE: " << FinalE << " InitialE: " << InitialE << " dx: " << dx <<  endl;*/
    
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

    cross = Cross[index-1];
    //Nout = (dx*cross)/(dx_max*maxCrossSection);

    Double_t maxNorm = Thickness(S);
    Nout = (dx*cross)/maxNorm;

 
    if(Nout >= 1.0)
    cout << " Nout: " << Nout << " maxNorm: " << maxNorm << " En: " << En << " Energy: " << Energy[index]
	 << " Cross[index]:  " << Cross[index] << " Cross[index-1]: " << cross << " dx: " << dx
	 << " maxCSEnergy: " << maxCSEnergy << " maxCS: " << maxCrossSection << " dx_max: " << dx_max << endl;
    
    return Nout;
  }

  
};

#endif
