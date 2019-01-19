#ifndef __STRAGGLING_CPP__
#define __STRAGGLING_CPP__

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

class Strag
{
public:

  Int_t val=0;
  vector <Double_t> Energy;
  vector <Double_t> de_dx_el;
  vector <Double_t> de_dx_nuc;
  vector <Double_t> range;
  vector <Double_t> long_strag;
  vector <Double_t> lateral_strag;

 
  
  Strag()
  {
    ;
  }


  Int_t ReadFile(const char* ELossFile)
  {
    ifstream eloss_file;
    string line;
    Int_t line_counter;
   

    eloss_file.open(ELossFile);
    if(!eloss_file.is_open())
      {
	cout << " eloss file: " << ELossFile << " is NOT open " << endl;
      }
    else
      {
	cout << " eloss file: " << ELossFile << " is being read " << endl;
	line_counter = 0;
	do
	  {
	    getline(eloss_file, line);
	    if(!line.empty())
	      line_counter++;
	  } while(!eloss_file.eof());
	eloss_file.close();

	val = line_counter-1;
	Energy.resize(val,0);
	de_dx_el.resize(val,0);
	de_dx_nuc.resize(val,0);
	range.resize(val,0);
	long_strag.resize(val,0);
	lateral_strag.resize(val,0);
	

	cout << val << endl;

	eloss_file.open(ELossFile);
	getline(eloss_file,line);
	cout << line << endl;

	for(Int_t i=0; i<val; i++)
	  {
	    eloss_file >> Energy[i] >> de_dx_el[i] >> de_dx_nuc[i]
		       >> range[i] >> long_strag[i] >> lateral_strag[i];
	    
	    //cout << Energy [i] << " " << lateral_strag[i] << endl;
	  }

	// std::vector<Double_t>::iterator it
	//auto it = std::max_element(Cross.begin(),Cross.end());
	//maxCrossSection = *it; 
	
      }

    return 1;
  }

  
  void Range_Straggling(Double_t dist_z, Double_t& y_strag, Double_t& z_strag)
  {    
    // convert distance input in cm from simulation to mm to compare to srim file's range
    dist_z = dist_z*10.0; 
    
    // returns an iterator pointing to the first element in the range
    // [first,last) which doesn't compare less than val
    auto it = std::lower_bound(range.begin(),range.end(),dist_z);

    if(dist_z > range[range.size()-1])
      cout << " distance out of range !! " <<  endl;
    
    Int_t index = std::distance(range.begin(),it);

    y_strag = lateral_strag[index];
    z_strag = long_strag[index];

    /*cout << " dist: " << dist_z << " mm " << " lateral_straggling: " << y_strag << " mm "
	 << " long_strag: " << z_strag << " mm " 
	 << " range from file: " << range[index] << " mm " << endl;*/
       
  }

  
};

#endif

