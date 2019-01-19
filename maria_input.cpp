#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

int main() {
  std::ofstream out("cross_high2.txt");

  out << "En (MeV)" << "\t" << "cross (mb)" << std::endl;

  int n=9150;
  vector <double> Energy(n,4.0);
  vector <double> cross(n,0);

  for(int i=0; i<n; i++)
    {
      Energy[i+1] = Energy[i]+0.001;
      cross[i] = 110;
	
      out << Energy[i] << "\t" << cross[i] << endl;     
    }

  cout << "ENERGY.SIZE(): " << Energy.size() << endl;

  
  return 0;
}
