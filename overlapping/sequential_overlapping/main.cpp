/*
  to compile the program
  g++ main.cpp -o 2DK_OMCS
  to execute it
  ./2DK_OMCS
*/

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include "2DK_OMCS.cpp"

using namespace std;

/*Implement*/
int main(){
  // Initialize n*n random matrix

  // Call function to compute K_max

  double size_rows = 5;
  double size_columns = 5;
  std::vector<double> M ;
  for (int i=0; i<size_rows*size_columns; i++) {
    M.push_back(i);
  }

  for (std::vector<double>::iterator it=M.begin(); it!=M.end(); ++it)
    std::cout << ' ' << *it;
  std::cout << '\n';

  std::cout<<"The solution"<<std::endl;

  std::vector<double> Sol = doublewoDKmax(M, 1, size_rows, size_columns);

  for (std::vector<double>::iterator it=Sol.begin(); it!=Sol.end(); ++it)
    std::cout << ' ' << *it;
  std::cout << '\n';

  return 0;
}
