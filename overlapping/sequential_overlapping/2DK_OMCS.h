#ifndef OMCS
#define OMCS
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>

using namespace std;


std::vector<double> TwoDKmax(vector<double> M, int K,int rows, int col);

/* Finds all possible F_w*/
std::vector<double> findConvex(std::vector<double> M, std::vector<int> sortedPairs, int K, int type,int rows,int col);

int numberOfCombinations(int m);

/**
 * Places all intervals of (s,t) in in a 2D matrix,
 * for all s<t
 * comb(i,j) = (s,t)
 **/
std::vector<int> findSortedCombinations(int m);

int summarize(vector<double> M,int k,int s,int t, int rows);

/*From K-tuple L_1,...,L_m max(L_1,...) returns a tuple with the K largest elements*/
vector<double> max(vector<double> L1, vector<double> L2, vector<double> L3);

void addInt(vector<double> &v, int a);

vector<double> add(vector<double> v1, vector<double> v2);



#endif
