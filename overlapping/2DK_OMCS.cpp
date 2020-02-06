//Include all the right things...
#include <algorithm> 
#include <vector>
#include <limits>
#include <iostream>
using namespace std;

void TwoDKmax(vector<double> M, int K,int rows, int col){
    // n is the number of columns??
    //int m = rows;
    //int n = col;
    // Set of all intervals (s,t) in increasing order
    vector<int> sortedPairs = findSortedCombinations(rows);
    // Last input == 1 for left, ==0 for right
    vector<double> F_w = findConvex(M,sortedPairs,K,1,rows,col);
    vector<double> F_n = findConvex(M,sortedPairs,K,0,rows,col);
    
    // Finilization: adding left and right convex shapes
    vector<double> F_wn;
    for (int k=1;k<col;k++){
        for(int s=0;s<rows;s++){
            for(int t=0;t<rows;t++){
                // set a range of values to a 
                int start = k*rows*rows*K + s*rows*K + t*K;
                int end = k*rows*rows*K+s*rows*K + t*K + K;

                //Herrejävlar, måste strukturera upp det här
                vector<double> sub1(&F_w[start],&F_w[end]);
                vector<double> sub2(&F_n[start],&F_n[end]);
                vector<double> temp = add(sub1,sub2);
                int tempint = summarize(M,k,s,t,rows);
                addInt(temp,summarize(M,k,s,t,rows));

                for(int j = 0;j<K;j++){
                    F_w[k*rows*rows*K + s*rows*K + t*K + j] = temp[j];
                }
                //std::fill(F_wn.begin() + start, F_wn.begin() + start + end,
                          //addInt(add(sub1,sub2),summarize(M,k,s,t)));
            }
        }
    }
}

/* Finds all possible F_w*/
std::vector<double> findConvex(std::vector<double> M, std::vector<int> sortedPairs, int K, int type,int rows,int col){
    //int m = sizeof(M,1); //N.of. rows
    //int n = sizeof(M,2); //N.of columns
    int num_of_comb = numberOfCombinations(rows);
    std::vector<double> F_w(col*rows*rows*K,-std::numeric_limits<double>::infinity());

    // set all k=0 values to 0 for sets s<= t
    for (int i = 0; i<num_of_comb;i++){
        int s = sortedPairs[i*2];
        int t = sortedPairs[i*2+1];
        F_w[0,s,t,0] = 0;
    }
    
    // TODO
    vector<int> k_range(K);
    if (type == 1){ //left
        int x = 0;
        std::generate(k_range.begin(), k_range.end(), [&]{ return x++; }); 
    }
    else if (type == 0){
        int x = K-1;
        std::generate(k_range.begin(), k_range.end(), [&]{ return x--; }); 
    }
    
    for(int k : k_range){
        for(int i=0;i<num_of_comb;i++){
            int s = sortedPairs[i*2];
            int t = sortedPairs[i*2+1];

            int start = k*rows*rows*K + s*rows*K + t*K;
            int end = k*rows*rows*K + s*rows*K + t*K + K;

            int prev = (type==1) ? (k-1) : (k+1);
            // Subvectors of F_w corresponding to the three cases
            std::vector<double> L1(&F_w[prev*rows*rows*K + s*rows*K + t*K],
                                   &F_w[prev*rows*rows*K + s*rows*K + t*K + K]);
            addInt(L1,summarize(M,k,s,t,rows));
            std::vector<double> L2(&F_w[k*rows*rows*K + (s+1)*rows*K + t*K],
                                   &F_w[k*rows*rows*K + (s+1)*rows*K + t*K + K]);
            addInt(L2,M[s,k]);
            std::vector<double> L3(&F_w[k*rows*rows*K + s*rows*K + (t+1)*K],
                                   &F_w[k*rows*rows*K + s*rows*K + (t+1)*K + K]);
            addInt(L3,M[t,k]);
            vector<double> L_max = max(L1,L2,L3);
            for(int j = 0;j<K;j++){
                F_w[k*rows*rows*K + s*rows*K + t*K + j] = L_max[j];
            }
            //std::fill(F_w.begin() + start, F_w.begin() + start + end, max(L1,L2,L3)); 
        }
    }
    return F_w;
}


int numberOfCombinations(int m){
    int num_of_comb = m*(m+1)/2;
}

/**
 * Places all intervals of (s,t) in in a 2D matrix,
 * for all s<t
 * comb(i,j) = (s,t)
 **/
std::vector<int> findSortedCombinations(int m){
    int num_of_comb = numberOfCombinations(m);
    // 2D vector of size num_of_comb*2 with initial value 0
    //std::vector<std::vector<int> > comb(num_of_comb,std::vector<int>(2, 0));
    // Each set of two represents a combination of (s,t) for a given iteration (counter)
    //int comb[num_of_comb*2];
    std::vector<int> comb(num_of_comb*2);
    int counter = 0;
    for(int diff=0;diff<m-1;diff++){
        for(int t=diff;t<m;t++){
            int s = t-diff;
            comb[counter*2] = s;
            comb[counter*2+1] = t;
            counter++;
        }
    }  
    return comb;
}

/**
 * Summarizes column k of matrix M,
 * from s to t.
 **/
int summarize(vector<double> M,int k,int s,int t, int rows){
    
    for(int i=s;i<=t;i++){
        double sum = sum + M[i*rows + k];
    }
}

/*From K-tuple L_1,...,L_m max(L_1,...) returns a tuple with the K largest elements*/
vector<double> max(vector<double> L1, vector<double> L2, vector<double> L3){
    // merge matrices and perform quicksort?
    return;
}

void addInt(vector<double> &v, int a){
    for (int i = 0;i<v.size();i++){
        v[i] = v[i]+a;
    }
}

vector<double> add(vector<double> v1, vector<double> v2){
    vector<double> v3;
    for(int i = 0; i<v1.size();i++){
        v3.push_back(v1[i] + v2[i]);
    }
    return v3;
}

/*Implement*/
int main(){
    // Initialize n*n random matrix
    // Call function to compute K_max

    double size_rows = 100;
    double size_columns = 100;
    vector<double> M;


    return 0;
}
