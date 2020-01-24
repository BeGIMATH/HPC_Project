//Include all the right things...
#include <algorithm> 
#include <Eigen/Core>
#include <vector>

using namespace Eigen;

/*Implement*/
void main(){
    // Initialize n*n random matrix
    // Call function to compute K_max
}

void 2DKmax(MatrixXd M, int K){
    // Set of all intervals (s,t) in increasing order
    MatrixXd sortedPairs = findSortedCombinations(n);
    MatrixXd F_w = findLeftConvex(M,sortedPairs);
    MatrixXd F_n = findRightConvex(M,sortedPairs);
    // Finilization: adding left and right convex shapes
    MatrixXd F_wn;
    for (int k=1;k<m;k++){
        for(int s=0;s<n;s++){
            for(int t=0;s<n;t++){
                F_wn(k,s,t) = F_w(k,s,t) + F_n(k,s,t) - summarize(M,k,s,t);
            }
        }
    }
}

/* Finds all possible F_w*/
MatrixXd findLeftConvex(MatrixXd M, MatrixXd sortedPairs){
    int m = sizeof(M,1);
    int n = sizeof(M,2);
    double F_w[m][n][n][K] = {0}; //=-infty
    for(int k=1;k<m;k++){
        // Find all intervals s-t in increasing order
        for(int i=0;i<num_of_comb;i++){
            int s = sortedPairs[i][1];
            int t = sortedPairs[i][2];
            F_w[k,s,t] = max(F_W(k-1,s,t)+summarize(M,k,s,t),
                             F_w(k,s+1,t) + M(s,k),
                             F_w(k,s,t+1) + a[t,k]);
        }
    }
}

/* Finds all possible F_n,
    might be unnecessary, could use the same function for F_n and F_w */
MatrixXd findRightConvex(MatrixXd M,MatrixXd sortedPairs){
    int m = sizeof(M,1);
    int n = sizeof(M,2)
    F_w = [m][n][n][K] = {0}; //=-infty
    for(int k=m-2;k>=0;k--){
        // Find all intervals s-t in increasing order
        for(int i=0;i<num_of_comb;i++){
            int s = sortedPairs[i][1];
            int t = sortedPairs[i][2];
            F_w[k,s,t] = max(F_W(k+1,s,t)+summarize(M,k,s,t),
                             F_w(k,s+1,t) + M(s,k),
                             F_w(k,s,t+1) + a[t,k]);
        }
    }
}

/**
 * Places all intervals of (s,t) in in a 2D matrix,
 * for all s<t
 * comb(i,j) = (s,t)
 **/
int** findSortedCombinations(n){
    int num_of_comb = n*(n+1)/2;
    MatrixXd comb(num_of_comb,2);
    int counter = 0;
    for(diff=0;diff<n-1;diff++){
        for(t=diff;t<n;t++){
            s = t-diff;
            comb(counter,1) = s;
            comb(counter,2) = t;
            counter++;

        }
    }  
    return comb;
}

/**
 * Summarizes column k of matrix M,
 * from s to t.
 **/
int summarize(MatrixXd M,int k,int s,int t){
    for(i=s;i<=t;i++){
        double sum = sum + M(i,k);
    }
}

/*From K-tuple L_1,...,L_m max(L_1,...) returns a tuple with the K largest elements*/
MatrixXd max(MatrixXd L1, MatrixXd L2, MatrixXd L3){
    
}


