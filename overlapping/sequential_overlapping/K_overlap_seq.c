// program to find the k overlapping maximum subarrays in sequential
/*
  * compile : gcc -o K_overlap_seq K_overlap_seq.c
  * execute : ./K_overlap_seq
*/
#include <stdio.h>
#include <math.h>



// void one_2D_subarray(int m, int n, float* a){
//   /*
//     Corresponds to Algorithm 8 Alternative 2D version of
//     prefix sum Algorithm 3
//   */
//   //initialization
//   float M = -INFINITY;
//   float min = -INFINITY;
//   float s[n];
//
//   for (int g = 0; g<m; g++){
//
//   }
//
//   for(int g = 0; g<n; g++){
//     printf("%f\n", s[n]);
//   }
//
// }


void main(){
  //initialization
  int K = 4; // number of subarrays
  int n = 5; // number of rows
  int m = 4; // number of lines

  float a[m][n] = {{1, 2, -1, -4, -20},
                     {-8, -3, 4, 2, 1},
                     {3, 8, 10, 1, 3},
                     {-4, -1, 1, 7, -6}
                    };
  float M[K]; //The resulting K sums

  for (int i = 0; i<K; i++){
    M[i] = -INFINITY;

  }

  //one_2D_subarray(m, n, a);

}
