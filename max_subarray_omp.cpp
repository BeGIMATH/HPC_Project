#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <omp.h>

using namespace Eigen;

void kadane(const VectorXd &array, double &maxSum, int &l, int &r)
{
  maxSum = -std::numeric_limits<double>::infinity();
  l = 0;
  r = 0;
  double sum = 0;
  int currentStartIndex = 0;
  for (int i = 0; i < array.size(); ++i)
  {
    sum += array(i);
    if (sum > maxSum)
    {
      maxSum = sum;
      l = currentStartIndex;
      r = i;
    }
    if (sum < 0)
    {
      sum = 0;
      currentStartIndex = i + 1;
    }
  }
}

void maxSubarray2D(const MatrixXd &array,
                   double &maxSum, int &left, int &right, int &top, int &bottom)
{

  maxSum = -std::numeric_limits<double>::infinity();
  left = -1;
  right = -1;
  top = -1;
  bottom = -1;

#pragma omp parallel
  {

    /// Create local variables for kadane 2D
    double l_maxSum = -std::numeric_limits<double>::infinity();
    int l_left = -1;
    int l_right = -1;
    int l_top = -1;
    int l_bottom = -1;
    double l_sum = 0;
    int l_start, l_finish;

    /// Numbers of elements to compute per precessors
    int size = (array.rows() + 1) * array.rows() / (2 * omp_get_num_threads());

    /// Define bounds of iteration

    int lo, hi;
    int id = omp_get_thread_num();

    for (int l = 0; l < 1; l++)
    {
      lo = (l == 0) ? id * size : array.rows() - (id + 1) * size;
      hi = (l == 0) ? (id + 1) * size : array.rows() - id * size;

      /// Loop over rows of 2D matrix
      for (int i = lo; i < hi; ++i)
      {
        /// Loop over column of 2D matrix
        VectorXd temp = VectorXd::Zero(array.cols());
        for (int j = i; j < array.rows(); ++j)
        {
          for (int k = 0; k < array.cols(); ++k)
          {
            temp(k) += array(j, k);
          }
          kadane(temp, l_sum, l_start, l_finish);
          if (l_sum > l_maxSum)
          {
            l_maxSum = l_sum;
            l_left = i;
            l_right = j;
            l_top = l_start;
            l_bottom = l_finish;
          }
        }
      }
    }

    /// Perform reduction of local information to global information
#pragma omp critical
    {
      if (l_maxSum > maxSum)
      {
        maxSum = l_maxSum;
        left = l_left;
        right = l_right;
        top = l_top;
        bottom = l_bottom;
      }
    }
  }
}

int main()
{
  /// Size of the matrix
  int n = 5;

  /// nxn Matrix filled with random numbers between (-1,1)
  MatrixXd m = MatrixXd::Random(n, n);
  auto start = std::chrono::system_clock::now();
  double maxSum;
  int left, right, top, bottom;

  maxSubarray2D(m, maxSum, left, right, top, bottom);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Final " << maxSum << " "
            << " " << left << " " << right << "  " << top << " " << bottom << std::endl;
}
