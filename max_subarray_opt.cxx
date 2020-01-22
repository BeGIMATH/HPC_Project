#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <omp.h>

using namespace Eigen;

void prefix_sum(VectorXd &in, int nt, VectorXd &buffd)
{
  {
    // Get current thread id
    int tid = omp_get_thread_num();

    // Compute workload for each thread
    int work = (in.size() + nt - 1) / nt;

    // Set upper and lower bounds for each threads
    int lo = work * tid;
    int hi = std::min(lo + work, int(in.size()));

    // Compute local prefix sum for each thread
    for (int i = lo + 1; i < hi; i++)
      in[i] = in[i] + in[i - 1];

    // Compute delta for each thread
    buffd[tid] = in[hi - 1];

    // Compute global prefix sum
#pragma omp barrier
    for (int j = 1; j < nt; j <<= 1)
    {
      if (tid >= j)
        buffd[tid] += buffd[tid - j];
#pragma omp barrier
    }

    // Reset prefix sum
    buffd[tid] -= in[hi - 1];

    // Update local prefix sum
    for (int i = lo; i < hi; i++)
      in[i] = in[i] + buffd[tid];
  }
}

void maxloc_omp(VectorXd &in, int nt, VectorXd &maxval, VectorXi &maxloc)
{
  int mloc, loc_m;
  double mval, val_m;
  int n = in.size();
  {
    int id = omp_get_thread_num();
    maxval[id] = -1.0e30;
#pragma omp for
    for (int i = 0; i < n; i++)
    {
      if ((in[i] - maxval[id]) > 1.e-14)
      {
        maxloc[id] = i;
        maxval[id] = in[i];
      }
    }
#pragma omp flush(maxloc, maxval)
#pragma omp master
    {
      int nt = omp_get_num_threads();
      mloc = maxloc[0];
      mval = maxval[0];
      for (int i = 1; i < nt; i++)
      {
        if (maxval[i] > mval)
        {
          mval = maxval[i];
          mloc = maxloc[i];
        }
      }
    }
    maxval[id] = -1.0e30;
    // Search for last position of max value
#pragma omp for
    for (int i = n - 1; i > -1; i--)
    {
      if ((in[i] - maxval[id]) > 1.e-14)
      {
        maxval[id] = in[i];
        maxloc[id] = i;
      }
    }
#pragma omp flush(maxloc, maxval)
#pragma omp single
    {
      loc_m = maxloc[nt - 1];
      for (int i = nt - 2; i > 0; i--)
      {
        if (maxval[i] > loc_m)
        {
          loc_m = maxloc[i];
        }
      }
    }
  }
}

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

void prefix_max(VectorXd &in, int nt, VectorXd &buffd)
{
  {
    // Get current thread id
    int tid = omp_get_thread_num();

    // Compute workload for each thread
    int work = (in.size() + nt - 1) / nt;

    // Set upper and lower bounds for each threads
    int lo = work * tid;
    int hi = std::min(lo + work, int(in.size()));

    // Compute local prefix sum for each thread
    for (int i = lo + 1; i < hi; i++)
      in[i] = std::max(in[i], in[i - 1]);

    // Compute delta for each thread
    buffd[tid] = in[std::max(lo - 1, 0)];

    // Compute global prefix sum
#pragma omp barrier
    for (int j = 1; j < nt; j <<= 1)
    {
      if (tid >= j)
        buffd[tid] = std::max(buffd[tid], buffd[tid - j]);
#pragma omp barrier
    }
    // Update local prefix sum
    for (int i = lo; i < hi; i++)
      in[i] = std::max(in[i], buffd[tid]);
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
  double sum = 0;

#pragma omp parallel
  {
    double l_maxSum = -std::numeric_limits<double>::infinity();
    int l_left = -1;
    int l_right = -1;
    int l_top = -1;
    int l_bottom = -1;
    double l_sum = 0;
    int l_start, l_finish;
    // Loop over rows of 2D matrix
#pragma omp for
    for (int i = 0; i < array.rows(); ++i)
    {
      // Loop over column of 2D matrix
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
  int n = 800;

  /// nxn Matrix filled with random numbers between (-1,1)
  MatrixXd m = MatrixXd::Random(n, n);
  for (int i = 0; i < m.rows(); i++)
  {
    for (int j = 0; j < m.cols(); j++)
    {
      m(i, j) = static_cast<int>(10.0 * m(i, j));
    }
  }
  double maxSum;
  int left, right, top, bottom;

  auto start = std::chrono::high_resolution_clock::now();
  maxSubarray2D(m, maxSum, left, right, top, bottom);
  auto stop = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration<double>(stop - start).count();
  std::cout << "Maxsum: " << maxSum << std::endl;
  std::cout << "Bounds: " << std::endl;
  std::cout << "Left: " << left << " Right: " << right << " Top: " << top << " Bottom: " << bottom << std::endl;
  std::cout << "--------------------------------" << std::endl;
  //std::cout << " [" << M << " ]" << std::endl;
  std::cout << elapsed << " seconds." << std::endl;
  VectorXd buffd;
  VectorXi buffi;

  VectorXd max = m.col(1);
  VectorXd temp = m.col(1);

  int nt;
#pragma omp parallel
  {
#pragma omp single
    {
      nt = omp_get_num_threads();
      buffd.resize(nt);
      buffi.resize(nt);
    }
#pragma omp barrier
    prefix_sum(temp, nt, buffd);
#pragma omp single
    {
      max -= temp;
      temp.reverseInPlace();
    }
#pragma omp barrier
    prefix_max(temp, nt, buffd);
#pragma omp single
    max += temp.reverse();
#pragma omp barrier

// Reverse prefix max sum
#pragma omp single
    temp = m.col(1).reverse();
#pragma omp barrier
    prefix_sum(temp, nt, buffd);
#pragma omp single
    {
      max -= temp.reverse();
      temp.reverseInPlace();
    }
#pragma omp barrier
    prefix_max(temp, nt, buffd);
#pragma omp single
    max += temp;
#pragma omp barrier

    // Search for max
    maxloc_omp(max, nt, buffd, buffi);
  }
}
