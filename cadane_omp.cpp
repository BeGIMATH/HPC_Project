#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <omp.h>

using namespace Eigen;

/*The function which performs the 1D Kadane algorithm*/
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
void maxSubarray2D(const MatrixXd &array, double &maxSum, int &left, int &right, int &top, int &bottom)
{
    maxSum = -std::numeric_limits<double>::infinity();
    left = -1;
    right = -1;
    top = -1;
    bottom = -1;
    /*Enter the parallel region*/

#pragma omp parallel
    {

        //Create local variables for kadane 2D
        double l_maxSum = -std::numeric_limits<double>::infinity();
        int l_left = -1;
        int l_right = -1;
        int l_top = -1;
        int l_bottom = -1;
        double l_sum = 0;
        int l_start, l_finish;
        //Why this way of parallelizing is slower ?
        //Then the sequential one
        //Number of elements to compute per processors
        int size = (array.rows() + 1) * array.rows() / (2 * omp_get_num_threads());
        //Define bounds of iteration

        int lo, hi;
        int id = omp_get_thread_num();
        for (int l = 0; l < 1; l++)
        {
            lo = (l == 0) ? id * size : array.rows() - (id + 1) * size;
            hi = (l == 0) ? (id + 1) * size : array.rows() - id * size;
            for (int i = 0; i < hi; ++i)
            {
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
                        l_left = l_start;
                        l_right = l_finish;
                        l_top = i;
                        l_bottom = j;
                    }
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
    //Size of teh matrix

    int n = 800;
    MatrixXd m = MatrixXd::Random(n, n);
    for (int i = 0; i < m.rows(); i++)
    {
        for (int j = 0; j < m.cols(); j++)
        {
            m(i, j) = static_cast<int>(10.0 * m(i, j));
        }
    }
    auto start = std::chrono::high_resolution_clock::now();
    double maxSum;
    int left, right, top, bottom;
    maxSubarray2D(m, maxSum, left, right, top, bottom);
    auto stop = std::chrono::high_resolution_clock::now();
    /*
    MatrixXd M = MatrixXd(bottom - top + 1, right - left + 1);
    
    for (int i = 0; i < M.rows(); i++)
    {
        for (int j = 0; j < M.cols(); j++)
        {
            M(i, j) = m(i + top, j + left);
        }
    }
    */
    auto elapsed = std::chrono::duration<double>(stop - start).count();
    std::cout << "Maxsum: " << maxSum << std::endl;
    std::cout << "Bounds: " << std::endl;
    std::cout << "Left: " << left << " Right: " << right << " Top: " << top << " Bottom: " << bottom << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << elapsed << " seconds." << std::endl;
    //std::cout << " [" << M << " ]" << std::endl;
}
