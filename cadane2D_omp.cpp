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
                    l_left = l_start;
                    l_right = l_finish;
                    l_top = i;
                    l_bottom = j;
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
    int n = 1000;

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
}
