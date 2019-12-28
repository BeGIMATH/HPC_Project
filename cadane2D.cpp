#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>

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
    int start, finish;
    for (int i = 0; i < array.rows(); ++i)
    {
        VectorXd temp = VectorXd::Zero(array.cols());
        for (int j = i; j < array.rows(); ++j)
        {
            for (int k = 0; k < array.cols(); ++k)
            {
                temp(k) += array(j, k);
            }
            kadane(temp, sum, start, finish);
            if (sum > maxSum)
            {
                maxSum = sum;
                left = i;
                right = j;
                top = start;
                bottom = finish;
            }
        }
    }
}

int main()
{
    /// Size of the matrix
    int n = 10;

    /// nxn Matrix filled with random numbers between (-1,1)
    MatrixXd m = MatrixXd::Random(n, n);
    double maxSum;
    int left, right, top, bottom;
    maxSubarray2D(m, maxSum, left, right, top, bottom);
    std::cout << m << std::endl;
    std::cout << (right - left) << (bottom - top) << std::endl;
    MatrixXd M = MatrixXd(bottom - top, right - left);
    for (int i = 0; i < M.rows(); i++)
    {
        for (int j = 0; j < M.cols(); j++)
        {
            M(i, j) = m(i + top, j + left);
        }
    }

    std::cout << maxSum << " "
              << " " << left << " " << right << "  " << top << " " << bottom << std::endl;
    std::cout << " [" << M << " ]" << std::endl;
}