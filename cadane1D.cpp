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

int main()
{
    //Size of the array
    int n = 8;
    //An array  o size n filled with random numbers between (-1,1)
    VectorXd v = VectorXd::Random(n);

    for (int i = 0; i < v.size(); i++)
    {
        v(i) = static_cast<int>(10.0 * v(i));
    }
    //std::cout << " [ " << v.transpose() << " ] " << std::endl;
    double maxSum;
    int left, right;
    kadane(v, maxSum, left, right);
    std::cout << " [ " << v.transpose() << " ] " << std::endl;
    //std::cout << left << " " << right << " " << std::endl;
    VectorXd V = VectorXd(right - left + 1);
    for (int i = 0; i < V.size(); i++)
    {
        V(i) = v(i + left);
    }
    std::cout << maxSum << "  " << std::endl;
    //std::cout << "  " << left << "  " << right << "  " << std::endl;
    std::cout << " [ " << V.transpose() << " ] " << std::endl;
}