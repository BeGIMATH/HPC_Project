#include <iostream>
#include <limits>
#include <eigen3/Eigen/Dense>
#include <chrono>
using namespace Eigen;
/*A function to to compute the maximum subarray in 1D*/

void Kadane(const VectorXd &array, double &M, int &x_1, int &x_2)
{
    /*Initialize the max sum to - infinity*/
    M = -std::numeric_limits<double>::infinity();
    x_1 = 0;
    x_2 = 0;
    double t = 0;
    int c_S_Ind = 0;
    for (int i = 0; i < array.size(); ++i)
    {
        t += array(i);
        if (t > M)
        {
            M = t;
            x_1 = c_S_Ind;
            x_2 = i;
        }
        if (t < 0)
        {
            t = 0;
            c_S_Ind = i + 1;
        }
    }
}
void MSP2D(const MatrixXd &array,
           double &M, int &c_1, int &c_2, int &r_1, int &r_2)
{

    M = -std::numeric_limits<double>::infinity();
    c_1 = -1;
    c_2 = -1;
    r_1 = -1;
    r_2 = -1;
    double t = 0;
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
            Kadane(temp, t, start, finish);
            if (t > M)
            {
                M = t;
                c_1 = start;
                c_2 = finish;
                r_1 = i;
                r_2 = j;
            }
        }
    }
}

int main()
{
    /// Size of the matrix
    int n = 500;

    /// nxn Matrix filled with random integers between (-10,10)
    MatrixXd m = MatrixXd::Random(n, n);
    for (int i = 0; i < m.rows(); i++)
    {
        for (int j = 0; j < m.cols(); j++)
        {
            m(i, j) = static_cast<int>(10.0 * m(i, j));
        }
    }
    double M;
    int left, right, top, bottom;
    auto start = std::chrono::high_resolution_clock::now();
    MSP2D(m, M, left, right, top, bottom);
    auto stop = std::chrono::high_resolution_clock::now();
    //std::cout << m << std::endl;
    //std::cout << "----------------------------" << std::endl;
    MatrixXd Mat = MatrixXd(bottom - top + 1, right - left + 1);
    for (int i = 0; i < Mat.rows(); i++)
    {
        for (int j = 0; j < Mat.cols(); j++)
        {
            Mat(i, j) = m(i + top, j + left);
        }
    }

    auto elapsed = std::chrono::duration<double>(stop - start).count();
    std::cout << "Maxsum: " << M << std::endl;
    std::cout << "Bounds: " << std::endl;
    std::cout << "Left: " << left << " Right: " << right << " Top: " << top << " Bottom: " << bottom << std::endl;
    std::cout << "--------------------------------" << std::endl;
    //std::cout << " [" << Mat << " ]" << std::endl;
    std::cout << elapsed << " seconds." << std::endl;
}