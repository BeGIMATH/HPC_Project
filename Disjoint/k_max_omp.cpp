
#include <iostream>
#include <limits>
#include "eigen/Eigen/Dense"
#include <algorithm>
#include <chrono>
#include <omp.h>

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
    double sum = 0;

#pragma omp parallel
    {
        double l_M = -std::numeric_limits<double>::infinity();
        int l_c_1 = -1;
        int l_c_2 = -1;
        int l_r_1 = -1;
        int l_r_2 = -1;
        double l_t = 0;
        int l_start, l_finish;
        // Loop over rows of 2D matrix

        int i, id, nthrds, nthreads;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if (id == 0)
            nthreads = nthrds;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();

        for (i = id; i < array.rows(); i = i + nthrds)
        {

            VectorXd temp = VectorXd::Zero(array.cols());
            for (int j = i; j < array.rows(); ++j)
            {

                for (int k = 0; k < array.cols(); ++k)
                {
                    temp(k) += array(j, k);
                }
                Kadane(temp, l_t, l_start, l_finish);
                if (l_t > l_M)
                {
                    l_M = l_t;
                    l_c_1 = l_start;
                    l_c_2 = l_finish;
                    l_r_1 = i;
                    l_r_2 = j;
                }
            }
        }
#pragma omp critical
        {
            if (l_M > M)
            {
                M = l_M;
                c_1 = l_c_1;
                c_2 = l_c_2;
                r_1 = l_r_1;
                r_2 = l_r_2;
            }
        }
    }
}

void update(MatrixXd &array, int &left, int &right, int &top, int &bottom)
{
    int i, j;
#pragma omp parallel for private(j)

    for (i = top; i <= bottom; i++)
    {
        for (j = left; j <= right; j++)
        {
            array(i, j) = -std::numeric_limits<double>::infinity();
        }
    }
}
int main()
{

    /*Size of the matrix*/
    int N[5] = {100, 200, 300, 400, 500};

    for (int l = 0; l < 5; l++)
    {
        MatrixXd m = MatrixXd::Random(N[l], N[l]);
        for (int i = 0; i < m.rows(); i++)
        {
            for (int j = 0; j < m.cols(); j++)
            {
                m(i, j) = static_cast<int>(10.0 * m(i, j));
            }
        }
        double maxSum;
        int left, right, top, bottom;
        int k = 5;
        //int k = 10;
        //int k = 20;

        double start_time = omp_get_wtime();
        for (int i = 0; i < k; i++)
        {
            MSP2D(m, maxSum, left, right, top, bottom);
            /*		
		    std::cout << "Maxsum: " << maxSum << std::endl;
    		std::cout << "Bounds: " << std::endl;
    		std::cout << "Left: " << left << " Right: " << right << " Top: " << top << " 			Bottom: " << 	  bottom << std::endl;
		*/
            update(m, left, right, top, bottom);
        }
        double stop = omp_get_wtime() - start_time;

        std::cout << "--------------------------------" << std::endl;
        std::cout << "Time for size" << N[l] << " " << stop << " seconds." << std::endl;
    }
}
