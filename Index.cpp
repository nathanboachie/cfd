#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * @ brief This script serves the purpose of trying to explain how to fill matrices based upon their diagonals
 * Demonstrates implementations for square and technically non-square matrices
 * Also provides codes for banded storage matrices and printing for visualisation sake
 * Creating the correct problem can be just as complex as solving it
**/


/**
 * @brief Populate General Matrix
 *
 * @param   n       Leading dimension of matrix, Assumes square
 * @param   Mat     Pointer to matrix storage of size n*n
 * @param   val1    Value 1
 * @param		val2		Value 2
 * @param   dx      Grid spacing
 */
void Matrix1(const int n,double* Mat,const int val1, const int val2) {
    for (int i = 1; i < n; ++i) {
        Mat[i*n + i - 1] = val1;
        Mat[i*n + i] = val2;
    }
}

/**
 * @brief Populate General Matrix
 *
 * @param   m				Number of rows
 * @param 	n 			Number of colummns
 * @param   Mat     Pointer to matrix storage of size n*n
 * @param   val1    Value 1
 * @param		val2		Value 2
 * @param   dx      Grid spacing
 */

void Matrix2(const int m, const int n, double* Mat, const int val1, const int val2) {
		for(int i=0; i<m; ++i)
		{
			for(int j=0;j<n;++j)
			{
				if(j==i)
				{
					Mat[i*n+j]=val1;
				}
				else if(j==i-1)
				{
					Mat[i*n+j]=val2;
				}
			}
		}
}

void PrintMatrix(const int m, const int n, double* Mat) {
		std::cout.precision(4);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
					std::cout << std::setw(6) << Mat[i*n+j] << " ";
        }
				std::cout << std::endl;
    }
		std::cout << std::endl;
}

int main()
{
  const int nx=20;
	const int ny=20;

  double* Mat = new double[nx*ny];// Helmholtz matrix storagea
	Matrix2(nx,ny,Mat,2,3);
	PrintMatrix(nx,ny,Mat);

	return 0;
}
