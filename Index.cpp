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


// Implementations for Banded Matrices //
/**
 * @brief Populate Banded Matrix
 * @param   ny     	Number of rows already known to be 2 as this is a banded matrix
 * @param   Mat       Matrix
 */
void BandedMatrix1(int ny, double* Mat,const double val1, const double val2) {
    const int ldh = 2;      // Diagonal and upper diagonal, No. of stored diagonals

    Mat[1] = val1;
    for (int i = 1; i < ny; ++i) {
        Mat[i*ldh    ] = val1;
        Mat[i*ldh + 1] = val2;
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
	const int ldh=2;

  double* Mat = new double[nx*ny]; // General matrix storage
	double*	BandedMat = new double[ldh*ny];

	BandedMatrix1(ny,BandedMat,2,3);	
	PrintMatrix(ldh,ny,BandedMat);

	return 0;
}
