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
 * @param   nsv     Leading dimension of matrix
 * @param   H       Pointer to matrix storage of size n*n
 * @param   lam     lambda value
 * @param   dx      Grid spacing
 */
/**
 * @brief Populate Helmholtz symmetric matrix (upper only)
 *
 * @param   nsv     Leading dimension of matrix
 * @param   H       Pointer to matrix storage of size nsv*nsv
 * @param   lam     Lambda coefficient
 * @param   dx      Grid spacing
 */
void FillHelmholtzMatrix1(int nsv, double* H, double lam, double dx) {
    const double oodx2 = 1.0/dx/dx;
    H[0] = -lam - 2.0*oodx2;
    for (int i = 1; i < nsv; ++i) {
        H[i*nsv + i - 1] = oodx2;
        H[i*nsv + i] = -lam - 2.0*oodx2;
    }
}

void FillHelmholtzMatrix2(int nsv, double* H, double lam, double dx) {
    const double oodx2 = 1.0/dx/dx;
		for(int i=0; i<nsv; ++i)
		{
			for(int j=0;j<nsv;++j)
			{
				if(j==i)
				{
					H[i*nsv+j]=-lam-2.0*oodx2;
				}
				else if(j==i-1)
				{
					H[i*nsv+j]=oodx2;
				}
			}
		}
}

void PrintMatrix(int nsv, double* H) {
		std::cout.precision(4);
    for (int i = 0; i < nsv; ++i) {
        for (int j = 0; j < nsv; ++j) {
					std::cout << std::setw(6) << H[i*nsv+j] << " ";
        }
				std::cout << std::endl;
    }
		std::cout << std::endl;
}

int main()
{
	const int    n   = 21;          // Number of grid-points
  const int    nsv = n - 2;       // Number of unknown DOFs
  const double lam = 1.0;         // Value of Lambda
  const double L   = 1.0;         // Length of domain
  const double dx  = L / (n - 1); // Grid-point spacing

  double* H = new double[nsv*nsv];// Helmholtz matrix storagea
	FillHelmholtzMatrix2(nsv,H,lam,dx);
	PrintMatrix(nsv,H);

	return 0;
}
