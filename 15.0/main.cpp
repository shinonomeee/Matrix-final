/*
 * @Description: main.cpp
 * @version: 10.0
 * @Author: shinonome
 * @Date: 2021-03-08 21:46:40
 * @LastEditor: shinonome
 * @LastEditTime: 2021-03-09 00:42:34
 */

#include "matrix.h"
using namespace Linear_algebra;

int main(int argc, const char *argv[]) {
	// matrix<double, 2, 3> Mymatrix;
	// cin >> Mymatrix;
	// cout << Mymatrix.I_Tran();
	// cout << Mymatrix << endl;
	matrix<double, 3, 3> Matrix(5);
	cin >> Matrix;
	cout << Matrix.Gaussian_elimination();
	cout << Matrix.det() << endl;

	return 0;
}
