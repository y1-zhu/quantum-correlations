#ifndef CEIGEN
#define CEIGEN
#include <complex>
//==================================================================================================
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>
using namespace Eigen;
typedef Matrix<std::complex<double>,2,2> Matrix2c;
typedef Matrix<std::complex<double>,Dynamic,1> VectorXc;
typedef Matrix<std::complex<double>,Dynamic,Dynamic> MatrixXc;
typedef Matrix<int,Dynamic,1> VectorXi;
typedef Matrix<int,Dynamic,Dynamic> MatrixXi;
typedef Matrix<bool,Dynamic,1> VectorXb;
typedef Matrix<bool,Dynamic,Dynamic> MatrixXb;
typedef Array<bool, Dynamic, 1> ArrayXb;
typedef Array<int, Dynamic, 1> ArrayXi;
typedef Array<double, Dynamic, 1> ArrayXd;
/*
void sapply(VectorXd & X, double (*f) (double const &))
{
	int N= X.size();
	for (int n= 0; n<N; n++) {X(n)= f(X(n));};
};
//
void sapply(VectorXd & Y, VectorXd const & X, double (*f) (double const &))
{
	int N= X.size(); Y.resize(X);
	for (int n= 0; n<N; n++) {Y(n)= f(X(n));};
};
*/
//==================================================================================================

#endif