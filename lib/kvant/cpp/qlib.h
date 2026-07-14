#ifndef QLIB
#define QLIB
//=============================================================================================
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "ccomplex.h"
#include "crandomx.h"
#include "ceigen.h"
//inline double const sqr(double const & x){return (x*x);};
//=============================================================================================
//Type of complex matrices and real matrix
//=============================================================================================
//Some generic class of quantum states (to be included in qlib.cpp):
//Type of complex matrices and real matrix
typedef Eigen::Matrix<std::complex<double>,2,1> Vector2c;
typedef Eigen::Matrix<std::complex<double>,2,2> Matrix2c;
typedef Eigen::Matrix<std::complex<double>,4,4> Matrix4c;

typedef Eigen::Matrix<double,4,1> Vector4d;
typedef Eigen::Matrix<double,4,4> Matrix4d;

typedef Eigen::Matrix<std::complex<double>,3,1> Vector3c;
typedef Eigen::Matrix<std::complex<double>,3,3> Matrix3c;
typedef Eigen::Matrix<std::complex<double>,9,9> Matrix9c;

typedef Eigen::Matrix<double,9,1> Vector9d;
typedef Eigen::Matrix<double,9,9> Matrix9d;
//=============================================================================================
//Qubit basis:
class ConverterQubit
{
	std::vector<Matrix2c> G;
public:
	ConverterQubit(): G (4)
	{
		//Set up the Gellmann matrices:
		G[0]<<	1.,0.,
				0.,1.;
		G[1]<<  0.,1.,
				1.,0.;
		G[2]<<  0.,-i,
				+i,0.;
		G[3]<<  1., 0.,
				0.,-1.;
	
		//for (int k=0; k<9; k++)
		//	std::cout<<(G[k]*G[k]).trace().real()<< std::endl;
	}
	//Return the gellmann matrix:
	Matrix2c const & operator () (int const & k) const {return G[k];}	 
	Matrix2c const & get(int const & k) const {return G[k];}	
	//Convert from vector to matrix:
	Matrix2c const toOperator(Vector4d const & x) const
	{
		
		Matrix2c X; X.setZero();
		for (int k=0; k<4; k++)
			X+=x(k)*G[k];
		X*=0.5;
		return (X);
	}
	
	
	//Convert from bloch tensor to joint operator:
	Matrix4c const toOperator(Matrix4d const & Theta) const
	{
		Matrix4c rho; rho.setZero();
		for (int k1=0; k1<4; k1++)
			for (int k2=0; k2<4; k2++)
				rho+=Theta(k1,k2)*Eigen::kroneckerProduct(G[k1],G[k2]);
		rho*=0.25;
		return rho;
	}

	
	//Convert from matrix to vector
	Vector4d const toBlochVector(Matrix2c const & X) const
	{
		Vector4d x; 
		for (int k=0; k<4; k++)
			x(k)= (X*G[k]).trace().real();
		return (x);
	}
	
	//Convert from joint operator to bloch tensor:
	Matrix4d const toBlochTensor(Matrix4c const & rho) const
	{
		Matrix4d Theta;
		for (int k1=0; k1<4; k1++)
			for (int k2=0; k2<4; k2++)
				Theta(k1,k2)=std::real((rho*Eigen::kroneckerProduct(G[k1],G[k2])).trace());
		return Theta;
	}
} Pauli;

//==============================================================================================
//Class to convert for operator (3x3) complex to real vector 9D
class ConverterQutrit
{
	std::vector<Matrix3c> G;
public:
	ConverterQutrit(): G (9)
	{
		//Set up the Gellmann matrices:
		G[0]<<	1.,0.,0.,
			0.,1.,0.,
			0.,0.,1.;
		G[0]*=std::sqrt(2./3.);
		G[1]<<  0.,1.,0.,
			1.,0.,0.,
			0.,0.,0.;
		G[2]<<  0.,-i,0.,
			+i,0.,0.,
			0.,0.,0.;
		G[3]<<  1., 0.,0.,
			0.,-1.,0.,
			0., 0.,0.;
		G[4]<<  0.,0.,1.,
			0.,0.,0.,
			1.,0.,0.;
		G[5]<<  0.,0.,-i,
			0.,0.,0.,
			+i,0.,0.;
		G[6]<<  0.,0.,0.,
			0.,0.,1.,
			0.,1.,0.;
		G[7]<<  0.,0.,0.,
			0.,0.,-i,
			0.,+i,0.;
		G[8]<<  1.,0., 0.,
			0.,1., 0.,
			0.,0.,-2.;
		G[8]*=std::sqrt(1./3.);
		//for (int k=0; k<9; k++)
		//	std::cout<<(G[k]*G[k]).trace().real()<< std::endl;
	}
	//Return the gellmann matrix:
	Matrix3c const & get(int const & k) const {return G[k];}	
	//Convert from vector to matrix:
	Matrix3c const toOperator(Vector9d const & x) const
	{
		Matrix3c X; X.setZero();
		for (int k=0; k<9; k++)
			X+=x(k)*G[k];
		X*=0.5;
		return (X);
	}
	
	//Convert from bloch tensor to joint operator:
	Matrix9c const toOperator(Matrix9d const & Theta) const
	{
		Matrix9c rho; rho.setZero();
		for (int k1=0; k1<9; k1++)
			for (int k2=0; k2<9; k2++)
				rho+=Theta(k1,k2)*Eigen::kroneckerProduct(G[k1],G[k2]);
		rho*=0.25;
		return rho;
	}

	//Convert from matrix to vector
	Vector9d const toBlochVector(Matrix3c const & X) const
	{
		Vector9d x; 
		for (int k=0; k<9; k++)
			x(k)= (X*G[k]).trace().real();
		return (x);
	}
	
	
	//Convert from joint operator to bloch tensor:
	Matrix9d const toBlochTensor(Matrix9c const & rho) const
	{
		Matrix9d Theta;
		for (int k1=0; k1<9; k1++)
			for (int k2=0; k2<9; k2++)
				Theta(k1,k2)=std::real((rho*Eigen::kroneckerProduct(G[k1],G[k2])).trace());
		return Theta;
	}
} GellMann;
//=============================================================================================
class OperatorBasis
{
	int D; //Dimension
	std::vector<MatrixXc> G;
public:	
	//-----------------------------------------------------------------------------------------
	OperatorBasis(int const & D_): D(D_), G(D_*D_)
	{
		if (D==2)
		{
			G.resize(4); for (int k=0; k<4; k++) G[k].resize(2,2);
			G[0]<<	1.,0.,
					0.,1.;
			G[1]<<  0.,1.,
					1.,0.;
			G[2]<<  0.,-i,
					+i,0.;
			G[3]<<  1., 0.,
					0.,-1.;
		};
		if (D==3)
		{
			G.resize(9); for (int k=0; k<9; k++) G[k].resize(3,3);
			
			G[0]<<	1.,0.,0.,
					0.,1.,0.,
					0.,0.,1.;
			G[0]*=std::sqrt(2./3.);
			G[1]<<  0.,1.,0.,
					1.,0.,0.,
					0.,0.,0.;
			G[2]<<  0.,-i,0.,
					+i,0.,0.,
					0.,0.,0.;
			G[3]<<  1., 0.,0.,
					0.,-1.,0.,
					0., 0.,0.;
			G[4]<<  0.,0.,1.,
					0.,0.,0.,
					1.,0.,0.;
			G[5]<<  0.,0.,-i,
					0.,0.,0.,
					+i,0.,0.;
			G[6]<<  0.,0.,0.,
					0.,0.,1.,
					0.,1.,0.;
			G[7]<<  0.,0.,0.,
					0.,0.,-i,
					0.,+i,0.;
			G[8]<<  1.,0., 0.,
					0.,1., 0.,
					0.,0.,-2.;
			G[8]*=std::sqrt(1./3.);
		};
		if (D>3)
		{
			std::cout<< "***Error: not yet implemented!!!\n";
		}
	};
	//-----------------------------------------------------------------------------------------
	//Get basis kth:
	MatrixXc &  operator () (int const & k) {return G[k];};
};
//=============================================================================================
#endif
