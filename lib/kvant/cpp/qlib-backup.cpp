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
class converter_qubit
{
	std::vector<Matrix2c> G;
public:
	converter_qubit(): G (4)
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
	Matrix2c const to_operator(Vector4d const & x) const
	{
		
		Matrix2c X; X.setZero();
		for (int k=0; k<4; k++)
			X+=x(k)*G[k];
		X*=0.5;
		return (X);
	}
	
	//Convert from matrix to vector
	Vector4d const to_vector(Matrix2c const & X) const
	{
		Vector4d x; 
		for (int k=0; k<4; k++)
			x(k)= (X*G[k]).trace().real();
		return (x);
	}
	
	//Convert from bloch tensor to joint operator:
	Matrix4c const to_joint_operator(Matrix4d const & Theta) const
	{
		Matrix4c rho; rho.setZero();
		for (int k1=0; k1<4; k1++)
			for (int k2=0; k2<4; k2++)
				rho+=Theta(k1,k2)*Eigen::kroneckerProduct(G[k1],G[k2]);
		rho*=0.25;
		return rho;
	}
	
	//Convert from joint operator to bloch tensor:
	Matrix4d const to_bloch_tensor(Matrix4c const & rho) const
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
class converter_qutrit
{
	std::vector<Matrix3c> G;
public:
	converter_qutrit(): G (9)
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
	Matrix3c const to_operator(Vector9d const & x) const
	{
		Matrix3c X; X.setZero();
		for (int k=0; k<9; k++)
			X+=x(k)*G[k];
		X*=0.5;
		return (X);
	}
	
	//Convert from matrix to vector
	Vector9d const to_vector(Matrix3c const & X) const
	{
		Vector9d x; 
		for (int k=0; k<9; k++)
			x(k)= (X*G[k]).trace().real();
		return (x);
	}
	
	//Convert from bloch tensor to joint operator:
	Matrix9c const to_joint_operator(Matrix9d const & Theta) const
	{
		Matrix9c rho; rho.setZero();
		for (int k1=0; k1<9; k1++)
			for (int k2=0; k2<9; k2++)
				rho+=Theta(k1,k2)*Eigen::kroneckerProduct(G[k1],G[k2]);
		rho*=0.25;
		return rho;
	}
	
	//Convert from joint operator to bloch tensor:
	Matrix9d const to_bloch_tensor(Matrix9c const & rho) const
	{
		Matrix9d Theta;
		for (int k1=0; k1<9; k1++)
			for (int k2=0; k2<9; k2++)
				Theta(k1,k2)=std::real((rho*Eigen::kroneckerProduct(G[k1],G[k2])).trace());
		return Theta;
	}
} GellMann;
//=============================================================================================
class operator_basis
{
	int D; //Dimension
	std::vector<MatrixXc> G;
public:	
	//-----------------------------------------------------------------------------------------
	operator_basis(int const & D_): D(D_), G(D_*D_)
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
class state: public mtrandom
{
protected:
	int d;
	MatrixXc rho; //The density matrix
public:
	//----------------------------------------------------------------------------------------
	//Constructor:
	state(){};
	state(int const & d_): d(d_), rho(d_,d_){};
	state(MatrixXc const & rho_): rho(rho_) {d=rho.rows();};
	//------------------------------------------------------------------------------------------
	//Get a matrix element:
	std::complex<double> &  operator () (int const & k1, int const & k2) {return rho(k1,k2);};
	//-----------------------------------------------------------------------------------------
	MatrixXc & matrix(){return rho;};
	MatrixXc const & matrix() const {return rho;};
	//-----------------------------------------------------------------------------------------
	//Get density operator:
	MatrixXc const & get_density_operator() const {return rho;};
	state & set_density_operator(MatrixXc const & rho_) {rho=rho_; d=rho.rows(); return (*this);};
	//Get bloch vector:
	VectorXd get_bloch_vector() const 
	{
		int D (d*d);
		VectorXd X(D); X.setZero();
		operator_basis Basis (d);
		for (int k=0; k<D; k++) X(k)=(Basis(k)*rho).trace().real();
		return X;
	};
	state & set_bloch_vector(VectorXd const & X) 
	{
	
		int D (X.size());
		if (D==4) {d=2;} else {d=3;};
		operator_basis Basis (d); 
		rho.resize(d,d);
		for (int k=0; k<D; k++) rho+= X(k)*Basis(k);
		rho*=0.5;
		return (*this);
	}; 
	//-----------------------------------------------------------------------------------------
	//Get dim:
	int const & get_dim() const {return d;};
	//-----------------------------------------------------------------------------------------
	//Get spectral:
	//VectorXd const get_spectral() const {return rho.eigenvalues();};
	//-----------------------------------------------------------------------------------------
	//Set maximally mixed:
	state & set_maximally_mixed() {rho.setIdentity(); rho/= std::complex<double>(d,0.); return (*this);};
	//-----------------------------------------------------------------------------------------
	//Set random:
	//It would be nice to be able to call the partial trace method below, but I do not know how at the moment
	//Technique of sampling is due to Kimo Luoma: double the system, sample a pure state, trace out the ancillar half
	state & set_random_hilbert_schmidt()
	{
	    VectorXc V(d*d); //pure state in dxd system
	    //haar distributed pure states:
	    for (int k=0; k<d*d; k++) V(k)= std::complex<double>(rnormd(mt),rnormd(mt)); //gaussian distributed
	    V/=std::sqrt(V.array().abs2().sum()); //normalise
	    MatrixXc rhoAB (V*V.adjoint());
	    //trace out ancillar:
	    for (int k1=0; k1<d; k1++)
	    {
			for (int k2=0; k2<d; k2++)
			{
				rho(k1,k2)=0.;
				for (int kc=0; kc<d; kc++) rho(k1,k2)+= rhoAB(kc+k1*d,kc+k2*d);
			}
	    }
	    return (*this);
	};
	//----------------------------------------------------------------------------------------
	//Embedd the system into higher dim:
	state & embed(int const & d_)
	{
		MatrixXc rho_old ( rho);
		rho.resize(d_,d_); rho.setZero();
		for (int k1=0; k1<d; k1++)
			for (int k2=0; k2<d; k2++)
				rho(k1,k2)= rho_old(k1,k2);
		d=d_;
		return (*this);
	};

	state & reset(int const & d_)
	{
		d=d_; rho.resize(d,d); rho.setZero();
		return (*this);
	}
};
//---------------------------------------------------------------------------------------------
class pure_state: public state
{
protected:
	VectorXc V;
public:
	//-------------------------------------------------------------------------------------
	//Constructor:
	pure_state(){};
	pure_state(int const & d_): state(d_), V(d_) {};
	pure_state(VectorXc const & V_): V(V_) {d=V.size();rho.resize(d,d);rho=V*V.adjoint();};
	//--------------------------------------------------------------------------------------
	//Set random:
	pure_state set_random_haar()
	{
		for (int k=0; k<d; k++) V(k)= std::complex<double>(rnormd(mt),rnormd(mt));
		V/=std::sqrt(V.array().abs2().sum());
		rho= V*V.adjoint();
		return (*this);
	};
	//-------------------------------------------------------------------------------------
	//Set random qubit:
	pure_state set_random_haar_qubit()
	{
		V.setRandom();
                for (int k=0; k<2; k++) V(k)= std::complex<double>(rnormd(mt),rnormd(mt));
                V/=std::sqrt(V.array().abs2().sum());
                rho= V*V.adjoint();
                return (*this);
	};
	//--------------------------------------------------------------------------------------
	//Embedd the system into higher dim:
	state & embed(int const & d_)
	{
		VectorXc V_old(V);
		V.resize(d_); V.setZero();
		for (int k=0; k<d; k++) V(k)=V_old(k);
		rho.resize(d_,d_); rho= V*V.adjoint();
		d=d_;
		return (*this);
	};
	//Reset:
	state & reset(int const & d_)
	{
		d=d_; V.resize(d); V.setZero(); rho.resize(d,d);
		return (*this);
	}
};
//---------------------------------------------------------------------------------------------
//Class local hidden state ensemble:
//Method of sampling is currently not very good, but should work somehow!
class ensemble: public state
{
protected:
	//int d; //Dimensions
	//MatrixXc rho; //State (center)
	int M; //Number of states
	std::vector<MatrixXc> mu;
public:	
	//-------------------------------------------------------------------------------------
	//Constructor:
	ensemble(int const & d_): state(d_) {};
	ensemble(MatrixXc const & rho_): state(rho_) {};
	ensemble(state const & state_): state(state_) {};
	//-------------------------------------------------------------------------------------
	MatrixXc const & get_weighted_states(int const & k) const {return mu[k];};
	MatrixXc const & operator() (int const & k) const {return mu[k];};
	//-------------------------------------------------------------------------------------
	//Sample method:
	ensemble & set_random(int const & M_)
	{
		M=M_;
		//create the weights:
		VectorXd pp (M);
		for (int k=0; k<M; k++) {pp[k]=(runifd(mt)+1.)/2.;};
		pp/=pp.array().sum();
		//create the states, weighted at the same time
		mu.resize(M);
		for (int k=0; k<M; k++) 
		{ 
			mu[k].resize(d,d);
			pure_state psi (d);
			mu[k]= pp[k]*psi.set_random_haar().get_density_operator();
			
		};
		//renormalise the center by local filtering:
		MatrixXc Center(d,d), RFilter(d,d), LFilter(d,d);
		Center.setZero(); for (int k=0; k<M; k++) {Center+= mu[k];}; //compute the center
		LFilter= Center.sqrt().inverse()*rho.sqrt(); //compute the filter
		RFilter= LFilter.adjoint();
		for (int k=0; k<M; k++) {mu[k]=RFilter*mu[k]*LFilter;};
		return (*this);
	};
	//--------------------------------------------------------------------------------------
	//Embed the whole structure into higher dim:
	ensemble & embed(int const & d_)
	{
		state::embed(d_);
		for (int k=0; k<M; k++)
		{
			state sigma (mu[k]);
			mu[k]= sigma.embed(d_).get_density_operator();
		};
		return (*this);
	};

	ensemble & reset(int const & d_, int const & M_)
	{
		d=d_;M=M_;
		rho.resize(d,d);
		mu.resize(M);
		for (int k=0; k<M; k++) mu[k].resize(d,d);
		return (*this);
	};
};
//---------------------------------------------------------------------------------------------
//Class of bipartite states: the basis is always canonical!
class bipartite_state: public state
{
protected:
	int Da, Db; //Dimension of parties A and B
	//------------------------------------------------------------------------------------
public:
	//Constructors:
	bipartite_state (int const & Da_, int const & Db_): Da(Da_), Db(Db_), state(Da_*Db_) {};
	bipartite_state (MatrixXc const & rho_, int const & Da_, int const & Db_): Da(Da_), Db(Db_), state(rho_) 
	{	
		if (d != Da*Db) {std::cout<<"***Error::bipartite_state(...): dimension mismatched!!! Exited!!! \n"; exit(0);};	
	};
	//----------------------------------------------------------------------------------
	MatrixXd get_bloch_tensor() const 
	{
		int La (Da*Da), Lb(Db*Db);
		MatrixXd Theta(La,Lb); Theta.setZero();
		operator_basis BasisA (Da); 
		operator_basis BasisB (Db); 
		for (int ka=0; ka<La; ka++) 
			for (int kb=0; kb<Lb; kb++)
				Theta(ka,kb)=(Eigen::kroneckerProduct(BasisA(ka),BasisB(kb))*rho).trace().real();
		return Theta;
	};
	bipartite_state & set_bloch_tensor(MatrixXd & Theta) //Unclear why constant is not possible for Theta??? 
	{	
		int La (Theta.rows()), Lb(Theta.cols());
		if (La==4) {Da=2;} else {Da=3;};
		if (Lb==4) {Db=2;} else {Db=3;};
		d= Da*Db;
		rho.resize(d,d);
		operator_basis BasisA(Da), BasisB(Db);
		for (int ka=0; ka<La; ka++)
			for (int kb=0; kb<Lb; kb++)
				rho+= (Theta(ka,kb)*Eigen::kroneckerProduct(BasisA(ka),BasisB(kb)));
		rho*=0.25;
		return (*this);
	}
	//------------------------------------------------------------------------------------
	//To write: set_werner() and set_isotropic()
	bipartite_state & set_isotropic()
	{
		VectorXc V(d); V.setZero();
		int Dmin (std::min(Da,Db));
		double s (1./std::sqrt(double(Dmin)));
		for (int k=0; k<Dmin; k++) V(tensor_id(k,k))=std::complex<double>(0.,s);
		rho= V*V.adjoint();
		return (*this);
	}
	//------------------------------------------------------------------------------------
	bipartite_state & set_two_qubit_isotropic()
	{
		VectorXc V(d); V.setZero();
		int Dmin (2);
		double s (1./std::sqrt(double(Dmin)));
		for (int k=0; k<Dmin; k++) V(tensor_id(k,k))=std::complex<double>(0.,s);
        rho= V*V.adjoint();
        return (*this);
	}
	//--------------------------------------------------------------------------------------
	bipartite_state & set_two_qubit_random()
    {
        //random 2x2 state:
        state rhoAB (4); rhoAB.set_random_hilbert_schmidt();
	    //std::cout<< rhoAB.get_density_operator()<< "\n\n\n"; //Looks good! 
	    //embed into dxd system:
	    rho.setZero();
	    for (int ka1=0; ka1<2; ka1++)
			for (int kb1=0; kb1<2; kb1++)
		    	for (int ka2=0; ka2<2; ka2++)
					for (int kb2=0;kb2<2;kb2++)
						rho(tensor_id(ka1,kb1),tensor_id(ka2,kb2))=rhoAB(ka1*2+kb1,ka2*2+kb2);
	    //std::cout<< rho<<"\n\n\n"; //also looks good!
		return (*this);
     }
	//------------------------------------------------------------------------------------
	//Partial traces:
	//TrA:
	state TrA() 
	{
		state rhoB (Db);
		for (int kb1=0; kb1<Db; kb1++)
		{
			for (int kb2=0; kb2<Db; kb2++)
			{
				rhoB(kb1,kb2)=0.;
				for (int ka=0; ka<Da; ka++)
					rhoB(kb1,kb2)+= rho(tensor_id(ka,kb1),tensor_id(ka,kb2));
			}
				
		}
		return rhoB;
	};
	//TrB:
	state TrB() 
	{
		state rhoA (Da);
                for (int ka1=0; ka1<Da; ka1++)
                {
                        for (int ka2=0; ka2<Da; ka2++)
                        {
                                rhoA(ka1,ka2)=0.;
                                for (int kb=0; kb<Db; kb++)
                                        rhoA(ka1,ka2)+= rho(tensor_id(ka1,kb),tensor_id(ka2,kb));
                        }

                }
                return rhoA;
	};
	//-------------------------------------------------------------------------------------
	/*bipartite_state filter_bob() //to be tested!!!!!!
	{
		MatrixXd BFilter ((TrA()).matrix().sqrt().inverse());
		MatrixXd BFilter (Eigen::kroneckerProduct(MatrixXc::Identity(Da,Da),(this*).TrA().matrix().sqrt().inverse()));
		rho=BFilter*rho*BFilter;
		rho/=rho.trace();
		return (*this);
	}*/
	//--------------------------------------------------------------------------------------
	//emdding into higher dim:
	bipartite_state & embed(int const & Da_, int const & Db_)
    {
		//copy the current state to a new one:
		MatrixXc rho_ (rho);  
	    //embed into dxd system:
	    rho.resize(Da_*Db_,Da_*Db_);
		rho.setZero();
	    for (int ka1=0; ka1<Da; ka1++)
			for (int kb1=0; kb1<Db; kb1++)
		    	for (int ka2=0; ka2<Da; ka2++)
					for (int kb2=0; kb2<Db;kb2++)
						rho(ka1*Db_+kb1,ka2*Db_+kb2)=rho_(ka1*Db+kb1,ka2*Db+kb2);
		Da=Da_; Db=Db_; d=Da*Db;
		return (*this);
     };

	bipartite_state & reset(int const & Da_, int const & Db_)
	{
		Da=Da_; Db=Db_; d= Da*Db;
		rho.resize(d,d);
		return (*this);
	};
private:
	inline int const tensor_id(int const & ka, int const & kb) {return (ka*Db+kb);};
};
//---------------------------------------------------------------------------------------------
//class bipartite_pure_state: public bipartite_state
//{
//	int Da, Db;
//	VectorXc Psi;
	
	
	//TrA
//	state TrA() {};
//	state TrB() {};
//	
//};
//=============================================================================================
#endif
