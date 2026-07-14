#ifndef CRANDOMHEADER
#define CRANDOMHEADER
#include <random>
#include <cmath>
//==================================================================================================
//std::random_device rd;
//std::mt19937 rgen(rd());
//std::normal_distribution<double> rnorm (0.,1.);
//std::uniform_real_distribution<double> runif (0.,1.);
//std::exponential_distribution<double> rexp (1.0);
//#define rnorm() rnorm(rgen)
//#define runif() runif(rgen)
//#define rexp() rexp(rgen)
//#define rint() rgen()
//=================================================================================================
std::random_device rdevice;
class mtrandom
{
protected:
	//Random number generators
	std::mt19937 mt;
	std::normal_distribution<double> rnormd;
	std::uniform_real_distribution<double> runifd;

public:
	//Constructor:
	mtrandom():  mt(rdevice()), rnormd (0.,1.), runifd (0.,1.) {};
	//Random generating methods:
	inline unsigned long int rint() {return mt();}
	inline double rnorm() {return rnormd(mt);}	
	inline double runif() {return runifd(mt);}	
} RGen;
//==================================================================================================
#endif
