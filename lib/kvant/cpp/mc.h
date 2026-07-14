#ifndef CANNEALHEADER 
#define CANNEALHEADER

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <fstream>
#include <string>
#include <cmath>
#include "crandomx.h"
#include "ceigen.h"
//==================================================================================================
//Virtual objective function:
class objective_function: public mtrandom
{
protected:	
	int const N; //Size of the problem (number of variables)
	double f; //Value of the objective function
public:
	//----------------------------------------------------------------------------------------------
	//Constructor
	objective_function(int const & N_): N(N_){};
	//----------------------------------------------------------------------------------------------
	inline int const & size() const {return N;};
	inline double const & value() const {return f;};
	//----------------------------------------------------------------------------------------------
	//Sample the range of the function
	double const sample_range(int const & n_sample=10000) 
	{
		double f_min(1.e20), f_max (-1.e20), f_sample (0.);
		for (int k_sample=0; k_sample<n_sample; k_sample++)
			for (int k=0; k<N; k++)
				{randomise(); evaluate(); f_min=std::min(f_min,f); f_max=std::max(f_max,f);};
		return (f_max-f_min);
	};
	//----------------------------------------------------------------------------------------------
	//Randomize the state:
	virtual inline void randomise()= 0;
	//Step the function forward:
	virtual inline double const step(double const & Delta=1.) = 0;
	//Step the function backward:
	virtual inline void reverse() = 0;
	//----------------------------------------------------------------------------------------------
	//Evaluate the function:
	virtual inline double const evaluate()  = 0;
};
//==================================================================================================
//Simulated annealer to minimise the objective function
class mc: public mtrandom
{
protected:
	double f_min;
	double T;
	//double Delta;
public:
	mc() {};
	//----------------------------------------------------------------------------------------------
	//Set and get functions
	double const & setT(double const & T_) {T= T_; f_min=1e99; return T;};
	double const & getT() const {return T;};
	double const & min() const {return f_min;};
	//----------------------------------------------------------------------------------------------
	//Sample the temperature from the object
	double const & sampleT(objective_function & objfunc, int const & n_sample=10000) {T= 10.*objfunc.sample_range(n_sample); return T;};
	//Do M00 monte carlo steps on SYS
	double const & mcrun(objective_function & objfunc, int const & n_step=1)
	{
		//Dry run n_step steps
		for (int k= 0; k<n_step; k++)
		{	
			double df= objfunc.step(sqrt(T));
			if (df>0)
				if (runifd(mt)>exp(-df/T))
					objfunc.reverse();
			f_min= std::min(f_min,objfunc.value());
		};
		return f_min;
	};
	//----------------------------------------------------------------------------------------------
	//Cool the system down until a minimal temperature is reached! This function is not verbose; use manually cooling for verbose access
	double const & cool(objective_function & objfunc, int const & n_step=1000, double const & gamma=0.95, double const & T_min=1e-8)
	{
		//double T_max(T);
		f_min= objfunc.evaluate();
		while (T>T_min)
		{
			for (int k= 0; k<n_step; k++)
			{	
				double df= objfunc.step(sqrt(T));
				if (df>0)
					if (runifd(mt)>exp(-df/T))
						objfunc.reverse();
				f_min= std::min(f_min,objfunc.value());
			}
			T=gamma*T; //exponential cooling
			//T-=(T_max-T_min)/6000.; //linear cooling
		};
		return f_min;
	};
	//-----------------------------------------------------------------------------------------------
	//Cool the system down until a minimal temperature is reached! This function is not verbose; use manually cooling for verbose access
	double const & cool(objective_function & objfunc, std::ostream & flog, int const & n_step=1000, double const & gamma=0.995, double const & T_min=1e-8)
	{
		f_min= objfunc.evaluate();
		flog<< T<< "\t"<< f_min<< "\t"<< f_min<< "\n";
		while (T>T_min)
		{
			for (int k= 0; k<n_step; k++)
			{	
				double df= objfunc.step(sqrt(T));
				//std::cout<< "Here df=" << df<< "\n";
				//exit(0);
				if (df>0)	
					//std::cout<< df<< "\n";
					if (runifd(mt)>exp(-df/T))
						objfunc.reverse();
				f_min= std::min(f_min,objfunc.value());
			}
			flog<< T<< "\t"<< objfunc.evaluate()<< "\t"<< f_min<< "\n";
			//anneal(objfunc,n_step);
			T=gamma*T;
		};
		return f_min;
	};
	//Cool and recool the system:
	double const & repeated_cool(objective_function & objfunc, double const & T_max=100., int const & n_repeat=16, int const & n_step= 1000, double const & gamma=0.95, double const & T_min=1e-8)
	{
		//objfunc.randomise();
		//std::cout<< "Coming in... \n";
		f_min=objfunc.evaluate();
		//std::cout<< "fmin= "<< fmin<< "\n";
		//if (n_step<=0){n_step = 1000*N;};
		for (int k=0; k<n_repeat; k++)
		{
			//std::cout<<"*Repeat= "<< k<<"\n";
			objfunc.randomise();
			T=T_max;
			while (T>T_min)
			{
				for (int k= 0; k<n_step; k++)
				{	
					double df= objfunc.step(sqrt(T));
					if (df>0)
						if (runifd(mt)>exp(-df/T))
							objfunc.reverse();
					f_min= std::min(f_min,objfunc.value());
				}
				//anneal(objfunc,n_step);
				T=gamma*T;
			};
		};
		
		return f_min;
	};
	
	//----------------------------------------------------------------------------------------------
	double const & dryrun(objective_function & objfunc, int const & n_step)
	{
		//Dry run n_step steps
		for (int k= 0; k<n_step; k++)
		{	
			double df= objfunc.step();
			if (df>0)
				objfunc.reverse();
			f_min= std::min(f_min,objfunc.value());
		};
		return f_min;
	};
};
//==================================================================================================
#endif
