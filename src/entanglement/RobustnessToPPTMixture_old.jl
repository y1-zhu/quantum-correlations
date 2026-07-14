using LinearAlgebra
using SCS;
using Convex;

function RobustnessToPPTMixture(rho::MultiState)::Float64 
    N= length(rho.dims);
	if (N>3)
		print("Error:: Current implementation is only for N<=3! (to be update ssoon) \n")
	end;
	t=Variable(1,Positive());
    X=Array{Any,1}(undef,N);
    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
    for n in 1:N
          X[n]= HermitianSemidefinite(prod(rho.dims),prod(rho.dims));
          S+= X[n];
    end;
    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I;
    problem= maximize(t);
    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S);
    for n in 1:N
        problem.constraints+= partialtranspose(X[n],n,rho.dims) in :SDP;
    end;
    solve!(problem,() -> SCS.Optimizer());
    print(problem.status,"\n");
    return problem.optval;
end;
