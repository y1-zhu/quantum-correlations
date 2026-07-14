using LinearAlgebra;
using SCS;

using Convex;

function RobustnessToBoundary(rho::MultiState)
	xx=1.0.-prod(rho.dims).*eigvals(Hermitian(rho.mat));
	xx= filter(x -> x>0, xx); #Can be reconsiderred!!!
	return 1.0/maximum(xx);
end;

function SepStateRobustnessToBoundary(rho::MultiState,sep_state::MultiState=TensorStates(RandomMultiState([rho.dims[1]]),RandomMultiState([rho.dims[2]])))
	t=Variable(1,Positive());
	S = ComplexVariable(prod(rho.dims),prod(rho.dims));
	problem= maximize(t);
	problem.constraints+= (t*rho.mat+(1.0-t)*sep_state.mat == S);
	problem.constraints+= S in :SDP; 
	solve!(problem,() -> SCS.Optimizer(LOG=0));
    	x = problem.optval;
	return x; 
end;
