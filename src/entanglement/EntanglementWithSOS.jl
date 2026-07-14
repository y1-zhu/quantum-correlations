using LinearAlgebra;
using PolyJuMP;
using SCS;
using DynamicPolynomials;
using SumOfSquares;

function RobustnessToEntanglementBySOS(rho::MultiState,deg::Int=0)
	if (length(rho.dims) != 2)
		print("***Error: this function works only for bipartite state...\n");
		return 0;
	end;

	dA= rho.dims[1];
	dB= rho.dims[2];
	d= prod(rho.dims);
	RRho=real.(rho.mat);
	IRho=imag.(rho.mat);

	#Polynomial variables:
	@polyvar RX[1:dA] IX[1:dA] RY[1:dB] IY[1:dB];
	RV=kron(RX,RY)-kron(IX,IY);
	IV=kron(RX,IY)+kron(IX,RY);

	solver = optimizer_with_attributes(SCS.Optimizer, MOI.Silent() => true);
	model= Model(solver);

	#Real part of W (decision variables)
	RW=@variable(model,[1:d,1:d],Symmetric);

	#Imaginary part of W (decision variables)
	IW=Array{GenericAffExpr{Float64,VariableRef},2}(UndefInitializer(),d,d);
	for k=1:d
		IW[k,k]=0;
	end;
	for k1=1:d
		for k2=(k1+1):d
			IW[k1,k2]=@variable(model);
			IW[k2,k1]=-IW[k1,k2];
		end;
	end;

	#Construct objective function Tr(W*rho.mat)
	obj=tr(RW*RRho-IW*IRho);
	@objective(model,Min,obj);

	#Construct the equality constraint:
	@constraint(model,tr(RW)==d*(1+obj));

	#Construct the separably positive constraint:
	p=RV'*RW*RV-RV'*IW*IV+IV'*RW*IV+IV'*IW*RV;
	pcoeff= sum([RX;IX;RY;IY].^2);
	@constraint(model,cp,pcoeff^deg*p in SOSCone());

	#Solving the model:
	optimize!(model);

	#Check if the model is solved:
	#if (termination_status(model))
	#	print("Error: The model is not solved!!!!\n");
	#end;
	@show termination_status(model);
	return (1.0+objective_value(model));
end;
