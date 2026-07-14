using LinearAlgebra;
using PolyJuMP;
using SCS;
using DynamicPolynomials;
using SumOfSquares;

function RobustnessToSN2bySOS(rho::MultiState,deg::Int=0)
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
	@polyvar RX1[1:rho.dA] IX1[1:dA] RY1[1:dB] IY1[1:dB] RX2[1:dA] IX2[1:dA] RY2[1:dB] IY2[1:dB];
	@polyvar c1 c2;
	#c1=1.0;
	#c2=1.0;
	RV=c1*(kron(RX1,RY1)-kron(IX1,IY1))+c2*(kron(RX2,RY2)-kron(IX2,IY2));
	IV=c1*(kron(RX1,IY1)+kron(IX1,RY1))+c2*(kron(RX2,RY2)+kron(IX2,IY2));

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
	M=Array{:Any,2,2};
	M[1,1]= differentiate(differentiate(p,c1),c1);
	M[1,2]= differentiate(differentiate(p,c1),c2)/2; #check commutativity assumption
	M[2,1]= M[1,2];
	M[2,2]= differentiate(differentiate(p,c2),c2);
	@constraint(model,cp, p in SOSCone());
	#@constraint(model,cp,M in SOSMatrixCone());
	#pcoeff= sum([RX;IX;RY;IY].^2); #precoefficient
	#@constraint(model,cp,pcoeff^deg*p in SOSCone());

	#Solving the model:
	optimize!(model);

	#Check if the model is solved:
	#if (termination_status(model))
	#	print("Error: The model is not solved!!!!\n");
	#end;
	@show termination_status(model);
	return (1.0+objective_value(model));
end;
