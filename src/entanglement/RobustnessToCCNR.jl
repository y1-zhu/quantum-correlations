
using LinearAlgebra;
using MultiStates;

function HermitianOpBasis(M::Int64=2)  
	A = Array{Any,1}(undef,0);
	for i in 1:M
		mat = zeros(Complex{Float64},M,M);
		mat[i,i] = 1.0;
		push!(A,mat);
		for j in 1:i-1
			mat = zeros(Complex{Float64},M,M);
			mat[j,i] = 1.0/sqrt(2);
			mat[i,j] = 1.0/sqrt(2);
			push!(A,mat);
			mat = zeros(Complex{Float64},M,M);
			mat[j,i] = -(im)*1.0/sqrt(2);
			mat[i,j] = (im)*1.0/sqrt(2);
			push!(A,mat);
		end;
		
	end;
	return A; 
end;


function CCNRvalue(rho::MultiState);
	N = length(rho.dims); 
	if (N > 2)
		print("The state is not bi-partite","\n");
		return;
	end;
	R = rho.dims;
	A = HermitianOpBasis(R[1]);
	B = HermitianOpBasis(R[2]);
	mat = zeros(Complex{Float64},(prod(rho.dims))^2,(prod(rho.dims))^2);
	for i in 1:length(A)
		for j in 1:length(B)
			mat[i,j] = tr(rho.mat*kron(A[i],B[j])); 
		end;
	end;
	Mat = mat*adjoint(mat);
	eigenvalues = eigvals(Mat);
	val = real(sum((eigenvalues.+0im).^(1/2)));
	return val; 		
end; 

function RobustnessByCCNR(rho::MultiState)
	R = rho.dims;
	t = 0.0
	II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I;
	C = 0.0; 
	while (C <= 1.0)
		state = t*rho.mat+(1.0-t)/prod(rho.dims)*II;
		C = CCNRvalue(MultiState(state,R));
		t += 0.01;
	end;
	t -= 0.01;
	return t; 
end; 
