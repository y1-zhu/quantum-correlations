using LinearAlgebra
using SCS;
using Convex;
#using JLD;
#using Polyhedra;
#using CDDLib;
using MathOptInterface;

function BlochPolytope(s::String="ipolytope_cvr_d2_L92_u092.jld")
    return load("../../data/"*s)["allProjections"];
end
LoadBlochPolytope=BlochPolytope;
export LoadBlochPolytope;

function RandomBlochPolytope(d::Int=2, L::Int=100)
    QPolytope= Array{Array{Complex{Float64},2},1}(undef,L);
    for k in 1:L
        V=randn(Complex{Float64},d);
        V/=norm(V);
        QPolytope[k]=V*conj(transpose(V));
    end
    return QPolytope;
end

function RobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope,num_iter::Int64=1,solver=SCS.Optimizer, convergence_recognition::Bool=true, convergence_accuracy::Float64=10^-4, silent::Bool=true)
	x = 0.0;
	c = 0;
	d = 1 
	Z=Vector{Matrix{ComplexF64}};
	ResPolytope=Array{Matrix{ComplexF64}}(undef,0);
	if (num_iter == 1)
		ResPolytope = QPolytope; 
	end; 
	while(c < num_iter) 
		U = FirstRobustnessToSeparabilityByBlochPolytope(rho,QPolytope,solver, silent);
		if (c>0 && convergence_recognition)
			if (abs(U[1] - x) <= convergence_accuracy)
				break
			end
		end
		x = U[1];
		print("after ",c," iterations: ",x,"\n")
		Z = U[2]; 
		if (num_iter > 1 && c < (num_iter - 1))
			A = PolytopeOptimizationForWhiteRobustness(rho,Z,solver, silent)
			QPolytope = A[1];
			ResPolytope = A[2]; 
		end; 
		d += 1
		c += 1
	end; 
	return x, Z, QPolytope, ResPolytope, d;
end; 

function FirstRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope,solver=SCS.Optimizer, silent::Bool=true)#::Float64
	    x = 0.0;
	    N=length(QPolytope);
	    Z=Vector{Matrix{ComplexF64}}(undef,N);
	    t=Variable(1,Positive());
	    X=Array{Any,1}(undef,N);
	    for k in 1:N
		X[k]= HermitianSemidefinite(rho.dims[2],rho.dims[2]);
	    end;
   	    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
	    for k in 1:N
		S+= kron(QPolytope[k],X[k]);
	    end;
	    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I; #Identity
	    problem= Convex.maximize(t);
	    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S);
	    solve!(problem, solver, silent_solver=silent);
	    x = problem.optval;
            for i in 1:N
                Z[i]= evaluate(X[i]);
		if (real(tr(Z[i])) < 1.0e-3)
			Z[i] = RandomBlochPolytope(rho.dims[2],1)[1];
		end;
	    end;
	    push!(Z,(1/rho.dims[2])*(zeros(Complex{Float64},rho.dims[2],rho.dims[2])+I));
    	    return x,Z;
end;

# optimizes the polytope (in bipartite systems)
function PolytopeOptimizationForWhiteRobustness(rho::MultiState,G::Vector{Matrix{ComplexF64}},solver=SCS.Optimizer, silent::Bool=true)
    t=Variable(1,Positive()); 
    rho_perm = PermuteSystems(rho,1,2);
    N=length(G);
    sigma=Array{Any,1}(undef,N);
    for i in 1:N
	sigma[i] = HermitianSemidefinite(rho.dims[1],rho.dims[1]);
    end; 
    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
    for i in 1:N
	S+= kron(G[i],sigma[i]);
    end;
    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I;
    problem= Convex.maximize(t);
    problem.constraints+= ((t*rho_perm.mat+(1.0-t)/prod(rho.dims)*II) == S);
    solve!(problem,solver, silent_solver=silent);
    x = problem.optval;
    #print("value: ",x,"\n")
    sol = x;
    QPolytope= Vector{Matrix{ComplexF64}}(undef,N);
    ResPolytope = Vector{Matrix{ComplexF64}}(undef,0);
    for i in 1:N
	QPolytope[i]= evaluate(sigma[i]);
	#print("trace: ",real(tr(kron(G[i],QPolytope[i])))); 
	if (real(tr(kron(G[i],QPolytope[i]))) > 0.01)
		u = copy(QPolytope[i]);
		u = (1/tr(u))*u;
		#print("hallo","\n"); 
		push!(ResPolytope,u);
	end;
	if (real(tr(QPolytope[i])) < 1.0e-3)
		QPolytope[i] = RandomBlochPolytope(rho.dims[1],1)[1];
	end;
        if (tr(QPolytope[i])==0.0)
            print("some trace is zero");
        end;
    end;
    return QPolytope,ResPolytope; 
end;

###polytope adaptions for four qubits -- to be checked (use 15-30 vertices)
## implementation maybe questionable
function FourQubitRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope_1, QPolytope_2,num_iter::Int64=1,solver=SCS.Optimizer,convergence_recognition::Bool=true, convergence_accuracy::Float64=10^-4, silent::Bool=true)
	x = 0.0;
	c = 0; 
	N = length(QPolytope_1); 
	M = length(QPolytope_2);
	if (N != M)
		println("Error: Polytopes do not have the same number of vertices")
	end
	mat = rho.mat;
	rho = MultiState(mat,[4,4]); 
	TensorPolytope = Array{Any,1}(undef,0)
	for i in 1:N
		push!(TensorPolytope,kron(QPolytope_1[i],QPolytope_2[i])) 
	end; 
	Z=Vector{Matrix{ComplexF64}};
	ResPolytope=Array{Matrix{ComplexF64}}(undef,0);
	if (num_iter == 1)
		ResPolytope = TensorPolytope; 
	end; 
	while(c < num_iter) 
		U = FourQubitFirstRobustnessToSeparabilityByBlochPolytope(rho,TensorPolytope,solver, silent)
		if (c>0 && convergence_recognition)
			if (abs(U[1] - x) <= convergence_accuracy)
				break
			end
		end
		x = U[1];
		print("after ",c," iterations: ",x,"\n")
		Z = U[2]; 
		if (num_iter > 1 && c < (num_iter - 1))
			A = FourQubitPolytopeOptimizationForWhiteRobustness(rho,Z,solver, silent)
			TensorPolytope = A[1];
			ResPolytope = A[2]; 
		end; 
		c+=1;
	end; 
	return x, Z, QPolytope_1, ResPolytope;
end;

function FourQubitFirstRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope,solver=SCS.Optimizer, silent::Bool=true)#::Float64
	    x = 0.0;
	    N=length(QPolytope);
	    Z=Vector{Matrix{ComplexF64}}(undef,N);
	    t=Variable(1,Positive());
	    X=Array{Any,1}(undef,N);
	    for k in 1:N
		X[k]= HermitianSemidefinite(rho.dims[2],rho.dims[2]);
	    end;
   	    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
	    for k in 1:N
		S+= kron(QPolytope[k],X[k]);
	    end;
	    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I; #Identity
	    problem= Convex.maximize(t);
	    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S);
	    for i in 1:N
	    	problem.constraints+= partialtranspose(X[i],1,[2,2]) in :SDP
	    end;
	    solve!(problem,solver, silent_solver = silent);
	    x = problem.optval;
	    i = 1;
	    #test = 0;
	    for i in 1:N
	    	Z[i]= evaluate(X[i]);
	    end; 
            while (i<N+1)
		if (real(tr(Z[i])) < 1.0e-3)
			Z[i] = kron(RandomBlochPolytope(2,1)[1],RandomBlochPolytope(2,1)[1]);
		end;
		i += 1; 
		#print("i = ",i,"\n");
	    end;
	    #print(test,"\n"); 
	    push!(Z,(1/rho.dims[2])*(zeros(Complex{Float64},rho.dims[2],rho.dims[2])+I));
    	    return x,Z;
end;


function FourQubitPolytopeOptimizationForWhiteRobustness(rho::MultiState,G::Vector{Matrix{ComplexF64}},solver=SCS.Optimizer, silent::Bool=true)
    t=Variable(1,Positive()); 
    rho_perm = PermuteSystems(rho,1,2);
    N=length(G);
    sigma=Array{Any,1}(undef,N);
    for i in 1:N
	sigma[i] = HermitianSemidefinite(rho.dims[1],rho.dims[1]);
    end; 
    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
    for i in 1:N
	S+= kron(G[i],sigma[i]);
    end;
    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I;
    problem= Convex.maximize(t);
    problem.constraints+= ((t*rho_perm.mat+(1.0-t)/prod(rho.dims)*II) == S);
    for i in 1:N
	problem.constraints+= partialtranspose(sigma[i],1,[2,2]) in :SDP
    end;
    solve!(problem,solver, silent_solver = silent)
    x = problem.optval;
    #print("value: ",x,"\n")
    sol = x;
    QPolytope= Vector{Matrix{ComplexF64}}(undef,N);
    ResPolytope = Vector{Matrix{ComplexF64}}(undef,0);
    i = 1;
    for i in 1:N
	QPolytope[i]= evaluate(sigma[i]);
    end; 
    while (i < N+1)
	#QPolytope[i]= evaluate(sigma[i]);
	#print("trace: ",real(tr(kron(G[i],QPolytope[i])))); 
	#if (real(tr(kron(G[i],QPolytope[i]))) > 0.01)
	#	u = copy(QPolytope[i]);
	#	u = (1/tr(u))*u;
	#	#print("hallo","\n"); 
	#	push!(ResPolytope,u);
	#end;
	if (real(tr(QPolytope[i])) < 1.0e-3)
		#deleteat!(QPolytope, i);
		#N -= 1;
		#i -= 1; 
		QPolytope[i] = kron(RandomBlochPolytope(2,1)[1],RandomBlochPolytope(2,1)[1]);
	end;
        if (tr(QPolytope[i])==0.0)
            print("some trace is zero");
        end;
	i +=1;
    end;
    return QPolytope,ResPolytope; 
end;

## takes 2 polytopes(should have the same length) 
function ThreeQutritRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope1,QPolytope2,num_iter::Int64=1, solver=SCS.Optimizer, convergence_recognition::Bool=true, convergence_accuracy::Float64=10^-4, silent::Bool=true)
	x = 0.0;
	c = 0; 
	N = length(QPolytope1); 
	TensorPolytope = Array{Any,1}(undef,0)
	for i in 1:N	
		push!(TensorPolytope,kron(QPolytope1[i],QPolytope2[i])); 	
	end; 
	Z=Vector{Matrix{ComplexF64}};
	D=TensorPolytope;
	#ResPolytope=Array{Matrix{ComplexF64}}(undef,0);
	#if (num_iter == 1)
	#	ResPolytope = TensorPolytope; 
	#end; 
	while(c < num_iter) 
		A = ThreeQutritFirstRobustnessToSeparabilityByBlochPolytope(rho,D,1,solver, silent);
		if (c>0 && convergence_recognition)
			if (abs(A[1] - x) <= convergence_accuracy)
				break
			end
		end
		x = A[1];
		print("after ",c," iterations: ",x,"\n")
		Z = A[2]; 
		part=Vector{Matrix{ComplexF64}}(undef,0);
		for i in 1:length(D)
			push!(part,partialtrace(D[i],2,[3,3]))
		end; 
		D=Vector{Matrix{ComplexF64}}(undef,0);
		for i in 1:length(Z)
			push!(D,kron(part[i],Z[i]));
		end; 
		push!(D,(1/9)*(zeros(Complex{Float64},9,9)+I)); 
		if (num_iter > 1 && c < (num_iter - 1))
			A = ThreeQutritFirstRobustnessToSeparabilityByBlochPolytope(rho,D,2,solver, silent);
			Z = A[2]; 	
			part=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(D)
				push!(part,partialtrace(D[i],1,[3,3]))
			end; 
			D=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(Z)
				push!(D,kron(part[i],Z[i]));
			end; 
			push!(D,(1/9)*(zeros(Complex{Float64},9,9)+I));
			#TensorPolytope = A[1];
			#ResPolytope = A[2]; 
			A = ThreeQutritFirstRobustnessToSeparabilityByBlochPolytope(rho,D,3,solver, silent);
			Z = A[2]; 	
			part=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(D)
				push!(part,partialtrace(D[i],1,[3,3]))
			end; 
			D=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(Z)
				push!(D,kron(Z[i],part[i]));
			end; 
			push!(D,(1/9)*(zeros(Complex{Float64},9,9)+I));
		end; 
		TensorPolytope = D; 
		c+=1;
	end; 
	#return x, Z, QPolytope, ResPolytope;
	return x, Z, D; 
end;

function ThreeQutritFirstRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope,mode::Int64,solver=SCS.Optimizer, silent::Bool=true)
	    x = 0.0;
            if (mode == 2)
		rho = PermuteSystems(rho,2,3);
	    end;
	    if (mode == 3)
		rho = PermuteSystems(rho,1,3);
	    end;
            mat = rho.mat; 
	    rho = MultiState(mat,[9,3]); 
	    N=length(QPolytope);
	    Z=Vector{Matrix{ComplexF64}}(undef,N);
	    t=Variable(1,Positive());
	    X=Array{Any,1}(undef,N);
	    for k in 1:N
		X[k]= HermitianSemidefinite(rho.dims[2],rho.dims[2]);
	    end;
   	    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
	    for k in 1:N
		S+= kron(QPolytope[k],X[k]);
	    end;
	    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I; #Identity
	    problem= Convex.maximize(t);
	    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S);
	    solve!(problem,solver, silent_solver=silent)
	    x = problem.optval;
            for i in 1:N
                Z[i]= evaluate(X[i]);
		if (real(tr(Z[i])) < 1.0e-3)
			Z[i] = RandomBlochPolytope(3,1)[1];
		end;
	    end;
	    #push!(Z,(1/rho.dims[2])*(zeros(Complex{Float64},rho.dims[2],rho.dims[2])+I));
    	    return x,Z;
end;

## 5 qubits (in progress)
function FiveQubitsRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope1,QPolytope2,QPolytope3,QPolytope4,num_iter::Int64=1,solver=SCS.Optimizer, convergence_recognition::Bool=true, convergence_accuracy::Float64=10^-4, silent::Bool=true)
	x = 0.0;
	c = 0; 
	N = length(QPolytope1); 
	TensorPolytope = Array{Any,1}(undef,0)
	for i in 1:N	
		push!(TensorPolytope,kron(kron(QPolytope1[i],QPolytope2[i]),kron(QPolytope3[i],QPolytope4[i]))); 	
	end; 
	Z=Vector{Matrix{ComplexF64}};
	D=TensorPolytope;
	push!(D,(1/16)*(zeros(Complex{Float64},16,16)+I));
	while(c < num_iter) 
		A = FiveQubitsFirstRobustnessToSeparabilityByBlochPolytope(rho,D,0,solver, silent);
		if (c>0 && convergence_recognition)
			if (abs(A[1] - x) <= convergence_accuracy)
				break
			end
		end
		x = A[1];
		print("after ",c," iterations: ",x,"\n")
		Z = A[2];
		part=Vector{Matrix{ComplexF64}}(undef,0);
		for i in 1:length(D)
			push!(part,partialtrace(D[i],4,[2,2,2,2]))
		end; 
		D=Vector{Matrix{ComplexF64}}(undef,0);
		for i in 1:length(Z)
			push!(D,kron(part[i],Z[i]));
		end; 
		push!(D,(1/16)*(zeros(Complex{Float64},16,16)+I));
		if (num_iter > 1 && c < (num_iter - 1))
			A = FiveQubitsFirstRobustnessToSeparabilityByBlochPolytope(rho,D,4,solver, silent);
			x = A[1];
			print(x,"\n");
			Z = A[2];
			part=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(D)
				push!(part,partialtrace(D[i],3,[2,2,2,2]))
			end; 
			D=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(Z) 
				push!(D,kron(part[i],Z[i]));
			end; 
			push!(D,(1/16)*(zeros(Complex{Float64},16,16)+I));		

			A = FiveQubitsFirstRobustnessToSeparabilityByBlochPolytope(rho,D,3,solver, silent);
			x = A[1];
			print(x,"\n");
			Z = A[2];
			part=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(D)
				push!(part,partialtrace(D[i],2,[2,2,2,2]))
			end; 
			D=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(Z)
				push!(D,kron(part[i],Z[i]));
			end; 
			push!(D,(1/16)*(zeros(Complex{Float64},16,16)+I));	

			A = FiveQubitsFirstRobustnessToSeparabilityByBlochPolytope(rho,D,2,solver, silent);
			x = A[1];
			print(x,"\n");
			Z = A[2];
			part=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(D)
				push!(part,partialtrace(D[i],1,[2,2,2,2]))
			end; 
			D=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(Z)
				v = kron(part[i],Z[i]);
				push!(D,v);
			end; 
			push!(D,(1/16)*(zeros(Complex{Float64},16,16)+I));

			A = FiveQubitsFirstRobustnessToSeparabilityByBlochPolytope(rho,D,1,solver, silent);
			x = A[1];
			print(x,"\n");
			Z = A[2];
			part=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(D)
				push!(part,partialtrace(D[i],1,[2,2,2,2]));
			end; 
			D=Vector{Matrix{ComplexF64}}(undef,0);
			for i in 1:length(Z)
				#fac_1 = partialtrace(part[i],1,[4,2]);
				#fac_3 = partialtrace(part[i],2,[2,4]);
				#fac_2 = partialtrace(partialtrace(part[i],1,[2,4]),2,[2,2]);
				#l = kron(fac_1,kron(fac_2,fac_3));
				v = kron(Z[i],part[i]); 
				#u = kron(part[i],Z[i]);
				#y = partialtrace(u,1,[4,4]);
				#fac_1 = kron(partialtrace(y,1,[2,2]),partialtrace(y,2,[2,2]));
				#z = partialtrace(u,2,[4,4]);
				#fac_2 = kron(partialtrace(z,1,[2,2]),partialtrace(z,2,[2,2]));
				#u = kron(fac_1,fac_2); 
				#v = kron(u,partialtrace(part[i],1,[2,4]))
				push!(D,v);
			end; 
			rho = PermuteSystems(rho,2,4); 
			push!(D,(1/16)*(zeros(Complex{Float64},16,16)+I));
		end;
		c +=1; 
	end;
	return x, Z;
end;

function FiveQubitsFirstRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope,mode::Int64,solver=SCS.Optimizer, silent::Bool=true)
	x = 0.0;
	if (mode == 1)
		rho = PermuteSystems(rho,1,5);
		rho = PermuteSystems(rho,2,4); 
	end;
	if (mode == 2)
		rho = PermuteSystems(rho,2,5);
		rho = PermuteSystems(rho,3,4);
	end;
	if (mode == 3)
		rho = PermuteSystems(rho,3,5);
	end;
	if (mode == 4)
		rho = PermuteSystems(rho,4,5);
	end;
	mat = rho.mat; 
	rho = MultiState(mat,[16,2]); 
	N=length(QPolytope);
	Z=Vector{Matrix{ComplexF64}}(undef,N);
	t=Variable(1,Positive());
	X=Array{Any,1}(undef,N);
	for k in 1:N
		X[k]= HermitianSemidefinite(rho.dims[2],rho.dims[2]);
	end;
   	S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
	for k in 1:N
		S+= kron(QPolytope[k],X[k]);
	end;
	II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I; #Identity
	problem= Convex.maximize(t);
	problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S);
	solve!(problem,solver, silent_solver = silent)
	x = problem.optval;
        for i in 1:N
                Z[i]= evaluate(X[i]);
		if (real(tr(Z[i])) < 1.0e-3)
			Z[i] = RandomBlochPolytope(2,1)[1];
		end;
	end;
    	return x,Z;
end;


function DualRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope,solver=SCS.Optimizer) 
	d = prod(rho.dims); 
	N=length(QPolytope);
	III= zeros(Complex{Float64},prod(rho.dims[2]),prod(rho.dims[2]))+I;
	Y_0 = ComplexVariable(d,d);
	t=Variable(1);
	problem= minimize(t);
	problem.constraints += (t == tr(rho.mat * Y_0));
	problem.constraints += (Y_0 == adjoint(Y_0)); 
	problem.constraints += (tr(Y_0) == d); # white noise robustness
	for i in 1:N
		problem.constraints += (partialtrace(Y_0*(kron(QPolytope[i],III)),1,[rho.dims[1],rho.dims[2]]) in :SDP);
	end;
	solve!(problem,() -> SCS.Optimizer(LOG=0));
    	x = problem.optval;
	x = -x;
	x = (1.0/(1.0+x)); 
	print("value (FSEP): ",x,"\n"); 
	U = evaluate(Y_0);
	return x, U; 		
end;

function StateRobustnessToSeparabilityByBlochPolytope(rho::MultiState,QPolytope,num_iter::Int64=1,sep_state::MultiState=TensorStates(RandomMultiState([rho.dims[1]]),RandomMultiState([rho.dims[2]])),solver=SCS.Optimizer, convergence_recognition::Bool=true, convergence_accuracy::Float64=10^-4, silent::Bool=true)#::Float64
    c=0;
    x = 0.0;
    red_state1 = ReducedState(sep_state,1).mat;
    push!(QPolytope,red_state1);
    red_state2 = ReducedState(sep_state,2).mat;
    N=length(QPolytope);
    Z=Vector{Matrix{ComplexF64}}(undef,N);
    while (c < num_iter)
	    t=Variable(1,Positive());
	    X=Array{Any,1}(undef,N);
	    for k in 1:N
		X[k]= HermitianSemidefinite(rho.dims[2],rho.dims[2]);
	    end;
	    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
	    for k in 1:N
		S+= kron(QPolytope[k],X[k]);
	    end;
	    problem= Convex.maximize(t);
	    problem.constraints+= ((t*rho.mat+(1.0-t)*sep_state.mat) == S);
	    solve!(problem,solver, silent_solver=silent);
            Z=Vector{Matrix{ComplexF64}}(undef,N);
            for i in 1:N
                Z[i]= evaluate(X[i]);
		if (real(tr(Z[i])) < 1.0e-5)
			Z[i] = RandomBlochPolytope(rho.dims[2],1)[1];
		end;
	    end;
	    push!(Z,red_state2);
	    if (num_iter > 1 && c < (num_iter-1))
            	QPolytope = PolytopeOptimizationForSepStateRobustness(rho,Z,sep_state, solver, silent)
				N = length(QPolytope); 
	    end; 
		if (c>0 && convergence_recognition)
			if ((problem.optval - x) <= convergence_accuracy)
				break
			end
		end
        x = problem.optval;
	    c += 1;
    end;
    return x, Z, QPolytope;
end;

function PolytopeOptimizationForSepStateRobustness(rho::MultiState,G::Vector{Matrix{ComplexF64}},sep_state::MultiState, solver=SCS.Optimizer, silent::Bool=true)
    t=Variable(1,Positive());
    N=length(G);
    sigma=Array{Any,1}(undef,N);
    for i in 1:N
	sigma[i] = HermitianSemidefinite(rho.dims[1],rho.dims[1]);
    end; 
    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
    for i in 1:N
        S+= kron(sigma[i],G[i]); 
    end;
    problem= Convex.maximize(t);
    problem.constraints+= ((t*rho.mat+(1.0-t)*sep_state.mat) == S);
    solve!(problem, solver, silent_solver=silent);
    QPolytope= Vector{Matrix{ComplexF64}}(undef,N);
    for i in 1:N
	QPolytope[i]= evaluate(sigma[i]);
	if (real(tr(QPolytope[i])) < 1.0e-5)
		QPolytope[i] = RandomBlochPolytope(rho.dims[1],1)[1];
	end;
        if (tr(QPolytope[i])==0.0)
            print("some trace is zero");
        end;
    end;
    return QPolytope; 
end;



# produces an equivalent state whose reduced state is equal to the average point of the polytope
function LocalFilterToCOM(rho::MultiState,QPolytope)
    isrhoA=inv(sqrt(ReducedState(rho,1).mat));
    N = length(QPolytope);
    sigma_bar = zeros(Complex{Float64},rho.dims[1],rho.dims[1]);
    for i in 1:N
        sigma_bar += QPolytope[i] 
    end;
    sigma_bar = (1/N)*sigma_bar;
    #print(MultiState(sigma_bar,[rho.dims[1]]));
    sigma_bar = sqrt(sigma_bar);
    II= zeros(rho.dims[2],rho.dims[2])+I;
    mat= kron(sigma_bar*isrhoA,II)*rho.mat*kron(isrhoA*sigma_bar,II);
    mat/=tr(mat);
    return MultiState(mat,deepcopy(rho.dims));
end;

function LocalFilterCOM_and_Mixed(rho::MultiState,QPolytope,iter::Int=20)
    i=0;
    N = length(QPolytope);
    sigma_bar = zeros(Complex{Float64},rho.dims[1],rho.dims[1]);
    for i in 1:N
        sigma_bar += QPolytope[i] 
    end;
    sigma_bar = (1/N)*sigma_bar;
    #print(eigvals(MultiState(sigma_bar,[rho.dims[1]]).mat),"\n");
    while (i < iter)
        rho = Canonicalise(rho,2);
        rho = LocalFilterToCOM(rho,QPolytope);
        i += 1;
    end;
    #print(eigvals(ReducedState(rho,1).mat),"\n");
    #print(eigvals(ReducedState(rho,2).mat),"\n");
    return rho; 
end; 

##robustness to genuine multipartite entanglement for 3 qubits 
function MRobustnessToGenuineEntanglementByBlochPolytope(rho::MultiState,QPolytope::Array{Array{Complex{Float64},2},1},solver=SCS.Optimizer(LOG=0), silent::Bool=true)::Float64 #Currently for three qubits
    #QPolytope= BlochPolytope(s);
    #QPolytope= RandomBlochPolytope(3,200);
    if ((rho.dims[1] != 2) | (rho.dims[2] != 2) | (rho.dims[3] != 2))
        print("***Error: The system are not three qubits!!!");
    end;

    t=Variable(1,Positive());
    N=length(QPolytope);
    X1=Array{Any,1}(undef,N);
    X2=Array{Any,1}(undef,N);
    X3=Array{Any,1}(undef,N);
    for k in 1:N
        X1[k]= HermitianSemidefinite(4,4);
        X2[k]= HermitianSemidefinite(4,4);
        X3[k]= HermitianSemidefinite(4,4);
    end;
    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
    for k in 1:N
        S+= kron(QPolytope[k],X1[k]);
        S+= kron(X3[k],QPolytope[k]);
        
        C= zeros(Complex{Float64},8,8);
	    for l in 1:8
			kk= LocalIndex(l,[2,2,2]); tmp=kk[1]; kk[1]=kk[2]; kk[2]= tmp;
			p=LinearIndex(kk,[2,2,2]);
            C[p,l]=1.0;
		end;

        S+=  C*kron(QPolytope[k],X2[k])*C; 

        
	end;
    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I;
    problem= maximize(t);
    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S);
    solve!(problem, solver, silent_solver=silent);
    print(problem.status,"\n");
    print("robustness to white noise: ",problem.optval, "\n"); 
    return problem.optval;
end;

## full separability
function MRobustnessToFullSeparabilityByBlochPolytope(rho::MultiState,QPolytope::Array{Array{Complex{Float64},2},1},solver=SCS.Optimizer, silent::Bool=true)::Float64 #Currently for three qubits
    #QPolytope= BlochPolytope(s);
    #QPolytope= RandomBlochPolytope(3,200);
    if ((rho.dims[1] != 2) | (rho.dims[2] != 2) | (rho.dims[3] != 2))
        print("***Error: The system are not three qubits!!!");
    end;
    t=Variable(1,Positive());
    N=length(QPolytope);
    X=Array{Any,1}(undef,N);
    for k in 1:N
        X[k]= HermitianSemidefinite(4,4);
    end;
    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
    for k in 1:N
        S+= kron(QPolytope[k],X[k]);        
	end;

    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I;
    problem= maximize(t);
    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S);
    for k in 1:N
        problem.constraints+= partialtranspose(X[k],1,[2,2]) in :SDP; 
    end 
    solve!(problem,solver, silent_solver=silent);
    print(problem.status,"\n");
    return problem.optval;
end;

function AbsoluteRobustnessOfEntanglementByBlochPolytope(rho::MultiState,QPolytope,num_iter::Int64=1,solver=SCS.Optimizer, convergence_recognition::Bool=true, convergence_accuracy::Float64=10^-4, silent::Bool=true)
	x = 0.0;
	d = 1
	c = 0; 
	N = length(QPolytope); 
	Z=Vector{Matrix{ComplexF64}};
	while(c < num_iter) 
		U = MRobustnessOfEntanglementByBlochPolytope(rho,QPolytope,solver, silent);
		if (c>0 && convergence_recognition)
			if (abs(U[1] - x) <= convergence_accuracy)
				break
			end
		end
		x = U[1];
		print("after ",c," iterations: ",x,"\n")
		Z = U[2]; 
		sep_state = U[3]; 
		if (num_iter > 1 && c < (num_iter - 1))
			QPolytope = MRobustnessOfEntanglementAdaption(rho,Z,sep_state,solver, silent);
		end; 
		c += 1;
		d += 1
	end; 
	return x, Z, QPolytope, d;
end

## to do: code adaption algorithm and test it 
## robustness of entanglement with respect to general separable noise, (maybe not perfectly working in current status) 
function MRobustnessOfEntanglementByBlochPolytope(rho::MultiState,QPolytope::Array{Array{Complex{Float64},2},1},solver=SCS.Optimizer, silent::Bool=true)
    t=Variable(1,Positive());
	N=length(QPolytope);
	X=Array{Any,1}(undef,N);
	Y=Array{Any,1}(undef,N);
	for k in 1:N
		X[k]= HermitianSemidefinite(rho.dims[2],rho.dims[2]);
		Y[k]= HermitianSemidefinite(rho.dims[2],rho.dims[2]);
	end;
	S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
	R=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
	for k in 1:N
		S+= kron(QPolytope[k],X[k]-Y[k]);
		R+= kron(QPolytope[k],Y[k]); 
	end;
	problem= Convex.maximize(t);
	problem.constraints+= (t*rho.mat == S);
	problem.constraints+= (t+tr(R) == 1.0);
	solve!(problem,solver, silent_solver=silent);
	print(problem.status,"\n");
	x = problem.optval;
	results1= Vector{Matrix{ComplexF64}}(undef,N);
	results2= Vector{Matrix{ComplexF64}}(undef,N);
	res= Vector{Matrix{ComplexF64}}(undef,N);
	for i in 1:N
		results1[i] = evaluate(X[i]);
		if (real(tr(results1[i])) < 1.0e-2)
			results1[i] = RandomBlochPolytope(rho.dims[2],1)[1];
		end;
		results2[i] = evaluate(Y[i]);
		#if ((real(tr(results1[i]))/(1.0-x)) < 1.0e-3)
	#		results2[i] = RandomBlochPolytope(rho.dims[2],1)[1];
		#end;
		res[i] = results1[i] - results2[i]
	end;
	## sepstate = sum(sigma_i \otimes y(i))  / (1 - x) 
	sep_state = zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
	for i in 1:N
		sep_state += kron(QPolytope[i],results2[i]);
	end;
	sep_state /= 1-x
	sep_state = MultiState(sep_state,rho.dims);
    return x, results1, sep_state;   
end; 



function MRobustnessOfEntanglementAdaption(rho::MultiState,G::Vector{Matrix{ComplexF64}},sep_state::MultiState,solver=SCS.Optimizer, silent::Bool=true)
    t=Variable(1,Positive()); 
    rho_perm = PermuteSystems(rho,1,2);
	sep_state_perm = PermuteSystems(sep_state,1,2)
    N=length(G);
    sigma=Array{Any,1}(undef,N);
    for i in 1:N
		sigma[i] = HermitianSemidefinite(rho.dims[1],rho.dims[1]);
    end; 
    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims));
    for i in 1:N
		S+= kron(G[i],sigma[i]);
    end;
    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I;
    problem= Convex.maximize(t);
    problem.constraints+= (t*rho_perm.mat + (1-t)*sep_state_perm.mat == S);
    solve!(problem,solver, silent_solver=silent);
    x = problem.optval;
    QPolytope= Vector{Matrix{ComplexF64}}(undef,N);
    for i in 1:N
		QPolytope[i]= evaluate(sigma[i]);
		if (real(tr(QPolytope[i])) < 10^(-2))
			QPolytope[i] = RandomBlochPolytope(rho.dims[1],1)[1]
		end; 
		#QPolytope[i] /= tr(QPolytope[i]);
	end;
    return QPolytope; 
end;

function GeneralRobustnessOfEntanglementByBlochPolytope(rho::MultiState,QPolytope::Array{Array{Complex{Float64},2},1},num_iter::Int64=1,solver=SCS.Optimizer, convergence_recognition::Bool=true, convergence_accuracy::Float64=10^-4, silent::Bool=true)
	x = 0.0
	c = 0
	d = 1 
	Z=Vector{Matrix{ComplexF64}};
	ResPolytope=Array{Matrix{ComplexF64}}(undef,0);
	if (num_iter == 1)
		ResPolytope = QPolytope; 
	end; 
	while(c < num_iter) 
		U = GeneralRobustnessSDP(rho,QPolytope,solver,1, silent)
		if (c>0 && convergence_recognition)
			if (abs(U[1] - x) <= convergence_accuracy)
				break
			end
		end
		x = U[1];
		print("after ",c," iterations: ",x,"\n")
		Z = U[2]; 
		if (num_iter > 1 && c < (num_iter - 1))
			A = GeneralRobustnessSDP(rho,Z,solver,2, silent)
			QPolytope = A[2];
		end; 
		d += 1
		c += 1
	end; 
	return x, Z, QPolytope, d
end; 





# todo: implement adaption
function GeneralRobustnessSDP(rho::MultiState, QPolytope::Array{Array{Complex{Float64},2},1},solver=SCS.Optimizer,mode::Int64=1, silent::Bool=true)
	if (mode==2)
		rho = PermuteSystems(rho,1,2)
	end
	t=Variable(1,Positive())
	Y = HermitianSemidefinite(prod(rho.dims),prod(rho.dims))
	N=length(QPolytope)
	X=Array{Any,1}(undef,N)
	for k in 1:N
		X[k]= HermitianSemidefinite(rho.dims[2],rho.dims[2])
	end
	S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))
	for k in 1:N
		S+= kron(QPolytope[k],X[k])
	end 
	problem= Convex.maximize(t)
	problem.constraints += (t*rho.mat + Y == S)
	problem.constraints += (t+tr(Y) == 1.0)
	solve!(problem, solver, silent_solver=silent)
	print(problem.status,"\n")
	x = problem.optval
	QPolytope= Vector{Matrix{ComplexF64}}(undef,N);
    for i in 1:N
		QPolytope[i]= evaluate(X[i])
		if (real(tr(QPolytope[i])) < 10^-3)
			QPolytope[i] = RandomBlochPolytope(rho.dims[1],1)[1]
		end; 
	end;
	return x, QPolytope
end


function FSEP_dualSDP_seesaw(rho::MultiState,QPolytope::Array{Array{Complex{Float64},2},1})
	d = prod(rho.dims); 
	N=length(QPolytope);
	II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I;
	III= zeros(Complex{Float64},prod(rho.dims[2:3]),prod(rho.dims[2:3]))+I;
	Y_0 = ComplexVariable(d,d);
	#U_0 = ComplexVariable(d,d);
	t=Variable(1);
	Y=Array{Any,1}(undef,N);
	W=Array{Any,1}(undef,N);
	for i in 1:N
		Y[i] = ComplexVariable(prod(rho.dims[2:3]),prod(rho.dims[2:3]));	
		##test for improvement of the approximation of FSEP
		W[i] = ComplexVariable(prod(rho.dims[2:3]),prod(rho.dims[2:3]));
	end;
	problem= minimize(t);
	#problem.constraints += (t == (1/d)*tr(Y_0)); 
	#problem.constraints += (t == 1.0 + tr(rho.mat * Y_0));
	problem.constraints += (t == tr(rho.mat * Y_0));
	#problem.constraints += (U_0 == partialtranspose(Y_0,1,[2,2,2]));
	#problem.constraints += (t == tr(rho.mat * U_0))
	#problem.constraints += (t == (1/d)*tr(Y_0));
	problem.constraints += (Y_0 == adjoint(Y_0)); 
	#problem.constraints += (-tr((rho.mat-(1/d)*II)*Y_0)-1 in :SDP);
	#problem.constraints += (-tr((rho.mat-(1/d)*II)*Y_0) == 1.0);
	#problem.constraints += ((II-Y_0) in :SDP); # for generalized robustness and negativity
	#problem.constraints += (Y_0 + II in :SDP); # best separable approximation
	#problem.constraints += (Y_0 in :SDP); # for negativity
	problem.constraints += (tr(Y_0) == d); # white noise robustness
	#problem.constraints += (tr(Y_0) == 0.0); 
	for i in 1:N
		problem.constraints += 	(Y[i] == adjoint(Y[i])); 	
		problem.constraints += 	(W[i] == adjoint(W[i]));
		problem.constraints += (-partialtranspose(Y[i],1,[2,2]) in :SDP);
		problem.constraints += (-partialtranspose(W[i],1,[2,2]) in :SDP);
		#problem.constraints += ((ReducedState(MultiState(Y_0*(kron(QPolytope[i],III)),[2,4]),2).mat + Y[i]) in :SDP);
		problem.constraints += (partialtrace(Y_0*(kron(QPolytope[i],III)),1,[2,4])+Y[i] in :SDP);
		problem.constraints += (partialtrace(Y_0*(kron(III,QPolytope[i])),2,[4,2])+W[i] in :SDP);
	end;
	solve!(problem,() -> SCS.Optimizer(LOG=0));
    	x = problem.optval;
	#print("first: ",x,"\n");
	x = -x;
	x = (1.0/(1.0+x)); 
	print("value (FSEP): ",x,"\n"); 
	U = evaluate(Y_0);
	return U,x; 			
end; 

function RobustnessToFullSeparability(rho::MultiState,QPolytope::Array{Array{Complex{Float64},2},1})
	if ((rho.dims[1] != 2) | (rho.dims[2] != 2) | (rho.dims[3] != 2))
        print("***Error: The system are not three qubits!!!")
    end;
	d = prod(rho.dims)
	N = length(QPolytope)
	II = zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I
	III = zeros(Complex{Float64},prod(rho.dims[2:3]),prod(rho.dims[2:3]))+I
	t=Variable(1)
	Y_0 = ComplexVariable(d,d)
	Y=Array{Any,1}(undef,N)
	W=Array{Any,1}(undef,N)
	for i in 1:N
		Y[i] = ComplexVariable(prod(rho.dims[2:3]),prod(rho.dims[2:3]))	
		W[i] = ComplexVariable(prod(rho.dims[2:3]),prod(rho.dims[2:3]))
	end;
	problem = minimize(t)
	problem.constraints += (t == tr(rho.mat * Y_0) + 1)
	problem.constraints += (Y_0 == adjoint(Y_0))
	problem.constraints += (tr(Y_0)/d - tr(rho.mat * Y_0) == 1.0)
	for i in 1:N
		problem.constraints += 	(Y[i] == adjoint(Y[i]))	
		problem.constraints += (-partialtranspose(Y[i],1,[2,2]) in :SDP)
		problem.constraints += (partialtrace(Y_0*(kron(QPolytope[i],III)),1,[2,4])+Y[i] in :SDP)
	end;
	solve!(problem,() -> SCS.Optimizer(LOG=0))
	x = problem.optval
	return x
end 
export RobustnessToFullSeparability

function RobustnessToBiSeparability_seesaw(QPolytope::Array{Array{Complex{Float64},2},1},U)
    N=length(QPolytope);
    #epsilon = 0.185;
    t=Variable(1);
    #s=Variable(1,Positive());
    rho = HermitianSemidefinite(8,8); 
    problem= minimize(t);
    X1=Array{Any,1}(undef,N);
    X2=Array{Any,1}(undef,N);
    X3=Array{Any,1}(undef,N);
    #X=Array{Any,1}(undef,N);
    #for k in 1:N
    #    X[k]= HermitianSemidefinite(4,4);
    ###end;
    #S=zeros(Complex{Float64},8,8);
    #for k in 1:N
    #    S+= kron(QPolytope[k],X[k]);
    #end;
    
    #constraints  += rho.mat in :SDP; 
    #for k in 1:N
    #    problem.constraints+= partialtranspose(X[k],1,[2,2]) in :SDP; 
    #end; 

    for k in 1:N
        X1[k]= HermitianSemidefinite(4,4);
        X2[k]= HermitianSemidefinite(4,4);
        X3[k]= HermitianSemidefinite(4,4);
    end;

    S1=zeros(Complex{Float64},8,8);
    S2=zeros(Complex{Float64},8,8);
    S3=zeros(Complex{Float64},8,8);
    
    for k in 1:N
	S1+= kron(QPolytope[k],X1[k]);
        S2+= kron(X2[k],QPolytope[k]);
        C= zeros(Complex{Float64},8,8);
	    for l in 1:8
			kk= LocalIndex(l,[2,2,2]); tmp=kk[1]; kk[1]=kk[2]; kk[2]= tmp;
			p=LinearIndex(kk,[2,2,2]);
                        C[p,l]=1.0;
	    end;
        S3+= C*kron(QPolytope[k],X3[k])*C;
    end;
    II= zeros(Complex{Float64},8,8)+I;
    #problem.constraints += ((rho + (1.0-t)*(1/8)*II) == S);
    problem.constraints+= (rho == S1);
    problem.constraints+= (rho == S2);
    problem.constraints+= (rho == S3);
    problem.constraints+= (tr(rho) == 1.0);
    #problem.constraints += (-tr((rho-(1.0/8.0)*II)*U) == 1.0);
    #problem.constraints+= ()
    #problem.constraints+= (-tr(rho)*1.0-tr(U*rho)+(1.0/8.0)*tr(rho)*tr(U)+epsilon in :SDP); 
    #problem.constraints+= (-tr(rho + (1.0-t)*(1/8)*II)*1.0-tr(U*(rho + (1.0-t)*(1/8)*II))+(1.0/8.0)*tr(rho + (1.0-t)*(1/8)*II)*tr(U)+epsilon in :SDP);
    #problem.constraints+= (-tr(U*(rho + (1.0-t)*(1/8)*II)) in :SDP)
    #problem.constraints+= (-tr(U*rho)-1.0+tr(rho)*(1/8)*tr(U) in :SDP);
    #problem.constraints+= (tr(rho) == t); 
    #problem.constraints += (t == 1.0 + tr(rho*U));
    problem.constraints += (t == tr(rho*U));
    #problem.constraints += (t == tr(rho*partialtranspose(U,1,[2,2,2])));
    solve!(problem,() -> SCS.Optimizer(LOG=0));
    x = problem.optval;
    #t = evaluate(t); 
    x = -x;
    x = (1.0/(1.0+x)); 
    rho_opt = evaluate(rho);
    rho_opt = (1/tr(rho_opt))*rho_opt; 
    print("value: (BSEP) ",x,"\n");
    #print(rho_opt,"\n");
    return rho_opt,x; 
    #print(a);
end; 

function RobustnessToFullBiseparability(rho::MultiState,QPolytope::Array{Array{Complex{Float64},2},1},solver=SCS.Optimizer, silent::Bool=true) #first basic for three qubits
	if ((rho.dims[1] != 2) | (rho.dims[2] != 2) | (rho.dims[3] != 2))
        print("***Error: The system are not three qubits!!!")
    end
    N=length(QPolytope)
    t=Variable(1)
	II= zeros(Complex{Float64},8,8)+I
    problem= maximize(t)
    X1=Array{Any,1}(undef,N)
    X2=Array{Any,1}(undef,N)
    X3=Array{Any,1}(undef,N)
    for k in 1:N
        X1[k]= HermitianSemidefinite(4,4)
        X2[k]= HermitianSemidefinite(4,4)
        X3[k]= HermitianSemidefinite(4,4)
    end
	S1=zeros(Complex{Float64},8,8)
    S2=zeros(Complex{Float64},8,8)
    S3=zeros(Complex{Float64},8,8)
	for k in 1:N
		S1+= kron(QPolytope[k],X1[k])
		S2+= kron(X2[k],QPolytope[k])
		C= zeros(Complex{Float64},8,8)
		for l in 1:8
			kk= LocalIndex(l,[2,2,2]); tmp=kk[1]; kk[1]=kk[2]; kk[2]= tmp;
			p=LinearIndex(kk,[2,2,2])
			C[p,l]=1.0
		end
		S3+= C*kron(QPolytope[k],X3[k])*C
	end
	problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S1)
    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S2)
    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S3)
	solve!(problem, solver, silent_solver=silent)
    x = problem.optval
	println(x)
	return x
end


function SeeSawTest(QPolytope::Array{Array{Complex{Float64},2},1}=BlochPolytope("ipolytope_cvr_d2_L92_u092.jld"),input::MultiState=RandomMultiState([2,2,2]),N_states::Int = 1,N_it::Int=10,N_vert::Int=250,N_poly=1)
	Z = Array{Any,1}(undef,0);
	k=0;
	X=Array{Any,1}(undef,2);
	#QPolytope = RandomBlochPolytope(2,50); 
	while (k<N_states)
	#Threads.@threads for i in 1:N_states
		#A = RandomMultiState([2,2,2]);
		A = input;
		#A = CandidateState();
		print(eigvals(A.mat)); 
		#A = BennettState();
		for i in 1:N_poly
			#QPolytope = RandomBlochPolytope(2,N_vert); 
			value = 1.0;
			matrix = zeros(Complex{Float64},8,8);
			while (true) 
				#A = GHZThreeQubitState();
				U=FSEP_dualSDP_seesaw(A,QPolytope)[1]; 
				#U=FSEP_dualSDP_seesaw(A,BlochPolytope("ipolytope_cvr_d2_L92_u092.jld"));
				#X=RobustnessToBiSeparability_seesaw(BlochPolytope("opolytope_cvr_d2_L92_u092.jld"),U);
				X=RobustnessToBiSeparability_seesaw(QPolytope,U);
				#X = Tensor_Polytope_test(QPolytope,U)
				if (X[2]<0.997)
					value = X[2]; 
					print(value,"\n"); 
					matrix = X[1];
					print(real(eigvals(matrix)),"\n"); 
					display(eigvecs(matrix)[:,1:4]);
					print("\n");
					#for i in 1:length(matrix[1,:])
					#	print(eigvals)
					#end;
					break;		
				else
					A = RandomMultiState([2,2,2]);	
				end; 
			end;
			c = 0;
			while (c<N_it)
				U=FSEP_dualSDP_seesaw(MultiState(X[1],[2,2,2]),QPolytope)[1];
				#U=FSEP_dualSDP_seesaw(MultiState(X[1],[2,2,2]),BlochPolytope("ipolytope_cvr_d2_L162_u162.jld"));
				#print(U[],"\n");
				X=RobustnessToBiSeparability_seesaw(QPolytope,U);
				#X=RobustnessToBiSeparability_seesaw(BlochPolytope("opolytope_cvr_d2_L92_u092.jld"),U);
				#X = Tensor_Polytope_test(QPolytope,U)
				value = X[2]; 
				#print(value,"\n");
				matrix = X[1]; 
				print(real(eigvals(matrix)),"\n");
				#display(EntryMinimization(matrix)); 
				#display(eigvecs(matrix)[:,1:4]);
				print("\n");
				#print(X[2],"\n");
				c += 1; 
			end;
			#label = [i,k];
			print("value for polytope ",i,": ",value,"\n");  
			push!(Z,[value,matrix]);  
		end; 
		k += 1;
	end;
	return Z; 	
end;


