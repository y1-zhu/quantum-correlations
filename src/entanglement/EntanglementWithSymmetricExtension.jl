using IterTools;
using Combinatorics;

using LinearAlgebra;
using SCS;
using Convex;

#Pay special attention: julia index starts from 1
function ToAbsoluteIndex(sss::Array{Int},d::Int)::Int
	n=length(sss);
	L=1;
	for k=1:n
		L+=(sss[k]-1)*d^(n-k);
	end;
	return L;
end;

function DecomposedIndex(L::Int,d::Int,n::Int)
	sss=zeros(Int,n);
	for k=1:n
		sss[n-k+1]=mod(L-1,d)+1;
		L=div(L-1,d)+1;
	end;
	return sss;
end;

#Quick test: seems to be correct
# extension always on the second system, it makes also sense when the dimension of the second system is smaller (use rho=PermuteSystems(rho,1,2), if necessary)
function SymmetricKernel(N::Int,d::Int,method="generic")
	if (N==2) #Bipartite method
		#The the symmetric space is simply v_mu= (|i,j>+|j,i>)/sqrt(2)
		M=Int64((d+1)*d/2);
		mu=0;
		V=zeros(Complex{Float64},M,d,d);
		for na=1:d
			for nb=na:d
				mu+=1;
				if (na==nb)
					V[mu,na,nb]=1;
				else
					V[mu,na,nb]=sqrt(0.5);
					V[mu,nb,na]=sqrt(0.5);
				end;
			end;
		end;	
		G=zeros(Complex{Float64},M,M,d,d);
		for mu=1:M
			for nu=1:M
				for j=1:d
					for k=1:d
						for l=1:d
							G[mu,nu,j,k]+= V[mu,j,l]*conj(V[nu,k,l]);
						end;
					end;
				end;
			end;
		end;
	elseif (d==2) #Qubit
		M=N+1; #Number of Diecke's states
		V=zeros(Complex{Float64},M,2^N);
		V[1,1]=1.0;
		for m=2:M
			snum=0;
			for jjj in subsets(1:N,m-1)
				sss= ones(Int,N);
				for j in jjj 
					sss[j]=2;
				end;
				snum+=1;
				V[m,ToAbsoluteIndex(sss,d)]=1.;
			end
			V[m,:]/=sqrt(snum);
		end;
		#print("select=2: V=\n", V, "\n")
		G=zeros(Complex{Float64},M,M,d,d);
		for mu=1:M
			for nu=1:M
				for j=1:d
					for k=1:d
						for L=1:(d^(N-1))
							G[mu,nu,j,k]+= V[mu,(j-1)*d^(N-1)+L]*conj(V[nu,(k-1)*d^(N-1)+L]);
						end;
					end;
				end;
			end;
		end;
	else #Generic method
		#This'd be a very time consuming version of the code (can be better), but lets go for it now
		#
		#The idea is to symmetrising the computational basis
		#
		#
		#
		L::Int64=d^N;
		mu::Int64=0;
		M::Int64=binomial(N+d-1,d-1);
		V=zeros(Float64,M,L); #Should change to the sparse format, lazy for now!!!!
		allStatesLinearlyIndexed=zeros(Int8,L);
		allNormalisations=zeros(Float64,M);
		for k=1:L
			if (allStatesLinearlyIndexed[k]==0) #Check if the state is not yet symmetrised
				mu+=1;
				sss=DecomposedIndex(k,d,N);#Corresponding kets
				SSS=unique(permutations(sss)); #All of its permutations
				for ttt in SSS
					l=ToAbsoluteIndex(ttt,d);
					V[mu,l]=1.;
					allStatesLinearlyIndexed[l]=mu;
				end;
				allNormalisations[mu]=length(SSS);
				V[mu,:]/=sqrt(allNormalisations[mu]);
			end;
		end;
		if (M!=mu) print("***Wrong! Try again!!!\n") end;
		G=zeros(Complex{Float64},M,M,d,d);
		for mu=1:M
			for nu=1:M
				for j=1:d
					for k=1:d
						for K=1:(d^(N-1))
							L1=(j-1)*d^(N-1)+K;
							L2=(k-1)*d^(N-1)+K;
							G[mu,nu,j,k]+= V[mu,L1]*conj(V[nu,L2]);
						end;
					end;
				end;
			end;
		end;
	end;
	return G;
end;


function RobustnessToSymmetricallyExtendible(rho::MultiState,G::Array{Complex{Float64},4},solver=SCS.Optimizer(LOG=0))
	t=Variable(1,Positive());
	M=size(G)[1];
	#X=HermitianSemidefinite(rho.dims[1]*M,rho.dims[2]*M);
	#X=HermitianSemidefinite(rho.dims[2]*M,rho.dims[2]*M);
	X=HermitianSemidefinite(rho.dims[1]*M,rho.dims[1]*M);
	#Construct the constraints corresponding to the equality for all matrix elements of rho	
	problem=maximize(t);
	L=rho.dims[1]*rho.dims[2];
	problem.constraints+= partialtranspose(X,1,[rho.dims[1],M]) in :SDP;
	for L1=1:L
		for L2=L1:L
			ki=LocalIndex(L1,rho.dims); #Which is [k i] in the note 
			lj=LocalIndex(L2,rho.dims); #Which is [l j] in the note
			s=t*rho.mat[L1,L2]; #The state with visibility t
			if ((ki[1]==lj[1]) && (ki[2]==lj[2]))
				s+=(1-t)/L; #The maximally mixed state with visibility (1-t)
			end;
			for mu=1:M
				for nu=1:M
					s-=X[LinearIndex([ki[1],mu],[rho.dims[1],M]),LinearIndex([lj[1],nu],[rho.dims[1],M])]*G[mu,nu,ki[2],lj[2]];
				end;
			end;	
			problem.constraints+= s==0;
		end;
	end;
	solve!(problem,solver);
	return problem.optval;
end;
