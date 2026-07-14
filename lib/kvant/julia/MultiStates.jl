#DESCRIPTION---------------------------------------------------------------------

#This module intend to implement common methods to manipulate multipartie states.
#The central type is "BipartiteState", which contains the density operator of a state
#over system A and B, together with the local dimensions dA and dB.
#Single partite states (which are obtained, such as, by tracing out over A or B)
#are simply matrices (no more structures!)

#ISSUES:
#TESTING PUSHING SUBMODULE, ANOTHER TEST
#BEGIN IMPLEMENTATION------------------------------------------------------------
module MultiStates

import Base.+, Base.-, Base.*

using LinearAlgebra
using Combinatorics
using Permutations

i=im;
export i;

#Basic datatype for bipartite state:
struct MultiState
	mat::Matrix{Complex{Float64}}
	dims::Array{Int,1}
end;
export MultiState

#Computed the accumulated size
function AccumSizes(dims::Array{Int,1})
    n= length(dims);
	sizes=Array{Int,1}(undef,n);
	s=1;
    for k in reverse(1:n)
		s*= dims[k];
		sizes[k]=s;
    end;
	return sizes;
end; #Tested

#Accessing the elements of the state
function LinearIndex(kk::Array{Int,1},dd::Array{Int,1})
	sizes= AccumSizes(dd);
	n= length(kk);
	s=kk[n]-1;
	for k in reverse(1:(n-1))
		s+= sizes[k+1]*(kk[k]-1);
	end;
	return s+1;
end;
export LinearIndex;

#There are also more generic function for multiparties cases
function LocalIndex(L::Int,dd::Array{Int,1})
	n=length(dd);
	sizes= AccumSizes(dd);
	kk=Array{Int,1}(undef,n);
	s=L-1;
	for k in reverse(1:n)
		kk[k]= mod(s,dd[k])+1;
		s= div(s,dd[k]);
	end;
	#kA=div(L-1,dd[2])+1;
	#kB=mod(L-1,dd[2])+1;
	return kk;
end;
export LocalIndex;

function ElementByLocalIndex(rho::MultiState,pp::Array{Int,1},qq::Array{Int,1})
	return rho.mat[LinearIndex(pp,rho.dims),LinearIndex(qq,rho.dims)];
end;
export ElementByLocalIndex;
#CREATE STATES-----------------------------------------------------------------
#Taking tensoring of two local states:
function TensorStates(rhoA::MultiState,rhoB::MultiState)
	mat=Matrix{Complex{Float64}}(kron(rhoA.mat,rhoB.mat));
	dims= vcat(rhoA.dims,rhoB.dims); #consider a column array
	return MultiState(mat,dims);
end; #To check a bit more
export TensorStates;

#Pure state tensoring
#function pure_tensor_AB(psiA,psiB)
#end;

#Random state (just uniform, no ensemble), new uniformly distributed version
function RandomMultiState(dd)
	#M=randn(Complex{Float64},prod(dd),prod(dd));
	#M= adjoint(M)*M;
	#M= M/tr(M);
	arr = deepcopy(dd)
	append!(arr,prod(dd))
	pure_state = RandomPureState(arr)
	n = length(arr)
	M = PartialTrace(pure_state,n)
	return M; 
	#return MultiState(M,dd);
end;
export RandomMultiState;

#Random pure state
function RandomPureState(dd)
	M=randn(Complex{Float64},prod(dd))
	M= M*adjoint(M);
	M= M/tr(M);
	return MultiState(M,dd)
end; 
export RandomPureState

#SIMPLE ADDITION AND MULTIPLICATION--------------------------------------------
#Addition: does not work, sadly!!! Still do not understand why
function +(rho1::MultiState,rho2::MultiState);
	return MultiState(Base.:+(rho1.mat,rho2.mat),deepcopy(rho1.dims));
end;
export +;


function -(rho1::MultiState,rho2::MultiState)
	return MultiState(Base.:-(rho1.mat,rho2.mat),deepcopy(rho1.dims));
end;
export -;

#Multiplication by a number: does not work, sadly!!! Still do not understand why
function *(x::Float64,rho::MultiState)
	return MultiState(x.*rho.mat,deepcopy(rho.dims));
end;
export *;

#FUNCTION FOR SUBSYSTEMS--------------------------------------------------------
#Merge systems:
function MergeSystemsLinearly(rho::MultiState,k1::Int,k2::Int)
	dims= deepcopy(rho.dims);
	dmerge= prod(dims[k1:k2]);
	for k in (k1+1):k2 #if k1=k2, nothing changes
		deleteat!(dims,k1+1);
	end;
	dims[k1]=dmerge;
	return MultiState(deepcopy(rho.mat),dims);
end;
export MergeSystemsLinearly;

#Permute system k1 and k2
function PermuteSystems(rho::MultiState,k1::Int,k2::Int)
	dims= deepcopy(rho.dims);
	tmp = dims[k1];
	dims[k1]= dims[k2];
	dims[k2]= tmp;

	n= prod(dims);
	rhoP= Array{Complex{Float64},2}(undef,n,n);
	for l1 in 1:n
		for l2 in 1:n
			kk1= LocalIndex(l1,dims); tmp=kk1[k1]; kk1[k1]=kk1[k2]; kk1[k2]= tmp;
			kk2= LocalIndex(l2,dims); tmp=kk2[k1]; kk2[k1]=kk2[k2]; kk2[k2]= tmp;
			rhoP[l1,l2]=ElementByLocalIndex(rho,kk1,kk2);
		end;
	end;
    return MultiState(rhoP,dims);
end;
export PermuteSystems;


# for dit to decimal conversion
function dit2int(b::Array{Int64, 1}, k::Int64)
    N = length(b)
    d = 0
    for i in 1:N
        d += b[i]*k^(i-1)
    end
    return d
end



#creates unitary representation of symmetric group as list of permutation operators

function UnitaryRepresentation_SymGroup(d::Int64, n::Int64)
	list = Array{Matrix{ComplexF64}, 1}(undef, 0)
	for i in PermGen(n)
		mat = Matrix(i)
		U = zeros(ComplexF64, d^n, d^n)
		for j in 0:(d^n) - 1
			string = digits(j, base=d, pad=n)
			trans_string = mat * string 
			k = dit2int(trans_string, d)
			U[k+1, j+1] = 1.0
		end
		push!(list, U)
	end 
	return list
end
export UnitaryRepresentation_SymGroup

#creates unitary swap operator
function swap(dim::Array{Int64,1},target_dim::Array{Int64,1})
    d = prod(dim)
    dim_new = deepcopy(dim)
    temp = dim_new[target_dim[1]]
    dim_new[target_dim[1]] = dim_new[target_dim[2]]
    dim_new[target_dim[2]] = temp
    C = zeros(Complex{Float64},d,d);
    for l in 1:d
        kk = LocalIndex(l,dim) 
        tmp = kk[target_dim[1]]
        kk[target_dim[1]]=kk[target_dim[2]]
        kk[target_dim[2]]= tmp
        p = LinearIndex(kk, dim_new)
        C[p,l]=1.0
    end
    return C
end 
export swap

#Trace over system k:
function PartialTrace(rho::MultiState,k::Int)
	n= length(rho.dims);
	ddLeft= deepcopy(rho.dims);
	deleteat!(ddLeft,k);
	m= prod(ddLeft);
	rhoLeft= Array{Complex{Float64},2}(undef,m,m);
	for l1 in 1:m
		for l2 in 1:m
			kk1= LocalIndex(l1,ddLeft);
			kk2= LocalIndex(l2,ddLeft);
			if (k==1)
				kkk1=[0;kk1];
				kkk2=[0;kk2];
			elseif (k==n)
				kkk1=[kk1;0];
				kkk2=[kk2;0];
			else
				kkk1= [kk1[1:(k-1)];0;kk1[(k+1):n]];
				kkk2= [kk2[1:(k-1)];0;kk2[(k+1):n]];
			end;
			rhoLeft[l1,l2]=0.;
			for p in 1:rho.dims[k]
				kkk1[k]= p;
				kkk2[k]= p;
				rhoLeft[l1,l2]+=ElementByLocalIndex(rho,kkk1,kkk2);
			end;
		end
	end
	return MultiState(rhoLeft,ddLeft);
end;
export PartialTrace;

#Partial transpose system k
function PartialTranspose(rho::MultiState,k::Int)
	n= prod(rho.dims);
	mat=Array{Complex{Float64},2}(undef,n,n);
	for l1 in 1:n
		for l2 in 1:n
			kk1= LocalIndex(l1,rho.dims);
			kk2= LocalIndex(l2,rho.dims);
			tmp= kk1[k];
			kk1[k]=kk2[k];
			kk2[k]=tmp;
			mat[l1,l2]= ElementByLocalIndex(rho,kk1,kk2);
		end
	end
	return MultiState(mat,deepcopy(rho.dims));
end;
export PartialTranspose;

function ReducedState(rho::MultiState,k::Int)
	n= length(rho.dims);
	if (k==1)
		if (n==1)
			return rho;
		else
			rho1=MergeSystemsLinearly(rho,2,n);
			return PartialTrace(rho1,2);
		end
	elseif (k==n)
		rho1=MergeSystemsLinearly(rho,1,n-1);
		return PartialTrace(rho1,1);
	else
		rho1= MergeSystemsLinearly(rho,(k+1),n);
		rho1= PartialTrace(rho1,k+1);
		rho1= MergeSystemsLinearly(rho1,1,(k-1));
		return PartialTrace(rho1,1);
	end;
end;
export ReducedState;

#CANONICALISATION-----------------------------------------------------------------
function Canonicalise1(rho::MultiState) #Currently only for two parties
	isrhoA=inv(sqrt(ReducedState(rho,1).mat));
	II= zeros(rho.dims[2],rho.dims[2])+I;
	mat= kron(isrhoA,II)*rho.mat*kron(isrhoA,II);
	mat/=tr(mat);
	return MultiState(mat,deepcopy(rho.dims));
end;
export Canonicalise1;

function Canonicalise2(rho::MultiState) #Currently only for two parties
	isrhoB=inv(sqrt(ReducedState(rho,2).mat));
	II= zeros(rho.dims[1],rho.dims[1])+I;
	mat= kron(II,isrhoB)*rho.mat*kron(II,isrhoB);
	mat/=tr(mat);
	return MultiState(mat,deepcopy(rho.dims));
end;
export Canonicalise2;

function Canonicalise(rho::MultiState,k::Int=1)
    if (k==1)
        return Canonicalise1(rho);
    else
        return PermuteSystems(Canonicalise1(PermuteSystems(rho,1,k)),1,k);
    end;
end;
export Canonicalise;

#SYMMETRIC STATES--------------------------------------------------------------
#Create a maximally entangled_state of dimension d:
function MaximallyEntangledBipartiteState(d::Int)
	V=zeros(Complex{Float64},d*d,1);
	for k=1:d
        V[LinearIndex([k,k],[d,d])]=1.;
    end;
	return MultiState((V*V')./d,[d,d]);
end;
export MaximallyEntangledBipartiteState;

#Create an isotropic state:
function IsotropicBipartiteState(d::Int,p::Real)
	V=zeros(Complex{Float64},d*d,1);
	for k=1:d
        V[LinearIndex([k,k],[d,d])]=1.;
    end;
	return MultiState((p/d).*(V*V')+((1.0-p)./(d*d))*I,[d,d]);
end;
export IsotropicBipartiteState;

#Create the fully antisymmetric projection:
function FlipBipartiteOperator(d::Int)::Matrix{Complex{Float64}}
    FF=zeros(Complex{Float64},d*d,d*d);
    for k1=1:d
        for k2=1:d
            FF[LinearIndex([k1,k2],[d,d]),LinearIndex([k2,k1],[d,d])]=1.0;
        end;
    end;
	return FF;
end;
export FlipBipartiteOperator;


function SymmetricBipartiteProjection(d::Int)::Matrix{Complex{Float64}}
	return (I+FlipBipartiteOperator(d))/2.0;
end;
export SymmetricBipartiteProjection;


function AntiSymmetricBipartiteProjection(d::Int)::Matrix{Complex{Float64}}
	return (I-FlipBipartiteOperator(d))/2.0;
end;
export AntiSymmetricBipartiteProjection;

#Create the Werner state:
function WernerBipartiteState(d::Int,p::Real)::MultiState
	return MultiState((2.0*p)/(d*(d-1)).*AntiSymmetricBipartiteProjection(d)+(1.0-p)/(d*d)*I,[d,d]);
end;
export WernerBipartiteState;

#QUBIT-------------------------------------------------------------------------
#Pauli matrices:
Pauli=[[1. 0.;0. 1.],[0. 1.; 1. 0.],[0. (-im); im 0.],[1. 0.; 0. (-1.)]];
export Pauli;

#Compute the Bloch vector:
function BlochVector2d(rho::Matrix{Complex{Float64}})
	return [real(tr(rho*Pauli[k])) for k in 1:4];
end;
export BlochVector2d;

#Compute the state from the Bloch vector:
function MatrixFromBlochVector2d(v::Vector{Float64})
	mat=zeros(Complex,2,2);
	for k=1:4
		mat+= v[k]*Pauli[k];
	end;
	mat*=0.5;
	return mat;
end;
export MatrixFromBlochVector2d;

#Compute the Bloch tensor for qubit:
function BlochTensor2d(rho::MultiState) #Only works for two parties
	ts=zeros(Float64,4,4);
	for kA=1:4
		for kB=1:4
			ts[kA,kB]=real(tr(rho.mat*kron(Pauli[kA],Pauli[kB])));
		end;
	end;
	return ts;
end;
export BlochTensor2d;

#Compute the state from Bloch tensor for qubit:
function StateFromBlochTensor2d(ts::Matrix{Float64})
	mat=zeros(Complex{Float64},4,4);
	for kA=1:4
		for kB=1:4
			mat+= ts[kA,kB]*kron(Pauli[kA],Pauli[kB]);
		end;
	end;
	mat*=0.25;
	return MultiState(mat,[2,2]);
end;
export StateFromBlochTensor2d;

#QUTRIT------------------------------------------------------------------------
#GellMann matrices:
GellMann=Array{Matrix{Complex{Float64}},1}(UndefInitializer(),9);

GellMann[1]= [1. 0. 0.;
              0. 1. 0.;
              0. 0. 1.];
GellMann[1]*=sqrt(2.0/3.0);

GellMann[2]= [0. 1. 0.;
        	  1. 0. 0.;
              0. 0. 0.];
GellMann[3]= [0. (-i) 0.;
             (+i) 0.  0.;
              0.  0.  0.];

GellMann[4]= [1. 0. 0.;
              0. (-1.) 0.;
              0. 0. 0.];

GellMann[5]= [0. 0. 1.;
              0. 0. 0.;
              1. 0. 0.];

GellMann[6]= [0. 0. (-i);
              0. 0. 0.;
              (+i) 0. 0.];

GellMann[7]= [0. 0. 0.;
              0. 0. 1.;
              0. 1. 0.];

GellMann[8]= [0. 0. 0.;
              0. 0. (-i);
              0. (+i) 0.]

GellMann[9]= [1. 0. 0.;
              0. 1. 0.;
              0. 0. (-2.)];
GellMann[9]*=sqrt(1.0/3.0);

export Gellmann;

#Compute the Bloch vector for qutrit:
function BlochVector3d(rho::Matrix{Complex{Float64}})
	return [real(tr(rho*GellMann[k])) for k in 1:9];
end;
export BlochVector3d;

#Compute the state from Bloch vector for qutrit
function MatrixFromBlochVector3d(v::Vector{Float64})
	mat=zeros(Complex{Float64},3,3);
	for k=1:9
		mat+= v[k]*GellMann[k];
	end;
	mat*=0.5; #Check the coefficient
	return mat;
end;
export MatrixFromBlochVector3d;

#Compute the Bloch tensor for qutrit:
function BlochTensor3d(rho::MultiState)
	ts=zeros(Float64,9,9);
	for kA=1:9
		for kB=1:9
			ts[kA,kB]=real(tr(rho.mat*kron(GellMann[kA],GellMann[kB])));
		end;
	end;
	return ts;
end;
export BlochTensor3d;

#Compute the state from Bloch tensor:
function StateFromBlochTensor3d(ts::Matrix{Float64})
	mat=zeros(Complex{Float64},9,9);
	for kA=1:9
		for kB=1:9
			mat+= ts[kA,kB]*kron(GellMann[kA],GellMann[kB]);
		end;
	end;
	mat*=0.25; #Check the coefficient
	return MultiState(mat,[3,3]);
end;
export StateFromBlochTensor3d;
#end;

#---------------------------------------------------------------------------
#Combine generic functions:
function MatrixFromBlochVector(v::Vector{Float64})
	if (length(v)==4)
		return MatrixFromBlochVector2d(v);
	elseif (length(v)==9)
		return MatrixFromBlochVector3d(v);
	else
		print("***Error::StateFromBlochVector: not yet implemented for this length!!!" );
	end;
end;
export MatrixFromBlochVector;

function BlochVector(mat::Matrix{Complex{Float64}})
	if (size(mat)==(2,2))
		return BlochVector2d(mat);
	elseif (size(mat)==(3,3))
		return BlochVector3d(mat);
	else 
		error("***Blochvector(mat::Matrix{Complex{Float64}}) is not implemented for dimension ",size(mat),"\n")
	end;
end;
export BlochVector;
#-------------------------------------------------------------------
#Some special states:
function HorodeckiTwoQutritState(a::Float64)::MultiState
	mat=zeros(Complex{Float64},9,9);
	mat[1,1]=a;
	mat[2,2]=a;
	mat[3,3]=a;
	mat[4,4]=a;
	mat[5,5]=a;
	mat[6,6]=a;
	mat[7,7]=(1.0+a)/2.0;
	mat[8,8]=a;
	mat[9,9]=(1.0+a)/2.0;
	mat[1,5]=a;mat[5,1]=a;
	mat[1,9]=a;mat[9,1]=a;
	mat[5,9]=a;mat[9,5]=a;
	mat[7,9]=sqrt(1.0-a^2)/2.0;mat[9,7]=sqrt(1.0-a^2)/2.0;
	mat/=(8.0*a+1.0);
	return MultiState(mat,[3,3]);
end;
export HorodeckiTwoQutritState;

function HorodeckiTwoQutritState2(q::Float64=0.5)::MultiState
	mat=zeros(Complex{Float64},9,9);
	plus = (5/2) + q;
	minus = (5/2) - q;
	mat[1,1]=2;
	mat[2,2]=minus;
	mat[3,3]=plus;
	mat[4,4]=plus;
	mat[5,5]=2;
	mat[6,6]=minus;
	mat[7,7]=minus;
	mat[8,8]=plus;
	mat[9,9]=2;
	mat[1,5]=2;mat[5,1]=2;
	mat[1,9]=2;mat[9,1]=2;
	mat[5,9]=2;mat[9,5]=2;
	mat = mat/21;
	return MultiState(mat,[3,3]);
end;
export HorodeckiTwoQutritState2;

function GHZThreeQubitState(p::Float64=1.0)::MultiState
    V=zeros(Complex{Float64},8);
    V[1]=1.0/sqrt(2.0);
    V[8]=1.0/sqrt(2.0);
    return MultiState(p*V*conj(transpose(V))+(1-p)/8.0*I,[2;2;2]);
end
export GHZThreeQubitState;

##general GHZ state: local dimension = d , parties = p
function GHZ(d::Int64=2,p::Int64=3)
	G = Array{Int64,1}(undef,0);
	for i in 1:p
		push!(G,d)	
	end;
	V=zeros(Complex{Float64},d^p);
	local_vector = zeros(Complex{Float64},d);
	for i in 1:d
		local_vector[i] = 1.0;
		vec = local_vector;
		for i in 2:p
			local_vector = kron(vec,local_vector);
		end;
		V += (1/sqrt(d))*local_vector;
		local_vector = zeros(Complex{Float64},d);
	end;
	return MultiState(V*conj(transpose(V)),G);
end
export GHZ;

##general W state: local dimension = d , parties = p
function Wstate(d::Int64=2,p::Int64=3)
	G = Array{Int64,1}(undef,0);
	for i in 1:p
		push!(G,d)	
	end;
	factor = (1/sqrt(p*(d-1.0)));
	V=zeros(Complex{Float64},d^p);
	for i in 2:d
		zero_in_vec = zeros(Complex{Float64},d);
		i_in_vec = zeros(Complex{Float64},d);
		zero_in_vec[1] = 1.0;
		i_in_vec[i] = 1.0;
		vec = zeros(Complex{Float64},d^p);
		for j in 1:p
			initial_vec = zero_in_vec
			if (j==1)
				initial_vec = i_in_vec;
			end; 
			for k in 2:p
				factor_vec = zero_in_vec;
				if (k == j)
					factor_vec = i_in_vec;
				end;
				initial_vec = kron(initial_vec,factor_vec);
				
			end;
			vec += initial_vec;	
		end;		
		V += vec;
	end;
	V = factor*V
	return MultiState(V*conj(transpose(V)),G);
end
export Wstate;

theta = acos(sqrt(1/2 + sqrt(1/12)));

function FBSEP_robust_state(X::Float64=theta)
	Y = zeros(Complex{Float64},4,8);
	Y[1,:] = (1/sqrt(2))*[cos(X) 0 0 -cos(X) sin(X) 0 0 -sin(X)];
	Y[2,:] = (1/sqrt(2))*[0 im*sin(X) -im*sin(X) 0 0 -cos(X) cos(X) 0];
	Y[3,:] = (1/sqrt(2))*[-cos(X) 0 0 -cos(X) sin(X) 0 0 sin(X)];
	Y[4,:] = (1/sqrt(2))*[0 -im*sin(X) -im*sin(X) 0 0 -cos(X) -cos(X) 0];
	Z = zeros(Complex{Float64},8,8);
	for i in 1:4
		Z += 0.25*(Y[i,:])*adjoint((Y[i,:]));
	end;
	Z = MultiState(Z,[2,2,2]); 
	return Z; 
end;
export FBSEP_robust_state


function WThreeQubitState(p::Float64=1.0)::MultiState
    V=zeros(Complex{Float64},8);
    V[2]=1.0/sqrt(3.0);
    V[3]=1.0/sqrt(3.0);
    V[5]=1.0/sqrt(3.0);
    return MultiState(p*V*conj(transpose(V))+(1-p)/8.0*I,[2;2;2]);
end
export WThreeQubitState;


function WFourQubitState(p::Float64=1.0)::MultiState
    V=zeros(Complex{Float64},16);
    V[2]=1.0/sqrt(4.0);
    V[3]=1.0/sqrt(4.0);
    V[5]=1.0/sqrt(4.0);
    V[9]=1.0/sqrt(4.0);	
    return MultiState(p*V*conj(transpose(V))+(1-p)/16.0*I,[2;2;2;2]);
end
export WFourQubitState;

function ClusterState()::MultiState
	V=zeros(Complex{Float64},16);
	V[1] = 1.0;
	V[16] = 1.0;
	V[4] = 1.0;
	V[13] = 1.0;
	V *= 1/4;
	return MultiState(V*conj(transpose(V)),[2,2,2,2]); 
end; 
export ClusterState;

function ClusterState_2()
	## qubit state 0,1,+,-
	zero = zeros(Complex{Float64},2);
	zero[1] = 1.0;
	one = zeros(Complex{Float64},2);
	one[2] = 1.0;
	plus = zeros(Complex{Float64},2);
	plus[1] = 1.0;
	plus[2] = 1.0;
	plus *= (1/sqrt(2));
	minus = zeros(Complex{Float64},2);
	minus[1] = 1.0;
	minus[2] = -1.0;
	minus *= (1/sqrt(2));
	vec = kron(kron(plus,zero),kron(plus,zero));
	vec += kron(kron(plus,zero),kron(minus,one));
	vec += kron(kron(minus,one),kron(minus,zero));
	vec += kron(kron(minus,one),kron(plus,one));
	vec /= 2;
	#print(conj(transpose(vec))*vec, "\n");
	return MultiState(vec*conj(transpose(vec)),[2,2,2,2]);
end;
export ClusterState_2

function Ket_to_Vector(ket::Array{Int64,1}=[1,0]) ## gives normed vector in computational basis, given a binary ket(could be extended to qudits)
	N = length(ket)
	n = 2^N
	V = zeros(Complex{Float64},n);
	s = "";
	for i in ket
		s = string(s,string(i));
	end;
	a = parse(Int64, string("0b",s))+1
	V[a] = 1.0;
	return V
end; 
export Ket_to_Vector;


function DickeState(d::Int64,ex::Int64) ##only qubit systems for now 
	dims = zeros(Int64,d).+2;
	N = 2^d;
	V=zeros(Complex{Float64},N);
	left_ex = zeros(Int64,d);
	for i in 1:ex
		left_ex[i] = 1;
	end;
	x = collect(permutations(left_ex));
	perm_num = length(x); 
	for i in x
		V += Ket_to_Vector(i);
	end;
	norm = sqrt(conj(transpose(V))*V);
	V /= norm
	M = MultiState(V*conj(transpose(V)),dims)
	return M; 
end; 
export DickeState;

function HyllusThreeQubitState(eta::Float64=sqrt(3.0/2.0))::MultiState
    mat=zeros(Complex{Float64},8,8);
    mat[1,1]=2*eta;
    mat[2,2]=1;
    mat[3,3]=1;
    mat[4,4]=1.0/eta;
    mat[5,5]=1;
    mat[6,6]=1.0/eta;
    mat[7,7]=1.0/eta;
    mat[8,8]=0.0;

    mat[2,3]=1; mat[3,2]=1;
    mat[2,5]=1; mat[5,2]=1;
    mat[3,5]=1; mat[5,3]=1;

    mat/= (3+2*eta+3.0/eta);
    return MultiState(mat,[2,2,2]);
end
export HyllusThreeQubitState;

function ABLSThreeQubitState(a::Float64=0.5,b::Float64=0.5,c::Float64=0.5);
    mat=zeros(Complex{Float64},8,8);
    mat[1,1]=1;
    mat[2,2]=a;
    mat[3,3]=b;
    mat[4,4]=c;
    mat[5,5]=1.0/c;
    mat[6,6]=1.0/b;
    mat[7,7]=1.0/a;
    mat[8,8]=1;
    mat[1,8]=1;
    mat[8,1]=1;
    mat/= (2+a+1.0/a+b+1/b+c+1/c);
    return MultiState(mat,[2,2,2]);
end;
export ABLSThreeQubitState;

function HorodeckiThreeQubitState(b::Float64=0.5);
    mat=zeros(Complex{Float64},8,8);
    mat[1,1]=b;
    mat[2,2]=b;
    mat[3,3]=b;
    mat[4,4]=b;
    mat[5,5]=(1.0+b)/2.0;
    mat[6,6]=b;
    mat[7,7]=b;
    mat[8,8]=(1.0+b)/2.0;

    mat[1,6]=b; mat[6,1]=b;
    mat[2,7]=b; mat[7,2]=b;
    mat[3,8]=b; mat[8,3]=b;

    mat[5,8]=sqrt(1.0-b^2)/2.0;
    mat[8,5]=sqrt(1.0-b^2)/2.0

    mat/= (1.0+7.0*b);
    return MultiState(mat,[2,2,2]);
end;
export HorodeckiThreeQubitState;

function HorodeckiQuquardQubitState(b::Float64=0.5);
    mat=zeros(Complex{Float64},8,8);
    mat[1,1]=b;
    mat[2,2]=b;
    mat[3,3]=b;
    mat[4,4]=b;
    mat[5,5]=(1.0+b)/2.0;
    mat[6,6]=b;
    mat[7,7]=b;
    mat[8,8]=(1.0+b)/2.0;

    mat[1,6]=b; mat[6,1]=b;
    mat[2,7]=b; mat[7,2]=b;
    mat[3,8]=b; mat[8,3]=b;

    mat[5,8]=sqrt(1.0-b^2)/2.0;
    mat[8,5]=sqrt(1.0-b^2)/2.0

    mat/= (1.0+7.0*b);
    return MultiState(mat,[2,4]);
end;
export HorodeckiQuquardQubitState;

function GHZState(d::Int=4); #only for d=4
    mat=zeros(Complex{Float64},d^2,d^2);
    mat[1,1]=1.0/2.0;
    mat[d^2,d^2]=1.0/2.0;
    mat[d^2,1]=1.0/2.0;
    mat[1,d^2]=1.0/2.0;
    return MultiState(mat,[2,2,2,2]);
end; 
export GHZState;

function IdTimesGHZ();
    id= zeros(Complex{Float64},2,2);
    id= id + (1/2)*I; 
    V=zeros(Complex{Float64},8);
    V[1]=1.0/sqrt(2.0);
    V[8]=1.0/sqrt(2.0);
    V=V*conj(transpose(V));
    mat= kron(id, V*conj(transpose(V))); 
    return MultiState(mat,[2,2,2,2]);
end; 
export IdTimesGHZ;

function BennettState();
    v = Array{Any,1}(undef,4);
    num= 1/sqrt(2);
    v[1]=kron(kron([1,0],[0,1]),[num,num]);
    v[2]=kron(kron([0,1],[num,num]),[1,0]);
    v[3]=kron(kron([num,num],[1,0]),[0,1]);
    v[4]=kron(kron([num,-num],[num,-num]),[num,-num]);
    #return MultiState(mat,[2,2,2])
    mat= (1/4)*(zeros(Complex{Float64},8,8)+I)
    for i in 1:4
        mat= mat - (1/4)*v[i]*conj(transpose(v[i]));  
    end;
    return MultiState(mat,[2,2,2]); 
    #return tr(mat);
    #return MultiState(mat,[2,2,2]); 
end;
export BennettState;

function LeonardoState()
	d = 8
	rho = zeros(Complex{Float64},8,8) + I 
	pauli_product = kron(kron(Pauli[4], Pauli[2]),Pauli[4])
	pauli_product *= 1/sqrt(2)
	rho = 1/8 * (rho - pauli_product)
	return MultiState(rho, [2,2,2])
end
export LeonardoState


##Bell states 
function psi_plus()
	mat= zeros(Complex{Float64},4,4); 
	mat[2,2] = 1.0;
 	mat[2,3] = 1.0;
	mat[3,2] = 1.0;
	mat[3,3] = 1.0;
	mat = (1/2)*mat;
	mat = MultiState(mat,[2,2])
	return mat;
end;

function psi_minus()
	mat= zeros(Complex{Float64},4,4); 
	mat[2,2] = 1.0;
 	mat[2,3] = -1.0;
	mat[3,2] = -1.0;
	mat[3,3] = 1.0;
	mat = (1/2)*mat;
	mat = MultiState(mat,[2,2])
	return mat;
end;

function phi_plus()
	mat= zeros(Complex{Float64},4,4); 
	mat[1,1] = 1.0;
 	mat[1,4] = 1.0;
	mat[4,1] = 1.0;
	mat[4,4] = 1.0;
	mat = (1/2)*mat;
	mat = MultiState(mat,[2,2])
	return mat;
end;

function phi_minus()
	mat= zeros(Complex{Float64},4,4); 
	mat[1,1] = 1.0;
 	mat[1,4] = -1.0;
	mat[4,1] = -1.0;
	mat[4,4] = 1.0;
	mat = (1/2)*mat;
	mat = MultiState(mat,[2,2])
	return mat;
end;

## state of the paper by Moroder and Mittsovich

function rho_BFP()
	mat= zeros(Complex{Float64},16,16);
	mat += kron(phi_plus().mat,psi_minus().mat);	
	mat += kron(psi_plus().mat,psi_plus().mat);
	mat += kron(psi_minus().mat,phi_minus().mat);
	mat += kron(phi_minus().mat,psi_plus().mat);
	mat += kron(phi_minus().mat,psi_minus().mat);
	mat += kron(phi_minus().mat,phi_minus().mat);
	mat = (1/6)*mat;
	mat = MultiState(mat,[2,2,2,2]);
	mat = PermuteSystems(mat,2,3); 
	mat = MergeSystemsLinearly(mat,3,4);
	mat = MergeSystemsLinearly(mat,1,2);
	return mat;
end; 

function Max_entangled_family(d::Int64,beta::Float64) 
    p_plus = MaximallyEntangledBipartiteState(d).mat
    one = zeros(ComplexF64, d^2, d^2) + I 
    Z = zeros(ComplexF64, d^2, d^2)
    for i in 1:d
        vec = zeros(ComplexF64,d)
        vec[i] = 1.0
        tensor_vec = kron(vec, vec)
        Z += tensor_vec*adjoint(tensor_vec)
    end
    rho_g = 1/(d^2 - d)*(one-Z)
    mat = (1/1+beta)*(beta*rho_g + p_plus)
    return MultiState(mat,[d, d])
end 
export Max_entangled_family

function IsotropicState(d::Int64, beta::Float64)
    p_plus = MaximallyEntangledBipartiteState(d).mat
    one = zeros(ComplexF64, d^2, d^2) + I
    mat = one + beta*p_plus
    mat /= tr(mat)
    return MultiState(mat,[d, d])
end 
export IsotropicState

## magic simplex states

W_alpha = exp(im *(2*pi/3));

W_X = [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]; 

W_Z = [1.0 0.0 0.0; 0.0 W_alpha 0.0; 0.0 0.0 W_alpha];

psi_00 = (1/sqrt(3))*adjoint([1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0])

function bell_3x3(x::Int64,y::Int64)
	II = zeros(ComplexF64, 3, 3) + I
	mat = (W_X^x) * (W_Z^y)
	vec = kron(II, mat)*psi_00
	return vec 
end
export bell_3x3

## families of two qutrit bell diagonal states 
function MSS_A(a::Float64=0.2,b::Float64=-0.005,c::Float64=0.0)
	II = 1/9 * (zeros(ComplexF64, 9, 9) + I)
	rho = zeros(ComplexF64, 9, 9)
	rho += (1-a-b-c) * II
	Psi_00 = bell_3x3(0,0)
	rho += a * Psi_00 * adjoint(Psi_00)
	Psi_01 = bell_3x3(0,1)
	rho += b * Psi_01 * adjoint(Psi_01)
	Psi_02 = bell_3x3(0,2)
	rho += c * Psi_02 * adjoint(Psi_02)
	return MultiState(rho, [3, 3])
end 
export MSS_A



#----------------------------------------------------------------------------------
#NEWLY ADDED
#Embeding Opetators
function EmbedOperator(X::Array{Complex{Float64},2},k::Int64,dims::Array{Int64,1})
	Y=1;
	for l in 1:length(dims)
		if (l==k)
			Y=kron(Y,X);
		#	print(l,"\n");
		else
			Y=kron(Y,Matrix{Complex{Float64}}(I,dims[l],dims[l]));
		#	print(-l,"\n");
		end
		#print(size(Y),"\n")
	end;
	return Y;
end
export EmbedOperator;

function ConditionalState(rho::MultiState,E::Array{Complex{Float64},2},k::Int64)
	rho1= MultiState(rho.mat*EmbedOperator(E,k,rho.dims),rho.dims);
	rho1= PartialTrace(rho1,k);
	return MultiState(rho1.mat/tr(rho1.mat),rho1.dims);
end
export ConditionalState;

function UnnormalisedConditionalState(rho::MultiState,E::Array{Complex{Float64},2},k::Int64)
	rho1= MultiState(rho.mat*EmbedOperator(E,k,rho.dims),rho.dims);
	rho1= PartialTrace(rho1,k);
	return MultiState(rho1.mat,rho1.dims);
end
export UnnormalisedConditionalState;

end
