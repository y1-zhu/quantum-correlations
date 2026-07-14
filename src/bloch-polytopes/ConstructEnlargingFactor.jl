include("../bipartite-states/BipartiteStates.jl");

using JLD;

using Polyhedra;
import CDDLib;
lib= CDDLib.Library();

using Convex;
using SCS;

using LinearAlgebra;

allFileNames=readdir("../../data/");
allFileNames=filter(x->occursin("ipolytope",x),allFileNames);
#Reorder the files

alld2FileNames=filter(x->occursin("d2",x),allFileNames);
alld3FileNames=filter(x->occursin("d3",x),allFileNames);

function ReorderByLength(allFileNames)
	nFiles= length(allFileNames);
	LL=Array{Int,1}(UndefInitializer(),nFiles);
	#print(LL,"\n")
	for k in 1:nFiles
		FileName=split(allFileNames[k],"_");
		#print(split(FileName[4],"L")[2],"\n");
		#print(parse(Int,split(FileName[4],"L")[2]),"\n");
		LL[k]=parse(Int,split(FileName[4],"L")[2]);
	end;
	#print(LL,"\n");
	#print(sort(LL),"\n\n");
	return allFileNames[sortperm(LL)];
end;

allFileNames=[ReorderByLength(alld2FileNames);ReorderByLength(alld3FileNames)];

#Starting constructing outer polytopes:

LL=Array{Int,1}();
tt=Array{Float64,1}();
for FileName in allFileNames
	print("*For "*FileName*"...\n")
	allProjections=load("../../data/"*FileName,"allProjections");
	L=length(allProjections);
	d=length(allProjections[1][1,:]);

	allVertices=Array{Float64,2}(UndefInitializer(),L,d^2-1);
	for k in 1:L
		allVertices[k,:]=BlochVector(allProjections[k])[2:(d^2)];
	end;

	print("-Computing the halfspaces...\n");
	V=polyhedron(vrep(allVertices),lib);
	HS=halfspaces(V);

	#Set up the program to compute the enlarging factors:
	print("-Computing the enlarging factors...\n");
	t::Float64=1;
	for H in HS
		#print(H,"\n");
		if (d==2)
			A=Hermitian(MatrixFromBlochVector([1;H.a]));
			#prob.constraints+= (2.0*real(tr(A*E))-1.0 <= t*H.β);
			if (H.β >= 0)
				t=max(t,(2.0*maximum(eigvals(A))-1.0)/H.β);
				#is the enlarging factor correct? <a,x> < b => <a,x> < t b => 2 max_X tr(AX) - 1 < tb
			else
				print("***Error:: H.β is negative!!!");
			end;
		end; #Implement d==3 later
		if (d==3)
			A=Hermitian(MatrixFromBlochVector([sqrt(2.0/3.0);H.a]));
			if (H.β >= 0) t=max(t,(2.0*maximum(eigvals(A))-2.0/3.0)/H.β);
			else print("***Error:: H.β is negative!!!");
			end;
		end;
		#print("\t ",t,"\n");
	end;
	#solve!(prob,SCS.Optimizer());
	if (d==2)
		allOuterProjections= [MatrixFromBlochVector([1.0;t*allVertices[k,:]]) for k in 1:L];
	end;

	if (d==3)
		allOuterProjections= [MatrixFromBlochVector([1.0/sqrt(6);t*allVertices[k,:]]) for k in 1:L];
	end;

	save("../../data/"*replace(FileName,"ipolytope" => "opolytope"),"allProjections",allOuterProjections);
	push!(LL,L);
	push!(tt,t);
end;
