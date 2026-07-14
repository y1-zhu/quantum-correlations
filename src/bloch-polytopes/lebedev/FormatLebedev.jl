using JLD
using LinearAlgebra
using DelimitedFiles
include("../../bipartite-states/BipartiteStates.jl");
function FormatLebedev(FileName,uid)
	mat=readdlm(FileName);
	L=length(mat[:,1]);
	allProjections=Array{Array{Complex{Float64},2},1}(UndefInitializer(),L);
	for k=1:L
			v=[1.0,cos(mat[k,1])*sin(mat[k,2]),sin(mat[k,1])*sin(mat[k,2]),cos(mat[k,2])];
		allProjections[k]=StateFromBlochVector2d(v);
	end;
	save("polytope_lbd_d2_L"*string(L)*"_u"*uid*".jld","allProjections",allProjections);
end;

IDs=["017","041","047","053","059","065","071","077","083","089","095","101","107","113","119","125","131"];

for id in IDs
	FormatLebedev("lebedev_"*string(id)*".txt",id);
end;


