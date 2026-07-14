using JLD
using LinearAlgebra
using DelimitedFiles
include("../../bipartite-states/BipartiteStates.jl");
function FormatCovering(FileName,uid)
	mat=readdlm(FileName);
	L=length(mat[:,1]);
	allProjections=Array{Array{Complex{Float64},2},1}(UndefInitializer(),L);
	for k=1:L
		v=[1.0;mat[k,:]];
		allProjections[k]=StateFromBlochVector2d(v);
	end;
	save("../../../data/ipolytope_cvr_d2_L"*string(L)*"_u"*uid*".jld","allProjections",allProjections);
end;

IDs=["092","162","252","362","492","642","812","1002","2252","4002"];

for id in IDs
	FormatCovering("covering_"*string(id)*".txt",id);
end;


