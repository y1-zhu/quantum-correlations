using LinearAlgebra
using JLD
function ReadComplexMatrix!(fin,d1,d2)
	mat=zeros(Complex{Float64},d1,d2);
	for k1=1:d1
		sline=[];
		line="\n";
		while (length(sline) == 0)
			line=readline(fin);
			#print("read the line:/n",line);
			sline=split(line);
		end;
		if (d2 != length(sline)) print("**Warning: length of line not match!", d2,"\t",length(sline),"\n"); end;
		for k2=1:d2
			mat[k1,k2]=parse(Complex{Float64},sline[k2]);
		end;
	end;
	return mat;
end;

function OrthogonaliseSTPolytopes(FileName)
	fin=open(FileName*".txt","r");
	line=readline(fin);
	print(line);
	line=readline(fin);
	print(line);
	d=parse(Int,split(readline(fin))[1]);
	L=parse(Int,split(readline(fin))[1]);
	print("d= ",d,"\t L= ",L);

	H=ReadComplexMatrix!(fin,d,d);

	allProjections=Array{Array{Complex{Float64}, 2},1}(UndefInitializer(),L);
	for k=1:L
		allProjections[k]=ReadComplexMatrix!(fin,d,d);
	end;

	#using LinearAlgebra;
	V=eigvecs(Hermitian(H));
	DH=V'*H*V;
	S=sqrt(DH)*V'; #H=S'*S
	iS=inv(S);
	for k=1:length(allProjections)
			allProjections[k]=S*allProjections[k]*iS;
	end;
	s=tr(allProjections[1]);
	for k=1:length(allProjections)
		allProjections[k]/=s;
	end;
	#using JLD;
	save(FileName*".jld","allProjections",allProjections);
end;


#FileName=["polytope_st27_d3_L360_u1","polytope_st27_d3_L360_u2","polytope_st27_d3_L360_u3","polytope_st27_d3_L360_u4","polytope_st27_d3_L360_u5","polytope_st27_d3_L360_u5"];
FileNames=["polytope_st27_d3_L180_u2","polytope_st27_d3_L36_u6","polytope_st27_d3_L60_u3","polytope_st27_d3_L360_u1","polytope_st27_d3_L45_u5","polytope_st27_d3_L60_u4"];

for FileName in FileNames
	OrthogonaliseSTPolytopes(FileName)
end;



