#This work for all subdiectories of quantum-correlations
cpath=pwd();
cpath=split(cpath,"quantum-correlations")[1];

push!(LOAD_PATH,cpath*"quantum-correlations/lib/kvant/julia/");
push!(LOAD_PATH,cpath*"quantum-correlations/src/entanglement/");
push!(LOAD_PATH,cpath*"quantum-correlations/src/steerability/");

