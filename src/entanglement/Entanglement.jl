#Main module which contains all relevant methods
#module Entanglement


#module Entanglement

using MultiStates; 


include("RobustnessToBoundary.jl");
#export RobustnessToBoundary;

include("RobustnessToPPT.jl");
#export RobustnessToPPT;


#include("EntanglementSymmetricExtension.jl");


#include("EntanglementWithSOS.jl");

include("EntanglementWithBlochPolytope.jl");
#export RandomBlochPolytope
#export BlochPolytope
#export RandomBlochPolytope
#export RobustnessToSeparabilityByBlochPolytope
#export DualRobustnessToSeparabilityByBlochPolytope
#export StateRobustnessToSeparabilityByBlochPolytope
#export AbsoluteRobustnessOfEntanglementByBlochPolytope
#export GeneralRobustnessOfEntanglementByBlochPolytope
#export FourQubitRobustnessToSeparabilityByBlochPolytope
#export MRobustnessToGenuineEntanglementByBlochPolytope
#export MRobustnessToFullSeparabilityByBlochPolytope
#export RobustnessToFullBiseparability

include("RobustnessToPPTMixture.jl");
#export PrimalRobustnessToPPTMixture

include("EntanglementWithSymmetricExtension.jl");
#export SymmetricKernel
#export RobustnessToSymmetricallyExtendible

#Wrapper for Entanglement Robustness
function EntanglementRobustness(
			rho::MultiState; #Good if dim[1]>dim[2]
			method::String="PPT", #BlochPolytope or BP, SymmetricExtension or SE, CCNR
			sep_type::String="FSEP",
			solver=SCS.Optimizer, #Generic SDP solver
			qpolytope=nothing,niter= 0, #Parameters for Bloch Polytope
			SymKern=nothing, next= 1, #Paramters for Symmetric Extension (number of extension and symmetric Kernel
			silent=true
			)
	#With PPT

	d = rho.dims[1]

	#check if the state is multipartite

	multiparticle_error_message = "Right now, multiparticle mode only supports 3 qubit FSEP, BSEP, FBSEP, 4 qubits FSEP, 5 qubits FSEP and 3 qutrits FSEP"

	if (length(rho.dims) > 2)
		if ((method=="BlochPolytope")|(method=="BP"))
			if (isnothing(qpolytope)) 
				qpolytope= RandomBlochPolytope(d, nvert) 
			end
			#3 qubits
			if (rho.dims == [2, 2, 2])
				if (sep_type == "FSEP")
					return MRobustnessToFullSeparabilityByBlochPolytope(rho, qpolytope, solver, silent)
				elseif (sep_type == "BSEP")
					return MRobustnessToGenuineEntanglementByBlochPolytope(rho, qpolytope, solver, silent)
				elseif (sep_type == "FBSEP") 
					return RobustnessToFullBiseparability(rho, qpolytope, solver, silent)
				else
					println(multiparticle_error_message)
					return
				end
			elseif (rho.dims == [2, 2, 2, 2]) 
				if (sep_type == "FSEP")
					qpolytope_2 = RandomBlochPolytope(d, nvert)
					return FourQubitRobustnessToSeparabilityByBlochPolytope(rho, qpolytope, qpolytope_2, niter, solver, convergence_recognition, convergence_accuracy, silent)[1]
				else
					println(multiparticle_error_message)
					return
				end
			elseif (rho.dims == [2, 2, 2, 2, 2])
				if (sep_type == "FSEP")
					qpolytope_2 = RandomBlochPolytope(d, nvert)
					qpolytope_3 = RandomBlochPolytope(d, nvert)
					qpolytope_4 = RandomBlochPolytope(d, nvert)
					return FiveQubitsRobustnessToSeparabilityByBlochPolytope(rho, qpolytope, qpolytope_2, qpolytope_3, qpolytope_4, niter, solver, convergence_recognition, convergence_accuracy, silent)[1]
				else
					println(multiparticle_error_message)
					return
				end
			elseif (rho.dims == [3, 3, 3])
				if (sep_type == "FSEP")
					qpolytope_2 = RandomBlochPolytope(d, nvert)
					return ThreeQutritRobustnessToSeparabilityByBlochPolytope(rho, qpolytope, qpolytope_2, niter, solver, convergence_recognition, convergence_accuracy, silent)[1]
				else
					println(multiparticle_error_message)
					return
				end
			else
				println(multiparticle_error_message)
				return
			end
		end
	end

	if (method=="PPT")
		return RobustnessToPPT(rho);
	end;

	# with PPT mixer
	if (method=="PPT_mixer")
		return PrimalRobustnessToPPTMixture(rho);
	end;
	
	#With Bloch polytope
	if ((method=="BlochPolytope")|(method=="BP"))
		if (isnothing(qpolytope)) 
			qpolytope= RandomBlochPolytope(d, nvert) 
		end
		if (robustness_type == "fixed")
			if (isnothing(noise_state))
				return RobustnessToSeparabilityByBlochPolytope(rho,qpolytope,niter,solver, convergence_recognition, convergence_accuracy, silent)[1]
			else 
				return StateRobustnessToSeparabilityByBlochPolytope(rho, qpolytope, niter, noise_state, solver, convergence_recognition, convergence_accuracy, silent)[1] 
			end 
		elseif (robustness_type == "absolute")
			return AbsoluteRobustnessOfEntanglementByBlochPolytope(rho, qpolytope, niter, solver, convergence_recognition, convergence_accuracy, silent)[1]
		elseif (robustness_type == "general")
			return  GeneralRobustnessOfEntanglementByBlochPolytope(rho, qpolytope, niter, solver, convergence_recognition, convergence_accuracy, silent)[1]
		end
	end

	#With Symmetric Extension
	if ((method=="SymmetricExtension")|(method=="SE"))
		rho1= PermuteSystems(rho,1,2); #Because currently Symmetric Extension is done on the second system
		if (isnothing(SymKern)) SymKern= SymmetricKernel(1+next,rho1.dims[2]); end;
		return RobustnessToSymmetricallyExtendible(rho,SymKern,solver);
	end
	
	#With Sum of Square
	if (method=="SOS")
	end;

	#With CCNR
	if (method=="CCNR")
	end;
end

export EntanglementRobustness

#end 


