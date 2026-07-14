include("startup.jl")
using MultiStates

include("Entanglement.jl")

rho = RandomMultiState([2,2,2]) # random three qubit states, replace this by your state here... 

bp = RandomBlochPolytope(2,400) #random qubit polytope with 400 vertices, the more vertices the better...

val = EntanglementRobustness(rho, method="BP", sep_type="BSEP", qpolytope=bp)

println("The visibility is: ", val,". If this is >= 1, we know that the state must be biseparable.")
