using LinearAlgebra
using SCS

using Convex
using IterTools


"""takes the partialtranspose of subsets labelled in tp_sys from a matrix
with local dimensions sys"""
function multi_partial_transpose(X, tp_sys::Array, sys::Array)
    Y = deepcopy(X)

    for i in tp_sys
        Y = partialtranspose(Y, i, sys)
    end 

    return Y 
end 

function DualRobustnessToPPTMixture(rho::MultiState)
    N = length(rho.dims)

    d = prod(rho.dims)

    n = Int64(floor(N/2))

    set = [collect(subsets(1:N, m)) for m in 1:n]

    X= [[HermitianSemidefinite(d, d) for i in 1:length(set[m])] for m in 1:n]

    Y= [[HermitianSemidefinite(d, d) for i in 1:length(set[m])] for m in 1:n]

    problem= minimize(real(tr(rho.mat * (X[1][1] + partialtranspose(Y[1][1], 1, rho.dims)))))

    problem.constraints += real(tr(X[1][1] + partialtranspose(Y[1][1], 1, rho.dims))) == prod(rho.dims)

    for m in 1:n 
        for i in 1:length(set[m])
            s = set[m][i]
            problem.constraints += X[1][1] + partialtranspose(Y[1][1],1,rho.dims)  == X[m][i] + multi_partial_transpose(Y[m][i], s, rho.dims)
        end 
    end 

    solve!(problem, SCS.Optimizer, silent_solver=true)
    
    print(problem.status,"\n")
    
    x = problem.optval
    
    t = 1/(1-x)
    
    V = evaluate(X[1][1])
    
    Z = partialtranspose(evaluate(Y[1][1]), 1, rho.dims)

    return t, V+Z
end 

function PrimalRobustnessToPPTMixture(rho::MultiState)
   
    N = length(rho.dims)

    d = prod(rho.dims)

    n = Int64(floor(N/2))

    set = [collect(subsets(1:N, m)) for m in 1:n]

    X= [[HermitianSemidefinite(d, d) for i in 1:length(set[m])] for m in 1:n]
    
    S=zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))
    
    for m in 1:n
        for i in 1:length(set[m])
             S += X[m][i]
        end 
    end

    II= zeros(Complex{Float64},prod(rho.dims),prod(rho.dims))+I

    t=Variable(1,Positive())

    problem= maximize(t)

    problem.constraints+= ((t*rho.mat+(1.0-t)/prod(rho.dims)*II) == S)

    for m in 1:n
        for i in 1:length(set[m])
             s = set[m][i]
             problem.constraints += multi_partial_transpose(X[m][i], s, rho.dims) in :SDP
        end 
    end

    solve!(problem, SCS.Optimizer, silent_solver=true)

    print(problem.status,"\n")

    return problem.optval
end