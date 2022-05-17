using PosDefManifold
using Plots
using SciPy
using LinearAlgebra
include("main.jl")

numQBits = 1
dimension = 2^numQBits
H=randP(ComplexF64, dimension)
hOpr = Opr(H)
# heatmap(broadcast(abs, hOpr.vals), yflip = true)

register = rBs0() ^ numQBits
# heatmap(broadcast(abs, kto(register).vals), yflip = true)

# U = randU(ComplexF64, dimension)
# uOpr = Opr(U)
# heatmap(broadcast(abs, uOpr.vals), yflip = true)

function uParam(a::Float64, b::Float64)::Opr
    return out = matmul(rotGate(paZ2(), b), rotGate(paX2(), a))
end

function cost(a::Float64, b::Float64, print::Bool)::Float64
    regCopy = register
    uOpr = uParam(a,b)
    regCopy = uOpr * regCopy
    eVal = oprExpectation(hOpr, regCopy)
    if print
        # println("Estimated eigenval:  ", eVal)
        println("Estimated eigenval:  ", eVal.re)
        println("Estimated eigenvec:  ", regCopy.vals[1], ", ", regCopy.vals[2])
        println("Eigenpair residual:  ", norm((hOpr * regCopy).vals - (eVal.re * regCopy.vals)))
    end
    return eVal.re
end 

function costVecIn(a::Array{Float64, 1})::Float64
    return cost(a[1], a[2], false)
end

x0 = [1.0, 0.0]

res = SciPy.optimize.minimize(costVecIn, x0, method="Nelder-Mead", tol=1e-6)

# Check
eValRe = cost(res["x"][1], res["x"][2], true)
trueEVal = eigvals(H, 1:1)[1]
println("Calculated eigenval: ", trueEVal)
println("Eigenval error:      ", abs(trueEVal - eValRe))