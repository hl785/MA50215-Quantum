using PosDefManifold
using Plots
using SciPy
using LinearAlgebra
include("main.jl")

numQBits = 1
initReg = rBs0() ^ numQBits
x0 = [1.0, 0.0]
mtd = "Nelder-Mead"

dimension = 2^numQBits
H=randP(ComplexF64, dimension)
hOpr = Opr(H)

function uParam(a::Array{Float64, 1})::Opr
    if numQBits == 1
        @assert length(x0) == 2 "Wrong number of initial parameters"
        return matmul(rotGate(paZ2(), a[2]), rotGate(paX2(), a[1]))
    else
        @assert false "Not coded this number of QBits"
    end
    return scalarArray(ComplexF64(1,0))
end

function cost(a::Array{Float64, 1}, print::Bool)::Float64
    register = initReg
    uOpr = uParam(a)
    register = uOpr * register
    eVal = oprExpectation(hOpr, register)
    if print
        # println("Estimated eigenval:  ", eVal)
        println("Estimated eigenval:  ", eVal.re)
        println("Estimated eigenvec:  ", register.vals[1], ", ", register.vals[2])  # TODO: longer evals print
        println("Eigenpair residual:  ", norm((hOpr * register).vals - (eVal.re * register.vals)))
    end
    return eVal.re
end 

function costNoOpt(a::Array{Float64, 1})::Float64
    return cost(a, false)
end

res = SciPy.optimize.minimize(costNoOpt, x0, method = mtd, tol=1e-6)

# Check
eValRe = cost(res["x"], true)
trueEVal = eigvals(H, 1:1)[1]
println("Calculated eigenval: ", trueEVal)
println("Eigenval error:      ", abs(trueEVal - eValRe))