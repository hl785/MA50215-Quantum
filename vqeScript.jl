using PosDefManifold
using Plots
using SciPy
using LinearAlgebra
include("main.jl")

numQBits = 3
initReg = rBs0()^numQBits
# x0 = [1.0, 0.0]
# x0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
x0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# mtd = "Nelder-Mead"
mtd = "COBYLA"

dimension = 2^numQBits
H = randP(ComplexF64, dimension)
hOpr = Opr(H)

function uParam1QB(a::Array{Float64,1})::Opr
    return matmul(rotGate(paZ2(), a[2]), rotGate(paX2(), a[1]))
end

function recurse(q::Int64, p1::Float64, p2::Float64, pc1::Array{Float64,1}, pc2::Array{Float64,1}, prevUParam)
    gate1::Opr = (rotGate(paZ2(), p1) * (eye2()^(q - 1))) * rotGate(paX2(), p2)
    subgate1::Opr = Opr(0.5 * (eye2().vals + paZ2().vals))
    subgate1 = (eye2()^q) * subgate1
    subgate2::Opr = Opr(0.5 * (eye2().vals - paZ2().vals))
    subgate3::Opr = prevUParam(pc1) * subgate2
    gate2::Opr = Opr(subgate1.vals + subgate3.vals)
    gate3::Opr = (eye2()^q) * paX2()
    subgate3 = prevUParam(pc2) * subgate2
    gate4::Opr = Opr(subgate1.vals + subgate3.vals)
    return matmul(matmul(gate4, gate3), matmul(gate2, gate1))
end

function uParam2QB(a::Array{Float64,1})::Opr
    return recurse(1, a[1], a[2], a[3:4], a[5:6], uParam1QB)
end

function uParam3QB(a::Array{Float64,1})::Opr
    return recurse(2, a[1], a[2], a[3:8], a[9:14], uParam2QB)
end

function uParam(a::Array{Float64,1})::Opr
    if numQBits == 1
        @assert length(x0) == 2 "Wrong number of initial parameters"
        return uParam1QB(a)
    elseif numQBits == 2
        @assert length(x0) == 6 "Wrong number of initial parameters"
        return uParam2QB(a)
    elseif numQBits == 3
        @assert length(x0) == 14 "Wrong number of initial parameters"
        return uParam3QB(a)
    else
        @assert false "Not coded this number of QBits"
    end
    return scalarArray(ComplexF64(1, 0))
end

function cost(a::Array{Float64,1}, print::Bool)::Float64
    register = initReg
    uOpr = uParam(a)
    register = uOpr * register
    eVal = oprExpectation(hOpr, register)
    if print
        println("Estimated eigenval:  ", eVal.re)
        println("Estimated eigenvec:  ", register.vals[1], ", ", register.vals[2])  # TODO: longer evals print
        println("Eigenpair residual:  ", norm((hOpr * register).vals - (eVal.re * register.vals)))
    end
    return eVal.re
end

function costNoOpt(a::Array{Float64,1})::Float64
    return cost(a, false)
end

res = SciPy.optimize.minimize(costNoOpt, x0, method = mtd, tol = 1e-24)

# Check
eValRe = cost(res["x"], true)
trueEVal = eigvals(H, 1:1)[1]
println("Calculated eigenval: ", trueEVal)
println("Eigenval error:      ", abs(trueEVal - eValRe))

println(res)