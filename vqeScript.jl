using PosDefManifold
using SciPy
using LinearAlgebra
include("main.jl")
include("temp.jl")

# -------------------------------
#       (Quantum) Functions      
# -------------------------------

# Base case: minimal (size), maximally expressive (up to global phase) circuit for 1 Qbit.
function uParam1QB(a::Array{Float64,1})::Opr
    return matmul(rotGate(paZ2(), a[2]), rotGate(paX2(), a[1]))
end

# Minimal, maximally expressive circuit for N+1 Qbits can be built up from two copies of the N Qbits circuit.
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

# Basic recursive glue code.
function uParam2QB(a::Array{Float64,1})::Opr
    return recurse(1, a[1], a[2], a[3:4], a[5:6], uParam1QB)
end

# Basic recursive glue code.
function uParam3QB(a::Array{Float64,1})::Opr
    return recurse(2, a[1], a[2], a[3:8], a[9:14], uParam2QB)
end

# Basic recursive glue code.
function uParam4QB(a::Array{Float64,1})::Opr
    return recurse(3, a[1], a[2], a[3:16], a[17:30], uParam3QB)
end

# Basic recursive glue code.
function uParam5QB(a::Array{Float64,1})::Opr
    return recurse(4, a[1], a[2], a[3:32], a[33:62], uParam4QB)
end

# Basic recursive glue code.
function uParam6QB(a::Array{Float64,1})::Opr
    return recurse(5, a[1], a[2], a[3:64], a[65:126], uParam5QB)
end

# Basic recursive glue code.
function uParam7QB(a::Array{Float64,1})::Opr
    return recurse(6, a[1], a[2], a[3:128], a[129:254], uParam6QB)
end

# Basic recursive glue code.
function uParam8QB(a::Array{Float64,1})::Opr
    return recurse(7, a[1], a[2], a[3:256], a[257:510], uParam7QB)
end

# Call different functions depending on variable.
function uParam(a::Array{Float64,1}, numQBits::Int64)::Opr
    @assert length(x0) == (2^(numQBits + 1) - 2) "Wrong number of initial parameters"
    if numQBits == 1
        return uParam1QB(a)
    elseif numQBits == 2
        return uParam2QB(a)
    elseif numQBits == 3
        return uParam3QB(a)
    elseif numQBits == 4
        return uParam4QB(a)
    elseif numQBits == 5
        return uParam5QB(a)
    elseif numQBits == 6
        return uParam6QB(a)
    elseif numQBits == 7
        return uParam7QB(a)
    elseif numQBits == 8
        return uParam8QB(a)
    else
        @assert false "Not coded this number of QBits"
    end
    return scalarArray(ComplexF64(1, 0)) # Dummy return hence choose simplest Opr.
end

# An actual run of the Quantum circuit (with optional printouts).
function cost(a::Array{Float64,1}, initReg::Ket, hOpr::Opr, numQBits::Int64, print::Bool)::Float64
    register = initReg                          # Set up register.
    uOpr = uParam(a, numQBits)                  # Get the appropriate (parameterised) circuit.
    register = uOpr * register                  # Apply the circuit to the register.
    eVal = oprExpectation(hOpr, register)       # Evaluate/take expectation with respect to H.
    if print
        println("Estimated eigenval:  ", eVal.re)
        # println("Estimated eigenvec:  ", register.vals[1], ", ", register.vals[2])  # TODO: longer evals print
        # println("Eigenpair residual:  ", norm((hOpr * register).vals - (eVal.re * register.vals)))
    end
    return eVal.re
end

# ------------------
#       Script      
# ------------------

# ----- Inputs -----
numQBits = 2

# ----- Setup -----
initReg = rBs0()^numQBits                   # Must be unit Ket vector 
x0 = zeros((2^(numQBits + 1) - 2))          # Params must be correct length
mtd = "L-BFGS-B"                              # Type of optimizer
dimension = 2^numQBits                      # <-|
# H = randP(ComplexF64, dimension)            #   |-> Generate Hermitian matrix (find its smallest eval)
H = Hermitian([1.84764+0.0im          0.481679-0.686225im  -0.0646636+0.376332im  -0.351826+0.00449551im;
   0.481679+0.686225im      1.89892+0.0im        -0.347868-0.406104im   0.270096-0.473325im;
 -0.0646636-0.376332im    -0.347868+0.406104im     1.38095+0.0im        0.104298+0.57518im;
  -0.351826-0.00449551im   0.270096+0.473325im    0.104298-0.57518im     1.18633+0.0im])
hOpr = Opr(H)                               # <-|

# ----- Capture variables in Lambda -----
function costLam(a::Array{Float64,1})::Float64
    return cost(a, initReg, hOpr, numQBits, false)
end

# res = SciPy.optimize.minimize(costLam, x0, method = mtd, jac = jacFn, tol = 1e-36, options=Dict("maxiter"=>1e6))

# for i = 1:10
i = 1
    # x0 = 2*pi*rand((2^(numQBits + 1) - 2))          # Params must be correct length

# ----- Classical Optimizer -----
println(length(x0))
res = SciPy.optimize.minimize(costLam, x0, method = mtd, jac = jacFn, tol = 1e-36, options=Dict("maxiter"=>1e6))
# res = SciPy.optimize.minimize(costLam, x0, method = mtd, tol = 1e-36, options=Dict("maxiter"=>1e6))

# ----- Error check -----
eValRe = cost(res["x"], initReg, hOpr, numQBits, true)
println(length(res["x"]))
trueEVals = eigvals(H, 1:2)
trueEVal1 = trueEVals[1]
trueEVal2 = trueEVals[2]
println("Calculated eigenval: ", trueEVal1)
println("Calculated eigenval2:", trueEVal2)
println("Eigenval error:      ", abs(trueEVal1 - eValRe))
println("Relative error:      ", abs(eValRe-trueEVal1)/abs(trueEVal1 - trueEVal2))
println()
println()

# ----- Classical Optimizer Debug Info -----
# println(res)

# end
println(length(res["x"]))