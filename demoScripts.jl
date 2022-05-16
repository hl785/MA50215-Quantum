include("main.jl")

#----------------------------------------------------------------------------------------------------
# If the register is in the state Rot(PauliX, theta) * |0>, then expectation wrt PauliZ is cos(theta)
theta = 1.0                                             # Test parameter
trueRes = cos(theta)                                    # Targer result

register = Ket([ComplexF64(1,0), ComplexF64(0,0)])      # Set up in initial condition |0>
opr = rotGate(paX2(), theta)                            # Create a gate
register = opr * register                               # Apply a gate to the register
testRes = oprExpectation(paZ2(), register)              # Measure / expectation

println("\nTest rotation gates and expectations:\n")    # Print
println("True result: ", trueRes)
println("Quan result: ", testRes)
println("Error: ", abs(trueRes - testRes))
#----------------------------------------------------------------------------------------------------
using Plots
println("\nTest CNot gate on non neighboring gates:\n") # Print
println("(Display graph if run through repl)")
println(broadcast(abs, CNot(2,1,3).vals))
heatmap(broadcast(abs, CNot(2,1,3).vals), yflip = true)
#----------------------------------------------------------------------------------------------------

# a = ComplexF64(1,2)
# b = ComplexF64(4,6)
# println(a+b)
# println(a*b)
#  x = [a, b]
# y = Bra(x)
# x = [a+b, a*b]
#  z = Ket(x)
# dotProd = y*z
# println(y)
# println(z)
# println(dotProd)
# elemWiseMult = a*z
# println(elemWiseMult)
# println(Ket([ComplexF64(0,1)])*z)
# z = Bra(x)
# println(Bra([ComplexF64(1,0),ComplexF64(0,1)])*z)
# z = Opr([a b; b a])
# println(z)
# z = Ket(x) * Bra(x)
# println(z)
# println(Ket([ComplexF64(1,0), ComplexF64(2,0), ComplexF64(3,0)]) * Bra([ComplexF64(4,0), ComplexF64(5,0)]))
# println(matmul(Opr([ComplexF64(1,0) ComplexF64(0,-1); ComplexF64(1,1) ComplexF64(4,-1)]), Opr([ComplexF64(0,1) ComplexF64(1,-1); ComplexF64(2,-3) ComplexF64(4,0)])))
# println(matmul(Opr([ComplexF64(1,0) ComplexF64(0,0) ComplexF64(0,0); ComplexF64(0,0) ComplexF64(1,0) ComplexF64(0,0)]), Opr([ComplexF64(0,1) ComplexF64(1,-1); ComplexF64(2,-3) ComplexF64(4,0); ComplexF64(1,-3) ComplexF64(6,7)])))
# println(Bra(x)*z)
# println(Opr([ComplexF64(1,0) ComplexF64(2,0); ComplexF64(3,0) ComplexF64(4,0); ComplexF64(1,0) ComplexF64(0,0)])*Opr([ComplexF64(0,0) ComplexF64(5,0) ComplexF64(2,0); ComplexF64(6,0) ComplexF64(7,0) ComplexF64(3,0)]))
# println(eye2()*eye2())
# println(had2()*Ket([ComplexF64(0,0), ComplexF64(1,0)]))
# x = Ket([ComplexF64(1,0), ComplexF64(0,0)])
# y = Ket([ComplexF64(0,0), ComplexF64(0,1)])
# z = Ket([ComplexF64(0,0), ComplexF64(1,0)])
# xyz = (x*y)*z
# println(xyz)
# xyz = x*(y*z)
# println(xyz)
# xyz = x*y*z
# println(xyz)
# println(oprExp(eye2()))
# println(ComplexF64(1,0)*ComplexF64(1,1))
# println(conj([ComplexF64(1,2), ComplexF64(0,1)]))
# println(ComplexF64(1,-2).*Ket([ComplexF64(1,2), ComplexF64(0,1)]).vals)
# println(size(scalarArray(ComplexF64(1,0)).vals))
# println(CNot(0,1,2))
# println(rotGate(eye2(), Float64(pi, RoundDown)))
# println(oprExpectation(paZ2(), rotGate(paX2(), 1.0)*Ket([ComplexF64(1,0), ComplexF64(0,0)])), cos(1.0))