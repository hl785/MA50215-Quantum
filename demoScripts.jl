using Plots
include("main.jl")

#----------------------------------------------------------------------------------------------------
# If the register is in the state Rot(PauliX, theta) * |0>, then expectation wrt PauliZ is cos(theta)
theta = 1.0                                             # Test parameter
trueRes = cos(theta)                                    # Target result

register = Ket([ComplexF64(1,0), ComplexF64(0,0)])      # Set up in initial condition |0>
opr = rotGate(paX2(), theta)                            # Create a gate
register = opr * register                               # Apply a gate to the register
testRes = oprExpectation(paZ2(), register)              # Measure / expectation

println("\nTest rotation gates and expectations:\n")    # Print
println("True result: ", trueRes)
println("Quan result: ", testRes)
println("Error: ", abs(trueRes - testRes))
#----------------------------------------------------------------------------------------------------
# Understand how CNot behaves for different sizes and configurations.
println("\nTest CNot gate on non neighboring gates:\n") # Print
println("(Display graph if run through repl)")
println(broadcast(abs, CNot(2,1,3).vals))
heatmap(broadcast(abs, CNot(2,1,3).vals), yflip = true)
#----------------------------------------------------------------------------------------------------