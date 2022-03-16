import Base: ==, +, -, *, /, \, inv, ^, <, <=, <<, >>, %, ÷

# row vector = bra
mutable struct Bra
    vals::Array{ComplexF64, 1}
end 

# col vector = ket
mutable struct Ket
    vals::Array{ComplexF64, 1}
end 

# bra to ket
function btk(a::Bra)::Ket
    return Ket(conj(a.vals))
end

# ket to bra
function ktb(a::Ket)::Bra
    return Bra(conj(a.vals))
end

# dot product 
function (*)(a::Bra, b::Ket)::ComplexF64
    @assert size(a.vals,1) == size(b.vals,1) "Dimensions of dot product disagree!"
    elemMult::Array{ComplexF64, 1} = a.vals .* b.vals
    sums::ComplexF64 = sum(elemMult)
    return sums
end 

# ket element wise mult 
function (*)(a::ComplexF64, b::Ket)::Ket
    return Ket(a .* b.vals)
end 

# bra element wise mult 
function (*)(a::ComplexF64, b::Bra)::Bra
    return Bra(a .* b.vals)
end 

# tensor kets 
function (*)(a::Ket, b::Ket)::Ket
    out::Array{ComplexF64, 1} = []
    for i in a.vals
        temp::Ket = i*b
        out = vcat(out, temp.vals)
    end
    return Ket(out)
end 

# tensor bras 
function (*)(a::Bra, b::Bra)::Bra
    return ktb(btk(a)*btk(b))
end 

# matrix = operator (TODO: error with 1x1 Opr)
mutable struct Opr
    vals::Array{ComplexF64, 2}
end 

# tensor ket and bra 
function (*)(a::Ket, b::Bra)::Opr
    dim1::Int64 = size(a.vals,1)
    dim2::Int64 = size(b.vals,1)
    arr::Array{ComplexF64, 2} = Array{ComplexF64, 2}(undef, dim1, dim2)
    for i = 1:dim1
        for j = 1:dim2
            arr[i,j] = a.vals[i] * b.vals[j]
        end 
    end
    return Opr(arr)
end 

# matrix mult
function matmul(a::Opr, b::Opr)::Opr
    @assert size(a.vals,2) == size(b.vals,1) "Dimensions of matrix product disagree!"
    return Opr(a.vals * b.vals)
end

# bra to operator
function bto(a::Bra)::Opr
    dim1::Int64 = size(a.vals,1)
    arr::Array{ComplexF64, 2} = Array{ComplexF64, 2}(undef, 1, dim1)
    for i = 1:dim1
        arr[1,i] = a.vals[i]
    end
    return Opr(arr)
end

# ket to operator
function kto(a::Ket)::Opr
    dim1::Int64 = size(a.vals,1)
    arr::Array{ComplexF64, 2} = Array{ComplexF64, 2}(undef, dim1, 1)
    for i = 1:dim1
        arr[i,1] = a.vals[i]
    end
    return Opr(arr)
end

# row matrix mult
function (*)(a::Bra, b::Opr)::Bra
    aOpr::Opr = bto(a)
    c::Opr = matmul(aOpr,b)
    return Bra(vec(c.vals))
end 

# matrix col mult
function (*)(a::Opr, b::Ket)::Ket
    bOpr::Opr = kto(b)
    c::Opr = matmul(a,bOpr)
    return Ket(vec(c.vals))
end 

# tensor operators
function (*)(a::Opr, b::Opr)::Opr
    dimA1::Int64 = size(a.vals,1)
    dimA2::Int64 = size(a.vals,2)
    dimB1::Int64 = size(b.vals,1)
    dimB2::Int64 = size(b.vals,2)
    arr::Array{ComplexF64, 2} = Array{ComplexF64, 2}(undef, dimA1*dimB1, dimA2*dimB2)

    for i = 1:dimA1
        for k = 1:dimB1
            for j = 1:dimA2
                for l = 1:dimB2
                    arr[(i-1)*dimB1+k, (j-1)*dimB2+l] = a.vals[i,j]*b.vals[k,l]
                end
            end
        end
    end

    return Opr(arr)
end

# matrix scalar mult
function (*)(a::ComplexF64, b::Opr)::Opr
    return Opr(a .* b.vals)

    ### Check code errors ###
    # scalarOpr::Opr = Opr([a])
    # return scalarOpr*b 
    ### Test with ###
    # println(Opr([ComplexF64(1,0)]))
    # println(Opr([ComplexF64(1,0) ComplexF64(2,0)]))
end 

# create 2x2 operator
function mat2(a::ComplexF64, b::ComplexF64, c::ComplexF64, d::ComplexF64)::Opr
    return Opr([a b; c d])
end

# create 2x2 identity
function eye2()::Opr
    a::ComplexF64 = ComplexF64(1,0)
    b::ComplexF64 = ComplexF64(0,0)
    return mat2(a, b, b, a)
end

# create hadamard gate 
function had2()::Opr
    nrm::ComplexF64 = ComplexF64(1/sqrt(2), 0)
    op::Opr = mat2(ComplexF64(1,0), ComplexF64(1,0), ComplexF64(1,0), ComplexF64(-1,0))
    return nrm*op
end

# create pauli X gate 
function paX2()::Opr
    a::ComplexF64 = ComplexF64(0,0)
    b::ComplexF64 = ComplexF64(1,0)
    return mat2(a, b, b, a)
end

# create pauli Y gate 
function paY2()::Opr
    a::ComplexF64 = ComplexF64(0,0)
    b::ComplexF64 = ComplexF64(0,1)
    return mat2(a, b, b, a)
end

# create pauli Z gate 
function paZ2()::Opr
    a::ComplexF64 = ComplexF64(1,0)
    b::ComplexF64 = ComplexF64(0,0)
    c::ComplexF64 = ComplexF64(-1,0)
    return mat2(a, b, b, c)
end

# prob measure Ket vector
function pMeas(a::Ket)::Array{Float64, 1}
    dim1::Int64 = size(a.vals,1)
    probs::Array{Float64, 1} = zeros(Float64, dim1)
    for i = 1:dim1
        val::ComplexF64 = a.vals[i]
        probs[i] = sqrt(val.real*val.real + val.imag*val.imag)
    end
    return probs
end

# prob measure Bra vector
function pMeas(a::Bra)::Array{Float64, 1}
    return pMeas(btk(a))
end

# single (most likely) measure Ket vector
function sMeas(a::Ket)::Int64
    dim1::Int64 = size(a.vals,1)
    currentProb::Float64 = 0.0
    highestProb::Float64 = 0.0
    maxI::Int64 = -1    # for error detection (-ve => error)
    for i = 1:dim1
        val::ComplexF64 = a.vals[i]
        currentProb = sqrt(val.real*val.real + val.imag*val.imag)
        if (currentProb > highestProb)
            maxI = i
            highestProb = currentProb
        end
    end
    return maxI
end

# single (most likely) measure Bra vector
function sMeas(a::Bra)::Int64
    return sMeas(btk(a))
end

function oprExp(a::Opr)::Opr
    return Opr(exp(a.vals))
end

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