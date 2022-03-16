import Base: ==, +, -, *, /, \, inv, ^, <, <=, <<, >>, %, ÷

mutable struct Complex
    real::Float64
    imag::Float64
    # way of getting +=, -=, *=, /= ???
end

function (+)(a::Complex, b::Complex)::Complex
    return Complex(a.real + b.real, a.imag + b.imag)
end 

function (*)(a::Complex, b::Complex)::Complex
    return Complex((a.real * b.real) - (a.imag * b.imag), (a.real * b.imag) + (a.imag * b.real))
end 

# row vector = bra
mutable struct Bra
    vals::Array{Complex, 1}
end 

# col vector = ket
mutable struct Ket
    vals::Array{Complex, 1}
end 

# bra to ket
function btk(a::Bra)::Ket
    newVals::Array{Complex, 1} = a.vals
    dim1 = size(newVals,1)
    for i = 1:dim1
        oldVal::Complex = newVals[i]
        newVals[i] = Complex(oldVal.real, -oldVal.imag)
    end
    return Ket(newVals)
end

# ket to bra
function ktb(a::Ket)::Bra
    newVals::Array{Complex, 1} = a.vals
    dim1 = size(newVals,1)
    for i = 1:dim1
        oldVal::Complex = newVals[i]
        newVals[i] = Complex(oldVal.real, -oldVal.imag)
    end
    return Bra(newVals)
end

# dot product 
function (*)(a::Bra, b::Ket)::Complex
    @assert size(a.vals,1) == size(b.vals,1) "Dimensions of dot product disagree!"
    elemMult::Array{Complex, 1} = a.vals .* b.vals
    sums::Complex = sum(elemMult)
    return sums
end 

# ket element wise mult 
function (*)(a::Complex, b::Ket)::Ket
    out::Array{Complex, 1} = []
    for i in b.vals
        out = vcat(out, a * i)
    end
    return Ket(out)
end 

# bra element wise mult 
function (*)(a::Complex, b::Bra)::Bra
    return ktb(a * btk(b))
end 

# tensor kets 
function (*)(a::Ket, b::Ket)::Ket
    out::Array{Complex, 1} = []
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
    vals::Array{Complex, 2}
end 

# tensor ket and bra 
function (*)(a::Ket, b::Bra)::Opr
    dim1::Int64 = size(a.vals,1)
    dim2::Int64 = size(b.vals,1)
    arr::Array{Complex, 2} = Array{Complex, 2}(undef, dim1, dim2)
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
    dim1::Int64 = size(a.vals,1)
    dim2::Int64 = size(b.vals,2)
    dimInt::Int64 = size(a.vals,2)
    arr::Array{Complex, 2} = Array{Complex, 2}(undef, dim1, dim2)
    for i = 1:dim1
        for j = 1:dim2
            sum::Complex = Complex(0,0)
            for k = 1:dimInt
                sum = sum + (a.vals[i,k]*b.vals[k,j])       # Use dot prod?
            end
            arr[i,j] = sum
        end 
    end
    return Opr(arr)
end

# bra to operator
function bto(a::Bra)::Opr
    dim1::Int64 = size(a.vals,1)
    arr::Array{Complex, 2} = Array{Complex, 2}(undef, 1, dim1)
    for i = 1:dim1
        arr[1,i] = a.vals[i]
    end
    return Opr(arr)
end

# ket to operator
function kto(a::Ket)::Opr
    dim1::Int64 = size(a.vals,1)
    arr::Array{Complex, 2} = Array{Complex, 2}(undef, dim1, 1)
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
    arr::Array{Complex, 2} = Array{Complex, 2}(undef, dimA1*dimB1, dimA2*dimB2)

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
function (*)(a::Complex, b::Opr)::Opr
    dim1::Int64 = size(b.vals,1)
    dim2::Int64 = size(b.vals,2)
    arr::Array{Complex, 2} = Array{Complex, 2}(undef, dim1, dim2)
    for i = 1:dim1
        for j = 1:dim2
            arr[i,j] = a * b.vals[i,j]
        end 
    end
    return Opr(arr)


    ### This code errors ###
    # scalarOpr::Opr = Opr([a])
    # return scalarOpr*b 
    ### Test with ###
    # println(Opr([Complex(1,0)]))
    # println(Opr([Complex(1,0) Complex(2,0)]))
end 

# create 2x2 operator
function mat2(a::Complex, b::Complex, c::Complex, d::Complex)::Opr
    return Opr([a b; c d])
end

# create 2x2 identity
function eye2()::Opr
    a::Complex = Complex(1,0)
    b::Complex = Complex(0,0)
    return mat2(a, b, b, a)
end

# create hadamard gate 
function had2()::Opr
    nrm::Complex = Complex(1/sqrt(2), 0)
    op::Opr = mat2(Complex(1,0), Complex(1,0), Complex(1,0), Complex(-1,0))
    return nrm*op
end

# # qubit = impure states are 
# mutable struct MonoQubit
#     numStates::Int64
#     ket::Array{Ket, numStates}
# end 

# TODO: bra/ket constructor

a = Complex(1,2)
b = Complex(4,6)
# println(a+b)
# println(a*b)
 x = [a, b]
# y = Bra(x)
# x = [a+b, a*b]
#  z = Ket(x)
# dotProd = y*z
# println(y)
# println(z)
# println(dotProd)
# elemWiseMult = a*z
# println(elemWiseMult)
# println(Ket([Complex(0,1)])*z)
# z = Bra(x)
# println(Bra([Complex(1,0),Complex(0,1)])*z)
# z = Opr([a b; b a])
# println(z)
# z = Ket(x) * Bra(x)
# println(z)
# println(Ket([Complex(1,0), Complex(2,0), Complex(3,0)]) * Bra([Complex(4,0), Complex(5,0)]))
# println(matmul(Opr([Complex(1,0) Complex(0,-1); Complex(1,1) Complex(4,-1)]), Opr([Complex(0,1) Complex(1,-1); Complex(2,-3) Complex(4,0)])))
# println(matmul(Opr([Complex(1,0) Complex(0,0) Complex(0,0); Complex(0,0) Complex(1,0) Complex(0,0)]), Opr([Complex(0,1) Complex(1,-1); Complex(2,-3) Complex(4,0); Complex(1,-3) Complex(6,7)])))
# println(Bra(x)*z)
# println(Opr([Complex(1,0) Complex(2,0); Complex(3,0) Complex(4,0); Complex(1,0) Complex(0,0)])*Opr([Complex(0,0) Complex(5,0) Complex(2,0); Complex(6,0) Complex(7,0) Complex(3,0)]))
# println(eye2()*eye2())
# println(had2()*Ket([Complex(0,0), Complex(1,0)]))
x = Ket([Complex(1,0), Complex(0,0)])
y = Ket([Complex(0,0), Complex(0,1)])
z = Ket([Complex(0,0), Complex(1,0)])
xyz = (x*y)*z
println(xyz)
xyz = x*(y*z)
println(xyz)
xyz = x*y*z
println(xyz)