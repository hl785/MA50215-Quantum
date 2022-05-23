function jacFn(x::Array{Float64,1})::Array{Float64,1}

out = zeros(length(x))
h = 1e-12

for i = 1:length(x)
    pert = zeros(length(x))
    pert[i] = h
    val1 = costLam(x)
    val2 = costLam(x+pert)
    out[i] = (val2 - val1)/h
    # println(abs(val1 - val2)/h)
    # println()

end
return out
end 