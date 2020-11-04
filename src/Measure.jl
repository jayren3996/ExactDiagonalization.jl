#-----------------------------------------------------------------------------------------------------
# Measurement
#-----------------------------------------------------------------------------------------------------
export measure
function measure(
    o::Operation, 
    v::AbstractVecOrMat
)
    dot(v, mul(o, v))
end
#-----------------------------------------------------------------------------------------------------
function measure(
    o1::Operation,
    o2::Operation,
    v::AbstractVecOrMat
)
    dot(mul(o1, v), mul(o2, v))
end
#-----------------------------------------------------------------------------------------------------
export covmat
function covmat(
    ol::Vector{<:Operation}, 
    v::AbstractVector
)
    n = length(ol)
    am = Vector{Float64}(undef, n)
    cm = Matrix{Float64}(undef, n, n)
    for i=1:n
        am[i] = real(measure(ol[i], v))
    end
    for i=1:n 
        for j=i:n 
            cm[i,j] = real(measure(ol[i], ol[j], v)) - am[i] * am[j]
        end
    end
    Hermitian(cm)
end
#-----------------------------------------------------------------------------------------------------
function covmat(
    ol::Vector{<:Operation}, 
    v::AbstractMatrix
)
    n = length(ol)
    am = Vector{Float64}(undef, n)
    cm = Matrix{Float64}(undef, n, n)
    num = size(v, 2)
    for i=1:n
        am[i] = real(measure(ol[i], v))/num
    end
    for i=1:n 
        for j=i:n 
            cm[i,j] = real(measure(ol[i], ol[j], v))/num - am[i] * am[j]
        end
    end
    Hermitian(cm)
end
#-----------------------------------------------------------------------------------------------------
export covmatc
function covmatc(
    ol::Vector{<:Operation}, 
    v::AbstractVector
)
    n = length(ol)
    am = Vector{ComplexF64}(undef, n)
    cm = Matrix{ComplexF64}(undef, n, n)
    for i=1:n
        am[i] = measure(ol[i], v)
    end
    for i=1:n 
        for j=i:n 
            cm[i,j] = measure(ol[i], ol[j], v) - conj(am[i]) * am[j]
        end
    end
    Hermitian(cm)
end
#-----------------------------------------------------------------------------------------------------
function covmatc(
    ol::Vector{<:Operation}, 
    v::AbstractMatrix
)
    n = length(ol)
    am = Vector{ComplexF64}(undef, n)
    cm = Matrix{ComplexF64}(undef, n, n)
    num = size(v, 2)
    for i=1:n
        am[i] = measure(ol[i], v)/num
    end
    for i=1:n
        for j=i:n
            cm[i,j] = measure(ol[i], ol[j], v)/num - conj(am[i]) * am[j]
        end
    end
    Hermitian(cm)
end
#-----------------------------------------------------------------------------------------------------
# Entropy
#-----------------------------------------------------------------------------------------------------
export entropy
function entropy(
    v::AbstractVector, 
    l::Integer,
    i::Integer
)
    base = Int(length(v)^(1/l))
    m = reshape(v, base^i, :)
    s = svdvals(m)
    st = s[s .> 1e-10]
    - dot(st, log.(st))
end