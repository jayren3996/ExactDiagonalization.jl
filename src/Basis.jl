"""
Functions for basis.

The core structure of ExactDiagonalization.jl is the mapping between different
Basis objects.
"""
#-----------------------------------------------------------------------------------------------------
# Type: Operator
#-----------------------------------------------------------------------------------------------------
struct Basis{
    T<:AbstractVector{<:Integer}
}
    bits::T
    base::Int64
    len::Int64
end
#-----------------------------------------------------------------------------------------------------
# Initiation
#-----------------------------------------------------------------------------------------------------
function basis(
    base::Integer, 
    len::Integer
)
    Basis(zeros(Int64, len), Int64(base), Int64(len))
end
#-----------------------------------------------------------------------------------------------------
function basis(
    bits::AbstractVector{<:Integer}, 
    base::Integer
)
    Basis(bits, Int64(base), length(bits))
end
#-----------------------------------------------------------------------------------------------------
# Convertion between index numbers and bits lists
#-----------------------------------------------------------------------------------------------------
function list2num(
    bits::AbstractVector, 
    base::Integer, 
    len::Integer
)
    sum(bits[i] * base^(len - i) for i = 1:len) + 1
end
#-----------------------------------------------------------------------------------------------------
function num2list!(
    bits::AbstractVector,
    index::Integer,
    base::Integer,
    len::Integer,
)
    i = index - 1
    for j = 1:len
        n, i = divrem(i, base^(len - j))
        bits[j] = n
    end
end
#-----------------------------------------------------------------------------------------------------
# Indexing
#-----------------------------------------------------------------------------------------------------
function index(b::Basis)
    list2num(b.bits, b.base, b.len)
end
#-----------------------------------------------------------------------------------------------------
# Change basis
#-----------------------------------------------------------------------------------------------------
function change!(
    b::Basis, 
    index::Integer
)
    num2list!(b.bits, index, b.base, b.len)
end
#-----------------------------------------------------------------------------------------------------
function change!(
    b::Basis, 
    bits::AbstractVector{<:Integer}
)
    b.bits .= bits
end
#-----------------------------------------------------------------------------------------------------
# Copy Basis
#-----------------------------------------------------------------------------------------------------
function copy(b::Basis)
    Basis(copy(b.bits), b.base, b.len)
end
#-----------------------------------------------------------------------------------------------------
# View Basis
#-----------------------------------------------------------------------------------------------------
function view(
    b::Basis, 
    inds::AbstractVector
)
    Basis(view(b.bits, inds), b.base, length(inds))
end
