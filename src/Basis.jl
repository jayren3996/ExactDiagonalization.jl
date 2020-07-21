"""
Functions for basis.

The core structure of ExactDiagonalization.jl is the mapping between different
Basis objects.
"""
#--- Type Basis
"""
    Basis{T<:AbstractVector{<:Integer}}

Basis object for exact diagonalization. Linear operator can be thought as a
map from one `Basis` object to another `Basis` object.
"""
struct Basis{T<:AbstractVector{<:Integer}}
    bits::T
    base::Int64
    len::Int64
end
#--- Initiation
"""
    basis(len::Integer, base::Integer)

Initiate a zero Basis object with length `len` and with local Hilbert space
dimension `base`.
"""
function basis(base::Integer, len::Integer)
    Basis(zeros(Int64, len), Int64(base), Int64(len))
end
"""
    basis(bits::AbstractVector{<:Integer}, base::Integer)

Initiate a Basis object with a bit list.
"""
function basis(bits::AbstractVector{<:Integer}, base::Integer)
    Basis(bits, Int64(base), length(bits))
end
#--- Convertion between index numbers and bits lists.
function list2num(bits::AbstractVector, base::Integer, len::Integer)
    sum(bits[i] * base^(len - i) for i = 1:len) + 1
end
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
#--- Indexing
"""
    index(b::Basis)

Return the index number which the `Basis` object represent.
"""
function index(b::Basis)
    list2num(b.bits, b.base, b.len)
end
#--- Change the Basis
function change!(b::Basis, index::Integer)
    num2list!(b.bits, index, b.base, b.len)
end
function change!(b::Basis, bits::AbstractVector{<:Integer})
    b.bits .= bits
end
#--- Copy Basis
function copy(b::Basis)
    Basis(copy(b.bits), b.base, b.len)
end
#--- View Basis
function view(b::Basis, inds::AbstractVector)
    Basis(view(b.bits, inds), b.base, length(inds))
end
