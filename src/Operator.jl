"""
Functions for Operators.

Operator is local operator represented by a matrix, together with a list of
indices.
"""
#--- Type Operator
struct Operator{Tm<:AbstractMatrix, Ti<:AbstractVector{<:Integer}}
    mat::Tm
    inds::Ti
end
#--- Filling vectors
function fillvec!(vec::AbstractVector, o::Operator, b::Basis; c::Number=1)
    vb = view(b, o.inds)
    oi = o.mat[:, index(vb)] * c
    for k = 1:size(o.mat, 1)
        change!(vb, k)
        i = index(b)
        vec[i] += oi[k]
    end
end
#--- Filling mat
function fillmat!(mat::AbstractMatrix, o::Operator, b::Basis)
    for j = 1:size(mat,1)
        change!(b,j)
        fillvec!(view(mat,:,j), b, ol)
    end
end
