"""
Functions for Operators.

Operator is local operator represented by a matrix, together with a list of
indices.

Operation is the collection of Operators, together with a Basis object that
store the base and length information. Functions on Operations are design to
mimic the behavior of many-body matrices.
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
#--- Type Operation
"""
    Operations{To<:AbstractVector{Operator},Tb<:Basis}

Object that store all information of a many-body operator.
"""
struct Operation{T1<:AbstractVector{<:Operator}, T2<:Basis}
    opts::T1
    basis::T2
end
#--- Operation initiate
"""
    operation(mats::AbstractVector{<:AbstractMatrix}, inds::AbstractVector{<:AbstractVector}, base::Integer, len::Integer)

Create an Operation object from a list of matrix `mats`, a list of indices
`inds`, a base number `base`, and system size `len`.
"""
function operation(
    mats::AbstractVector{<:AbstractMatrix},
    inds::AbstractVector{<:AbstractVector},
    base::Integer,
    len::Integer
)
    b = basis(base, len)
    opts = [Operator(mats[i],inds[i]) for i=1:length(mats)]
    Operation(opts, b)
end
#--- Apply to vec
function fillvec!(vec::AbstractVector, op::Operation, j::Integer; c::Number=1)
    len = length(op.opts)
    b = op.basis
    ol = op.opts
    for i = 1:len
        change!(b, j)
        fillvec!(vec, ol[i], b, c=c)
    end
end
function mul!(vec::AbstractVector, op::Operation, state::AbstractVector)
    len = length(op.opts)
    b = op.basis
    ol = op.opts
    for j = 1:length(state)
        for i = 1:len
            change!(b, j)
            fillvec!(vec, ol[i], b, c=state[j])
        end
    end
end
function mul!(mat::AbstractMatrix, op::Operation, states::AbstractMatrix)
    len = length(op.opts)
    b = op.basis
    ol = op.opts
    for j = 1:size(states,1)
        for i = 1:len
            change!(b, j)
            o = ol[i]
            vb = view(b, o.inds)
            oi = o.mat[:, index(vb)]
            for k = 1:size(o.mat, 1)
                change!(vb, k)
                m = index(b)
                mat[m, :] .+= oi[k] * states[j,:]
            end
        end
    end
end
function mul(op::Operation, vom::AbstractVecOrMat{T}) where T <: Number
    out = zeros(T, size(vom))
    mul!(out, op, vom)
    out
end
#--- fill matrix
"""
    fillmat!(mat::AbstractMatrix, op::Operation)

Add to the matrix `mat` with operation `op`.
"""
function fillmat!(mat::AbstractMatrix, op::Operation)
    for j = 1:size(mat,1)
        fillvec!(view(mat, :, j), op, j)
    end
end
