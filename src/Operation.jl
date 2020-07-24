"""
Functions for Operations.

Operation is the collection of Operators, together with a Basis object that
store the base and length information. Functions on Operations are design to
mimic the behavior of many-body matrices.
"""
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
function operation(
    mats::AbstractVector{<:AbstractMatrix},
    inds::AbstractVector{<:Integer},
    base::Integer,
    len::Integer
)
    ind = [[i] for i in inds]
    operation(mats,ind,base,len)
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
    To = promote_type([eltype(o.mat) for o in op.opts]...)
    out = zeros(promote_type(To, T), size(vom))
    mul!(out, op, vom)
    out
end
*(op::Operation, vom::AbstractVecOrMat) = mul(op, vom)
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
