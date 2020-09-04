"""
Functions for Operations.

Operation is the collection of Operators, together with a Basis object that
store the base and length information. Functions on Operations are design to
mimic the behavior of many-body matrices.
"""
#-----------------------------------------------------------------------------------------------------
# Type: Operator
#-----------------------------------------------------------------------------------------------------
struct Operator{
    Tm<:AbstractMatrix, 
    Ti<:AbstractVector{<:Integer}
}
    mat::Tm
    inds::Ti
end
#-----------------------------------------------------------------------------------------------------
# Basis functions for type Operator
#-----------------------------------------------------------------------------------------------------
function *(
    c::Number, 
    o::Operator
)
    Operator(c * o.mat, o.inds)
end
#-----------------------------------------------------------------------------------------------------
function fillvec!(
    vec::AbstractVector, 
    o::Operator, 
    b::Basis; 
    c::Number=1
)
    vb = view(b, o.inds)
    oi = o.mat[:, index(vb)] * c
    for k = 1:size(o.mat, 1)
        change!(vb, k)
        i = index(b)
        vec[i] += oi[k]
    end
end
#-----------------------------------------------------------------------------------------------------
# Type: Operation
#-----------------------------------------------------------------------------------------------------
export Operation
struct Operation{
    T1<:AbstractVector{<:Operator}, 
    T2<:Basis
}
    opts::T1
    basis::T2
end
#-----------------------------------------------------------------------------------------------------
# Operation initiate
#-----------------------------------------------------------------------------------------------------
export operation
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
#-----------------------------------------------------------------------------------------------------
export onsite_operation
function onsite_operation(
    mats::AbstractVector{<:AbstractMatrix},
    inds::AbstractVector{<:Integer},
    base::Integer,
    len::Integer
)
    ind = [[i] for i in inds]
    operation(mats,ind,base,len)
end
function operation(
    mats::AbstractVector{<:AbstractMatrix},
    inds::AbstractVector{<:Integer},
    base::Integer,
    len::Integer
)
    onsite_operation(mats, inds, base, len)
end
#-----------------------------------------------------------------------------------------------------
export trans_inv_operation
function trans_inv_operation(
    mat::AbstractMatrix,
    ind::AbstractVector{<:Integer},
    base::Integer,
    len::Integer
)
    mats = fill(mat, len)
    inds = [mod.(ind .+ (i-1), len) .+ 1 for i=0:len-1]
    operation(mats,inds,base,len)
end
#-----------------------------------------------------------------------------------------------------
# Basic Functions
#-----------------------------------------------------------------------------------------------------
function *(
    c::Number, 
    o::Operation
)
    Operation(c .* o.opts, o.basis)
end
#-----------------------------------------------------------------------------------------------------
function +(
    o1::Operation, 
    o2::Operation
)
    Operation(vcat(o1.opts, o2.opts), o1.basis)
end
#-----------------------------------------------------------------------------------------------------
function sum(
    ol::AbstractVector{<:Operation}
)
    basis = ol[1].basis
    opts = vcat([oi.opts for oi in ol]...)
    Operation(opts, basis)
end
#-----------------------------------------------------------------------------------------------------
# Apply to vector
#-----------------------------------------------------------------------------------------------------
function fillvec!(
    vec::AbstractVector, 
    op::Operation, 
    j::Integer; 
    c::Number=1
)
    len = length(op.opts)
    b = op.basis
    ol = op.opts
    for i = 1:len
        change!(b, j)
        fillvec!(vec, ol[i], b, c=c)
    end
end
#-----------------------------------------------------------------------------------------------------
function mul!(
    vec::AbstractVector, 
    op::Operation, 
    state::AbstractVector
)
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
#-----------------------------------------------------------------------------------------------------
# Apply to matrix
#-----------------------------------------------------------------------------------------------------
export mul!
function mul!(
    mat::AbstractMatrix, 
    op::Operation, 
    states::AbstractMatrix
)
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
#-----------------------------------------------------------------------------------------------------
# General multiplying
#-----------------------------------------------------------------------------------------------------
export mul
function mul(
    op::Operation, 
    vom::AbstractVecOrMat{T}
) where T <: Number

    To = promote_type([eltype(o.mat) for o in op.opts]...)
    out = zeros(promote_type(To, T), size(vom))
    mul!(out, op, vom)
    out
end
#-----------------------------------------------------------------------------------------------------
function *(
    op::Operation, 
    vom::AbstractVecOrMat
)
    mul(op, vom)
end
#-----------------------------------------------------------------------------------------------------
# Fill matrix
#-----------------------------------------------------------------------------------------------------
export fillmat!
function fillmat!(
    mat::AbstractMatrix, 
    op::Operation
)
    for j = 1:size(mat,1)
        fillvec!(view(mat, :, j), op, j)
    end
end
