"""
Functions for Operations.

Operation is the collection of Operators, together with a Basis object that
store the base and length information. Functions on Operations are design to
mimic the behavior of many-body matrices.
"""
#-----------------------------------------------------------------------------------------------------
# Type: Operator & Operation
#-----------------------------------------------------------------------------------------------------
struct Operator{
    MatType <: AbstractMatrix, 
    IndType <: AbstractVector{<:Integer}
}
    mat::MatType
    inds::IndType
end
#-----------------------------------------------------------------------------------------------------
# Basic functions for type Operator
#-----------------------------------------------------------------------------------------------------
function *(
    c::Number, 
    o::Operator
)
    Operator(c * o.mat, o.inds)
end
#-----------------------------------------------------------------------------------------------------
function /( 
    o::Operator,
    c::Number
)
    Operator(o.mat / c, o.inds)
end
#-----------------------------------------------------------------------------------------------------
# Operator to fill vector(s)
#-----------------------------------------------------------------------------------------------------
function addtovec!(
    vec::AbstractVector, 
    o::Operator, 
    b::Basis,
    c::Number = 1
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
function addtovecs!(
    vecs::AbstractMatrix, 
    o::Operator, 
    b::Basis,
    c::AbstractVector{<:Number}
)
    vb = view(b, o.inds)
    oi = o.mat[:, index(vb)]
    for k = 1:size(o.mat, 1)
        change!(vb, k)
        i = index(b)
        vecs[i, :] += oi[k] * c
    end
end
#-----------------------------------------------------------------------------------------------------
# Type: Operation
#-----------------------------------------------------------------------------------------------------
export Operation
struct Operation{
    OptType <: AbstractVector{<:Operator}, 
    BasType <: Basis
}
    opts::OptType
    basis::BasType
end
#-----------------------------------------------------------------------------------------------------
# Operation initiate
#-----------------------------------------------------------------------------------------------------
export operation
function operation(
    mats::AbstractVector{<:AbstractMatrix},
    inds::AbstractVector{<:AbstractVector},
    len::Integer;
    base::Int64 = 0
)
    if base == 0
        base = Int(size(mats[1], 1)^(1/length(inds[1])))
    end
    b = basis(base, len)
    opts = [Operator(mats[i],inds[i]) for i=1:length(mats)]
    Operation(opts, b)
end
#-----------------------------------------------------------------------------------------------------
export onsite_operation
function onsite_operation(
    mats::AbstractVector{<:AbstractMatrix},
    ind::AbstractVector{<:Integer},
    len::Integer
)
    base = size(mats[1], 1)
    inds = [[i] for i in ind]
    operation(mats, inds, len, base=base)
end
#-----------------------------------------------------------------------------------------------------
function onsite_operation(
    mats::AbstractVector{<:AbstractMatrix}
)
    base = size(mats[1], 1)
    len = length(mats)
    inds = [[i] for i in 1:len]
    operation(mats, inds, len, base=base)
end
#-----------------------------------------------------------------------------------------------------
function onsite_operation(
    mat::AbstractMatrix,
    len::Integer
)
    mats = fill(mat, len)
    base = size(mat, 1)
    inds = [[i] for i in 1:len]
    operation(mats, inds, len, base=base)
end
#-----------------------------------------------------------------------------------------------------
export trans_inv_operation
function trans_inv_operation(
    mat::AbstractMatrix,
    ind::AbstractVector{<:Integer},
    len::Integer
)
    mats = fill(mat, len)
    inds = [mod.(ind .+ (i-1), len) .+ 1 for i=0:len-1]
    operation(mats, inds, len)
end
#-----------------------------------------------------------------------------------------------------
# Basic Functions for type Operation
#-----------------------------------------------------------------------------------------------------
function *(
    c::Number, 
    o::Operation
)
    Operation(c .* o.opts, o.basis)
end
#-----------------------------------------------------------------------------------------------------
function /( 
    o::Operation,
    c::Number
)
    Operation(o.opts ./ c, o.basis)
end
#-----------------------------------------------------------------------------------------------------
function +(
    o1::Operation, 
    o2::Operation
)
    Operation(vcat(o1.opts, o2.opts), o1.basis)
end
#-----------------------------------------------------------------------------------------------------
function -(
    o1::Operation, 
    o2::Operation
)
    Operation(vcat(o1.opts, (-1) * o2.opts), o1.basis)
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
# Apply to vector(s)
#-----------------------------------------------------------------------------------------------------
function addtovec!(
    vec::AbstractVector, 
    op::Operation, 
    j::Integer,
    c::Number=1
)
    len = length(op.opts)
    b = op.basis
    ol = op.opts
    for i = 1:len
        change!(b, j)
        addtovec!(vec, ol[i], b, c)
    end
end
#-----------------------------------------------------------------------------------------------------
function addtovecs!(
    vecs::AbstractMatrix, 
    op::Operation, 
    j::Integer,
    c::AbstractVector{<:Number}
)
    len = length(op.opts)
    b = op.basis
    ol = op.opts
    for i = 1:len
        change!(b, j)
        addtovecs!(vecs, ol[i], b, c)
    end
end
#-----------------------------------------------------------------------------------------------------
export mul!
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
            addtovec!(vec, ol[i], b, state[j])
        end
    end
end
#-----------------------------------------------------------------------------------------------------
function mul!(
    mat::AbstractMatrix, 
    op::Operation, 
    states::AbstractMatrix
)
    len = length(op.opts)
    b = op.basis
    ol = op.opts
    for j = 1:size(states, 1)
        for i = 1:len
            change!(b, j)
            addtovecs!(mat, ol[i], b, states[j, :])
        end
    end
end
#-----------------------------------------------------------------------------------------------------
# General multiplication
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
        addtovec!(view(mat, :, j), op, j)
    end
end
