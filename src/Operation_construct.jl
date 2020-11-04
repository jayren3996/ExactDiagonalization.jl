#-----------------------------------------------------------------------------------------------------
# More construction method for Operation
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
