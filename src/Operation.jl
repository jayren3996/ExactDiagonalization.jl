#-----------------------------------------------------------------------------------------------------
# Functions for Operations.

# Operation is the collection of "Operator"s, together with a "Basis".
# Functions on "Operation"s are design to mimic the behavior of matrices.
#-----------------------------------------------------------------------------------------------------
"""
    Operator{MatType<:AbstractMatrix, IndType<:AbstractVector{<:Integer}}

Type that contains 2 fields:
1. mat : Matrix representation of local operator.
2. inds: Index of sites the operator acts on.
"""
struct Operator{
    MatType <: AbstractMatrix, 
    IndType <: AbstractVector{<:Integer}
}
    mat::MatType
    inds::IndType
end
#-----------------------------------------------------------------------------------------------------
# Basic functions on "Operator":

# 1. *: Multiplication by a number.
# 2. /: Division by a number.
#-----------------------------------------------------------------------------------------------------
*(c::Number, o::Operator) = Operator(c * o.mat, o.inds)
/(o::Operator, c::Number) = Operator(o.mat / c, o.inds)
#-----------------------------------------------------------------------------------------------------
# Operator to fill vector/matrix.

# To finally get the multiplication function, here we focus on a single row manipulation:

# 1. Start with a given product state(digits), we cut a segment from it;
# 2. Iterate all possible digits in the segment, and read the matrix element of operator;
# 3. Add to the vector/matrix element.
#-----------------------------------------------------------------------------------------------------
"""
    addtovec!(vec::AbstractVector, opt::Operator, basis::Basis, coeff::Number=1)

Elementary multiplication step. Calculate single vector element of an operator * vector.

vec  : Vector to fill
opt  : Operator
basis: Indicate the specific row
coeff: The vector element of input state
"""
function addtovec!(
    vec::AbstractVector, 
    opt::Operator, 
    basis::Basis,
    coeff::Number=1
)
    # Create a view on the segment of the digits
    basis_view = view(basis, opt.inds)
    # Get the initial index of the viewed segment
    index_view = index(basis_view)
    # Get the elements in row of index_view
    opt_column = opt.mat[:, index_view] * coeff
    row_number = length(opt_column)
    # Fill the vector
    for k = 1:row_number
        change!(basis_view, k)
        i = index(basis)
        vec[i] += opt_column[k]
    end
    # Reset the segment
    change!(basis_view, index_view)
    return nothing
end
#-----------------------------------------------------------------------------------------------------
"""
    addtovecs!(vecs::AbstractMatrix, opt::Operator, basis::Basis, coeff::AbstractVector)

Elementary multiplication step. Calculate single row of matrix element of an operator * matrix.

vecs : Vectors to fill
opt  : Operator
basis: Indicate the specific row
coeff: The row of matrix elements of input states
"""
function addtovecs!(
    vecs::AbstractMatrix, 
    opt::Operator, 
    basis::Basis,
    coeff::AbstractVector{<:Number}
)
    # Create a view on the segment of the digits
    basis_view = view(basis, opt.inds)
    # Get the initial index of the viewed segment
    index_view = index(basis_view)
    # Get the elements in row of index_view
    opt_column = opt.mat[:, index_view]
    row_number = length(opt_column)
    # Fill the matrix
    for k = 1:row_number
        change!(basis_view, k)
        i = index(basis)
        vecs[i, :] += opt_column[k] * coeff
    end
    # Reset the segment
    change!(basis_view, index_view)
    return nothing
end
#-----------------------------------------------------------------------------------------------------
# Type: Operation
#-----------------------------------------------------------------------------------------------------
export Operation
"""
    Operation{OptType<:AbstractVector{<:Operator}, BasType<:Basis}

Type that contains 2 fields:
opts : List of operators
basis: Basis for the system
"""
struct Operation{
    OptType <: AbstractVector{<:Operator}, 
    BasType <: Basis
}
    opts::OptType
    basis::BasType
end
#-----------------------------------------------------------------------------------------------------
# Basis Operation initiation
#-----------------------------------------------------------------------------------------------------
export operation
function operation(
    mats::AbstractVector{<:AbstractMatrix},
    inds::AbstractVector{<:AbstractVector},
    len::Integer=0;
    base::Int64=0
)
    if base == 0
        base = Int(size(mats[1], 1)^(1/length(inds[1])))
    end
    if len == 0
        len = maximum(maximum.(inds))
    end
    b = basis(base, len)
    opts = [Operator(mats[i],inds[i]) for i=1:length(mats)]
    Operation(opts, b)
end
#-----------------------------------------------------------------------------------------------------
# Basic Functions for type Operation
#-----------------------------------------------------------------------------------------------------
*(c::Number, o::Operation) = Operation(c .* o.opts, o.basis)
/(o::Operation, c::Number) = Operation(o.opts ./ c, o.basis)
+(o1::Operation, o2::Operation) = Operation(vcat(o1.opts, o2.opts), o1.basis)
-(o1::Operation, o2::Operation) = Operation(vcat(o1.opts, (-1) * o2.opts), o1.basis)
sum(ol::AbstractVector{<:Operation}) = begin
    basis = ol[1].basis
    opts = vcat([oi.opts for oi in ol]...)
    Operation(opts, basis)
end
#-----------------------------------------------------------------------------------------------------
# Operation to fill vector/matrix:

# The method is basically iterate each operator to the given digits-represented basis
#-----------------------------------------------------------------------------------------------------
"""
    addtovec!(vec::AbstractVector, opt::Operation, coeff::Number=1)

Elementary multiplication step. Calculate single vector element of an operation * vector.
The row information is stored in the basis of operation.

vec  : Vector to fill
opt  : Operation
coeff: The vector element of input state
"""
function addtovec!(
    vec::AbstractVector, 
    opt::Operation, 
    coeff::Number=1
)
    # Get basis and operators
    basis = opt.basis
    opts = opt.opts
    num_of_opts = length(opts)
    # Multiply by each operator and add them all to vector
    for i = 1:num_of_opts
        addtovec!(vec, opts[i], basis, coeff)
    end
    return nothing
end
#-----------------------------------------------------------------------------------------------------
"""
    addtovecs!(vecs::AbstractMatrix, opt::Operation, coeff::AbstractVector)

Elementary multiplication step. Calculate single row of matrix element of an operation * matrix.
The row information is stored in the basis of operation.

vecs : Vectors to fill
opt  : Operation
coeff: The row of matrix elements of input states
"""
function addtovecs!(
    vecs::AbstractMatrix, 
    opt::Operation, 
    coeff::AbstractVector{<:Number}
)
    # Get basis and operators
    basis = opt.basis
    opts = opt.opts
    num_of_opts = length(opts)
    # Multiply by each operator and add them all to vectors
    for i = 1:num_of_opts
        addtovecs!(vecs, opts[i], basis, coeff)
    end
    return nothing
end
#-----------------------------------------------------------------------------------------------------
# Multiplication

# The idea is to iterate all product state basis, and get all the vector/matrix elements
#-----------------------------------------------------------------------------------------------------
export mul!
"""
    mul!(vec::AbstractVector, opt::Operation, state::AbstractVector)

Full multiplication for operation and vector.

vec  : Vector to fill 
opt  : Operation
state: Input vector
"""
function mul!(
    vec::AbstractVector, 
    opt::Operation, 
    state::AbstractVector
)
    # Get basis
    basis = opt.basis
    for j = 1:length(state)
        change!(basis, j)
        addtovec!(vec, opt, state[j])
    end
    return nothing
end
#-----------------------------------------------------------------------------------------------------
"""
    mul!(mat::AbstractVector, opt::Operation, state::AbstractVector)

Full multiplication for operation and matrix.

mat   : Matrix to fill 
opt   : Operation
states: Input states(matrix)
"""
function mul!(
    mat::AbstractMatrix, 
    opt::Operation, 
    states::AbstractMatrix
)
    # Get basis
    basis = opt.basis
    for j = 1:size(states, 1)
        change!(basis, j)
        addtovecs!(mat, opt, states[j, :])
    end
    return nothing
end
#-----------------------------------------------------------------------------------------------------
export mul
"""
    mul(opt::Operation, vec_or_mat::AbstractVecOrMat)

General multiplication for operation and vector/matrix.

opt       : Operation
vec_or_mat: Input vector/matrix
"""
function mul(
    opt::Operation, 
    vec_or_mat::AbstractVecOrMat{T}
) where T <: Number
    # Find the promoted type
    shared_type = promote_type([eltype(o.mat) for o in opt.opts]...)
    result_type = promote_type(T, shared_type)
    # Initiate zero vector/matrix
    out = zeros(result_type, size(vec_or_mat))
    mul!(out, opt, vec_or_mat)
    return out
end
*(opt::Operation, vom::AbstractVecOrMat) = mul(opt, vom)
#-----------------------------------------------------------------------------------------------------
# Fill matrix
#-----------------------------------------------------------------------------------------------------
export fillmat!
"""
    fillmat!(mat::AbstractMatrix, opt::Operation)

Fill the zero matrix with matrix element from operation.

mat : Matrix to fill
opt : Operation
"""
function fillmat!(
    mat::AbstractMatrix, 
    opt::Operation
)
    basis = opt.basis
    row_number = size(mat, 1)
    for j = 1:row_number
        change!(basis, j)
        addtovec!(view(mat, :, j), opt)
    end
    return nothing
end
