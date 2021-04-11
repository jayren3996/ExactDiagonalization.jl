# ExactDiagonalization.jl

 Julia package for exact diagonalization.

## Installation

Run the following script in the ```Pkg REPL``` :

```julia
pkg> add https://github.com/jayren3996/ExactDiagonalization.jl
```

## Representation of basis

In this package, the many-body basis is the tensor-product basis, which can be represented by a vector of integer. We define type `Basis` :

```julia
struct Basis{
    T<:AbstractVector{<:Integer}
}
    bits::T
    base::Int64
    len::Int64
end
```

In the definition, the field `bits` is the vector that represents the basis state, the field `base` is the nunmber of the local degrees of freedom (for spin-1/2, `base=2` , for spin-1, `base=3` , etc.), and field `len` is the size of the basis state.

The `Basis` type is an inner type of the package. In most of the case we do not need to work with it directly.

## Representation of operation

In `ExactDiagonalization` , a many-body operator is represented by the type `Operator`:

```julia
struct Operator{T1<:Number, T2<:Integer}
    mat::Matrix{T1}
    inds::Vector{T2}
end
```

In this definition,  the field `mat` is the matrix representation of the local operator, and the field `inds` is the indices of sites it acts on. In many case, the operator we care about is the sum of the local operator. To better deal with the sum of operator, we define another type `Operation`:

```julia
struct Operation{T1 <: Number, T2 <: Integer, T3}
    opts::Vector{Operator{T1, T2}}
    basis::Basis{T3}
end
```

In this definition, the field `opts` is the list of `Operator` that represent the summation of local operators, the field `basis` is the `Basis` object.

## Initialization of operation

The general initialization for an `Operation` object is by using the function `operation`:

```julia
function operation(
    mats::AbstractVector{<:AbstractMatrix},
    inds::AbstractVector{<:AbstractVector},
    len::Integer=0;
    base::Int64=0
)
```

The basis input is the list of matrices `mat` (which represent the local operators), the list of indices `inds` (which indicate the indices of site on which the operators act), and the size of the sistem `len`. The optional input is the `base` number.

For on-site operation (which means each operator acts on a single site), there is a special initialization function:

```julia
function onsite_operation(
    mats::AbstractVector{<:AbstractMatrix},
    ind::AbstractVector{<:Integer},
    len::Integer
)
```

This function is basically the same as `operation` but the input `ind` is now a list of integers, rather than list of vectors. If the on-site operation acts on every site, there is another method:

```julia
function onsite_operation(
    mats::AbstractVector{<:AbstractMatrix}
)
```

If the on-site operator are translational invariant, we can use the method:

```julia
function onsite_operation(
    mat::AbstractMatrix,
    len::Integer
)
```

which takes a single matrix, rather than a list of matrices.

Also, for a general translational invariant operation, there is another special initialization function:

```julia
function trans_inv_operation(
    mat::AbstractMatrix,
    ind::AbstractVector{<:Integer},
    len::Integer
)
```

The input `ind` is the indices of the first operator. 

## Functions for Operations

### Convert to matrix

In `ExactDiagonalization` , an `Operation` object is basically like a matrix. If we want a full matrix representation of the many-body operator. We first initiate an zero matrix with the correct dimension, and then use the function:

```julia
function fillmat!(
    mat::AbstractMatrix, 
    op::Operation
)
```

The `fillmat!` function fill the matrix `mat` using `op` . Note that `fillmat!` is an add-to function, so the `mat` should be a zero matrix rather than a random one.

### Multiply to vector or matrix

We can directly using `*` to do the multiplycation. Or, we could use the function:

```julia
function mul!(mat::AbstractMatrix, opt::Operation, states::AbstractMatrix)
```

## Spin tools

We also provide a helper function to create spin-s operators (represented by matrices):

```julia
function spin(spins...; D::Integer=2)
```

In the definition, `spins` are arbituary number of tuples such as `(1.0, "xzx")`. The supported characters are `x`, `y`, `z`, `1`, `+`, `-`, `Y`. Note that `Y=iSʸ`. The other input `D` is the dimension of the matrix (`D = 2s+1`).

It return an operation which is the sum of the same spin operators.

## Measurement

The function `measure` returns 1-point or 2-point expectation value:

```julia
# return ⟨v|O₁|v⟩
function measure(
    o::Operation, 
    v::AbstractVecOrMat
)
# return  ⟨v|(O₁')O₂|v⟩
function measure(
    o1::Operation, 
    o2::Operation,
    v::AbstractVecOrMat
)
```

The function `covmat` returns the covariant matrix of a given list of (Hermitian) operations:

```julia
function covmat(
    ol::Vector{<:Operation}, 
    v::AbstractVecOrMat
)
```

In the definition, `v` could be a vector, or a set of orthogonal vecotrs represented by a matrix. For non-Hermitian operations, use the alternative function ```covmatc```

```julia
function covmatc(
    ol::Vector{<:Operation}, 
    v::AbstractVecOrMat
)
```

The function `entropy` returns the entanglement entropy of the given state:

```julia
function entropy(
    v::AbstractVector, 
    l::Integer,
    i::Integer
)
```

The input `v` is the many-body vecotr, `l` is the size of the system, and `i` is the index of the site at which the system is cut.

## Miscellaneous

`ExactDiagonalization` provides a function ```krylovspace``` for computing Krylov space from an operator and a linear space (represented by a vector or a matrix):

```julia
function krylovspace(
    H, 
    vs::AbstractVecOrMat; 
    rtol::Real=1e-3
)
```

