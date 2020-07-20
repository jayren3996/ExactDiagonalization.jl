module ExactDiagonalization

export Operator, basis, fillmat!, chain
export spinop, spinop!

using LinearAlgebra
using SparseArrays
import Base: view, copy
include("Basis.jl")
include("Spin.jl")

end # module
