module ExactDiagonalization

export Operator, basis, fillmat!
export spinop, spinop!

using LinearAlgebra
using SparseArrays
import Base: view, copy
include("Basis.jl")
include("Spin.jl")

end # module
