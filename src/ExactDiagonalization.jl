module ExactDiagonalization

export operation, fillmat!, mul, mul!, spin

using LinearAlgebra
using SparseArrays
import Base: view, copy
import LinearAlgebra: mul!
include("Basis.jl")
include("Operator.jl")
include("Spin.jl")

end # module
