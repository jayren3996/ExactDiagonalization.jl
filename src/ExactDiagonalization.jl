module ExactDiagonalization

export Operation, operation, fillmat!, mul, mul!, spin

using LinearAlgebra
using SparseArrays
import Base: view, copy
import Base.:*
import Base.:+
import Base: sum
import LinearAlgebra: mul!
include("Basis.jl")
#include("Operator.jl")
include("Operation.jl")
include("Correlation.jl")
include("Spin.jl")

end # module
