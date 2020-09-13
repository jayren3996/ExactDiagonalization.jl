module ExactDiagonalization

using LinearAlgebra
using SparseArrays
import Base: view, copy
import LinearAlgebra: mul!
import Base: sum
import Base.:+
import Base.:-
import Base.:*
import Base.:/


include("Basis.jl")
include("Operation.jl")
include("Measure.jl")
include("Spin.jl")
#-----------------------------------------------------------------------------------------------------
# Krylov space
#-----------------------------------------------------------------------------------------------------
function reducespace(
    vs::AbstractMatrix; 
    rtol::Real=1e-3
)
    mat = vs' * vs
    e,v = eigen(Hermitian(mat))
    pos = e .> rtol
    ep = e[pos]
    vp = v[:,pos]
    nvs = vs * vp
    for i=1:length(ep)
        nvs[:,i] ./= sqrt(ep[i])
    end
    nvs
end
#-----------------------------------------------------------------------------------------------------
function reducespace(
    v::AbstractVector; 
    rtol::Real=1e-3
)
    vs = reshape(v,:,1)
    reducespace(vs, rtol=rtol)
end
#-----------------------------------------------------------------------------------------------------
export krylovspace
function krylovspace(
    H, 
    vs::AbstractMatrix; 
    rtol::Real=1e-3
)
    r = rank(vs)
    v = reducespace(hcat(vs, H * vs), rtol=rtol)
    nr = size(v,2)
    while nr > r
        r = nr
        v = reducespace(hcat(v, H * v), rtol=rtol)
        nr = size(v,2)
    end
    v
end
#-----------------------------------------------------------------------------------------------------
function krylovspace(
    H, 
    v::AbstractVector; 
    rtol::Real=1e-3
)
    vs = reshape(v,:,1)
    krylovspace(H, vs, rtol=rtol)
end


end # module ExactDiagonalization
