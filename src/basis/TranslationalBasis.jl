#-----------------------------------------------------------------------------------------------------
# Translational Basis
#-----------------------------------------------------------------------------------------------------
struct TranslationalBasis <: AbstractBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    C::ComplexF64
    B::Int
end
#-----------------------------------------------------------------------------------------------------
function change!(b::TranslationalBasis, i::Integer)::Float64
    change!(b.dgt, i, base=b.B)
    b.R[i]
end
#-----------------------------------------------------------------------------------------------------
function index(b::TranslationalBasis)::Tuple{ComplexF64, Int}
    Im, T = indmin(b.dgt, b.B)
    N = b.C ^ T * b.R[Im]
    ind = binary_search(b.I, Im)
    N, ind
end
#-----------------------------------------------------------------------------------------------------
size(b::TranslationalBasis, i::Integer) = (i == 2 || i == 1) ? length(b.I) : 1
size(b::TranslationalBasis) = (l=length(b.I); (l,l))
#-----------------------------------------------------------------------------------------------------
movebits!(dgt::AbstractVector{<:Integer}) = vcat(dgt[1:1], dgt[2:end])
function indmin(dgt::AbstractVector{<:Integer}, base::Integer)
    Im, T = index(dgt, base=base), 0
    for i=1:dgt.L-1
        dgt = movebits!(dgt)
        In = index(dgt, base=base)
        if In < Im
            Im, T = In, i
        end
    end
    Im, T
end
