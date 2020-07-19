#--- product state basis
struct Basis{T<:AbstractVector{<:Integer}}
    vec::T
    d::Int64
    l::Int64
end
basis(l::Int64, d::Int64) = Basis(zeros(Int64, l),d,l)

list2num(v::AbstractVector,d::Integer,l::Integer)=sum(v[i]*d^(l-i) for i=1:l)+1
function num2list!(v::AbstractVector,i::Integer,d::Integer,l::Integer)
    i -= 1
    for j=1:l
        n,i = divrem(i, d^(l-j))
        v[j] = n
    end
end

index(s::Basis) = list2num(s.vec, s.d, s.l)
change!(s::Basis, i::Integer) = num2list!(s.vec, i, s.d, s.l)
change!(s::Basis, vec::AbstractVector{<:Integer}) = (s.vec .= vec)
copy(s::Basis) = Basis(copy(s.vec), s.d, s.l)
view(s::Basis, inds) = Basis(view(s.vec, inds), s.d, length(inds))

#--- Filling matrix
struct Operator{T<:AbstractMatrix}
    mat::T
    inds::AbstractVector{Int64}
end

function fillvec!(mj::AbstractVector, s::Basis, o::Operator)
    vs = view(s, o.inds)
    hi = o.mat[index(vs),:]
    for k = 1:size(o.mat,1)
        change!(vs,k)
        i = index(s)
        mj[i] += hi[k]
    end
end

function fillvec!(mj::AbstractVector, s::Basis, ol::AbstractVector{<:Operator})
    for i = 1:length(ol)
        cs = copy(s)
        fillvec!(mj,cs,ol[i])
    end
end

function fillmat!(m::AbstractMatrix, s::Basis, ol::AbstractVector{<:Operator})
    for j = 1:size(m,1)
        change!(s,j)
        fillvec!(view(m,:,j),s,ol)
    end
end
