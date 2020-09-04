export covmat
function covmat(ol::Vector{<:Operation}, v::AbstractVector)
    n = length(ol)
    #cm = Array{ComplexF64}(undef, n, n)
    cm = Array{Float64}(undef, n, n)
    for i=1:n 
        v1 = ol[i] * v
        for j=i:n 
            v2 = ol[j] * v
            #cm[i,j] = dot(v1, v2) - dot(v, v1) * dot(v, v2)
            cm[i,j] = real(dot(v1, v2) - dot(v, v1) * dot(v, v2))
        end
    end
    Hermitian(cm)
end

function covmat(ol::Vector{<:Operation}, v::AbstractMatrix)
    n = length(ol)
    cm = zeros(Float64, n, n)
    num = size(v, 2)
    for i=1:n 
        v1 = ol[i] * v
        for j=i:n 
            v2 = ol[j] * v
            ele = real(sum(conj.(v1) .* v2))
            ele -= real( sum(conj.(v) .* v1) * sum(conj.(v) .* v2) )
            cm[i,j] = ele
        end
    end
    Hermitian(cm) ./ num
end