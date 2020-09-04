export covmat
function covmat(ol::Vector{<:Operation}, v::AbstractVector)
    n = length(ol)
    cm = Array{ComplexF64}(undef, n, n)
    for i=1:n 
        v1 = ol[i] * v
        for j=i:n 
            v2 = ol[j] * v
            cm[i,j] = dot(v1, v2) - dot(v, v1) * dot(v, v2)
        end
    end
    Hermitian(cm)
end

function covmat(ol::Vector{<:Operation}, v::AbstractMatrix)
    n = length(ol)
    cm = zeros(ComplexF64, n, n)
    num = size(v, 2)
    for i=1:n 
        v1 = ol[i] * v
        for j=i:n 
            v2 = ol[j] * v
            for k=1:num 
                cm[i,j] += dot(v1[:,k], v2[:,k]) - dot(v[:,k], v1[:,k]) * dot(v[:,k], v2[:,k])
            end
        end
    end
    Hermitian(cm)
end