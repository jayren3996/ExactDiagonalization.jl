using LinearAlgebra
#--- Import test functions
include("../src/ExactDiagonalization.jl")
import .ExactDiagonalization: spin, Operation, fillmat!, operation, covmat
l = 10
ol = Vector{Operation}(undef, 12*l)
rl = randn(12*l)
for i=1:l
    j = (i-1)*12
    ol[j+ 1] = operation([spin("1x",2)],[[i,i%l+1]],2,10)
    ol[j+ 2] = operation([spin("1y",2)],[[i,i%l+1]],2,10)
    ol[j+ 3] = operation([spin("1z",2)],[[i,i%l+1]],2,10)
    ol[j+ 4] = operation([spin("xx",2)],[[i,i%l+1]],2,10)
    ol[j+ 5] = operation([spin("xy",2)],[[i,i%l+1]],2,10)
    ol[j+ 6] = operation([spin("xz",2)],[[i,i%l+1]],2,10)
    ol[j+ 7] = operation([spin("yx",2)],[[i,i%l+1]],2,10)
    ol[j+ 8] = operation([spin("yy",2)],[[i,i%l+1]],2,10)
    ol[j+ 9] = operation([spin("yz",2)],[[i,i%l+1]],2,10)
    ol[j+10] = operation([spin("zx",2)],[[i,i%l+1]],2,10)
    ol[j+11] = operation([spin("zy",2)],[[i,i%l+1]],2,10)
    ol[j+12] = operation([spin("zz",2)],[[i,i%l+1]],2,10)
end
function ham(cl,ol)
    nol = .*(cl,ol)
    M = zeros(ComplexF64, 2^l,2^l)
    for i=1:12 * l
        fillmat!(M, nol[i])
    end
    M
end
H0 = ham(rl, ol)
es,vs = eigen(Hermitian(H0))
num = 2^(l-1)
v = vs[:, num]
cm = covmat(ol, v)
ce,cv = eigen(cm)
println(ce[1:2])
cv1 = cv[:,1]
rt = rl ./ cv1
err = sum(abs2.(rt .- sum(rt)/length(rt)))
println(err)
