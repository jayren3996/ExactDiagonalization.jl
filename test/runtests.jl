using Test
using LinearAlgebra
using SparseArrays
#--- Import test functions
include("../src/ExactDiagonalization.jl")
import .ExactDiagonalization: basis, index, change!, copy, view
import .ExactDiagonalization: Operator, operation, fillvec!, fillmat!, mul!, mul
import .ExactDiagonalization: spin
import .ExactDiagonalization: Operation, covmat
#--- Test Basis
@testset "Basis" begin
    b1 = basis(3, 4)
    b2 = basis([1,2,0,1], 3)
    # test index:
    @test index(b2) == 47
    # test function change! with integer:
    change!(b1, 47)
    @test b1.bits == b2.bits
    # test copy and view:
    b3 = copy(b1)
    b4 = view(b3, [2,3])
    @test b4.len == 2
    # test function change! with vector:
    change!(b4, [0,1])
    @test index(b4) == 2
    @test index(b3) == 32
    @test index(b1) == 47
end
#--- Test Operators
@testset "Operator" begin
    m1 = rand(3,3)
    m2 = rand(3,3)
    mc = rand(ComplexF32,3,3)
    k13 = kron(m1, I(3), m2, I(3))       # answer
    ktt = k13 + kron(I(3), m1, m2, I(3)) # answer
    ktc = k13 + kron(I(3), m2, mc, I(3)) # answer
    op13 = Operator(kron(m1,m2), [1,3])
    ot13 = operation([kron(m1,m2)],[[1,3]], 3, 4)
    ott = operation([kron(m1,m2), kron(m1,m2)],[[1,3],[2,3]], 3, 4)
    otc = operation([kron(m1,m2), kron(m2,mc)],[[1,3],[2,3]], 3, 4)
    # test fillvec! for o13
    vec = zeros(81)
    b = basis(3, 4)
    ind = rand(1:81)
    change!(b, ind)
    fillvec!(vec, op13, b)
    @test vec ≈ k13[:, ind] atol=1e-5
    # test fillvec! for o13
    vec = zeros(81)
    ind = rand(1:81)
    fillvec!(vec,ot13,ind)
    @test vec ≈ k13[:, ind] atol=1e-5
    # test fillvec! for ott
    vec = zeros(81)
    ind = rand(1:81)
    fillvec!(vec,ott,ind)
    @test vec ≈ ktt[:, ind] atol=1e-5
    # test fillmat! for ott
    mat = zeros(81,81)
    fillmat!(mat,ott)
    @test mat ≈ ktt atol=1e-5
    # test mul! for vec
    vec = rand(81)
    zvec = zeros(81)
    mul!(zvec, ott, vec)
    @test zvec ≈ ktt * vec atol=1e-5
    # test mul! for mat
    mat = rand(81,3)
    zmat = zeros(81,3)
    mul!(zmat, ott, mat)
    @test zmat ≈ ktt * mat atol=1e-5
    # test mul
    vec = rand(81)
    mat = rand(81,3)
    @test mul(ott, vec) ≈ ktt * vec atol=1e-5
    @test mul(ott, mat) ≈ ktt * mat atol=1e-5
    @test mul(otc, mat) ≈ ktc * mat atol=1e-5
    # test *
    @test ott * mat ≈ ktt * mat atol=1e-5
    @test otc * mat ≈ ktc * mat atol=1e-5
end

@testset "Spin" begin
    sx = spin('x', 3)
    sy = spin('Y', 3)
    sz = spin('z', 3)
    # check single spin operators
    @test sx ≈ [0 sqrt(2) 0;sqrt(2) 0 sqrt(2);0 sqrt(2) 0] ./ 2
    @test sy ≈ [0 sqrt(2) 0;-sqrt(2) 0 sqrt(2);0 -sqrt(2) 0] ./ 2
    @test sz ≈ Diagonal([1,0,-1])
    # 3-site spin
    mx = kron(sx,I(9)) + kron(I(3), sx, I(3)) + kron(I(9), sx)
    my = kron(sy,I(9)) + kron(I(3), sy, I(3)) + kron(I(9), sy)
    mz = kron(sz,I(9)) + kron(I(3), sz, I(3)) + kron(I(9), sz)
    opx = spin('x', 3, 3)
    opy = spin('Y', 3, 3)
    opz = spin('z', 3, 3)
    ex = Diagonal(ones(27))
    ey = Diagonal(ones(27))
    ez = Diagonal(ones(27))
    # test mul
    rmat = rand(ComplexF32, 27,27)
    @test mx ≈ mul(opx,ex) atol=1e-5
    @test my ≈ mul(opy,ey) atol=1e-5
    @test mz ≈ mul(opz,ez) atol=1e-5
    @test mx*rmat ≈ mul(opx,rmat) atol=1e-5
    @test my*rmat ≈ mul(opy,rmat) atol=1e-5
    @test mz*rmat ≈ mul(opz,rmat) atol=1e-5
    # test mulmul
    @test mx^2 ≈ mul(opx,mul(opx,ex)) atol=1e-5
    @test my^2 ≈ mul(opy,mul(opy,ey)) atol=1e-5
    @test mz^2 ≈ mul(opz,mul(opz,ez)) atol=1e-5
    @test mx^2*rmat ≈ mul(opx,mul(opx,rmat)) atol=1e-5
    @test my^2*rmat ≈ mul(opy,mul(opy,rmat)) atol=1e-5
    @test mz^2*rmat ≈ mul(opz,mul(opz,rmat)) atol=1e-5
    # test xyz
    m1 = kron(sx,sy,sz)
    m2 = spin("xYz",3)
    @test m1 ≈ m2 atol=1e-5
end

@testset "correlation" begin
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
    @test ce[1] < 1e-7
    @test ce[2] > 1e-3
    cv1 = cv[:,1]
    rt = rl ./ cv1
    err = sum(abs2.(rt .- sum(rt)/length(rt)))
    @test err < 1e-7
end