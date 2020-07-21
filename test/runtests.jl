using Test
using LinearAlgebra
using SparseArrays
#--- Import test functions
include("../src/ExactDiagonalization.jl")
import .ExactDiagonalization: basis, index, change!, copy, view
import .ExactDiagonalization: Operator, operation, fillvec!, fillmat!, mul!, mul
import .ExactDiagonalization: spin
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
    k13 = kron(m1, I(3), m2, I(3))       # answer
    ktt = k13 + kron(I(3), m1, m2, I(3)) # answer
    op13 = Operator(kron(m1,m2), [1,3])
    ot13 = operation([kron(m1,m2)],[[1,3]], 3, 4)
    ott = operation([kron(m1,m2), kron(m1,m2)],[[1,3],[2,3]], 3, 4)
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
end

@testset "Spin" begin
    sx = spin('x', 3)
    sy = spin('Y', 3)
    sz = spin('z', 3)
    # test xxx
    m1 = kron(sx,I(9)) + kron(I(3), sx, I(3)) + kron(I(9), sx)
    op = spin('x', 3, 3)
    m2 = zeros(27,27)
    fillmat!(m2, op)
    @test m1 ≈ m2 atol=1e-5
    # test xyz
    m1 = kron(sx,sy,sz)
    m2 = spin("xYz",3)
    @test m1 ≈ m2 atol=1e-5
end
