using LinearAlgebra
using SparseArrays
include("../src/ExactDiagonalization.jl")
include("deprecate/ExactDiagonalization.jl")
using .ExactDiagonalization
import .OldED as OED
using Test
#-------------------------------------------------------------------------------------------------------------------------
# Pauli matrices
#-------------------------------------------------------------------------------------------------------------------------
const X = [0 1; 1 0]
const Y = [0 -1; 1 0] * 1im
const Z = [1 0; 0 -1]
@testset "PauliKron" begin
    B = tensorbasis(3, base=2)
    XII = operator([X], [[1]], B) |> Array
    @test XII ≈ kron(X, I(4))
    IYI = operator([Y], [[2]], B) |> Array
    @test IYI ≈ kron(I(2), Y, I(2))
    IIZ = operator([Z], [[3]], B) |> Array
    @test IIZ ≈ kron(I(4), Z)
end

#-------------------------------------------------------------------------------------------------------------------------
# Randon Matrices
#-------------------------------------------------------------------------------------------------------------------------
@testset "Real Random Matrices" begin
    L = 12
    println("\nL = $L System:")
    mats = [rand(4, 4) |> Hermitian |> Array for i=1:L]
    inds = [mod.([i-1, i], L) .+ 1 for i=1:L]
    B = tensorbasis(L, base=2)
    print("EDKit                :")
    @time op1 = operator(mats, inds, B) |> Array
    print("ExactDiagonalization :")
    @time op2 = OED.operation(mats, inds, L) |> Array
    @test op1 ≈ op2
end

@testset "ComplexF64 Random Matrices" begin
    L = 13
    println("\nL = $L System:")
    mats = [rand(ComplexF64, 4, 4) |> Hermitian |> Array for i=1:L]
    inds = [mod.([i-1, i], L) .+ 1 for i=1:L]
    B = tensorbasis(L, base=2)
    print("EDKit                :")
    @time op1 = operator(mats, inds, B) |> Array
    print("ExactDiagonalization :")
    @time op2 = OED.operation(mats, inds, L) |> Array
    @test op1 ≈ op2
end

@testset "Benchmark Sparse Matrix" begin
    L = 15
    println("\nL = $L System:")
    mat = kron(X, Y, Z)
    B = tensorbasis(L, base=2)
    print("EDKit                :")
    @time op1 = trans_inv_operator(mat, 3, B) |> Array
    print("ExactDiagonalization :")
    @time op2 = OED.trans_inv_operation(mat, 1:3, L) |> Array
end
