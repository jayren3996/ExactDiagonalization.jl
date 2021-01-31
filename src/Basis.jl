#-----------------------------------------------------------------------------------------------------
# Functions for "Basis"

# The "Basis" object is used to label product states:

# 1. A single digit starts from 0, and ends with (base-1)
# 2. The index number starts from 1, and ends with (base^len)
# 3. Digits are read from left to right, i.e., the left digits have greater weights
#-----------------------------------------------------------------------------------------------------
"""
    Basis{T<:AbstractVector{<:Integer}} 

Type that contains 3 fields:
1. bits: vector of integer represent the product states
2. base: dimension of local Hilbert space
3. len : size of the product state
"""
struct Basis{
    T<:AbstractVector{<:Integer}
}
    bits::T
    base::Int64
    len::Int64
end
#-----------------------------------------------------------------------------------------------------
# Initiations of "Basis":

# 1. Initiate with zero-bits.
# 2. Initiate with a given bits.
#-----------------------------------------------------------------------------------------------------
# Initiate zero basis
function basis(
    base::Integer, 
    len::Integer
)
    return Basis(zeros(Int64, len), Int64(base), Int64(len))
end

# Initiate with given bits
function basis(
    bits::AbstractVector{<:Integer}, 
    base::Integer
)
    return Basis(bits, Int64(base), length(bits))
end
#-----------------------------------------------------------------------------------------------------
# Convertion between index numbers and bits lists:

# 1. List of bits -> index number
# 2. Index number -> list of bits
#-----------------------------------------------------------------------------------------------------
# List of bits -> index number
function list2num(
    bits::AbstractVector{<:Integer}, 
    base::Integer, 
    len::Integer
)
    # Left ⟶ right
    # 1st digit : base^(len-1)
    # ith digit: base^(len - i)
    # last digit: 1
    # Since (0...0) ⟶ 1, we add 1 to the result
    num = bits[len]
    basen = 1
    for i = 1:len-1
        basen *= base
        num += bits[len-i] * basen
    end
    return num + 1
end

# Index number -> list of bits
function num2list!(
    bits::AbstractVector{<:Integer},
    index::Integer,
    base::Integer,
    len::Integer,
)
    # Since 1 ⟶ (0...0), we substract 1 from the input
    # The result is remaining_mumber
    remaining_number = index - 1
    # Iterate from left ⟶ right, labeled by i
    for i = 1:len
        # remaining_mumber = ith_digit * base^(len-i) + next_remaining_number
        # Use divrem function to get ith_digit & next_remaining_number
        ith_digit, remaining_number = divrem(remaining_number, base^(len - i))
        bits[i] = ith_digit
    end
    return nothing
end
#-----------------------------------------------------------------------------------------------------
# Functions defined on "Basis":

# 1. index  : Get the index of a given "Basis"
# 2. change!: Change the digits of a "Basis", given:
#    a. list of digits
#    b. index number
# 3. copy   : Copy a "Basis"
# 4. view   : Get a view on part of the "Basis", the view is also a "Basis".
#-----------------------------------------------------------------------------------------------------
function index(
    b::Basis
)
    return list2num(b.bits, b.base, b.len)
end

function change!(
    b::Basis, 
    index::Integer
)
    num2list!(b.bits, index, b.base, b.len)
    return nothing
end

function change!(
    b::Basis, 
    bits::AbstractVector{<:Integer}
)
    b.bits .= bits
    return nothing
end

function copy(
    b::Basis
)
    return Basis(copy(b.bits), b.base, b.len)
end

function view(
    b::Basis, 
    inds::AbstractVector
)
    basis_view = Basis(view(b.bits, inds), b.base, length(inds))
    return basis_view
end
