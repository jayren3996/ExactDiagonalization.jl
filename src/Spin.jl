#-----------------------------------------------------------------------------------------------------
# Basic spin operators
#-----------------------------------------------------------------------------------------------------
coefficient(D::Integer) = [sqrt(i*(D-i)) for i = 1:D-1]
Sp(D::Integer) = sparse(1:D-1, 2:D,coefficient(D), D,D)
Sm(D::Integer) = Sp(D)'
Sx(D::Integer) = (sp=Sp(D); (sp+sp')/2)
iSy(D::Integer) = (sp=Sp(D); (sp-sp')/2)
Sy(D::Integer) = (sp=Sp(D); 0.5im*(sp'-sp))
Sz(D::Integer) = (J=(D-1)/2; sparse(1:D,1:D, J:-1:-J))
#-----------------------------------------------------------------------------------------------------
# Combination
#-----------------------------------------------------------------------------------------------------
export spin
function spin(
    s::Char, 
    D::Integer
)
    if s == '+' return Sp(D)
    elseif s == '-' return Sm(D)
    elseif s == 'x' return Sx(D)
    elseif s == 'Y' return iSy(D)
    elseif s == 'y' return Sy(D)
    elseif s == 'z' return Sz(D)
    elseif s == '1' return sparse(I,D,D)
    end
end
#-----------------------------------------------------------------------------------------------------
function spin(
    s::Char, 
    D::Integer, 
    L::Integer
)
    sop = spin(s, D)
    operation(fill(sop, L), 1:L, D, L)
end
#-----------------------------------------------------------------------------------------------------
function spin(
    s::String, 
    D::Integer
) 
    kron([spin(si,D) for si in s]...)
end
