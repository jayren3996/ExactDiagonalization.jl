#--- Basic
coefficient(D::Integer) = [sqrt(i*(D-i)) for i = 1:D-1]
Sp(D::Integer) = sparse(1:D-1, 2:D,coefficient(D), D,D)
Sm(D::Integer) = Sp(D)'
Sx(D::Integer) = (sp=Sp(D); (sp+sp')/2)
iSy(D::Integer) = (sp=Sp(D); (sp-sp')/2)
Sy(D::Integer) = (sp=Sp(D); 0.5im*(sp'-sp))
Sz(D::Integer) = (J=(D-1)/2; sparse(1:D,1:D, J:-1:-J))
function spinop(s::Char,D::Integer)
    if s == '+' return Sp(D)
    elseif s == '-' return Sm(D)
    elseif s == 'x' return Sx(D)
    elseif s == 'Y' return iSy(D)
    elseif s == 'y' return Sy(D)
    elseif s == 'z' return Sz(D)
    elseif s == '1' return sparse(I,D,D)
    end
end
#--- large
spinop(s::String, D::Integer) = kron([spinop(si,D) for si in s]...)
function spinop!(mat::AbstractMatrix,
                 s::Char,
                 D::Integer,
                 L::Integer)
    m = spinop(s,D)
    #println(typeof(m))
    op = chain(m,1,L)
    b = basis(L,D)
    fillmat!(mat,b,op)
end
