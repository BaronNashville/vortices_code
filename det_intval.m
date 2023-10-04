function d = det_intval(A)
%Returns the determinant of A using PLU decompositon but accepts interval inputs

N = size(A,1);
M = size(A,2);

if N ~= M
    d = NaN;
else
    [L,U,p] = luwpp(A);

    d = parity(p) * prod(diag(U));
end

end