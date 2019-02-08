function ndof = polyDim(k, dim)
% Computes the dimension of the space of polynomials of degree k or less
    if k < 0
        ndof = 0;
    else
        ndof = nchoosek(k+dim, k);
    end
end