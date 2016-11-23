function nk = polyDim(k, dim)
%   Calculates the dimension of the space of polynomials of degree less
%   than or equal to k in dim dimensions.

    if k == -1
        nk = 0;
    else
    nk = nchoosek(k+dim,k);
    end
    
end