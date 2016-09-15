function state = calculateCellPressure(state, G, S)

    k = 1;
    
    alpha = [0 1 0 0]; beta = [0 0 1 0]; gamma = [0 0 0 1];
    
    S.PiNstar;

end
    
function coeff = trinomialExpansion(a, b, c, n)
    
    if n == 0
        alpha = 0; beta = 0; gamma = 0;
    elseif n == 1
        alpha = [1,0,0]; beta = [0,1,0]; gamma = [0,0,1];
    elseif n == 2
        alpha = [2,1,1,0,0,0]; beta = [0,1,0,2,1,0]; gamma = [0,0,1,0,1,2];
    else
        alpha = [3 2 2 1 1 1 0 0 0 0];
        beta  = [0 1 0 2 0 1 3 2 1 0];
        gamma = [0 0 1 0 2 1 0 1 2 3];
    end
    
    r = size(a,1);     
    coeff = repmat(factorial(n)./(factorial(alpha).*factorial(beta).*factorial(gamma)), r, 1);
    coeff = coeff.*bsxfun(@power, a,repmat(alpha,r,1))...
         .*bsxfun(@power, b,repmat(beta,r,1))...
         .*bsxfun(@power, c,repmat(gamma,r,1));
    
end