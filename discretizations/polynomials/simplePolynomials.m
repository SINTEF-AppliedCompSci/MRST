function p = simplePolynomials(degree)
    
    n = polyDim(degree, 1);
    p = cell(n,1);
    
    for k = 1:n
        p{k} = Polynomial(k-1,1);
    end
    
end