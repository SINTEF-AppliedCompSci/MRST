function t = chebyshevPolynomials(degree)
    
    n = polyDim(degree, 1);
    t = cell(n,1);
    
    t{1} = Polynomial(0, 1);
    x    = Polynomial(1, 1);
    if degree > 0
        t{2} = Polynomial(1,1);
        for k = 1:n-2
            t{k+2} = 2*x*t{k+1} - t{k};
        end
    end
    
end