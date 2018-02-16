function l = legendrePolynomials(degree)
    
    n = polyDim(degree, 1);
    l = cell(n,1);
    
    l{1} = Polynomial(0, 1);
    if degree > 0
        l{2} = Polynomial(1,1);
        for k = 1:n-2
            l{k+2} = ((2*k+1)*l{2}*l{k+1} - k*l{k})./(k+1);
        end
    end
    
end