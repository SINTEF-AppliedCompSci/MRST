function res = W(v) 
% well function    
    res = -0.5772 - log(v) + v;
    for i = 2:11
        res = res - (-1)^i * (v.^i)/(i*factorial(i));
    end
end
    
