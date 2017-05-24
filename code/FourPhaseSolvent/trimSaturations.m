function s = trimSaturations(s)

    tol = 1e-15;
    s = s - s.*(s < tol);
    
    
end