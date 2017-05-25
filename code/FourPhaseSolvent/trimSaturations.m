function s = trimSaturations(s)

%     tol = 0;
%     s = s - s.*(s < tol);
%     
    tol = 1e-15;
    s = max(s, tol);
    
    
end