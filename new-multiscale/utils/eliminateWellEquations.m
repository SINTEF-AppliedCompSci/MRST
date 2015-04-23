function [A_pp, q_p, A_ww, A_wp, q_w] = eliminateWellEquations(A, q, nc)
    iscell = false(size(q));
    iscell(1:nc) = true;
    
    A_pp = A(iscell, iscell);
    A_pw = A(iscell, ~iscell);
    q_p = q(iscell);
    
    A_ww = A(~iscell, ~iscell);
    A_wp = A(~iscell, iscell);
    q_w = q(~iscell);
    
    A_pp = A_pp - A_pw*(A_ww\A_wp);
    q_p = q_p - A_pw*(A_ww\q_w);
end
