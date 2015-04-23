function p = recoverWellSolution(A_ww, A_wp, q_w, p)
    w = A_ww\(q_w - A_wp*p);
    p = [p; w];
end
