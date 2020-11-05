function [gangiFnOut] = gangiFn(p, Pc, alpha, Pmax, m, k0, matrix_perm)
%GANGIFN Summary of this function goes here
%   Detailed explanation goes here
    gangiFnOut = (k0./matrix_perm).*power((1-power(((Pc - alpha.*p)./Pmax),m)),3);
end

