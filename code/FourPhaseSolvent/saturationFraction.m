function f = saturationFraction(sa, sb, tol)

    f = (sa./(sa + sb + tol)).*(1+tol);
    f(isnan(double(f))) = 0;

%     tol = eps;
%     noa = sa < tol;
%     f = sa./(sa + sb).*(~noa) + 0.*noa;
%     f(isnan(double(f))) = 0;
    
    
end