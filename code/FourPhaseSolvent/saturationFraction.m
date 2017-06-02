function f = saturationFraction(sa, sb)

    tol = 1e-4;
%     tol = 0;
    sb = sb+tol;
    c = 1+tol;  
    f = sa./(sa + sb).*c;
%     f(isnan(double(f))) = 0;

end