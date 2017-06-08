function f = saturationFraction(sa, sb, tol)

    c = 1 + tol;
    f = (sa./(sa + sb + tol)).*c;
    f(isnan(double(f))) = 0;

end