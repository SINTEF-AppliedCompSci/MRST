function f = saturationFraction(sa, sb)

%     tol = 10*eps;
%     zz = num < tol & den < tol;
% %     num = num - (num < tol).*num;
% %     den = den - (den < tol).*den;
%     f = (num./(den + (den < tol).*tol)).*(~zz) + 0.*zz;
    
%     tol = 10*eps;
%     num = num - (num < tol).*num;
%     den = den - (den < tol).*den;

%     num = num - (num < 0).*num;
%     den = den - (den < 0).*den;

%     f = num./den;
%     f(isnan(double(f))) = 0;

    tol = 1e-2;
    sb(sb < tol) = tol;
    f = sa./(sa + sb);
    
end