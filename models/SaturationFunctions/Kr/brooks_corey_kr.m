function [sw, krw, kro] = brooks_corey_kr(swc, sor, krwSor, kroSwc, lambda)
    assert(swc >= 0 & swc < 1,"Wrong input for Swc")
    assert(sor >= 0 & sor < 1-swc,"Wrong input for Sor")
    nPoints = 50;
    sw  = round(linspace(swc, 1 - sor, nPoints)',6);
    Se  = round((sw - swc) ./ (1 - swc - sor),6);
    Se(1) = 0; Se(end) = 1;
    krw = krwSor .* (Se .^ ((2 + 3 .* lambda) / lambda));
    kro = kroSwc .* (((1 - Se) .^ 2) .* (1 - Se .^ ...
        ((2 + lambda) / lambda)));
    krw(1) = 0; kro(end) = 0;
    [sw, krw, kro] = kr_correction(sw, krw, kro);
end
