function [sw, krw, kro] = LET_kr(swc, sor, krwSor, kroSwc, Lw, Ew, Tw, Lnw, Enw, Tnw)
    assert(swc >= 0 & swc < 1,"Wrong input for Swc")
    assert(sor >= 0 & sor < 1-swc,"Wrong input for Sor")
    nPoints = 50;
    sw  = round(linspace(swc, 1 - sor, nPoints)',6);
    Se  = round((sw - swc) ./ (1 - swc - sor),6);
    Se(1) = 0; Se(end) = 1;
    krw = krwSor * ( (Se .^ Lw) ./ (Se .^ Lw + Ew .*...
        (1 - Se) .^ Tw));
    kro = kroSwc * ( ( (1 - Se) .^ Lnw ) ./ ( (1 - Se) .^ Lnw...
        + Enw .* (Se .^ Tnw) ) );
    krw(1) = 0; kro(end) = 0;
    [sw, krw, kro] = kr_correction(sw, krw, kro);
end
