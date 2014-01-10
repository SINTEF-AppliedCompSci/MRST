function f = assignDENSITY(f, density, reg)
% dens of size ntpvtx3
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    f.rhoOS = density(1, 1);
    f.rhoWS = density(1, 2);
    f.rhoGS = density(1, 3);
else
    f.rhoOS = density(reg.PVTNUM, 1);
    f.rhoWS = density(reg.PVTNUM, 2);
    f.rhoGS = density(reg.PVTNUM, 3);
end
