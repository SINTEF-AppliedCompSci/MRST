function f = assignPLYSHEAR(f, plyshear, reg)
% Polymer shear thinning/thickening
    f.plyshearMult = @(VW, varargin) plyshearMult(VW, plyshear, reg, ...
        varargin{:});
end

function v = plyshearMult(VW, plyshear, reg, varargin)
    satinx = getRegMap(VW, reg.PVTNUM, reg.PVTINX, varargin{:});
    plyshear = extendTab(plyshear);
    v = interpReg(plyshear, VW, satinx);
end
