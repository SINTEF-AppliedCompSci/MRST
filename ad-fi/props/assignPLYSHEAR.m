function f = assignPLYSHEAR(f, plyshear, reg)
    f.plyshearMult = @(c, varargin)plyshearMult(c, plyshear, reg, varargin{:});
end

function v = plyshearMult(c, plyshear, reg, varargin)
    satinx = getRegMap(c, reg.PVTNUM, reg.PVTINX, varargin{:});
    plyshear = extendTab(plyshear);
    v = interpReg(plyshear, c, satinx);
end
