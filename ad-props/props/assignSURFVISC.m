function f = assignSURFVISC(f, surfvisc, reg)
f.muWSftMult = @(c, varargin) muWSftMult(c, surfvisc, reg, varargin{:});
end

function v = muWSftMult(c, surfvisc, reg, varargin)
satinx = getRegMap(c, reg.PVTNUM, reg.PVTINX, varargin{:});
surfvisc = extendTab(surfvisc);
v = interpReg(surfvisc, c, satinx);
end
