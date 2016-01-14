function f = assignSURFVISC(f, surfvisc, reg)
f.muWSurfMult = @(c, varargin) muWSurfMult(c, surfvisc, reg, varargin{:});
end

function v = muWSurfMult(c, surfvisc, reg, varargin)
satinx = getRegMap(c, reg.PVTNUM, reg.PVTINX, varargin{:});
surfvisc = extendTab(surfvisc);
v = interpReg(surfvisc, c, satinx);
end
