function f = assignSURFVISC(f, surfvisc, reg)
f.muWSft = @(c, varargin) muWSft(c, surfvisc, reg, varargin{:});
end

function v = muWSft(c, surfvisc, reg, varargin)
satinx = getRegMap(c, reg.PVTNUM, reg.PVTINX, varargin{:});
surfvisc = extendTab(surfvisc);
v = interpReg(surfvisc, c, satinx);
end
