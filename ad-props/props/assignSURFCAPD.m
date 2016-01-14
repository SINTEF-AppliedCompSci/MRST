function f = assignSURFCAPD(f, surfcapd, reg)
f.surfads = @(c, varargin) surfads(c, surfcapd, reg, varargin{:});
end

function v = surfads(c, surfcapd, reg, varargin)
satinx = getRegMap(c, reg.SATNUM, reg.SATINX, varargin{:});
surfcapd = extendTab(surfcapd);
v = interpReg(surfcapd, c, satinx);
end
