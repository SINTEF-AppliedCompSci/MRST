function f = assignSURFCAPD(f, surfcapd, reg)
   f.surfcapd = @(Nc, varargin) surfcapd(Nc, surfcapd, reg, varargin{:});
end

function m = surfcapd(Nc, surfcapd, reg, varargin)
   satinx = getRegMap(Nc, reg.SATNUM, reg.SATINX, varargin{:});
   surfcapd = extendTab(surfcapd);
   m = interpReg(surfcapd, Nc, satinx);
end
