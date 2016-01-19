function f = assignSURFCAPD(f, surfcapd, reg)
   f.miscfact = @(Nc, varargin) miscfact(Nc, surfcapd, reg, varargin{:});
end

function m = miscfact(Nc, surfcapd, reg, varargin)
   satinx = getRegMap(Nc, reg.SATNUM, reg.SATINX, varargin{:});
   surfcapd = extendTab(surfcapd);
   m = interpReg(surfcapd, Nc, satinx);
end
