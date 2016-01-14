function f = assignSURFST(f, surfst, reg)
f.ift = @(c, varargin) ift(c, surfst, reg, varargin{:});
end

function v = ift(c, surfst, reg, varargin)
satinx = getRegMap(c, reg.PVTNUM, reg.PVTINX, varargin{:});
surfst = extendTab(surfst);
v = interpReg(surfst, c, satinx);
end
