function f = assignSWFN(f, swfn, reg)
f.krW  = @(sw, varargin)krW(sw, swfn, reg, varargin{:});
f.pcOW = @(sw, varargin)pcOW(sw, swfn, reg, varargin{:});
swcon = cellfun(@(x)x(1,1), swfn);
ntsat = numel(reg.SATINX);
if ntsat == 1
    f.sWcon = swcon(1);
else
    f.sWcon = swcon(reg.SATNUM);
end
end

function v = krW(sw, swfn, reg, varargin)
satinx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), swfn, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, satinx);
end

function v = pcOW(sw, swfn, reg, varargin)
satinx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), swfn, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, satinx);
end

