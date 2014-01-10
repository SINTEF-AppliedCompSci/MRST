function f = assignSOF2(f, sof2, reg)
f.krO  = @(so, varargin)krO(so, sof2, reg, varargin{:});
end

function v = krO(so, sof2, reg, varargin)
satinx = getRegMap(so, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), sof2, 'UniformOutput', false);
v = interpReg(T, so, satinx);
end
