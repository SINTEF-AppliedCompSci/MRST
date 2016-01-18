function f = assignSOF2(f, sof2, reg)
f.krO  = @(so, varargin)krO(so, sof2, reg, varargin{:});

if isfield(reg, 'SURFNUM')
   % Assign miscible relperm for surfactant
   f.krOSft  = @(so, varargin)krOSft(so, sof2, reg, varargin{:});
end

end

function v = krO(so, sof2, reg, varargin)
satinx = getRegMap(so, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), sof2, 'UniformOutput', false);
v = interpReg(T, so, satinx);
end

function v = krOSft(so, sof2, reg, varargin)
surfinx = getRegMap(so, reg.SURFNUM, reg.SURFINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), sof2, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, so, surfinx);
end

