function f = assignSWOF(f, swof, reg)
f.krW  = @(sw, varargin)krW(sw, swof, reg, varargin{:});
f.krOW = @(so, varargin)krOW(so, swof, reg, varargin{:});
f.pcOW = @(sw, varargin)pcOW(sw, swof, reg, varargin{:});
swcon  = cellfun(@(x)x(1,1), swof);
ntsat = numel(reg.SATINX);
if ntsat == 1
    f.sWcon = swcon(1);
else
    f.sWcon = swcon(reg.SATNUM);
end

if isfield(reg, 'SURFNUM')
   % Assign miscible relperm for surfactant
   f.krWSft  = @(sw, varargin)krWSft(sw, swof, reg, varargin{:});
   f.krOWSft  = @(so, varargin)krOWSft(so, swof, reg, varargin{:});
end

end

function v = krW(sw, swof, reg, varargin)
satinx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, satinx);
end

function v = krOW(so, swof, reg, varargin)
satinx = getRegMap(so, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, 1-so, satinx);
end

function v = pcOW(sw, swof, reg, varargin)
satinx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = cellfun(@(x)x(:,[1,4]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, satinx);
end

function v = krWSft(sw, swof, reg, varargin)
surfinx = getRegMap(sw, reg.SURFNUM, reg.SURFINX, varargin{:});
T = cellfun(@(x)x(:,[1,2]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, sw, surfinx);
end

function v = krOWSft(so, swof, reg, varargin)
surfinx = getRegMap(so, reg.SURFNUM, reg.SURFINX, varargin{:});
T = cellfun(@(x)x(:,[1,3]), swof, 'UniformOutput', false);
T = extendTab(T);
v = interpReg(T, so, surfinx);
end

