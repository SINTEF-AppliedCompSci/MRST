function f = assignPVTG(f, pvtg, reg)
f.BG  = @(pg, rv, flag, varargin)BG(pg, rv, pvtg, flag, reg, varargin{:});
f.bG  = @(pg, rv, flag, varargin)bG(pg, rv, pvtg, flag, reg, varargin{:});
f.muG = @(pg, rv, flag, varargin)muG(pg, rv, pvtg, flag, reg, varargin{:});
f.rvSat = @(pg, varargin)rvSat(pg, pvtg, reg, varargin{:});
end

function v = BG(pg, rv, pvtg, flag, reg, varargin)
pvtinx = getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvtg;
for k = 1:numel(T), T{k}.data = T{k}.data(:,1:2); end
v = interpRegPVT(T, rv, pg, flag, pvtinx);
end

function v = bG(pg, rv, pvtg, flag, reg, varargin)
pvtinx = getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvtg;
for k = 1:numel(T), T{k}.data = [T{k}.data(:,1), 1./T{k}.data(:,2)]; end
v = interpRegPVT(T, rv, pg, flag, pvtinx);
end

function v = muG(pg, rv, pvtg, flag, reg, varargin)
pvtinx = getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvtg;
for k = 1:numel(T), T{k}.data = T{k}.data(:,[1 3]); end
v = interpRegPVT(T, rv, pg, flag, pvtinx);
end

function v = rvSat(pg, pvtg, reg, varargin)
pvtinx = getRegMap(pg, reg.PVTNUM, reg.PVTINX, varargin{:});
T = cellfun(@(x)[x.key x.data(x.pos(1:end-1),1)], pvtg, 'UniformOutput', false);
v = interpReg(T, pg, pvtinx);
end
