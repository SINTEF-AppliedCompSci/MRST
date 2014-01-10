function f = assignPVTO(f, pvto, reg)
f.BO  = @(po, rs, flag, varargin)BO(po, rs, pvto, flag, reg, varargin{:});
f.bO  = @(po, rs, flag, varargin)bO(po, rs, pvto, flag, reg, varargin{:});
f.muO = @(po, rs, flag, varargin)muO(po, rs, pvto, flag, reg, varargin{:});
f.rsSat = @(po, varargin)rsSat(po, pvto, reg, varargin{:});
end

function v = BO(po, rs, pvto, flag, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvto;
for k = 1:numel(T), T{k}.data = T{k}.data(:,1:2); end
v = interpRegPVT(T, po, rs, flag, pvtinx);
end

function v = bO(po, rs, pvto, flag, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvto;
for k = 1:numel(T), T{k}.data = [T{k}.data(:,1), 1./T{k}.data(:,2)]; end
v = interpRegPVT(T, po, rs, flag, pvtinx);
end

function v = muO(po, rs, pvto, flag, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = pvto;
for k = 1:numel(T), T{k}.data = T{k}.data(:,[1 3]); end
v = interpRegPVT(T, po, rs, flag, pvtinx);
end

function v = rsSat(po, pvto, reg, varargin)
pvtinx = getRegMap(po, reg.PVTNUM, reg.PVTINX, varargin{:});
T = cellfun(@(x)[x.data(x.pos(1:end-1),1) x.key], pvto, 'UniformOutput', false);
v = interpReg(T, po, pvtinx);
end
