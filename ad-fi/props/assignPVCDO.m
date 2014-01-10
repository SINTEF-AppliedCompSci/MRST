function f = assignPVCDO(f, pvcdo, reg)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    f.cO  = pvcdo(1, 3);
else
    f.cO  = pvcdo(reg.PVTNUM, 3);
end
f.BO     = @(po, varargin)BO(po, pvcdo, reg, varargin{:});
f.bO     = @(po, varargin)bO(po, pvcdo, reg, varargin{:});
f.BOxmuO = @(po, varargin)BOxmuO(po, pvcdo, reg, varargin{:});
end

function v = BO(po, pvcdo, reg, inx)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    pvtnum = 1;
elseif nargin < 4
    pvtnum = reg.PVTNUM;
    assert(numel(pvtnum)==size(po,1));
else
    pvtnum = reg.PVTNUM(inx);
end
por  = pvcdo(pvtnum,1); % ref pres
bor  = pvcdo(pvtnum,2); % ref fvf
co   = pvcdo(pvtnum,3); % compress
X = co.*(po-por);
v = bor.*exp(-X);
end

function v = bO(po, pvcdo, reg, inx)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    pvtnum = 1;
elseif nargin < 4
    pvtnum = reg.PVTNUM;
    assert(numel(pvtnum)==size(po,1));
else
    pvtnum = reg.PVTNUM(inx);
end
por  = pvcdo(pvtnum,1); % ref pres
bor  = pvcdo(pvtnum,2); % ref fvf
co   = pvcdo(pvtnum,3); % compress
X = co.*(po-por);
v = exp(X)./bor;
end

function v = BOxmuO(po, pvcdo, reg, inx)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    pvtnum = 1;
elseif nargin < 4
    pvtnum = reg.PVTNUM;
    assert(numel(pvtnum)==size(po,1));
else
    pvtnum = reg.PVTNUM(inx);
end
por  = pvcdo(pvtnum,1); % ref pres
bor  = pvcdo(pvtnum,2); % ref fvf
co   = pvcdo(pvtnum,3); % compress
muor = pvcdo(pvtnum,4); % ref visc
vbo  = pvcdo(pvtnum,5); % viscosibility
Y = (co-vbo).*(po-por);
v = bor.*muor.*exp(-Y);
end
