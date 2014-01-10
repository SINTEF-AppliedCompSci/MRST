function f = assignRelPerm(f,opt)
if ~isfield(f, 'krOG')      % two-phase water/oil
    f.relPerm = @(sw, varargin)relPermWO(sw, f, varargin{:});
elseif ~isfield(f, 'krOW')  % two-phase oil/gas
    f.relPerm = @(sg, varargin)relPermOG(sg, f, varargin{:});
else                        % three-phase
    f.relPerm = @(sw, sg, sgmax, varargin)relPermWOG(sw, sg, f, opt, varargin{:});
end
end

function [krW, krO] = relPermWO(sw, f, varargin)
krW = f.krW(sw, varargin{:});
if isfield(f, 'krO')
    krO = f.krO(1-sw, varargin{:});
else
    krO = f.krOW(1-sw, varargin{:});
end
end

function [krO, krG] = relPermOG(sg, f, varargin)
krG = f.krG(sg, varargin{:});
if isfield(f, 'krO')
    krO = f.krO(1-sg, varargin{:});
else
    krO = f.krOG(1-sg, varargin{:});
end
end

function [krW, krO, krG] = relPermWOG(sw, sg, sgmax, f, opt, varargin)
swcon = min(opt.swcr, double(sw)-1e-5);


d  = (sg+sw-opt.swcr);
ww = (sw-opt.swcr)./d;


swm = (sw-opt.swcr)./(opt.swu-opt.sowcr-opt.swcr);
krW = f.krW(swm, varargin{:});

wg = (1-ww);

sg_m =@(sg) (sg-opt.sgcr)./(opt.sgu-opt.sogcr-opt.sgcr);
%sgm=(sg-opt.sgcr)./(opt.sgu-opt.sogcr-opt.sgcr);
sgm = sg_m(sg);
krG = f.krG(sgm, varargin{:});

% hysteresis
sg_eff_max=sg_m(sgmax)
%krGmax=f.krG(sg_eff_max,varargin{:})
sg_m_in =@(sg) (sg-opt.isgcr)./(opt.isgu-opt.isogcr-opt.isgcr);
% since the two curves are equal the crossing in the scaled variable, so
% calculating the corrosing in real variables
sg_m_in_max=sg_eff_max.*(opt.isgu-opt.isogcr-opt.isgcr)+opt.isgcr;
krG_inem=max(f.krG(sg_m_in(sg_m_in_max-(sg_max-sg)),varargin{:}),0);
krG(sg<sgmax)=krG_inem;

so = 1-sw-sg;
som=(so-opt.sowcr)./(1-opt.swcr-opt.sgl-opt.sowcr);
krow = f.krOW(som, varargin{:});
som=(so-opt.sogcr)./(1-opt.isgl-opt.swcr-opt.sgl-opt.sogcr);
krog = f.krOG(som,  varargin{:});
krO  = wg.*krog + ww.*krow;
end
