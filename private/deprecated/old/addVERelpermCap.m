function fluid = addVERelpermCap(fluid,varargin)   
    opt=struct('alpha',2,'beta',1,'H',[],'cap_scale',[],'kr_pressure',false);
    opt = merge_options(opt,varargin{:});
    assert(opt.alpha>1);    
    assert(~isempty(opt.cap_scale))
    opt=merge_options(opt, varargin{:});        
    if(~opt.kr_pressure)
        fake_pressure=200*barsa;
        fluid.krG=@(sg,varargin) krG(sg, fake_pressure, fluid, opt, varargin{:});
        fluid.krWG=@(so,varargin) krWG(so, fake_pressure, fluid, opt, varargin{:});
        fluid.pcWG=@(sg, p, varargin) pcWG(sg,p ,fluid,opt, varargin{:});
        fluid.cutValues=@(state,varargin) cutValues(state,opt);
        fluid.S3D=@(SVE, samples, H) S3D(SVE,fake_pressure, samples, H, fluid, opt);
        
    else
        fluid.krG=@(sg,varargin) krG(sg, p, opt, fluid, varargin{:});
        fluid.krWG=@(so,varargin) krWG(so, p,opt, fluid, varargin{:});
        fluid.pcWG=@(sg, p, varargin) pcWG(sg, p,fluid,opt,varargin{:});
        fluid.cutValues=@(state,varargin) cutValues(state,opt);

    end
end
function kr= krG(sg,p, fluid, opt,varargin)
    beta=1;
    krGH=S_beta(opt.H,p,opt.beta, fluid, opt);
    %krGH=1;
    sGH=S_beta(opt.H, p, beta, fluid, opt);
    h = invS_beta(sg, p, beta, fluid, opt);
    kr = S_beta(h,p, opt.beta, fluid, opt);
    bind=sg>sGH;
    kr(bind) = krGH(bind)+(1-krGH(bind)).*(sg(bind)-sGH(bind));
end
function kr= krWG(so,opt,varargin)%#ok
    %beta=1;
    % fo now linear
    kr = so;
end
function pc= pcWG(sg, p, fluid,opt,varargin)
    beta=1;
    h = invS_beta(sg, p, beta, fluid, opt);
    pc = norm(gravity)*(fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p)).*h;
end
function state = cutValues(state,opt)
    sg=state.s(:,2);
    sGmax=state.smax(:,2);
    sg=max(sGmax*opt.res_gas,sg);
    sg=min(sg,1);
    state.s=[1-sg,sg];
    state.rs=max(state.rs,0);
    %state.rs=min(state.rs,max_rs);
end
function S_out = S_beta(h, p, beta, fluid, opt)
% integral of power of Saturation
    alpha = opt.alpha;
    %beta  = opt.beta;
    C     = opt.cap_scale;
    p1=C;
    drho  =  (fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p));
    fac =drho*norm(gravity)/C;
    n = -beta/alpha;
    S_out = (fac.*h + (p1/C)).^(n+1) - (p1/C ).^(n+1);
    S_out = (1/(fac*(n+1))).*S_out;
    S_out =S_out./opt.H;
end
function h = invS_beta(S, p, beta, fluid, opt)
    alpha = opt.alpha;
    %beta  = opt.beta;
    C     = opt.cap_scale;
    p1=C;
    drho  =  (fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p));
    fac =drho*norm(gravity)/C;
    n = -beta/alpha;
    h=((((S.*opt.H).*(fac*(n+1))+(p1/C)).^(n+1)).^(1/(n+1)) -(p1/C))./fac;
end
function [s, z_out] = S3D(SVE,p, samples, H, fluid, opt)
    % utility function for plotting saturation frofile and definition of
    % the saturation profile the other functions are obtained from
    alpha = opt.alpha;
    %beta  = opt.beta;
    C     = opt.cap_scale;
    p1=C;
    drho  =  (fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p));
    fac =drho*norm(gravity)/C;
    n = -1/alpha;
    opt.H=H;
    h = invS_beta(SVE, p, 1, fluid, opt);
    z= linspace(0,h,samples);
    z_out= h - z;
    s = (fac.*z+(p1/C)).^(n);
end
    
    