function fluid = addVERelpermCapLinear(fluid,varargin)   
    opt=struct('res_oil',0,'res_gas',0,'beta',1,'H',[],'cap_scale',[],'kr_pressure',false);
    opt = merge_options(opt,varargin{:});    
    % should also include endpoint scaling
    assert(~isempty(opt.cap_scale))
    opt=merge_options(opt, varargin{:});
    end_scale=(1-opt.res_oil).^opt.beta;
    if(~opt.kr_pressure)
        fake_pressure=200*barsa;
        fluid.krG=@(sg, p, varargin) end_scale.*krG(sg, fake_pressure, fluid, opt, varargin{:});
        fluid.krOG=@(so, p, varargin) krOG(so,  opt, varargin{:});
        fluid.pcOG=@(sg, p, varargin) pcOG(sg,p ,fluid,opt, varargin{:});
        fluid.cutValues=@(state,varargin) cutValues(state,opt);
        fluid.S3D=@(SVE, samples, H) S3D(SVE,fake_pressure, samples, H, fluid, opt);
        
    else
        fluid.krG=@(sg, p, varargin) end_scale.*krG(sg, p, fluid, opt, varargin{:});
        fluid.krOG=@(so, p, varargin) krOG(so,  opt, varargin{:});
        fluid.pcOG=@(sg, p, varargin) pcOG(sg, p,fluid,opt,varargin{:});
        fluid.cutValues=@(state,varargin) cutValues(state,opt);
    end
    fluid.res_oil=opt.res_oil;
    fluid.res_gas=opt.res_gas;    
    fluid.invPc3D =@(p) invPc(p,opt);
    fluid.kr3D =@(s) end_scale.*(s./(1-opt.res_oil)).^opt.beta;
end

function s= invPc(p,opt)
    s=(p/opt.cap_scale).*(1-opt.res_oil);
    s(s>(1-opt.res_oil))=1-opt.res_oil;
    s(s<0)=0;
    s=1-s;
end

function kr= krG(sg, p, fluid,opt,varargin)
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))
        sg_free = free_sg(sg,loc_opt.sGmax,opt);
        %ineb=sg>loc_opt.sGmax;
        %sg_free=sg-(loc_opt.sGmax-sg)*opt.res_gas;
        %sg_free(ineb)=sg(ineb);
        %sg_free(sg_free<0)=0.0*sg_free(sg_free<0);
        h = invS(sg_free, p, fluid, opt);
        kr = S_beta(h, p, opt.beta, fluid, opt);
    else
        h = invS(sg, p, fluid, opt);
        kr = S_beta(h,p, opt.beta, fluid, opt);        
    end
end

function kr= krOG(so,opt,varargin)
    %beta=1;
    % fo now linear
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))
        sg=1-so;
        ineb=(sg)>loc_opt.sGmax;
        sg_res=(loc_opt.sGmax-sg);
        %assert(all(sg_res>=0));
        so_free=1-loc_opt.sGmax;     
        kr=so_free+(1-opt.res_gas)*sg_res;
        % this to avoid errors in ADI derivative
        kr(~ineb)=so(~ineb);
        kr(kr<0)=0.0*kr(kr<0);
        assert(all(kr>=0));
    else
        kr = so;
    end
end
function pc= pcOG(sg, p, fluid,opt,varargin)
    %h = invS(sg, p, fluid, opt);    
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))
        % could been put in separate function
        sg_free = free_sg(sg,loc_opt.sGmax,opt);
        %ineb=(sg)>loc_opt.sGmax;
        %sg_free=sg-(loc_opt.sGmax-sg)*opt.res_gas;        
        %sg_free(ineb)=sg(ineb);
        %sg_free(sg_free<0)=0.0*sg_free(sg_free<0);
        h = invS(sg_free, p, fluid, opt);
        assert(all(sg_free>=0))
        pc = norm(gravity)*(fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p)).*h;
    else
       h = invS(sg, p, fluid, opt); 
       pc = norm(gravity)*(fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p)).*h;
    end
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
    C=opt.cap_scale;  
    drho  =  (fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p));
    cap_norm=(C./(drho*norm(gravity)))./opt.H;
    %V_cap=cap_norm.^(beta+1)/(beta+1);
    %V_cap=cap_norm/(beta+1);
    h_norm=h./opt.H;
    S_out=(h_norm-cap_norm)+cap_norm./(beta+1);
    ind=h_norm<cap_norm;
    S_out(ind)=(h_norm(ind)./cap_norm(ind)).^(beta+1)/(beta+1).*cap_norm(ind); 
    ind=h_norm>1;
    if(any(ind))
        S_out(ind)=S_out(ind)-((h_norm(ind)-1)./cap_norm(ind)).^(beta+1)/(beta+1).*cap_norm(ind);
    end
end
function h = invS(S, p, fluid, opt)
    % scale availabe gas volume
    % hack to avoid inf
    %{
    ind_zero=double(S)==0;    
    if(isa(S,'ADI'))
        S.val(ind_zero)=sqrt(eps);
    else
        S(ind_zero)=sqrt(eps);
    end
    %}
    
    S=S/(1-opt.res_oil);
    C=opt.cap_scale;  
    drho  =  (fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p));
    %drho  =  (fluid.rhoOS-fluid.rhoGS);
    cap_norm=(C./(drho*norm(gravity)))./opt.H;
    %vol_cap=cap_norm.^2/2;
    vol_cap=(1/2)*cap_norm;
    assert(all(cap_norm<1));
    vol_capb=vol_cap+1-cap_norm;
    h=(2*(S)./cap_norm).^(1/2).*cap_norm;
    
    % reguralize sqare root function near zero
    h_tol=1e-5;S_tol=0.5 *(h_tol./double(cap_norm)).^2.*double(cap_norm);
    h_tmp=h_tol*S./S_tol;
    ind=h<h_tol;
    h(ind)=h_tmp(ind);
         
    
    ind=(S>vol_cap) & (S<vol_capb);
    if(any(ind))
        h(ind)=(S(ind)-vol_cap(ind))+cap_norm(ind);
    end
    ind=(S>=vol_capb);
    if(any(ind))
        %
        %h(ind)=1+S-vol_capb;
        a=-1./(2*cap_norm(ind));b=1;c=-S(ind)+(vol_cap(ind)-cap_norm(ind)+1);
        res=(b).^2-4*a.*c;
        res(res<0)=0;
        h(ind)=(-b+res.^(1/2))./(2*a)+1;

        %h(ind)
    end
    %plot(S,(2*S).^(1/2),S,(S-vol_cap)+cap_norm,S,(vol_cap-0*S),S,(vol_capb-0*S),S(ind),h(ind))
    h=h.*opt.H;
    %{
    if(isa(S,'ADI'))
        h.val(ind_zero)=0;
    else
       h(ind_zero)=0; 
    end
    %}
end
function [s, z_out] = S3D(SVE,p, samples, H, fluid, opt)
    % utility function for plotting saturation frofile and definition of
    % the saturation profile the other functions are obtained from    
    C=opt.cap_scale;  
    drho  =  (fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p));
    cap_norm=(C./(drho*norm(gravity)))./H;
    %vol_cap=cap_norm.^2/2;
    vol_cap=1/2*cap_norm;
    vol_capb=vol_cap+1-cap_norm;
    h=invS(S,p,fluid,opt);
    z=linspace(0,h,samples);
    if(SVE<vol_cap)
       s=z*(h./(opt.H*cap_norm));
    else (SVE<vol_capb)
       s=1; 
       s(z<cap_norm*opt.H)=z./(cap_norm*opt.H); 
    end
    z_out=h-z;
    s=s(end:-1:1);
end
 