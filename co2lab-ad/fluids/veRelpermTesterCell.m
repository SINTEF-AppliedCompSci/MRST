function [s,pc,kr, s_max, s_free,fval]= veRelpermTesterCell(hs, drho, fluid, H, varargin)
    opt=struct('samples',100,'hs_max',[],'kscale',[]);
    opt=merge_options(opt, varargin{:});
    if(isfield(fluid,'is_kscaled'))
        if(fluid.is_kscaled)
             invPc3D =@(p) fluid.invPc3D(p, opt.kscale);
        else
            invPc3D =@(p) fluid.invPc3D(p);
        end
    else
       invPc3D =@(p) fluid.invPc3D(p);        
    end
    assert(all(hs<H));
    h=linspace(0,hs,opt.samples)';
    dh=h(2)-h(1);
    %drho=(fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p))*norm(gravity);
    drho=drho*norm(gravity);
    s_h = invPc3D(h.*drho); 
    sg_h=1-s_h;
    kr_h=fluid.kr3D(sg_h);
    pc=hs.*drho;
    SH=(sum(sg_h)-sg_h(end)/2)*dh;
    krH=(sum(kr_h)-kr_h(end)/2)*dh;
    if(~isempty(opt.hs_max))
     hs_max=opt.hs_max;
      h_max=linspace(0,hs_max,opt.samples)';
      dh_max=h_max(2)-h_max(1);
      s_hmax= invPc3D(h_max.*drho);
      sg_hmax=1-s_hmax;
      SH_max=(sum(sg_hmax)-sg_hmax(end)/2)*dh_max;
      s_max=SH_max./H;      
      s=SH./H;
      s_free=s;
      %SH=((SH_max-SH)/(1-fluid.res_water))*fluid.res_gas+SH;
      SH=((SH_max-SH)/(1-fluid.res_water))*fluid.res_gas+SH;
      %SH=(SH_max-SH)*fluid.res_gas+SH;
      %s = free_sg(s,s_max,fluid);
      s=SH./H;
    else
      s=SH./H;  
    end
    kr=krH./H;
    fval=struct('h',h,'s_h',s_h,'kr_h',kr_h,'s_hmax',s_hmax,'h_max',h_max);
end
%function fsg = free_sg(sg,sGmax,opt)
%        ineb=sg>sGmax;
%        %sg_free=sg-(loc_opt.sGmax-sg)*opt.res_gas;                
%        fsg=((1-opt.res_water)*sg-(sGmax*opt.res_gas))./(1-opt.res_gas-opt.res_water);
%        %assert(all(fsg>=-sqrt(eps)));
%        fsg(ineb)=sg(ineb);
%        %fsg(fsg<0)=0.0*fsg(fsg<0);
%end    