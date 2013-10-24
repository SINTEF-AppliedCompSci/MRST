function fluid = addVERelperm1DTables(fluid,varargin)

%{
#COPYRIGHT#
%}

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

opt = struct('height',[],'table_co2',[],'table_water',[],'res_gas',0,'res_oil',0);
opt = merge_options(opt, varargin{:});


%prop = @(  varargin) properties(opt, varargin{:});
fluid.krG   = @(sg, p, varargin) krG(sg, p, opt.height, fluid, opt,varargin{:});
fluid.krOG   = @(so, p, varargin) krOG(so, p, opt.height,  fluid, opt,varargin{:});
fluid.pcOG   = @(sg, p, varargin) cap_press(sg, p, opt.height, fluid, opt, varargin{:});
if(opt.table_co2.is_kscaled)
    kscale=sqrt(opt.rock.perm/opt.rock.poro)*fluid.surf_tension;
    fluid.invPc3D  = @(p) opt.table_co2.invPc3D(p./kscale);
else
    fluid.invPc3D  = @(p) opt.table_co2.invPc3D(p);
end
fluid.kr3D   = @(s) opt.table_co2.kr3D(s);
fluid.res_gas=opt.res_gas;
fluid.res_oil=opt.res_oil;

end

%---------------------------------------------------------------------
function varargout = cap_press(sg, p, H, fluid, opt, varargin)
   % this trasformation has to be doen twise as long as
   % pc aand relperm i separate functions
   loc_opt=struct('sGmax',[]);
   loc_opt=merge_options(loc_opt,varargin{:});
   sg = free_sg(sg,loc_opt.sGmax,opt);
   SH= sg.*H;
   h=interpTable(opt.table_co2.SH,opt.table_co2.h,SH);
   %dh=dinterpTable(opt.table_co2.SH,opt.table_co2.h,SH).*H;
   if(any(h>H))
       disp('Some hight lager than H')
   end
   drho= (fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p))*norm(gravity);   
   varargout{1} = h.*drho;%(bsxfun(@times,h,drho));
   %if nargout==2
   %  varargout{2} = drho*dh;
   %end
   if(nargout>1)
     error('3 output arguments not implemented');
   end
   
end
function varargout = krG(sg, p, H, fluid, opt, varargin)
   loc_opt=struct('sGmax',[]);
   loc_opt=merge_options(loc_opt,varargin{:});
   sg = free_sg(sg,loc_opt.sGmax,opt);
   SH= sg.*H;   
   gh=interpTable(opt.table_co2.SH,opt.table_co2.h,SH);
   %dgh=dinterpTable(opt.table_co2.SH,opt.table_co2.h,SH).*H;         
   kr=interpTable(opt.table_co2.h,opt.table_co2.krH, gh)./H;
   %dkr=dinterpTable(opt.table_co2.h,opt.table_co2.krH, gh)./H;
   if(any(gh>H))
       disp(['Capillary is larger than H for ', num2str(sum(gh>H)),'values'])
       sH=interpTable(opt.table_co2.h,opt.table_co2.SH,H(gh>H));
       kr_H=interpTable(opt.table_co2.h,opt.table_co2.SH,sH)./H(gh>H);
       kr(gh>H)=kr_H+((sg(gh>H)-sH./H(gh>H))./(1-sH./H(gh>H))).*(1-kr_H);
   %    dkr(gh>H)=(1./(1-(sH./H(gh>H)))).*(1-kr_H);
   end
   varargout{1} = kr;
   %varargout{2} = dkr.*dgh;
end
function varargout = krOG(so, p, H, fluid, opt, varargin)
   kr=interpTable(opt.table_water.h,opt.table_water.krH, so);
   %dkr=dinterpTable(opt.table_water.S,opt.table_water.kr, so); 
   varargout{1} = kr;
   %varargout{2} = dkr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%

