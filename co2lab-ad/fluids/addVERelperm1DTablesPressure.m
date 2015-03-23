function fluid = addVERelperm1DTablesPressure(fluid,varargin)

%{
#COPYRIGHT#
%}

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

opt = struct('height',[],...
             'table_co2',[],...
             'table_water',[],...
             'kr_pressure',true,...
             'rock',[],...
             'res_water',0,...
             'res_gas',0);
opt = merge_options(opt, varargin{:});



if(~opt.kr_pressure)
        fake_pressure=200*barsa;
        fluid.krG=@(sg,varargin) krG(sg, fake_pressure, fluid, opt, varargin{:});
        fluid.krW=@(so,varargin) krW(so, fake_pressure, fluid, opt, varargin{:});
        fluid.pcWG=@(sg, p, varargin) pcWG(sg,p ,fluid,opt, varargin{:});
        fluid.cutValues=@(state,varargin) cutValues(state,opt);
        %fluid.S3D=@(SVE, samples, H) S3D(SVE,fake_pressure, samples, H, fluid, opt);
        
else
        fluid.krG=@(sg, p, varargin) krG(sg, p, fluid, opt, varargin{:});
        fluid.krW=@(so, p,varargin) krW(so, p, fluid,  opt, varargin{:});
        fluid.pcWG=@(sg, p, varargin) pcWG(sg, p,fluid,opt,varargin{:});
        fluid.cutValues=@(state,varargin) cutValues(state,opt);

end
if(opt.table_co2.is_kscaled)
    %kscale=sqrt(opt.rock.perm./opt.rock.poro)*fluid.surface_tension;
    fluid.invPc3D  = @(p, kscale) opt.table_co2.invPc3D(p./kscale);
else
    fluid.invPc3D  = @(p) opt.table_co2.invPc3D(p);
end
fluid.is_kscaled=opt.table_co2.is_kscaled;
fluid.kr3D  =@(s) opt.table_co2.kr3D(s);
fluid.res_gas=opt.res_gas;
fluid.res_water=opt.res_water;
end

%---------------------------------------------------------------------
function varargout = pcWG(sg, p, fluid, opt, varargin)
   % this trasformation has to be doen twise as long as
   % pc aand relperm i separate functions
   loc_opt=struct('sGmax',[]);
   loc_opt=merge_options(loc_opt,varargin{:});
   sg = free_sg(sg,loc_opt.sGmax,opt);
   drho=((fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS*fluid.bG(p))*norm(gravity));
   H=opt.height;   
   dSP=H.*drho;
   SP= sg.*dSP;
   if(opt.table_co2.is_kscaled)
      kscale=sqrt(opt.rock.poro./opt.rock.perm).*fluid.surface_tension; 
      SP = SP ./kscale;
      dSP =dSP ./kscale;
   end
   pc=interpTable(opt.table_co2.SP,opt.table_co2.p,SP);      
   %dp=dinterpTable(opt.table_co2.SP,opt.table_co2.p,SP).*dSP;
   if(opt.table_co2.is_kscaled)     
      pc=pc.*(kscale);
   %   dp=p.*(kscale);
   end
   varargout{1}=pc;
   if(any(SP>dSP))
       disp('Some hight lager than H')
   end
   %if nargout==2
   %  varargout{2} = dp;
   %end
   
   
   if(nargout>1)
     error('3 output arguments not implemented');
   end
end
function varargout = krG(sg, p, fluid, opt, varargin)
   loc_opt=struct('sGmax',[]);
   loc_opt=merge_options(loc_opt,varargin{:});
   sg = free_sg(sg,loc_opt.sGmax,opt);
   H=opt.height;
   drho=((fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS*fluid.bG(p))*norm(gravity));
   dSP=H.*drho;
   SP= sg.*dSP;   
   if(opt.table_co2.is_kscaled)
      kscale=sqrt(opt.rock.poro./opt.rock.perm).*fluid.surface_tension; 
      SP = SP ./kscale;
      dSP = dSP./kscale;   
   end
   pc=interpTable(opt.table_co2.SP,opt.table_co2.p,SP);      
   %dp=dinterpTable(opt.table_co2.SP,opt.table_co2.p, SP).*dSP;
   h_int=pc./(drho);
   if(opt.table_co2.is_kscaled)
       h_int=h_int.*kscale;
   end
   
   kr=interpTable(opt.table_co2.p,opt.table_co2.krP, pc)./dSP;
   %dkr=dinterpTable(opt.table_co2.p,opt.table_co2.krP, p).*dp./(dSP);
   
   
   %if(any(pc>dpH))
   h_ind=h_int>opt.height;
   
   if(any(h_ind))
       dpH=opt.height(h_ind).*drho(h_ind);%*(1-opt.res_water);
       disp(['Capillary is larger than H for ', num2str(sum(h_ind)),'values'])       
       if(opt.table_co2.is_kscaled)
           dpH=dpH./kscale(h_ind);           
       end
       spH=interpTable(opt.table_co2.p,opt.table_co2.SP,dpH);
       sH=spH./drho(h_ind);
       if(opt.table_co2.is_kscaled)
          sH=sH.*kscale(h_ind); 
       end
       
       kr_H=interpTable(opt.table_co2.p,opt.table_co2.krP,dpH)./dpH;
       kr_end=1;%opt.table_co2.kr3D(1-opt.res_water);
       kr(h_ind)=kr_H+((sg(h_ind)-sH./H(h_ind))./(1-sH./H(h_ind))).*(kr_end-kr_H);
        %kr(h_ind)=kr_H+((sg(h_ind)-sH./H(h_ind))).*(1-kr_H);
    %   dkr(p>dpH)=(1./(1-(sH./H(p>dpH)))).*(1-kr_H);
   end
   
   varargout{1} = kr;%./H;
   %varargout{2} = dkr;%./H;%.*dh;
end
function kr = krW(so, p, fluid, opt, varargin)
   kr=interpTable(opt.table_water.S,opt.table_water.kr, 1-so);
   %dkr=dinterpTable(opt.table_water.S,opt.table_water.kr, so); 
   %varargout{1} = kr;
   %varargout{2} = dkr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55


