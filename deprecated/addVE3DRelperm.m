function fluid = addVE3DRelperm(fluid,varargin)   
    opt=struct('res_water',0,'res_gas',0,'Gt',[]);
    opt=merge_options(opt, varargin{:});        
    fluid.krG=@(sg,varargin) krG(sg,opt,varargin{:});
    fluid.krW=@(so,varargin) krW(so,opt,varargin{:});
    fluid.pcWG=@(sg, varargin) pcWG(sg,fluid,opt,varargin{:});
    fluid.cutValues=@(state,varargin) cutValues(state,opt);
    fluid.invPc3D = @(p) (sign(p)+1)/2;
    fluid.kr3D =@(s) s;
end
function kr= krG(sg,opt,varargin)
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))
        
        ineb=sg>loc_opt.sGmax;
        sg_free=sg-(loc_opt.sGmax-sg)*opt.res_gas;                
        sg_free(ineb)=sg(ineb);
        sg_free(sg_free<0)=0.0*sg_free(sg_free<0);
        
        kr=sg_free;
        kr(kr<0)=0.0.*sg_free(kr<0);
        assert(all(kr>=0));
    else
        kr = sg;
    end
end
function kr= krW(so,opt,varargin)
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
function pc= pcWG(sg, fluid,opt,varargin)
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))
        % could been put in separate function
        ineb=(sg)>loc_opt.sGmax;
        sg_free=sg-(loc_opt.sGmax-sg)*opt.res_gas;        
        sg_free(ineb)=sg(ineb);
        sg_free(sg_free<0)=0.0*sg_free(sg_free<0);
        
        assert(all(sg_free>=0))
        pc = norm(gravity)*(fluid.rhoWS-fluid.rhoGS).*sg_free.*opt.Gt.cells.H;
    else
       pc = norm(gravity)*(fluid.rhoWS-fluid.rhoGS).*sg.*opt.Gt.cells.H;
    end
end
function state = cutValues(state,opt);
    sg=state.s(:,2);
    sGmax=state.smax(:,2);
    sg=max(sGmax*opt.res_gas,sg);
    sg=min(sg,1);
    state.s=[1-sg,sg];
    state.rs=max(state.rs,0);
    %state.rs=min(state.rs,max_rs);
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
