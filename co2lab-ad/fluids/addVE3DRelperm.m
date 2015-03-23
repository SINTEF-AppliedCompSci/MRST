function fluid = addVE3DRelperm(fluid,varargin)   
    opt=struct('res_water',0,'res_gas',0,'Gt',[]);
    opt=merge_options(opt, varargin{:});        
    fluid.krG=@(sg,varargin) krG(sg,opt,varargin{:});
    fluid.krOG=@(so,varargin) krOG(so,opt,varargin{:});
    fluid.pcOG=@(sg, varargin) pcOG(sg,fluid,opt,varargin{:});
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
function kr= krOG(so,opt,varargin)
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
function pc= pcOG(sg, fluid,opt,varargin)
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))
        % could been put in separate function
        ineb=(sg)>loc_opt.sGmax;
        sg_free=sg-(loc_opt.sGmax-sg)*opt.res_gas;        
        sg_free(ineb)=sg(ineb);
        sg_free(sg_free<0)=0.0*sg_free(sg_free<0);
        
        assert(all(sg_free>=0))
        pc = norm(gravity)*(fluid.rhoOS-fluid.rhoGS).*sg_free.*opt.Gt.cells.H;
    else
       pc = norm(gravity)*(fluid.rhoOS-fluid.rhoGS).*sg.*opt.Gt.cells.H;
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