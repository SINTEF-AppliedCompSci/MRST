function fluid = addVERelperm(fluid,varargin)   
    opt=struct('res_oil',       0,...
               'res_gas',       0,...
               'Gt',            [],...
               'kro',           1,...
               'krg',           1,...
               'top_trap',      [],...
               'surf_topo',     'smooth');
   
    opt = merge_options(opt, varargin{:});        
    fluid.krG=@(sg, p, varargin) krG(sg,opt,varargin{:});
    fluid.krOG=@(so, p, varargin) krOG(so,opt,varargin{:});
    fluid.pcOG=@(sg, p, varargin) pcOG(sg,p ,fluid,opt,varargin{:});
    fluid.cutValues=@(state,varargin) cutValues(state,opt);
    fluid.invPc3D = @(p) invPc3D(p,opt);
    fluid.kr3D =@(s) s;
    fluid.res_gas = opt.res_gas;
    fluid.res_oil =opt.res_oil;
    %fluid=rmfield(fluid,'relPerm');
end

function s = invPc3D(p,opt)
         s=(sign(p+eps)+1)/2*(1-opt.res_oil);
         s=1-s;
end

function kr= krG(sg,opt,varargin)
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))                
        sg_free = free_sg(sg,loc_opt.sGmax,opt);
    else
        sg_free = sg;
    end
    switch opt.surf_topo
        case 'inf_rough'
            % for infinite rough upper surface
            kr = (sg_free.*opt.Gt.cells.H-opt.top_trap)./opt.Gt.cells.H;% simple model
        case 'sinus'
            kr2 = ((sg_free.*opt.Gt.cells.H).^2-opt.top_trap.^2)./(opt.Gt.cells.H.^2-opt.top_trap.^2);
            %sigma=1e4;
            kr2=kr2+1e-4*sg_free;
            kr2(kr2<0)=0*sg_free(kr2<0);
            kr=(kr2).^(0.5);
            %{
            kr(kr2<=0)=0*kr2(kr<=0);
            kr(kr2>0)=kr2(kr2>0).^(0.5);
            %}
        case 'square'
           kr_s=(sg_free.^2-(opt.top_trap./opt.Gt.cells.H).^2)./(sg_free.*(1-(opt.top_trap./opt.Gt.cells.H).^2));
           %kr=sg_free;
           kr=kr_s;
           kr((opt.top_trap./opt.Gt.cells.H)>sg_free)=0*sg_free((opt.top_trap./opt.Gt.cells.H)>sg_free);
           %kr(opt.top_trap>0*sg_free)=kr_s(opt.top_trap>0*sg_free);
           %kr(opt.top)=kr_s(opt.top_trap>0);
           %kr(sg_free<=0)=0.0*sg_free(sg_free<=0);
           %kr(kr<0)=kr(kr<0).*0.0;
        case 'smooth'
           kr = sg_free;
        otherwise
           error('Unknown surface topology')
    end
    if(any(double(kr)<0))
        kr(kr<0)=0.0.*sg_free(kr<0);
        kr(kr<0)=0.0.*sg_free(kr<0);
    end
    assert(all(kr>=0));
    kr=kr.*opt.krg;
end

function kr= krOG(so,opt,varargin)
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))
        sg=1-so;
        %sg_free = free_sg(sg,sGmax,opt); 
        ineb=(sg)>loc_opt.sGmax;
        sg_res=(loc_opt.sGmax-sg);
        %assert(all(sg_res>=0));
        so_free=1-(loc_opt.sGmax/(1-opt.res_oil));     
        kr=so_free+(1-opt.res_gas)*sg_res;
        % this to avoid errors in ADI derivative
        %kr(~ineb)=so(~ineb)
        kr(ineb)=(1-sg(ineb)/(1-opt.res_oil));
        kr(kr<0)=0.0*kr(kr<0);
        assert(all(kr>=0));
    else
        kr = so;
    end
    kr=kr.*opt.krg;
    assert(all(double(kr(so<=opt.res_oil))==0));
end
function pc= pcOG(sg, p, fluid,opt,varargin)
    loc_opt=struct('sGmax',[]);    
    loc_opt=merge_options(loc_opt,varargin{:});
    if(~isempty(loc_opt.sGmax))
        % could been put in separate function
        sg_free = free_sg(sg,loc_opt.sGmax,opt);
        assert(all(sg_free>=0))
        pc = norm(gravity)*(fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p)).*sg_free.*opt.Gt.cells.H;
    else
       pc = norm(gravity)*(fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p)).*sg.*opt.Gt.cells.H;
    end
    pc=pc/(1-opt.res_oil);
end
function state = cutValues(state,opt)
    sg=state.s(:,2);
    %sGmax=state.smax(:,2);
    %sg=max(sGmax*opt.res_gas,sg);
    %sg=min(sg,1);
    %sg=min(sg,1-opt.res_oil);
    sg=max(sg,0);
    state.s=[1-sg,sg];
    state.rs=max(state.rs,0);
    %state.sGmax=min(1,state.sGmax);
    %state.sGmax=max(0,state.sGmax);
    %state.rs=min(state.rs,max_rs);
end
