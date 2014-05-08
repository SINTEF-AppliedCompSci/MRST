function fluid = addVERelperm(fluid, Gt, varargin)
% Add VE-upscaled rel.perm. (and related functions) for a two-phase fluid object
%
% SYNOPSIS:
%   function fluid = addVERelperm(fluid, varargin)
%
% DESCRIPTION:
% Add VE-upscaled rel.perm. (and related functions) for a two-phase fluid
% object.  The two phases are described as oil and gas, but could also be
% interpreted as water and gas for a CO2 system.
%
% PARAMETERS:
%   fluid    - Fluid object to modify
%   Gt       - Top surface grid
%   varargin - Option/value pairs, where the following options are available:
%              res_oil   - residual oil saturation (scalar)
%              res_gas   - residual gas saturation (scalar)
%              kro       - rel. perm of oil at full flowing saturation
%              krg       - rel. perm of gas at full flowing saturation
%              top_trap  - Thickness of sub-resolution caprock rugosity
%              surf_topo - Sub-resolution rugosity geometry type.  Can be
%                          'smooth', 'square', 'sinus' or 'inf_rough'.
%
% RETURNS:
%   fluid - The modified fluid object, endowed with the additional
%   functions/fields:
%   res_gas   - residual gas saturation (scalar)
%   res_oil   - residual oil saturation (scalar)
%   krG       - upscaled rel.perm. of gas as a function of (gas) saturation
%   krOG      - upscaled rel.perm. of oil as a function of (oil) saturation
%   pcOG      - upscaled 'capillary pressure' as a function of gas saturation
%   invPc3D   - Fine-scale oil saturation as function of cap. pressure
%   cutValues - Function to ensure consistency and nonzero values in
%               saturation fields of the 'state' variable ('state.s' and
%               'state.rs')
%   kr3D      - Dummy function, returning a rel.perm. value that is simply
%               equal to the input saturation.


    opt=struct('res_oil',       0,...
               'res_gas',       0,...
               'kro',           1,...
               'krg',           1,...
               'top_trap',      [],...
               'surf_topo',     'smooth');

    opt = merge_options(opt, varargin{:});

    fluid.krG=@(sg, p, varargin) krG(sg, Gt, opt, varargin{:});
    fluid.krOG=@(so, p, varargin) krOG(so,opt,varargin{:});

    fluid.pcOG=@(sg, p, varargin) pcOG(sg, p ,fluid, Gt, opt, varargin{:});

    fluid.cutValues = @(state,varargin) cutValues(state,opt);
    fluid.invPc3D   = @(p) invPc3D(p,opt);
    fluid.kr3D      = @(s) s;
    fluid.res_gas   = opt.res_gas;
    fluid.res_oil   = opt.res_oil;

end

% ============================================================================
function s = invPc3D(p, opt)
% Fine-scale oil saturation, considered equal to residual saturation
% ('res_oil') in the gas zone and 1 in the oil zone.  @@ It doesn't take
% hysteresis into account).
         s=(sign(p+eps)+1)/2*(1-opt.res_oil);
         s=1-s;
end

% ----------------------------------------------------------------------------
function kr= krG(sg, Gt, opt,varargin)

    % Check for records of residual saturation
    loc_opt=struct('sGmax',[]);
    loc_opt=merge_options(loc_opt,varargin{:});

    % Determine how much of the gas saturation that  can be considered 'free'
    % (and how much is locked up as residual saturation)
    if(~isempty(loc_opt.sGmax))
        sg_free = free_sg(sg,loc_opt.sGmax,opt);
    else
        sg_free = sg;
    end

    switch opt.surf_topo
        case 'inf_rough'
            % for infinite rough upper surface: kr = (h-dh)/H
            kr = (sg_free .* Gt.cells.H - opt.top_trap) ./ Gt.cells.H;

      case 'sinus'
            kr2 = ((sg_free .* Gt.cells.H).^2 - opt.top_trap.^2) ./ ...
                  (Gt.cells.H.^2 - opt.top_trap.^2);
            factor=1e-4;% @@ Adding small 'fudge factor' to avoid singularity
            %kr2=kr2+1e-4*sg_free;  
            kr2(kr2<0)=0*sg_free(kr2<0);% @@ Really necessary?? 
            kr=(kr2).^(0.5);            
            kr(kr2<factor)=(kr2(kr2<factor)/factor)*(factor^(0.5));

      case 'square'
           kr_s=(sg_free.^2 - (opt.top_trap ./ Gt.cells.H).^2) ./ ...
                (sg_free .* (1 - (opt.top_trap ./ Gt.cells.H).^2));

           kr=kr_s;
           kr((opt.top_trap ./ Gt.cells.H) > sg_free) = ...
               0 * sg_free((opt.top_trap ./ Gt.cells.H) > sg_free);

      case 'smooth'
           kr = sg_free;

      otherwise
           error('Unknown surface topology')
    end

    if(any(double(kr)<0))
        kr(kr<0)=0.0.*sg_free(kr<0);
    end
    assert(all(kr>=0));
    kr=kr .* opt.krg;
end

% ----------------------------------------------------------------------------
function kr= krOG(so,opt,varargin)

    loc_opt=struct('sGmax',[]);
    loc_opt=merge_options(loc_opt,varargin{:});

    if(~isempty(loc_opt.sGmax))
        sg = 1 - so;

        ineb = (sg) > loc_opt.sGmax;
        sg_res = (loc_opt.sGmax - sg);

        so_free = 1 - (loc_opt.sGmax / (1 - opt.res_oil));
        kr=so_free+(1-opt.res_gas)*sg_res;
        % this to avoid errors in ADI derivative

        if any(ineb) % test necessary since otherwise we risk subtracting an
                     % array of size 0 from a scalar, which will crash
            kr(ineb)=(1-sg(ineb)/(1-opt.res_oil));
        end
        kr(kr<0)=0.0*kr(kr<0);
        assert(all(kr>=0));
    else
        kr = so;
    end
    %kr=kr.*opt.krg;
    kr=kr.*opt.kro;  % @@ Odd: the above line most likely a typo...?
    assert(all(double(kr(so<=opt.res_oil))==0));
end

% ----------------------------------------------------------------------------
function pc = pcOG(sg, p, fluid, Gt, opt, varargin)
    loc_opt = struct('sGmax',[]);
    loc_opt = merge_options(loc_opt, varargin{:});
    if(~isempty(loc_opt.sGmax))
        % could been put in separate function
        sg_free = free_sg(sg, loc_opt.sGmax, opt);
        assert(all(sg_free>=0));
        pc = (fluid.rhoOS .* fluid.bO(p) - fluid.rhoGS .* fluid.bG(p)) * ...
             norm(gravity) .* sg_free .* Gt.cells.H;
    else
       pc = (fluid.rhoOS.*fluid.bO(p)-fluid.rhoGS.*fluid.bG(p)) * ...
            norm(gravity) .*sg.*Gt.cells.H;
    end
    pc = pc / (1-opt.res_oil);
end

% ----------------------------------------------------------------------------
function state = cutValues(state, opt)
    sg = state.s(:,2);
    %sGmax=state.smax(:,2);
    %sg=max(sGmax*opt.res_gas,sg);
    %sg=min(sg,1);
    %sg=min(sg,1-opt.res_oil);
    sg = max(sg,0);
    state.s = [1-sg, sg];
    state.rs = max(state.rs, 0);  % @@ Anything missing?  state.rs capped,
                                  %    but not otherwise updated
    %state.sGmax=min(1,state.sGmax);
    %state.sGmax=max(0,state.sGmax);
    %state.rs=min(state.rs,max_rs);
end
