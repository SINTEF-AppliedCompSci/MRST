function fluid = addVERelperm(fluid, Gt, varargin)
% Add VE-upscaled rel.perm. (and related functions) for a two-phase fluid
% object and the sharp-interface assumption.
%
% SYNOPSIS:
%   function fluid = addVERelperm(fluid, varargin)
%
% DESCRIPTION:
% Add VE-upscaled relative permeability functions and capillary pressure
% function for a two-phase fluid object (water and CO2).
%
% PARAMETERS:
%   fluid    - Fluid object to modify
%   Gt       - Top surface grid
%   varargin - Option/value pairs, where the following options are available:
%              res_water - residual oil saturation (scalar)
%              res_gas   - residual gas saturation (scalar)
%              krw       - rel. perm of water at full flowing saturation
%              krg       - rel. perm of gas at full flowing saturation
%              top_trap  - Thickness of sub-resolution caprock rugosity
%              surf_topo - Sub-resolution rugosity geometry type.  Can be
%                          'smooth', 'square', 'sinus' or 'inf_rough'.
%
% RETURNS:
%   fluid - The modified fluid object, endowed with the additional
%   functions/fields:
%   res_gas   - residual gas saturation (scalar)
%   res_water   - residual oil saturation (scalar)
%   krG       - upscaled rel.perm. of gas as a function of (gas) saturation
%   krW      - upscaled rel.perm. of water as a function of (water) saturation
%   pcWG      - upscaled 'capillary pressure' as a function of gas saturation
%   invPc3D   - Fine-scale water saturation as function of cap. pressure
%   kr3D      - Dummy function, returning a rel.perm. value that is simply
%               equal to the input saturation.
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


    opt=struct('res_water',     0,...
               'res_gas',       0,...
               'krw',           1,...
               'krg',           1,...
               'top_trap',      [],...
               'surf_topo',     'smooth');

    opt = merge_options(opt, varargin{:});

    fluid.krG=@(sg, p, varargin) krG(sg, Gt, opt, varargin{:});
    fluid.krW=@(sw, p, varargin) krW(sw,opt,varargin{:});

    fluid.pcWG=@(sg, p, varargin) pcWG(sg, p ,fluid, Gt, opt, varargin{:});

    fluid.invPc3D   = @(p) invPc3D(p,opt);
    fluid.kr3D      = @(s) s; % @@ should we rather return nothing here,
                              % since the underlying 3D relperm does not
                              % necessarily have to be linear?
    fluid.res_gas   = opt.res_gas;
    fluid.res_water = opt.res_water;

end

% ============================================================================

function s = invPc3D(p, opt)
% Fine-scale water saturation, considered equal to residual saturation
% ('res_water') in the gas zone and 1 in the water zone.  @@ It doesn't take
% hysteresis into account).
   s = (sign(p + eps) + 1) / 2 * (1 - opt.res_water); 
   s = 1-s;
end

% ----------------------------------------------------------------------------

function kr= krG(sg, Gt, opt,varargin)

    % Check for records of residual saturation
    loc_opt=struct('sGmax',[]);
    loc_opt=merge_options(loc_opt,varargin{:});

    % Determine how much of the gas saturation that  can be considered 'free'
    % (and how much is locked up as residual saturation)
    if(~isempty(loc_opt.sGmax))
        sg_free = free_sg(sg, loc_opt.sGmax, opt);
    else
        sg_free = sg;
    end

    switch opt.surf_topo
      case 'inf_rough'
        % for infinite rough upper surface: kr = (h - dh) / H
        kr = (sg_free .* Gt.cells.H - opt.top_trap) ./ Gt.cells.H; 

      case 'sinus'
        kr2 = ((sg_free .* Gt.cells.H).^2 - opt.top_trap.^2) ./...
              (Gt.cells.H.^2 - opt.top_trap.^2); 
        factor = 1e-4; % Adding small 'fudge factor' to avoid singularity
                       % kr2 = kr2 + 1e-4 * sg_free; 
        kr2(kr2<0) = 0 * sg_free(kr2<0); % @@ Really necessary?
        kr = (kr2).^(0.5); 
        kr(kr2<factor) = (kr2(kr2<factor) / factor) * (factor^(0.5)); 

      case 'square'
        kr_s = (sg_free.^2 - (opt.top_trap ./ Gt.cells.H).^2) ./...
               (sg_free .* (1 - (opt.top_trap ./ Gt.cells.H).^2)); 

        kr = kr_s; 
        kr((opt.top_trap ./ Gt.cells.H) > sg_free) = ...
            0 * sg_free((opt.top_trap ./ Gt.cells.H) > sg_free); 

      case 'smooth'
        kr = sg_free; 

      otherwise
        error('Unknown surface topology')
    end

    if(any(value(kr)<0))
       kr(kr<0) = 0.0 .* sg_free(kr<0); 
    end
    assert(all(kr >= 0)); 
    kr = kr .* opt.krg;
 end

% ----------------------------------------------------------------------------
function kr= krW(sw,opt,varargin)

   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 

   if(~isempty(loc_opt.sGmax))
      sg = 1 - sw; 

      ineb = (sg) > loc_opt.sGmax; 
      
      % compute fraction of aquifer thickness where CO2 saturation is
      % residual ( equivalent to (h_max - h)/H in the height formulation)
      sg_res = (loc_opt.sGmax - sg) / (1 - opt.res_water - opt.res_gas);  

      sw_free = 1 - (loc_opt.sGmax / (1 - opt.res_water)); 
      kr = sw_free + (1 - opt.res_gas) * sg_res; 
      % this to avoid errors in ADI derivative

      if any(ineb)% test necessary since otherwise we risk subtracting an
                  % array of size 0 from a scalar, which will crash
                  % kr = kr .* double(~ineb) + double(~ineb) .* (1 - sg / (1 - opt.res_water)); 
                  % kr(ineb) = (1 - sg(ineb) / (1 - opt.res_water)); 
         kr = ifcond(kr, 1 - sg / (1 - opt.res_water), ~ineb); 
         % kr = min(kr, (1 - sg / (1 - opt.res_water))
      end
      % kr(kr<0) = 0.0 * kr(kr<0); 
      kr = max(kr, 0.0); 
      assert(all(kr >= 0)); 
   else
      kr = sw; 
   end
   kr = kr .* opt.krw; 
   
   kr(sw <= opt.res_water) = 0; % normally it is the case, but slight inaccuracies may
                                % invalidate this assumption.
   % assert(all(double(kr) == 0 | double(sw)>opt.res_water));
   % assert(all(double(kr(double(sw) <= opt.res_water)) == 0));
end

% ----------------------------------------------------------------------------
function pc = pcWG(sg, p, fluid, Gt, opt, varargin)
    loc_opt = struct('sGmax',[], 'T', []);
    loc_opt = merge_options(loc_opt, varargin{:});
    if(~isempty(loc_opt.sGmax))
       % Adjusting the gas saturation to be used for computing cap. press
        sg_free = free_sg(sg, loc_opt.sGmax, opt);
        assert(all(sg_free>=0));
        sg = sg_free; % to be used below
    end

    if isempty(loc_opt.T)
       pc = (fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p)) *...
            norm(gravity) .* sg .* Gt.cells.H;
    else
       % temperature-dependent formation-volume factors
       pc = (fluid.rhoWS .* fluid.bW(p, loc_opt.T) -...
             fluid.rhoGS .* fluid.bG(p, loc_opt.T)) *... 
            norm(gravity) .* sg .* Gt.cells.H; 
    end
    
    pc = pc / (1-opt.res_water);
end

