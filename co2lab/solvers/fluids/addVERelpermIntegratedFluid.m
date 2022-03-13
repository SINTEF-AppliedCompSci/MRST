function fluid = addVERelpermIntegratedFluid(fluid, varargin)
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

  opt = struct('res_water'   , 0     , ...  
                'res_gas'     , 0     , ...
                'kr_pressure' , false , ...
                'Gt'          , []    , ...
                'int_poro'    , false , ...
                'rock'        , []); 

   opt = merge_options(opt, varargin{:}); 
   % should also include endpoint scaling    

   assert(~isempty(opt.Gt)); 
   assert(~isempty(opt.rock)); 

   % precalculate the complete perm and pore volume
   kr_H = integrateVertically(opt.rock.parent.perm(:, 1), opt.Gt.cells.H, opt.Gt); 
   opt.perm2D = kr_H ./ opt.Gt.cells.H; 
   assert(norm(opt.perm2D - opt.rock.perm) ./ norm(opt.rock.perm)<1e-6); 
   opt.kr_H = kr_H; 
   pv_3D(opt.Gt.columns.cells) = opt.rock.parent.poro(opt.Gt.columns.cells)...
       .* rldecode(opt.Gt.cells.volumes, diff(opt.Gt.cells.columnPos)); 
   opt.volumes = integrateVertically(pv_3D', inf, opt.Gt) .* opt.Gt.cells.volumes; 
   
   if(~opt.kr_pressure)
      fluid.krG     = @(sg, varargin) krG(sg, opt, varargin{:}); 
      fluid.krW    = @(so, varargin) krW(so, opt, varargin{:}); 
      fluid.pcWG    = @(sg, p, varargin) pcWG(sg, p, fluid, opt, varargin{:}); 
      %fake_pressure = 200 * barsa;  % @@ make the value of this constant optional?
      %fluid.S3D     = @(SVE, samples, H) S3D(SVE, fake_pressure, samples, H, fluid, opt); 
   else
      fluid.krG  = @(sg, p, varargin) krG(sg, opt, varargin{:}); 
      fluid.krW = @(so, p, varargin) krW(so, opt, varargin{:}); 
      fluid.pcWG = @(sg, p, varargin) pcWG(sg, p, fluid, opt, varargin{:}); 
   end 
   fluid.invPc3D   = @(p) invPc3D(p, opt); 
   fluid.kr3D      = @(s) s; 
   fluid.res_gas   = opt.res_gas;
   fluid.res_water = opt.res_water;
end

% ============================================================================

function s = invPc3D(p,opt)
         s=(sign(p+eps)+1)/2*(1-opt.res_water);
         s=1-s;
end

% ----------------------------------------------------------------------------

function pc = pcWG(sg, p, fluid, opt, varargin)
% this transformation has to be done twice as long as pc and relperm are
% separate functions
   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 
   if opt.int_poro
      error('Int_poro: not implemented jet!!')
      [h, h_max] = saturation2HeightIntporo(sg, opt, loc_opt); 
   else
      [h, h_max] = saturation2Height(sg, opt, loc_opt); 
   end   
   assert(all(h >= 0)); 
   drho = ((fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS * fluid.bG(p)) * norm(gravity)); 
   pc = drho .* h; 
end

% ----------------------------------------------------------------------------

function kr = krG(sg, opt, varargin)
   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 
   if opt.int_poro
      error('Int_poro: not implemented yet!!')
      [h, h_max] = saturation2HeightIntporo(sg, opt, loc_opt); 
   else
      [h, h_max] = saturation2Height(sg, opt, loc_opt); 
   end
   assert(all(h >= 0)); 
   if(isa(h, 'ADI'))
      [kr_tmp, dkr_tmp] = integrateVertically(opt.rock.parent.perm(:, 1), h.val, opt.Gt); 
      kr = ADI(kr_tmp, lMultDiag(dkr_tmp, h.jac)); 
   else
      kr = integrateVertically(opt.rock.parent.perm(:, 1), h, opt.Gt); 
   end
   kr = kr * (1 - opt.res_water); 
   assert(all(kr >= 0)); 
   kr = kr ./ opt.kr_H; % bsxfun(@rdivide, kr, opt.kr_H);
end

% ----------------------------------------------------------------------------

function kr = krW(sg, opt, varargin)
   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 
   if opt.int_poro
      error('Int_poro: not implemented yet!!')
      [h, h_max] = saturation2HeightIntporo(sg, opt, loc_opt); 
   else
      [h, h_max] = saturation2Height(sg, opt, loc_opt); 
   end
   assert(all(h >= 0)); 
   if(isa(h, 'ADI'))
      [kr, dkr] = integrateVertically(opt.rock.parent.perm(:, 1), h.val, opt.Gt); 
      kr = ADI(kr, lMultDiag(dkr, h.val)); 
   else
      kr = integrateVertically(opt.rock.parent.perm(:, 1), h, opt.Gt); 
   end
   if(isa(h_max, 'ADI'))
      [kr_max, dkr_max] = integrateVertically(opt.rock.parent.perm(:, 1), h_max.val, opt.Gt); 
      kr_max = ADI(kr_max, lMultDiag(dkr_max, h_max.val)); 
   else
      kr_max = integrateVertically(opt.rock.parent.perm(:, 1), h_max, opt.Gt); 
   end
   
   kr = bsxfun(@rdivide, (opt.kr_H - kr_max) + (1 - opt.res_gas(1)) .* (kr_max - kr), opt.kr_H); 
   kr = ((opt.kr_H - kr_max) + (1 - opt.res_gas(1)) .* (kr_max - kr)) ./ opt.kr_H;

end

% ----------------------------------------------------------------------------

function [h h_max] = saturation2Height(sg, opt, loc_opt)
% this transformation is based on the simple transformation
% s * H = h * (1 - sr(2)) + (h_max - h) * sr(1)
% s_max * H = h_max * (1 - sr(2))
   
   s = free_sg(sg, loc_opt.sGmax, opt); 
   h_max = loc_opt.sGmax .* opt.Gt.cells.H ./ (1 - opt.res_water); 
   h = s .* (opt.Gt.cells.H) ./ (1 - opt.res_water); 
   if(any(h<0))
      h(h<0) = 0; 
   end
   if(any(h>opt.Gt.cells.H))
      h(h>opt.Gt.cells.H) = opt.Gt.cells.H(h>opt.Gt.cells.H);
   end
end

% ----------------------------------------------------------------------------

function [h h_max] = saturation2HeightIntPoro(sg, opt, loc_opt)
% this transformation is based on the simple tranformation
% s * H = V(h) * (1 - sr(2)) + V(h_max) - V(h)) * sr(1)
% s_max * V(H) = V(h_max) * (1 - sr(2))
   s_max = loc_opt.sGmax; 
   s = sg; 
   assert(numel(s) == numel(opt.H)); 
   Vh_max = bsxfun(@rdivide, s_max .* opt.volumes, (1 - opt.res_water)); 
   V_h = s .* opt.volumes - bsxfun(@times, Vh_max, opt.res_gas); 
   V_h = bsxfun(@rdivide, V_h, (1 - opt.res_water - opt.res_gas)); 
   h = opt.Vinv(V_h); 
   h_max = opt.Vinv(Vh_max); 

end

% ----------------------------------------------------------------------------

function J = lMultDiag(d, J1)
   n = numel(d); 
   D = sparse((1:n)', (1:n)', d, n, n); 
   J = cell(1, numel(J1)); 
   for k = 1:numel(J)
      J{k} = D * J1{k}; 
   end
end