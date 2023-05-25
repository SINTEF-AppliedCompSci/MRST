function fluid = addVERelpermIntegratedFluid(fluid, varargin)
%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
                'int_poro'    , true , ...
                'krw'         , [], ...
                'krg'         , [], ...
                'rock'        , []); 

   opt = merge_options(opt, varargin{:}); 
   
   assert(~isempty(opt.Gt)); 
   assert(~isempty(opt.rock)); 
   
   if isempty(opt.krw)
       opt.krw = 1 - opt.res_gas;
   end 
   if isempty(opt.krg)
       opt.krg = 1 - opt.res_water;
   end
   
   % precalculate the complete perm
   kr_H = integrateVertically(opt.rock.parent.perm(:, end), opt.Gt.cells.H, opt.Gt); 
   opt.perm2D = kr_H ./ opt.Gt.cells.H; % arithmetic mean of permeability along column
   opt.kr_H = kr_H; % integrated permeability along column
   
   if opt.int_poro
       opt.poro3D = opt.rock.parent.poro;
   else
       opt.poro3D = [];
   end
   
   if(~opt.kr_pressure)
      fluid.krG     = @(sg, varargin) krG(sg, opt, varargin{:}); 
      fluid.krW    = @(sw, varargin) krW(sw, opt, varargin{:}); 
   else
      fluid.krG  = @(sg, p, varargin) krG(sg, opt, varargin{:}); 
      fluid.krW = @(sw, p, varargin) krW(sw, opt, varargin{:}); 
   end
   
   fluid.pcWG      = @(sg, p, varargin) pcWG(sg, p, fluid, opt, varargin{:}); 
   fluid.invPc3D   = @(p) invPc3D(p, opt.res_water); 
   %fluid.kr3D      = @(s) s;  % @@ needed?
   fluid.res_gas   = opt.res_gas;
   fluid.res_water = opt.res_water;
end

% ============================================================================

function s = invPc3D(p, res_water)
% return water saturation that is either 1 (if below sharp interface) or
% res_water (if above it)
    s = (sign(p + eps) + 1) / 2 * (1 - res_water);  % CO2 saturation
    s = 1 - s; % converting to water saturation
end

% ----------------------------------------------------------------------------

function pc = pcWG(sg, p, fluid, opt, varargin)
% this transformation has to be done twice as long as pc and relperm are
% separate functions
   loc_opt = struct('sGmax', sg); % default is current saturation
   loc_opt = merge_options(loc_opt, varargin{:}); 
   [h, h_max] = saturation2Height(sg, loc_opt.sGmax, opt.Gt, ...
                                  opt.res_gas, opt.res_water, 'poro', opt.poro3D); 
   if opt.int_poro
      error('Int_poro: not implemented yet!!')
      [h, h_max] = saturation2HeightIntporo(sg, opt, loc_opt); 
   else
      [h, h_max] = saturation2Height(sg, opt, loc_opt.sGmax); 
   end   
   assert(all(h >= 0)); 
   drho = ((fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS * fluid.bG(p)) * norm(gravity)); 
   pc = drho .* h; 
end

% ----------------------------------------------------------------------------

function kr = krG(sg, opt, varargin)

   loc_opt = struct('sGmax', sg);  % default is current saturation
   loc_opt = merge_options(loc_opt, varargin{:}); 

   [h, h_max] = saturation2Height(sg, loc_opt.sGmax, opt.Gt, ...
                                  opt.res_gas, opt.res_water, 'poro', opt.poro3D); 
   if opt.int_poro
      error('Int_poro: not implemented yet!!')
      [h, h_max] = saturation2HeightIntporo(sg, opt, loc_opt); 
   else
      [h, h_max] = saturation2Height(sg, opt, loc_opt.sGmax); 
   end
   assert(all(h >= 0)); 
   if(isa(h, 'ADI'))
      [kr_tmp, dkr_tmp] = integrateVertically(opt.rock.parent.perm(:, 1), h.val, opt.Gt); 
      kr = ADI(kr_tmp, lMultDiag(dkr_tmp, h.jac)); 
   else
      kr = integrateVertically(opt.rock.parent.perm(:, 1), h, opt.Gt); 
   end
   assert(all(kr >= 0)); 
   kr = kr ./ opt.kr_H; % compute fraction of full permeablity
   kr = kr * opt.krg; % end point scaling

end

% ----------------------------------------------------------------------------

function kr = krW(sw, opt, varargin)
   sg = 1 - sw;
   loc_opt = struct('sGmax', sg);  % default is current saturation
   loc_opt = merge_options(loc_opt, varargin{:}); 
   [h, h_max] = saturation2Height(sg, loc_opt.sGmax, opt.Gt, ...
                                  opt.res_gas, opt.res_water, 'poro', opt.poro3D); 
   if opt.int_poro
      error('Int_poro: not implemented yet!!')
      [h, h_max] = saturation2HeightIntporo(sg, opt, loc_opt); 
   else
      [h, h_max] = saturation2Height(sg, opt, loc_opt.sGmax); 
   end
   assert(all(h >= 0)); 
   if(isa(h, 'ADI'))
      [kr, dkr] = integrateVertically(opt.rock.parent.perm(:, 1), h.val, opt.Gt); 
      kr = ADI(kr, lMultDiag(dkr, h.jac)); 
   else
      kr = integrateVertically(opt.rock.parent.perm(:, 1), h, opt.Gt); 
   end
   if(isa(h_max, 'ADI'))
      [kr_max, dkr_max] = integrateVertically(opt.rock.parent.perm(:, 1), h_max.val, opt.Gt); 
      kr_max = ADI(kr_max, lMultDiag(dkr_max, h_max.jac)); 
   else
      kr_max = integrateVertically(opt.rock.parent.perm(:, 1), h_max, opt.Gt); 
   end
   
   kr = ((opt.kr_H - kr_max) + opt.krw .* (kr_max - kr)) ./ opt.kr_H;

   %kr = bsxfun(@rdivide, (opt.kr_H - kr_max) + (1 - opt.res_gas(1)) .* (kr_max - kr), opt.kr_H); 
   kr = ((opt.kr_H - kr_max) + (1 - opt.res_gas(1)) .* (kr_max - kr)) ./ opt.kr_H;
end

% ----------------------------------------------------------------------------

function [h h_max] = saturation2Height(sg, opt, sGmax)
% this transformation is based on the simple transformation
% s * H = h * (1 - sr(2)) + (h_max - h) * sr(1)
% s_max * H = h_max * (1 - sr(2))
   
   s = free_sg(sg, sGmax, opt); 
   h_max = sGmax .* opt.Gt.cells.H ./ (1 - opt.res_water); 
   h = s .* (opt.Gt.cells.H) ./ (1 - opt.res_water); 
   if(any(h<0))
      h(h<0) = 0; 
   end
   if(any(h>opt.Gt.cells.H))
      h(h>opt.Gt.cells.H) = opt.Gt.cells.H(h>opt.Gt.cells.H);
   end
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