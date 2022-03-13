function fluid = addVERelperm1DTables(fluid, varargin)

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

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

opt = struct('height'      , [] , ...
             'table_co2'   , [] , ...
             'table_water' , [] , ...
             'res_gas'     , 0  , ...
             'res_water'   , 0); 
opt = merge_options(opt, varargin{:}); 


% prop = @(  varargin) properties(opt, varargin{:}); 
fluid.krG = @(sg, p, varargin) krG(sg, p, opt.height, fluid, opt, varargin{:}); 
fluid.krW = @(sw, p, varargin) krW(sw, p, opt.height, fluid, opt, varargin{:}); 
fluid.pcWG = @(sg, p, varargin) cap_press(sg, p, opt.height, fluid, opt, varargin{:}); 

if(opt.table_co2.is_kscaled)
   kscale = sqrt(opt.rock.perm / opt.rock.poro) * fluid.surf_tension; 
   fluid.invPc3D = @(p) opt.table_co2.invPc3D(p ./ kscale); 
else
   fluid.invPc3D = @(p) opt.table_co2.invPc3D(p); 
end

fluid.kr3D = @(s) opt.table_co2.kr3D(s); 
fluid.res_gas = opt.res_gas; 
fluid.res_water = opt.res_water;

end

%---------------------------------------------------------------------
function varargout = cap_press(sg, p, H, fluid, opt, varargin)
% this transformation has to be done twice as long as
% pc and relperm are separate functions
   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 
   sg = free_sg(sg, loc_opt.sGmax, opt); 
   SH = sg .* H; 
   h = interpTable(opt.table_co2.SH, opt.table_co2.h, SH); 
   % dh = dinterpTable(opt.table_co2.SH, opt.table_co2.h, SH) .* H; 
   if(any(h>H))
      disp('Some height are larger than H')
   end
   drho = (fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p)) * norm(gravity); 
   varargout{1} = h .* drho; % (bsxfun(@times, h, drho)); 
                             % if nargout == 2
                             % varargout{2} = drho * dh; 
                             % end
   if(nargout>1)
      error('3 output arguments not implemented');
   end
   
end

% ----------------------------------------------------------------------------

function varargout = krG(sg, p, H, fluid, opt, varargin)
   loc_opt = struct('sGmax', []); 
   loc_opt = merge_options(loc_opt, varargin{:}); 
   sg = free_sg(sg, loc_opt.sGmax, opt); 
   SH = sg .* H; 
   gh = interpTable(opt.table_co2.SH, opt.table_co2.h, SH); 
   kr = interpTable(opt.table_co2.h, opt.table_co2.krH, gh) ./ H; 
   assert(all(SH ./ H <= 1)); 
   varargout{1} = kr;
end

% ----------------------------------------------------------------------------

function varargout = krW(sw, p, H, fluid, opt, varargin)

   varargout{1} = sw; 

end

