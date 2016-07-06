function fluidimp = veFluid2Fluid(G,fluid,H)   
% Function to convert ve-fluid (made with initVEFluid to fluid that can be
% used with s-formulation, at the moment hack and hardcoded geometry

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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

      fluidimp.pc=@(rsol) pc(rsol, fluid, H);
      fluidimp.relperm =@(s, varargin) relperm(s);
      
      fluidimp.saturation = @(rsol) rsol.s;  
      fluidimp.properties = @(varargin) properties(fluid,varargin{:});

end

function varargout = pc(rsol, fluid, H)   
   varargout{1}    =  norm(gravity)*(fluid.rho(1)-fluid.rho(2))*fluid.pc(rsol.s.*H);

   if nargout > 1,
      varargout{2} = norm(gravity)*(fluid.rho(1)-fluid.rho(2))*H.*ones(numel(rsol.s),1);
   end
end


function [mu, rho] = properties(fluid,varargin)
   mu = fluid.mu; 
   rho =fluid.rho;
end
function [kr, dkr] = relperm(s)
   kr = [s,1-s]; 
   dkr = repmat([1,-1],numel(s),1);
end
