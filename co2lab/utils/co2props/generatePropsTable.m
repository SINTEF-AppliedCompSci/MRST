function generatePropsTable(savedir, fluidname, pname, P_range, T_range, P_num, T_num)
% Generates and saves a sampled table of fluid properties, using 'coolprops'.
%  
% NB: 'coolprops' needs to be available, which means that the directory
% containing the function 'propsSI' needs to be in the search path, and the
% corresponding mex-files must be compiled and available
% 
% SYNOPSIS:
%   function tables = generatePropsTable(fluidname, pname, P_range, T_range, P_num, T_num)
%
% DESCRIPTION:
%
% PARAMETERS:
%   fluidname - name of fluid.  should be 'Carbon dioxide' or 'Water'
%   pname     - property name. 'D' for density; 'V' for viscosity or 'H' for enthalpy
%   P_range   - [pmin, pmax] : lower and upper end of sampled pressure range (Pascal)
%   T_range   - [tmin, tmax] : lower and upper end of sampled temp. range (Kelvin)
%   P_num     - number of (evenly spaced) pressure samples within 'P_range'
%   T_num     - number of (evenly spaced) temperature samples within 'T_range'
%
% RETURNS:
%   Nothing.  The generated table is saved directly, along with the
%             information on pressure and temperature ranges. 
%
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

   [P, T, vals]  = make_prop_table(P_range, T_range, P_num, T_num, fluidname, pname);%#ok

   save([savedir propFilename(P_range, T_range, P_num, T_num, fluidname, pname)], 'P', 'T' , 'vals');

end

% ----------------------------------------------------------------------------

function [P, T, vals] = make_prop_table(P_range, T_range, P_num, T_num, fluid, prop)
   
    pvals = linspace(P_range(1), P_range(2), P_num);
    tvals = linspace(T_range(1), T_range(2), T_num);
    
    P.span = P_range;
    P.stepsize = (P_range(2)-P_range(1))/(P_num-1);
    P.num = P_num;
    
    T.span = T_range;
    T.stepsize = (T_range(2)-T_range(1))/(T_num-1);
    T.num = T_num;
    
    vals = zeros(P_num, T_num);
    for c = 1:T_num
        for r = 1:P_num
           vals(r, c) = PropsSI(prop, 'P', pvals(r), 'T', tvals(c), fluid);
        end
    end
end
