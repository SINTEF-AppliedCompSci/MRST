function generatePropsTable(fluidname, pname, P_range, T_range, P_num, T_num)
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
   [P, T, vals]  = make_prop_table(P_range, T_range, P_num, T_num, fluidname, pname);

   save(propFilename(P_range, T_range, P_num, T_num, fluidname, pname), 'P', 'T' , 'vals');

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
