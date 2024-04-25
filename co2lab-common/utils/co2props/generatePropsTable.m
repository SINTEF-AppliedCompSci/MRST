function generatePropsTable(savedir, fluidname, pname, P_range, T_range, P_num, T_num)
% Generates and saves a sampled table of fluid properties, using `CoolProp`.
%  
% SYNOPSIS:
%   function tables = generatePropsTable(fluidname, pname, P_range, T_range, P_num, T_num)
%
% DESCRIPTION:
%   
% Generate and save a sampled table of fluid properties, in a format that can be used 
% by makeVEFluid.  Note that this function requires `CoolProp` (see below).
% 
% PARAMETERS:
%   fluidname - name of fluid.  should be 'CarbonDioxide' or 'Water'
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
% NOTE: `CoolProp` is available at www.coolprop.org, but is also directly
% accessible to MATLAB through its Python interface.  To run this function, you
% need to have Python available for your MATLAB installation.  If Python is
% present but the `CoolProp` package is not installed, the user will be
% given the option to install it using `pip` when running this function.
% 
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
    ensure_coolprops_available();
    
    [P, T, vals]  = make_prop_table(P_range, T_range, P_num, T_num, fluidname, pname);%#ok

    save(fullfile(savedir, ...
                  propFilename(P_range, T_range, P_num, T_num, fluidname, pname)), ...
         'P', 'T' , 'vals');
end

% ----------------------------------------------------------------------------
function res = coolprop_available()
    res = ~(isa(py.importlib.util.find_spec("CoolProp"), 'py.NoneType'));
end

% ----------------------------------------------------------------------------
function ensure_coolprops_available()
    
    % check if python is available
    [ver, exe, loaded] = pyversion();
    if ~loaded && isempty(ver)
        error('Python not available.  Required for running CoolProp.');
    end
    
    % check if CoolProp package is installed
    if ~coolprop_available()
        try_install = userConsent(['The Python `CoolProp package was not found. ' ...
                                   'Try to install with `pip`']);
        if try_install
            % try to install
            pe = pyenv();
            system([pe.Executable{:}, ' -m pip install --user -U CoolProp']);
        
            if ~coolprop_available()
                error('Did not succeed in installing `CoolProp`.');
            end
        else
            error('Cannot continue.  `CoolProp` not available.');
        end
    end
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
            vals(r, c) = py.CoolProp.CoolProp.PropsSI(prop, 'P', pvals(r), ...
                                                      'T', tvals(c), fluid);
             %vals(r, c) = PropsSI(prop, 'P', pvals(r), 'T', tvals(c), fluid);
        end
    end
end
