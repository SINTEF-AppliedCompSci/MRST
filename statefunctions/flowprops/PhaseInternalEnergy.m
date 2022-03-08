classdef PhaseInternalEnergy < StateFunction
% State function for internal enery in each fluid phase
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = PhaseInternalEnergy(model, varargin)
            gp@StateFunction(model, varargin{:});
            switch model.thermalFormulation
                case 'temperature'
                    % Temperature formulation - use fluid function handles
                    % to compute internal energy
                    gp = gp.dependsOn({'PhasePressures'              , ...
                                       'ComponentPhaseMassFractions'}, ...
                                       'PVTPropertyFunctions'        );
                    gp = gp.dependsOn({'temperature'}, 'state');
                case 'enthalpy'
                    % Enthalpy formulation - compute internal energy from
                    % enthalpy, pressure and density
                    gp = gp.dependsOn({'PhasePressures'      , ...
                                       'PhaseEnthalpy'       , ...
                                       'Density'             }, ...
                                       'PVTPropertyFunctions');
            end
            gp.label = 'u_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function u = evaluateOnDomain(prop,model, state)
            switch model.thermalFormulation
                case 'enthalpy'
                    % Get pressure, enthalpy and density
                    [p, h, rho] = model.getProps(state, 'PhasePressures', ...
                                                        'PhaseEnthalpy' , ...
                                                        'Density'       );
                    % Compute u from enthalpy, pressure and density
                    u = cellfun(@(p,h,rho) h - p./rho, p, h, rho, 'UniformOutput', false);
                case 'temperature'
                    % Get pressure, temperature and mass fraction
                    [p, T, X] = model.getProps(state, 'PhasePressures'             , ...
                                                      'temperature'                , ...
                                                      'ComponentPhaseMassFractions');
                    % Identyfy NaCl (if present)
                    cnames = model.getComponentNames();
                    ix     = strcmpi(cnames, 'NaCl');
                    if any(ix), X = X{ix}; else, X = []; end
                    phases = model.getPhaseNames();
                    nph    = model.getNumberOfPhases();
                    u      = cell(1, nph);
                    % Compute u using fluid function handles
                    for i = 1:nph
                        ix    = model.getPhaseIndex(phases(i));
                        u{ix} = prop.evaluateFluid(model, ['u', phases(i)], p{ix}, T, X);
                    end
            end
        end       
    end
end

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