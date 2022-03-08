classdef ThermalDensity < StateFunction
%State function for computing pressure- and temperature-dependent density
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = ThermalDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures'             , ...
                               'ComponentPhaseMassFractions'});
            gp = gp.dependsOn({'temperature'}, 'state');
            if model.getNumberOfPhases() == 2
                gp = gp.dependsOn({'pressure', 'enthalpy'}, 'state');
            end
            gp.label = '\rho_\alpha';
        end
        
        %-----------------------------------------------------------------%
        function rho = evaluateOnDomain(prop, model, state)
            nph = model.getNumberOfPhases();
            rho = cell(1, nph);
            if nph == 2
                [p, h] = model.getProps(state, 'pressure', 'enthalpy');
                rhoMix = prop.evaluateFluid(model, 'rho', p, h);
                [rho{:}] = deal(rhoMix);
                twoPhase = state.flag == 3;
                if any(twoPhase)
                    rhoPh = prop.computePhaseDensites(model, state);
                    for i = 1:nph
                        rho{i}(twoPhase) = rhoPh{i}(twoPhase);
                    end
                end
            else
                rho = prop.computePhaseDensites(model, state);
            end
        end
        
        %-----------------------------------------------------------------%
        function rho = computePhaseDensites(prop, model, state)
            % Get phase pressures, mass fractions and temperature
            [p, X] = prop.getEvaluatedDependencies(state, 'PhasePressures'             , ...
                                                          'ComponentPhaseMassFractions');
            T = model.getProp(state, 'temperature');
            % Identify NaCl (if present)
            cnames = model.getComponentNames();
            ix     = strcmpi(cnames, 'NaCl');
            if any(ix), X = X{ix}; else, X = []; end
            phases = model.getPhaseNames();
            rho    = cell(1, numel(phases));
            % Compute density using fluid function handle
            for i = 1:numel(phases)
                rho{i} = prop.evaluateFluid(model, ['rho', phases(i)], p{i}, T, X);
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