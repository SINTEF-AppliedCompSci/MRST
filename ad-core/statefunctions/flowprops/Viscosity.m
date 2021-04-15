classdef Viscosity < StateFunction
    
    methods
        function mu = Viscosity(model, varargin)
            mu@StateFunction(model, varargin{:});
            mu = mu.dependsOn('PhasePressures');
            mu.label = '\mu_\alpha';
            mu.outputRange = [0, inf];
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            names = model.getPhaseNames();
            nph = numel(names);
            mu = cell(1, nph);
            p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
            [sample, isAD] = getSampleAD(p_phase{:});
            for ph = 1:nph
                mu{ph} = prop.evaluatePhaseViscosity(model, state, names(ph), p_phase{ph});
            end
            if isAD
                for i = 1:numel(mu)
                    if ~isa(mu{i}, 'ADI')
                        mu{i} = model.AutoDiffBackend.convertToAD(mu{i}, sample);
                    end
                end
            end
        end
        
        function mu = evaluatePhaseViscosity(prop, model, state, name, p) %#ok
            mu = prop.evaluateFluid(model, ['mu', name], p);
        end
        
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
