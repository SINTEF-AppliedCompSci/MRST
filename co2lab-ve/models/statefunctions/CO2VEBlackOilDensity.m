classdef CO2VEBlackOilDensity < StateFunction

    properties
    end
    
    methods
        function gp = CO2VEBlackOilDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'ShrinkageFactors', 'SurfaceDensity'});
            if model.disgas
                gp = gp.dependsOn({'rs'}, 'state');
            end
            gp.label = '\rho_\alpha';
            gp.outputRange = [0, inf];
        end
        
        function rho = evaluateOnDomain(prop, model, state)
            [b, rhoS] = prop.getEvaluatedDependencies(state, 'ShrinkageFactors', ...
                                                             'SurfaceDensity');
            nph = size(rhoS, 2);
            rho = cell(1, nph);

            % First, include the pure phase densities
            for i = 1:nph
                rho{i} = rhoS{i}.*b{i};
            end
            
            % Include dissolved CO2 density in brine if applicable
            if model.disgas
                [wix, gix] = model.getPhaseIndex('W', 'G');
                rs = model.getProp(state, 'rs');
                rho{wix} = rho{wix} + rs.*b{wix}.*rhoS{gix};
            end
        end
    end
end

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
