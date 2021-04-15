classdef PhaseFluxFixedTotalVelocity < StateFunction
    properties
    end
    
    methods
        function pf = PhaseFluxFixedTotalVelocity(varargin)
            pf@StateFunction(varargin{:});
            pf = pf.dependsOn({'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceMobility', 'FaceTotalMobility'});
        end
        function f = evaluateOnDomain(prop, model, state)
            [vT, G, T, fmob, mobT] = prop.getEvaluatedDependencies(state,...
                 'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceMobility', 'FaceTotalMobility');
            nph = numel(fmob);
            f = cell(1, nph);
            for i = 1:nph
                mobG = 0;
                for j = 1:nph
                    if i ~= j
                        mobG = mobG + fmob{j}.*(G{i} - G{j});
                    end
                end
                f{i} = (fmob{i}./mobT).*(vT + T.*mobG);
            end
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
