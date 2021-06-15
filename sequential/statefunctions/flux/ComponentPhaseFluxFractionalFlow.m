classdef ComponentPhaseFluxFractionalFlow < StateFunction
    properties
    end
    
    methods
        function cf = ComponentPhaseFluxFractionalFlow(backend, upwinding)
            cf@StateFunction(backend);
            cf = cf.dependsOn({'FaceComponentMobility', 'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceMobility', 'FaceTotalMobility'});
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [compMob, vT, G, T, fmob, mobT] = prop.getEvaluatedDependencies(state,...
                 'FaceComponentMobility', 'TotalFlux', 'PhaseInterfacePressureDifferences', 'Transmissibility', 'FaceMobility', 'FaceTotalMobility');
            w = 1./mobT;
            nf = numelValue(vT);
            kgrad = cell(1, nph);
            for i = 1:nph
                mobG = zeros(nf, 1);
                for j = 1:nph
                    if i ~= j
                        mobG = mobG + fmob{j}.*(G{i} - G{j});
                    end
                end
                kgrad{i} = w.*(vT + T.*mobG);
            end
            v = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    mob = compMob{c, ph};
                    if ~isempty(mob)
                        v{c, ph} = mob.*kgrad{ph};
                    end
                end
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
