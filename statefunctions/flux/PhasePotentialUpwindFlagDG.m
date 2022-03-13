classdef PhasePotentialUpwindFlagDG < PhasePotentialUpwindFlag
    properties
    end
    
    methods
        function gp = PhasePotentialUpwindFlagDG(varargin)
            gp@PhasePotentialUpwindFlag(varargin{:});
        end

        function flags = evaluateOnDomain(prop, model, state)
            mob = model.getProp(state, 'Mobility');
            nph = numel(mob);
            flags = cell(1, nph);
            if strcmpi(state.type, 'face')
                [G, T, vT] = prop.getEvaluatedDependencies(state, ...
                                'PhaseInterfacePressureDifferences', 'Transmissibility', 'TotalFlux');
                if prop.includeTotalVelocity
                    v = vT;
                else
                    v = zeros(size(value(vT)));
                end
                flag_array = multiphaseUpwindIndices(G, v, T, mob, model.operators.faceUpstr);
                for i = 1:nph
                    flags{i} = flag_array(:, i);
                end
            else
                [flags{:}] = deal(true);
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
