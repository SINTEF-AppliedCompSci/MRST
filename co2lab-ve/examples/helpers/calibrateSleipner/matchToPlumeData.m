function obj = matchToPlumeData(model, states, plumes)
%Undocumented Utility Function

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

    numSteps = length(states); % in practice: either 1 or all the states

    obj = repmat({[]}, numSteps, 1);
    for step = 1:numSteps

        if numSteps == 1
            assert(~iscell(states)); % should just be a single state
            state = states; 
            plume = plumes{1};
        else
            state = states{step};
            plume = plumes{step};
        end

        sG = state.s(:,2);
        if iscell(sG)
            assert(length(sG) == 1)
            sG = sG{1};
        end
        
        if ~isempty(plume)
            o = matchToCo2Surface(sG, plumes{step}, model.G, model.fluid);
            alpha = 1/9; % tested alpha = 0, 0.1, 1 (a scaling value for dz)
            obj{step} = o + sum(alpha * (model.G.dz).^2 .* model.G.cells.volumes);
        else
            obj{step} = double2ADI(0, sG);
        end
    end
end
