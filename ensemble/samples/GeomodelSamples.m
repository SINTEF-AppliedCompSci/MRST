classdef GeomodelSamples < RockSamples
    
    methods
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            model  = problem.SimulatorSetup.model;
            model.rock = sampleData.rock;
            model.G = sampleData.G;
            model = model.setupOperators();
            problem.SimulatorSetup.model = model;
            W = problem.SimulatorSetup.schedule.control(1).W;
            for i = 1:max(model.G.cells.tag)
                W(i).cells = find(model.G.cells.tag == i);
            end
            problem.SimulatorSetup.schedule.control(1).W = W;
            if samples.updateWells
                % Update well indices
                schedule = problem.SimulatorSetup.schedule;
                schedule = samples.updateWellIndices(model, schedule);
                problem.SimulatorSetup.schedule = schedule;
            end
            state0 = problem.SimulatorSetup.state0;
            state0 = initResSol(model.G, state0.pressure(1), state0.s(1,:));
            problem.SimulatorSetup.state0 = state0;
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
