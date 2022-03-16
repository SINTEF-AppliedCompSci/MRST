classdef NetworkRockSamples < RockSamples
    % Class that maps RockSamples to a network-type reduced reservoir model
    % (GPSNet, INSIM, etc). It has the same properties as RockSamples,
    % along with an additional connectionIndices property that tells us how
    % to map the values from samples.data to the different connections.
    %

    
    properties
        % Inherited properties only
        connectionIndices
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            
            assert(~isempty(samples.connectionIndices), ...
                'Not able to map sampleData to problem since no connection indices given');
            
            numConnections = numel(samples.connectionIndices.cells);
            assert(numel(sampleData.poro) == numConnections, ...
                'mismatch between number of connections and sample data size');
            
            sampleRock = problem.SimulatorSetup.model.rock;
            
            for c = 1:numConnections
                sampleRock.poro(samples.connectionIndices.cells{c}) = sampleData.poro(c);
                sampleRock.perm(samples.connectionIndices.cells{c}) = sampleData.perm(c);
            end
            
            % Update the rock of the model in the problem.
            problem = setSample@RockSamples(samples, sampleRock, problem);
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
