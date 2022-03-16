classdef NetworkState0Samples < State0Samples
    % Class that maps State0Samples to a network-type reduced reservoir 
    % model (GPSNet, INSIM, etc). It has the same properties as 
    % State0Samples, along with an additional connectionIndices property 
    % that tells us how to map the values from samples.data to the 
    % different connections.
    %
    % The connectionIndices should be a struct containing (either one or
    % both fields) 
    %   - faces: Cell array with the face indices for which to map the
    %            transmissibility values in data{i}.T
    %   - cells: Cell array with the cell indices for which to map the
    %            porevolumes values in data{i}.pv
    properties
        % Inherited properties only
        connectionIndices
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Applies the sample realization of the initial states to a
            % problem.
            %
            % SYNOPSIS:
            %   problem = sample.setSample(sampleData, problem)
            %
            % PARAMETERS:
            %   sampleData - The data for the specific sample realization.
            %
            %   problem - The problem which the sampleData will be applied
            %             to.
            % RETURNS:
            %   problem - A problem representing a single ensemble member
            
            % For now, only transmissibility and porevolume is supported
                     
            
            initSwIsField = isfield(sampleData, 'initSw');
            initSoIsField = isfield(sampleData, 'initSo');
            initSgIsField = isfield(sampleData, 'initSg');
            pressureIsField = isfield(sampleData, 'pressure');
            fluxIsField = isfield(sampleData, 'flux');
            
            assert(initSwIsField + initSoIsField + initSgIsField + pressureIsField + ...
                fluxIsField == numel(fields(sampleData)), ...
                'Invalid fields in sample data for State0Samples');
            
            % Initial saturation
            if (initSwIsField || initSoIsField || initSgIsField)
                numPhases = size(problem.SimulatorSetup.state0.s, 2);       
                sampleData = samples.fillInitialSaturation(sampleData, numPhases, ...
                                                           initSwIsField, initSoIsField, initSgIsField);
            end
            
            for conn = 1:numel(samples.connectionIndices.cells)
                if (initSwIsField || initSoIsField || initSgIsField)
                    problem.SimulatorSetup.state0.s(samples.connectionIndices.cells{conn}, 1) = sampleData.initSw(conn);
                    problem.SimulatorSetup.state0.s(samples.connectionIndices.cells{conn}, 2) = sampleData.initSo(conn);
                    if (numPhases == 3)
                        problem.SimulatorSetup.state0.s(samples.connectionIndices.cells{conn}, 3) = sampleData.initSg(conn);
                    end
                end % if initial saturation

                % Pressure
                if pressureIsField
                    problem.SimulatorSetup.state0.pressure(samples.connectionIndices.cells{conn}) = sampleData.pressure(conn);
                end

                % flux
                if fluxIsField
                    problem.SimulatorSetup.state0.flux(samples.connectionIndices.faces{conn}) = sampleData.flux(conn);
                end
            end
                
        end % function setSample
        
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
