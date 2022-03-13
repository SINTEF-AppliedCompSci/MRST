classdef OperatorSamples < BaseSamples
    % Class for holding stochastic samples that represent uncertain
    % operator properties, such as transmissibility and porevolume, for an
    % MRSTEnsemble.
    %
    % DESCRIPTION:
    %   This class is an extention of `BaseSamples` and is used within a 
    %   MRSTEnsemble to organize stochastic operator properties. These  
    %   properties can either be pre-computed or generated on the fly.
    % 
    % SYNOPSIS
    %   samples = OperatorSamples('data', data);
    %   samples = OperatorSamples('generatorFn', generatorFn);
    %
    % OPTIONAL PARAMETERS
    %   'data' - precomputed sample data. Supported formats:
    %               * Cell array of data samples
    %               * Instance of ResultHandler class with
    %                 information about storage location and names. See
    %                 ResultHandler class for details
    %
    %   'generatorFn' - Function for generating a stochastic sample
    % 
    % NOTE:
    %   Either 'data' or 'generatorFn' must be provided.
    %   Each sample should consist of a struct with the fields also found
    %   in SimulatorSetup.model.operators
    %
    % SEE ALSO:
    %   `RockSamples`, `DeckSamples`, `BaseSamples`, `MRSTExample`, `BaseQoI`
    
    methods
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Applies the sample realization of the operator properties to  
            % a problem.
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
            
            transmissibilityIsField = isfield(sampleData, 'T');
            porevolumeIsField = isfield(sampleData, 'pv');
            assert(transmissibilityIsField + porevolumeIsField  == numel(fields(sampleData)), ...
                'Only T and pv are currently supported in OperatorSamples');
            
            if transmissibilityIsField
                assert(all(size(sampleData.T) == size(problem.SimulatorSetup.model.operators.T)), ...
                    'sampleData.T and problem...operators.T does not match');
                problem.SimulatorSetup.model.operators.T = sampleData.T;
            end
            
            if porevolumeIsField
                assert(all(size(sampleData.pv) == size(problem.SimulatorSetup.model.operators.pv)), ...
                    'sampleData.pv and problem...operators.pv does not match');
                problem.SimulatorSetup.model.operators.pv = sampleData.pv;
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
