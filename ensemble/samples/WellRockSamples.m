classdef WellRockSamples < WellSamples & RockSamples
    % Class that combines ensemble configurations for both rock samples and
    % well samples.
    % For use with an MRSTEnsemble.
    %
    % DESCRIPTION:
    %   This class combines `RockSamples` and `WellSamples` so that
    %   stochastic rock and well properties can be applied within the
    %   same MRSTEnsemble. These well properties can either be pre-computed
    %   or generated on the fly.
    %
    %   The data property for this class is expected to be a cell array with
    %   structs containing the fields 'rock' and 'well', and each of them
    %   should be compatible with the data properties expected in RockSamples
    %   and WellSamples, respectively.
    % 
    % SYNOPSIS
    %   samples = WelRocklSamples('data', data);
    %   samples = WellRockSamples('generatorFn', generatorFn);
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
    %   'wells' - Names or indices of wells that the samples should be
    %             applied to.
    % 
    % NOTE:
    %   Either 'data' or 'generatorFn' must be provided.
    %
    % SEE ALSO:
    %   `RockSamples`, `WellSamples`, `BaseSamples`, `MRSTExample`, `BaseQoI`
    
    properties
        % Inherited properties only
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Applies the sample realization of the rock and well 
            % properties to a problem.
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
            
            % Update the well settings for the problem
            problem = setSample@WellSamples(samples, sampleData.well, problem);
            
            % Update the rock of the model in the problem.
            problem = setSample@RockSamples(samples, sampleData.rock, problem);
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