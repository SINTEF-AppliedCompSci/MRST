classdef WellSamples < BaseSamples
    % Class that holds an ensemble of sampled well-related properties and
    % methods to apply such samples to a given problem. 
    % For use with an MRSTEnsemble.
    %
    % DESCRIPTION:
    %   This class is an extention of `BaseSamples` and is used within a 
    %   MRSTEnsemble to organize stochastic well properties. These well 
    %   properties can either be pre-computed or generated on the fly.
    % 
    % SYNOPSIS
    %   samples = WellSamples('data', data);
    %   samples = WellSamples('generatorFn', generatorFn);
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
    %   Each sample should consist of a struct with the fields that are
    %   also found in the typical well schedule struct in MRST.
    %
    % SEE ALSO:
    %   `RockSamples`, `DeckSamples`, `BaseSamples`, `MRSTExample`, `BaseQoI`
    
    properties
        % Inherits properties from BaseSamples only.
        wells;
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Applies the sample realization of the well properties to a 
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
            
            % We assume that all parameters in the samplData are found in
            % the well control of the problem
            sampleFields = fieldnames(sampleData);
            numFields = numel(sampleFields);
            
            numSchedules = numel(problem.SimulatorSetup.schedule.control);
            for s = 1:numSchedules
                W = problem.SimulatorSetup.schedule.control(s).W;

                % If the wells property is not defined, we assume that there
                % will be sampled values for all wells.
                if isempty(samples.wells)
                    samples.wells = 1:numel(W);
                end

                % Map field from sampleData to well
                for i = 1:numFields
                    assert(isfield(W, sampleFields{i}), ...
                        "Sample contained invalid well field");

                    % Map values to the correct well
                    for w = 1:numel(samples.wells)
                        W(samples.wells(w)).(sampleFields{i}) = sampleData.(sampleFields{i})(w);
                    end
                end

                % TODO:
                % Do some sanity checking. E.g., if W.val is changed, do we
                % also specify the same W.type?
                %
                % What if we only have samples for selected Wells? E.g., we
                % want to have an uncertain injection rate, but always produce
                % at the same pressure. Makes sense or not?

                problem.SimulatorSetup.schedule.control(s).W = W;
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
