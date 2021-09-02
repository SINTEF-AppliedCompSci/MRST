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
#COPYRIGHT#
%}