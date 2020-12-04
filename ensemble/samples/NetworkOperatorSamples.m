classdef NetworkOperatorSamples < OperatorSamples
    % Class that maps OperatorSamples to a network-type reduced reservoir 
    % model (GPSNet, INSIM, etc). It has the same properties as 
    % OperatorSamples, along with an additional connectionIndices property 
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
            
            assert(~isempty(samples.connectionIndices), ...
                'Not able to map sampleData to problem since no connection indices given');
     
            transmissibilityIsField = isfield(sampleData, 'T');
            porevolumeIsField = isfield(sampleData, 'pv');
            assert(transmissibilityIsField + porevolumeIsField  == numel(fields(sampleData)), ...
                'Only T and pv are currently supported in OperatorSamples');
            
            if transmissibilityIsField
                assert(numel(samples.connectionIndices.faces) == numel(sampleData.T), ...
                    'mismatch between number of connections and sample data size for T');
                for conn = 1:numel(samples.connectionIndices.faces)
                    problem.SimulatorSetup.model.operators.T(samples.connectionIndices.faces{conn}) = sampleData.T(conn);
                end
            end
            
            if porevolumeIsField
                assert(numel(samples.connectionIndices.cells) == numel(sampleData.pv), ...
                    'mismatch between number of connections and sample data size for pv');
                for conn = 1:numel(samples.connectionIndices.faces)
                    problem.SimulatorSetup.model.operators.pv(samples.connectionIndices.cells{conn}) = sampleData.pv(conn);
                end
            end
           
        end
        
        
    end
end