classdef WellSamples < BaseSamples
    % Class that holds an ensemble of sampled well-related properties and
    % methods to apply such samples to a given problem.
    
    properties
        % Inherits properties from BaseSamples only.
        wells;
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Returns a new problem that combines the input sampleData with
            % the input problem.
            
            % We assume that all parameters in the samplData are found in
            % the well control of the problem
            sampleFields = fieldnames(sampleData);
            numFields = numel(sampleFields);
            
            W = problem.SimulatorSetup.schedule.control.W;
            
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
            
            problem.SimulatorSetup.schedule.control.W = W;
            
        end
        
        
    end
end

