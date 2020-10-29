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
        
        
             
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function sampleVectors = getSampleVectors(samples)

            assert(~isempty(samples.data), ...
                'This function only work with pre-generated samples');
            
            
            wellIsField = isfield(samples.data{1}, 'well');
            if wellIsField
                sampleFields = fieldnames(samples.data{1}.well);
                numWells = numel(samples.data{1}.well.(sampleFields{1}));
            else
                sampleFields = fieldnames(samples.data{1});
                numWells = numel(samples.data{1}.(sampleFields{1}));
            end
            numFields = numel(sampleFields);
            
            sampleVectors = zeros(numFields*numWells, samples.num);
            
            for i = 1:samples.num
                for f = 1:numFields
                    field = sampleFields{f};
                    startIndex = (f-1)*numWells + 1;
                    endIndex   = f*numWells;
                    
                    if wellIsField
                        fieldData = samples.data{i}.well.(field);
                    else
                        fieldData = samples.data{i}.(field);
                    end
                    
                    if strcmp(field, 'WI') && samples.transformSampleVectors
                        % Logarithmic transform of well production indices
                         sampleVectors(startIndex:endIndex, i)= log(fieldData);
                    else
                        sampleVectors(startIndex:endIndex, i) = fieldData;
                    end
                end
            end
        end
        
        function samples = setSampleVectors(samples, newSampleVectors)
            
            assert(size(newSampleVectors, 2) == samples.num, ...
                'number of columns of new samples does not match ensemble size');
            
            wellIsField = isfield(samples.data{1}, 'well');
            if wellIsField
                sampleFields = fieldnames(samples.data{1}.well);
                numWells = numel(samples.data{1}.well.(sampleFields{1}));
            else
                sampleFields = fieldnames(samples.data{1});
                numWells = numel(samples.data{1}.(sampleFields{1}));
            end
            numFields = numel(sampleFields);
            
            assert(size(newSampleVectors, 1) == numWells*numFields, ...
                'number of rows of new sample does not match old sample size');
            
            for i = 1:samples.num
                for f = 1:numFields
                    field = sampleFields{f};
                    startIndex = (f-1)*numWells + 1;
                    endIndex   = f*numWells;

                    if strcmp(field, 'WI') && samples.transformSampleVectors
                        % Logarithmic transform of well production indices
                        fieldData = exp(newSampleVectors(startIndex:endIndex, i));
                    else
                        fieldData = newSampleVectors(startIndex:endIndex, i);
                    end
                    
                    if wellIsField
                        samples.data{i}.well.(field)(:) = fieldData;
                    else
                        samples.data{i}.(field)(:) = fieldData;
                    end
                end
            end           
        end
        
        
        
    end
end

