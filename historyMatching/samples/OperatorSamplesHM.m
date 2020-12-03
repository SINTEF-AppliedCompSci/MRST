classdef OperatorSamplesHM < BaseSamplesHM & OperatorSamples
    % Class that combines ensemble configurations for both rock samples and
    % well samples.
    %
    % The data property for this class is expected to be a cell array with
    % structs containing the fields 'rock' and 'well', and each of them
    % should be compatible with the data properties expected in RockSamples
    % and WellSamples, respectively.

    
    methods
        
        
        %-----------------------------------------------------------------%
        function samples = OperatorSamplesHM(varargin)
            samples = samples@OperatorSamples(varargin{:});
        end
        
        
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function sampleVectors = getSampleVectors(samples)
            assert(~isempty(samples.data), ...
                'This function only work with pre-generated samples');
            
            transmissibilityIsField = isfield(samples.data{1}, 'T');
            porevolumeIsField = isfield(samples.data{1}, 'pv');
            assert(transmissibilityIsField + porevolumeIsField  == numel(fields(samples.data{1})), ...
                'Only T and pv are currently supported in OperatorSamples');
            
            vectorSizeT = 0;
            vectorSizePv= 0;
            if transmissibilityIsField
                vectorSizeT = numel(samples.data{1}.T);
            end
            if porevolumeIsField
                vectorSizePv = numel(samples.data{1}.pv);
            end 
            
            sampleVectors = zeros(vectorSizeT+vectorSizePv, samples.num);
            
            for i = 1:samples.num
                
                if transmissibilityIsField
                    if samples.transformSampleVectors
                        % Logarithmic transform of transmissibility 
                        sampleVectors(1:vectorSizeT, i) = log(samples.data{i}.T(:));
                    else
                        sampleVectors(1:vectorSizeT, i) = samples.data{i}.T(:);
                    end
                end
                if porevolumeIsField
                    sampleVectors(vectorSizeT+1:vectorSizeT+vectorSizePv, i) = ...
                        samples.data{i}.pv;
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function samples = setSampleVectors(samples, newSampleVectors)
            
            transmissibilityIsField = isfield(samples.data{1}, 'T');
            porevolumeIsField = isfield(samples.data{1}, 'pv');
            assert(transmissibilityIsField + porevolumeIsField  == numel(fields(samples.data{1})), ...
                'Only T and pv are currently supported in OperatorSamples');
            
            vectorSizeT = 0;
            vectorSizePv= 0;
            if transmissibilityIsField
                vectorSizeT = numel(samples.data{1}.T);
            end
            if porevolumeIsField
                vectorSizePv = numel(samples.data{1}.pv);
            end 
            
            assert(size(newSampleVectors, 1) == vectorSizeT+vectorSizePv, ...
                'number of rows of new sample does not match old sample size');
            
            for i = 1:samples.num
                
                if transmissibilityIsField
                    if samples.transformSampleVectors
                        samples.data{i}.T(:) = exp(newSampleVectors(1:vectorSizeT, i));
                    else
                        samples.data{i}.T(:) = newSampleVectors(1:vectorSizeT, i);
                    end
                end
                if porevolumeIsField
                    samples.data{i}.pv(:) = abs(newSampleVectors(vectorSizeT+1:vectorSizeT+vectorSizePv, i));
                end
            end
        end
        
        
    end
end