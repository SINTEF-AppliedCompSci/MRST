classdef State0SamplesHM < BaseSamplesHM & State0Samples
    % Class that combines ensemble configurations for both initial states
    %
    % The data property for this class is expected to be a cell array with
    % structs containing the fields, and each of them
    % should be compatible with the data properties expected in RockSamples
    % and WellSamples, respectively.

    properties
        pvScale = 100;
    end
    
    methods
        
        
        %-----------------------------------------------------------------%
        function samples = State0SamplesHM(varargin)
            samples = samples@State0Samples(varargin{:});
        end
        
        
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function sampleVectors = getSampleVectors(samples)
            assert(~isempty(samples.data), ...
                'This function only work with pre-generated samples');
                        
            sizes = samples.getSizes();
            % sizes = [sizeInitSw, sizeInitSo, sizeInitSg, sizePressure, sizeFlux]
            
            totalVectorSize = sum(sizes);
            cumSizes = [0 cumsum(sizes)];
            fieldNames = {'initSw', 'initSo', 'initSg', 'pressure', 'flux'};
            
            sampleVectors = zeros(totalVectorSize, samples.num);
            
            for i = 1:samples.num
                
                for f = 1:5
                    if cumSizes(f) ~= cumSizes(f+1)
                        sampleVectors(cumSizes(f)+1:cumSizes(f+1), i ) = ...
                            samples.data{i}.(fieldNames{f})(:);
                    end
                end
            end
             
        end
        
        %-----------------------------------------------------------------%
        function samples = setSampleVectors(samples, newSampleVectors)
            
           sizes = samples.getSizes();
            % sizes = [sizeInitSw, sizeInitSo, sizeInitSg, sizePressure, sizeFlux]
            
            cumSizes = [0 cumsum(sizes)];
            fieldNames = {'initSw', 'initSo', 'initSg', 'pressure', 'flux'};
                        
            for i = 1:samples.num
                for f = 1:5
                    if cumSizes(f) ~= cumSizes(f+1)
                        samples.data{i}.(fieldNames{f})(:) = newSampleVectors(cumSizes(f)+1:cumSizes(f+1), i );
                        if (f < 5) % initSw, initSo, initSg, pressure
                            samples.data{i}.(fieldNames{f})(samples.data{i}.(fieldNames{f}) < 0) = 0;
                        end
                        if (f < 4) % initSw, initSo, initSg
                            samples.data{i}.(fieldNames{f})(samples.data{i}.(fieldNames{f}) > 1) = 1;
                        end
                    end
                end
            end
        end
        
        
    end
    
    methods (Access = protected)
        function sizes = getSizes(samples)
            sizeInitSw = 0;
            sizeInitSo = 0;
            sizeInitSg = 0;
            sizePressure = 0;
            sizeFlux = 0;
            
            if isfield(samples.data{1}, 'initSw')
                sizeInitSw = numel(samples.data{1}.initSw);
            end
            if isfield(samples.data{1}, 'initSo')
                sizeInitSo = numel(samples.data{1}.initSo);
            end
            if isfield(samples.data{1}, 'initSg')
                sizeInitSg = numel(samples.data{1}.initSg);
            end
            if isfield(samples.data{1}, 'pressure')
                sizePressure = numel(samples.data{1}.pressure);
            end
            if isfield(samples.data{1}, 'flux')
                sizeFlux = numel(samples.data{1}.flux);
            end
            
            sizes = [sizeInitSw, sizeInitSo, sizeInitSg, sizePressure, sizeFlux];
        end

    end
end