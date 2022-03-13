classdef OperatorSamplesHM < BaseSamplesHM & OperatorSamples
    % Class that combines ensemble configurations for both rock samples and
    % well samples.
    %
    % The data property for this class is expected to be a cell array with
    % structs containing the fields 'rock' and 'well', and each of them
    % should be compatible with the data properties expected in RockSamples
    % and WellSamples, respectively.

    properties
        pvScale = 100;
        TScale = 0; % If TScale is zero and transformSampleVectors is true,
                    % we transform the transmissibilities using log
                    % instead.
                    
        minPvValue = 0;
        minTValue = 0;
        maxPvValue = inf;
        maxTValue = inf;
    end
    
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
                        if samples.TScale > 0
                            sampleVectors(1:vectorSizeT, i) = ...
                                samples.data{i}.T(:)./samples.TScale;
                        else
                            % Logarithmic transform of transmissibility 
                            sampleVectors(1:vectorSizeT, i) = log(samples.data{i}.T(:));
                        end
                    else
                        sampleVectors(1:vectorSizeT, i) = samples.data{i}.T(:);
                    end
                end
                if porevolumeIsField
                    if samples.transformSampleVectors
                        sampleVectors(vectorSizeT+1:vectorSizeT+vectorSizePv, i) = ...
                            samples.data{i}.pv./samples.pvScale;
                    else
                        sampleVectors(vectorSizeT+1:vectorSizeT+vectorSizePv, i) = ...
                            samples.data{i}.pv;
                    end
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
                        if samples.TScale > 0
                            samples.data{i}.T(:) = ...
                                newSampleVectors(1:vectorSizeT, i).*samples.TScale;
                        else                     
                            samples.data{i}.T(:) = exp(newSampleVectors(1:vectorSizeT, i));
                        end
                    else
                        samples.data{i}.T(:) = newSampleVectors(1:vectorSizeT, i);
                    end
                    
                    % Check/ensure that the new values are acceptable
                    if any(samples.data{i}.T < samples.minTValue)
                        samples.data{i}.T(samples.data{i}.T < samples.minTValue) = samples.minTValue;
                    end
                    if any(samples.data{i}.T > samples.maxTValue)
                        samples.data{i}.T(samples.data{i}.T > samples.maxTValue) = samples.maxTValue;
                    end
                    if any(isnan(samples.data{i}.T))
                        error('Got NaN in transmissibility');
                    end
                end
                if porevolumeIsField
                    if samples.transformSampleVectors
                        samples.data{i}.pv(:) = ...
                            abs(newSampleVectors(vectorSizeT+1:vectorSizeT+vectorSizePv, i)).*samples.pvScale;
                    else
                        samples.data{i}.pv(:) = ...
                            abs(newSampleVectors(vectorSizeT+1:vectorSizeT+vectorSizePv, i));
                    end
                    
                    % Check/ensure that the new values are acceptable
                    if any(samples.data{i}.pv < samples.minPvValue)
                        samples.data{i}.pv(samples.data{i}.pv < samples.minPvValue) = samples.minPvValue;
                    end
                    if any(samples.data{i}.pv > samples.maxPvValue)
                        samples.data{i}.pv(samples.data{i}.pv > samples.maxPvValue) = samples.maxPvValue;
                    end
                    if any(isnan(samples.data{i}.pv))
                        error('Got NaN in pore volume');
                    end
                end
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
