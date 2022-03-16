classdef BaseSamplesHM
    % Template class for holding the stochastic samples relevant for all
    % ensemble members in an MRSTEnsemble
    %
    % NOTE:
    %   Not intended for direct use.
    %
    % DESCRIPTION:
    %   This class (and its super classes) is used within a MRSTEnsemble
    %   to organize the stochastic component/variables/parameters that
    %   separates the ensemble members.
    %   This base class defines the main API for interacting with such
    %   samples.
    %
    % SEE ALSO:
    %   `RockSamples`, `WellSamples`, `DeckSamples`, `MRSTExample`, `BaseQoI`
    
    properties

        %% Properties related to history matching
        transformSampleVectors = true % Boolean. Indicates whether the sample vectors
                                      % should be transformed to make them more
                                      % suitable for history matching
        
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function sampleSize = getSampleSize(samples)
            % Returns the size of a sample realization.
            % Note that this function does not work properly if a
            % generatorFn is used for generating samples on the fly, and
            % will in those cases return -1.
            
            if isempty(samples.data)
                sampleSize = -1;
                return;
            end
            
            function ssize = getElementSize(c)
                ssize = 0;
                if isnumeric(c)
                    ssize = numel(c);
                elseif iscell(c)
                    for i=1:numel(c)
                        ssize = ssize + getElementSize(c{i});
                    end
                elseif isstruct(c)
                    fn = fieldnames(c);
                    for i=1:numel(fn)
                        ssize = ssize + getElementSize(c.(fn{i}));
                    end
                else
                    assert(false, 'Found unexpected element');
                end
            end
            
            sampleSize = getElementSize(samples.data{1});  
        end
        
        
        function sampleVectors = getSampleVectors(samples)
            % Structure the samples in a matrix so that each column
            % consists of the sampled parameters of its corresponding
            % ensemble member.
            
            error('Template class not meant for direct use!');        
        end
        
        function samples = setSampleVectors(samples, newSampleVectors)
            % Replaces the samples based on the provided sample vectors.
            % The columns in newSampleVectors should represent parameter
            % values for individual ensemble members, and should be on the
            % same structure as what is provided by 'getSampleVectors()'
            
            error('Template class not meant for direct use!');
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
