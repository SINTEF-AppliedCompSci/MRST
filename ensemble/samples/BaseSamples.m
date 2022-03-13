classdef BaseSamples
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

        num = inf        % inf means that we can sample new ensemble
                         % members on the fly
        data             % precomputed sample data. Supported formats:
                         %      * Cell array of data samples
                         %      * Instance of ResultHandler class with
                         %        information about storage location and
                         %        names. See ResultHandler class and
                         %        <EXAMPLE> for details
        generatorFn      % Function for generating a stochastic sample
        processProblemFn % Function handle to for processing problem after
                         % sample has been set. Default: empty
    
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function samples = BaseSamples(varargin)
            % Class constructor.
            if nargin == 0, return; end
            samples = merge_options(samples, varargin{:});
            assert(xor(isempty(samples.data), isempty(samples.generatorFn))        , ...
                   'Please provide either a generatorFn or precomputed sample data');
            if isempty(samples.data)
                % Samples will be computed on the fly using generatorFn
            elseif iscell(samples.data)
                % We have a cell array of data samples
                samples.num = numel(samples.data);
            elseif isa(samples.data, 'ResultHandler')
                % Samples are stored to file - we are given a ResultHandler
                % for loading sample data
                samples.num = samples.data.numelData();
            else
                % Input is not on the correct format - throw an error
                error(['data must either be a cell array of samples ' , ...
                       'or an instance of the ''ResultHandler'' class']);
            end
        end
        
        %-----------------------------------------------------------------%
        function problem = getSampleProblem(samples, baseProblem, seed, sampleData)
            % Takes a realization of the sample according to the seed and
            % applies it to the baseProblem, thus creating a self standing
            % problem definition of an ensemble member.
            %
            % SYNOPSIS:
            %   problem = sample.getSampleProblem(baseProblem, seed)
            %
            % PARAMETERS:
            %   baseProblem - An MRST problem containing all but the
            %                 stochastic component
            %   seed        - Uniquely identify the sample realization that
            %                 will be applied. Is either the index of the
            %                 data (ensemble member id), or the seed for
            %                 the random generator.
            %
            % RETURNS:
            %   problem - A problem representing a single ensemble member
            
            % Based on the baseProblem, create a problem with sample
            % given by the seed
            if nargin < 4 || isempty(sampleData)
                sampleData = samples.getSample(seed, baseProblem);
            end
            problem    = samples.setSample(sampleData, baseProblem);
            if ~isempty(samples.processProblemFn)
                % Process problem
                problem = samples.processProblemFn(problem);
            end
            % Set output directory for the sample
            problem.OutputHandlers.wellSols.dataFolder = num2str(seed);
            problem.OutputHandlers.states.dataFolder   = num2str(seed);
            problem.OutputHandlers.reports.dataFolder  = num2str(seed);
            % Check if data directory exists, and make it if it does'nt
            dataDir = problem.OutputHandlers.states.getDataPath();
            if ~exist(dataDir, 'dir')
                mkdir(dataDir);
            end
        end
        
        %-----------------------------------------------------------------%
        function sampleData = getSample(samples, seed, problem)
            % Get a single sample realization based on the seed.
            %
            % SYNOPSIS:
            %   sampleData = sample.getSample(seed, problem)
            %
            % PARAMETERS:
            %   seed    - Uniquely identify the sample realization that
            %             will be applied. Is either the index of the data
            %             (ensemble member id), or the seed for the random
            %             generator.
            %   problem - The problem for which to apply the sample. Only
            %             used for generating data through a function on
            %             the fly.
            %
            % RETURNS:
            %   sampleData - The data of the sample realization 
            
            % Get a simple sample based on the seed, which is the seed for
            % the random generator or an index of data. Second input
            % argument ''problem'' is only needed if samples are computed
            % on the fly with generatorFn
            if ~isempty(samples.data)
                assert(seed >= 1 && seed <= samples.num, ...
                    ['seed must be in the range ' , ...
                     '[1, %d] = [1, samples.num]'], samples.num);
                sampleData = samples.data{seed};
            else
                sampleData = samples.generatorFn(problem, seed);
            end
        end
        
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem) %#ok
            % Applies the sample data realization to a problem.
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
            %
            % NOTE:
            %   This function must be implemented by each specific sample
            %   type.
            
            % Set sampleData to the problem. The specific implementation of
            % this function depends on the type of sample
            error('Template class not meant for direct use!');
        end
        
        function meanSample = getMeanSample(sample)
            % Computes the mean of the sample distributions and structure
            % the output in the same form as the relevant sample object.
            %
            % SYNOPSIS:
            %   meanSamples = sample.getMeanSamples()
            %
            % RETURNS:
            %   meanSamples - same object type as sample but with one num
            %   only.
            %
            % NOTE:
            %   This function ust be implemented by each specific sample
            %   type.
                        
            if isinf(sample.num)
                error('Function cannot be used when generating samples on the fly');
            elseif ~iscell(sample.data)
                error('Function currently only implemented for cell array samples');
            end

            fields = fieldnames(sample.data{1});
            meanData = sample.data{1};
            for f = 1:numel(fields)
                for i = 2:sample.num
                    meanData.(fields{f}) = meanData.(fields{f})+sample.data{i}.(fields{f});
                end
                meanData.(fields{f}) = meanData.(fields{f})/sample.num;
            end
            
            meanSample = sample;
            meanSample.data = {meanData};
            meanSample.num = 1;
            
        end
        
          function varSample = getVarianceSample(sample)
            % Computes the variance of the sample distributions and structure
            % the output in the same form as the relevant sample object.
            %
            % SYNOPSIS:
            %   meanSamples = sample.getMeanSamples()
            %
            % RETURNS:
            %   meanSamples - same object type as sample but with one num
            %   only.
            %
            % NOTE:
            %   This function ust be implemented by each specific sample
            %   type.
                        
            if isinf(sample.num)
                error('Function cannot be used when generating samples on the fly');
            elseif ~iscell(sample.data)
                error('Function currently only implemented for cell array samples');
            end
            
            meanSample = sample.getMeanSample();

            fields = fieldnames(sample.data{1});
            varData = sample.data{1};
            for f = 1:numel(fields)
                varData.(fields{f}) = (varData.(fields{f}) - meanSample.data{1}.(fields{f})).^2;
                for i = 2:sample.num
                    varData.(fields{f}) = varData.(fields{f}) + ...
                        (sample.data{i}.(fields{f}) - meanSample.data{1}.(fields{f})).^2;
                end
                varData.(fields{f}) = varData.(fields{f})./(sample.num-1);
            end
            
            varSample = sample;
            varSample.data = {varData};
            varSample.num = 1;
            
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
