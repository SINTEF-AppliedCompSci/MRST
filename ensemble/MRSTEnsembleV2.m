classdef MRSTEnsembleV2 < BaseEnsemble
    % Class that facilitates ensembles in MRST.
    %
    % SYNOPSIS:
    %   ensembles = Ensemble(mrstExample, samples, qoi)
    %
    % DESCRIPTION:
    %   This class is used to organize, set-up, run and post-process 
    %   ensemble simulations in MRST. It can be used for a varity of 
    %   problems, such as uncertainty quantification, optimization under 
    %   uncertainty, and history matching.
    % 
    % REQUIRED PARAMETERS:
    %   mrstExample - Instance of MRSTExample, or a function name that 
    %                 generates a MRSTExample, that serves as a base 
    %                 problem for the ensemble. If will contain all aspects
    %                 and configurations that will be common for all
    %                 ensemble members
    %
    %   samples - Sample object, typically of a superclass of BaseSamples.
    %             It represents stochastic parameters or configurations
    %             that are specific for each individual ensemble member.
    %             This object will also have functionality to map sample
    %             realizations to the base problem and thereby create the
    %             ensemble member simulations.
    %
    %   qoi - Quantity of interest object, typically a superclass of
    %         BaseQoI. It contains a post-processing mapping of a simulated
    %         problem (ensemble member) to a quantity of interest, and is
    %         therefore used to store and read relevant simulation results
    %         for the ensemble.
    %
    % OPTIONAL PARAMETERS:
    %   'directory' - Path to where we store results when simulating
    %                 the ensemble. Automatically generated if not provided
    %                 Default:
    %                 mrstOutputDirectory()/ensemble/mrstExample.name
    %   'reset' - Boolean (default: false). Will delete any old results
    %             existing in the 'directory' so that any simulation
    %             results related to the new ensemble will have to be
    %             regenerated.
    %   'storeOutput' - Boolean. If false, the simulation results (states,
    %                   wellSols, etc) of individual ensemble members will 
    %                   be deleted after the quantity of interest is 
    %                   computed and stored.
    %                   Default: false
    %   'solve' - Function handle for simulating/running an ensemble
    %             member.
    %             Default: @simulatePackedProblem(problem)
    %   'simulationStrategy' - String to control how to run the ensemble.
    %                          Acceptable values:
    %                       'serial'     - No parallelization (default)
    %                       'parallel'   - Run ensemble in parallel using 
    %                           the Parallel Computing Toolbox,
    %                       'background' - Run ensemble members by spawning
    %                           new matlab sessions in the background.
    %   'maxWorkers' - Maximum number of parallel workers to use for
    %                  processing the ensemble members (default is system
    %                  dependent)
    %   'verbose' - Boolean (default: true). Print some extra info
    %   'verboseSimulation' - Boolean (default: false). Print the output of
    %                         each ensemble member simulation   
    %   'matlabBinary' - If 'simulationStrategy' is set to 'background',
    %                    it is possible to provide a specific location of
    %                    the matlab binary (typically used if you have
    %                    several matlab installations on your system).
    %   'backgroundEvalFn' - Function name for running ensemble members in
    %                        the background. 
    %                        Default: @simulateEnsembleMembersStandalone
    %
    %   If the 'mrstExample' parameter is a string, it is possible to add
    %   extra parameters which will then be passed on to the MRSTExample
    %   class.
    %
    % RETURNS:
    %   Class instance.
    %
    % SEE ALSO:
    %   `MRSTExample`, `BaseSamples`, `BaseQoI`
    
    properties
        qoi
        storeOutput
    end
    
    methods
        %-----------------------------------------------------------------%
        function ensemble = MRSTEnsembleV2(mrstExample, samples, qoi, varargin)
            % handle local options
            opt  = struct('reset',             false, ...
                          'prepareSimulation', true,  ...
                          'storeOutput',       false, ...
                          'exampleArgs',       {{}});
            [opt, other] = merge_options(opt, varargin{:});
            % Set example. This defines the base problem
            if isa(mrstExample, 'MRSTExample')
                % Example given
                setup = mrstExample;
            else
                % Example name given - set up example
                setup = MRSTExample(mrstExample, opt.exampleArgs{:});
            end
            
            % Call BaseEnsemble
            ensemble = ensemble@BaseEnsemble(samples, other{:}, ...
                                        'reset',             false, ...  
                                        'prepareSimulation', false, ...
                                        'setup',             setup);
            
            ensemble.storeOutput = opt.storeOutput;
            
            % Validate qoi
            baseProblem  = ensemble.getBaseProblem();
            ensemble.qoi = qoi.validateQoI(baseProblem);
            % Delete existing results if requested
            if opt.reset
                ensemble.reset('prompt', false, ...
                               'prepareSimulation', false);
            end
            
            % Prepare ensemble
            if opt.prepareSimulation
                ensemble.prepareEnsembleSimulation();
            end
        end
    
        %-----------------------------------------------------------------%
        function reset(ensemble, varargin)
            % Deletes any old results so that we can start the ensemble
            % simulation fra scratch.
            opt = struct('prepareSimulation', true);
            [opt, extra] = merge_options(opt, varargin{:});
            % call reset at BaseEnsemble (don't do prepareSimulation)
            flag = ensemble.reset@BaseEnsemble(extra{:}, 'prepareSimulation', false);
            if flag
                % Delete QoIs
                ensemble.qoi.ResultHandler.resetData();
                % Prepare ensemble
                if opt.prepareSimulation
                    ensemble.prepareEnsembleSimulation();
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMember(ensemble, seed, varargin)
            % Run simulation for a specific ensemble member.
            %
            % SYNOPSIS:
            %   ensemble.simulateEnsembleMember(seed);
            %
            % PARAMETERS:
            %   seed - Integer specifying which ensemble member to run.
            %            
            if ensemble.qoi.isComputed(seed)
                % QoI is already computed - nothing to do here!
                if ensemble.verbose
                    disp(strcat("Simulation for ", num2str(seed), " found on disk"));
                end
                return;
            end
            problem = ensemble.simulateEnsembleMember@BaseEnsemble(seed, varargin{:});
            if ensemble.getSimulationStatus(seed) > 0
                % Compute QoI
                ensemble.qoi.getQoI(problem);
            else
                warning('Could not compute QOI for seed %d due to failed simulation', ...
                        seed);
            end
            % Clear output if requested
            if ~ensemble.storeOutput
                clearPackedSimulatorOutput(problem, 'prompt', false);
                ensemble.simulationStatus.resetData(seed);
            end
        end
        
        %-----------------------------------------------------------------%
        function h = plotQoI(ensemble, varargin)
            % Creates plot(s) of the quantity of interest for all simulated
            % ensemble members
            %
            % SYNOPSIS:
            %   h = ensemble.plotQoI();
            %
            % OPTIONAL PARAMETERS:
            %   'h' - Figure handle 
            opt = struct('h', []);
            [opt, extra] = merge_options(opt, varargin{:});
            ids    = ensemble.qoi.ResultHandler.getValidIds();
            sample = ensemble.qoi.ResultHandler{ids(1)};
            if isscalar(sample{1})
                if ~isempty(opt.h), clf(opt.h); end
                n = min(ceil(numel(ids)/3), 10);
                h = ensemble.qoi.plotQoIHistogram(opt.h, ...
                                                  'edges'      , n   , ...
                                                  'includeMean', true, ...
                                                  extra{:}           );
            else
                h = ensemble.qoi.plotEnsembleQoI(ensemble, opt.h, extra{:});
            end
            drawnow
        end
        
        %-----------------------------------------------------------------%
        function plotProgress(ensemble, range)
            % Utility function for showing the progress of simulating a
            % range of ensemble members. Only available for
            % 'simulationStrategy' = 'background'.
            [h_progress, h_qoi] = deal([]);
            n = 0;
            while true
                pause(0.1);
                progress = ensemble.getEnsembleMemberProgress(range);
                h_progress = plotEnsembleProgress(ensemble, progress, range, h_progress);
                if ensemble.qoi.ResultHandler.numelData > n
                    h_qoi = ensemble.plotQoI('h', h_qoi);
                    n = ensemble.qoi.ResultHandler.numelData;
                end
                drawnow
                if all(isinf(progress)), break; end
            end
        end
        
    end % methods
    
    methods (Access = protected)
        
        %-----------------------------------------------------------------%
        function progress = getEnsembleMemberProgress(ensemble, range)
            % Utility function for monitoring the progression of an
            % ensemble member that is being run right now.
            if nargin < 2, range = ensemble.num; end
            progress = zeros(numel(range),1);
            nsteps   = numel(ensemble.setup.schedule.step.val);
            for i = 1:numel(range)
                if exist(fullfile(ensemble.directory(), ...
                        [ensemble.qoi.ResultHandler.dataPrefix, num2str(range(i)), '.mat']), 'file')
                    progress(i) = inf;
                    continue
                end
                dataDir = fullfile(ensemble.directory(), num2str(range(i)));
                if ~exist(dataDir, 'dir')
                    continue;
                end
                files = ls(dataDir);
                progress(i) = numel(folderRegexp(files, 'state\d+\.mat', 'match'))/nsteps;
            end
        end
        
        
                
    end
end

%% Helpers
function matches = folderRegexp(list, expression, outputFormat)
    if size(list, 1) > 1
        % Windows behavior
        % Pad with a space at the end and reshape into a single
        % long string.
        pad = repmat(' ', size(list, 1), 1);
        list = reshape([list, pad]', 1, []);
    end
    matches = regexp(list, expression, outputFormat);
end

function name = getExampleName(example)
% get ensemble name proir to calling BaseEnsemble constructor
if isa(example, 'MRSTExample')
    name = example.name;
elseif ischar(example)
    name = example;
else
    error('Unexepected class ''%s'' of input ''mrstExample''.', class(example));
end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
