classdef MRSTEnsemble < BaseEnsemble
    % Class that facilitates ensembles in MRST.
    %
    % SYNOPSIS:
    %   ensembles = Ensemble(baseCase, samples, qoi)
    %
    % DESCRIPTION:
    %   This class is used to organize, set-up, run and post-process 
    %   ensemble simulations in MRST. It can be used for a varity of 
    %   problems, such as uncertainty quantification, optimization under 
    %   uncertainty, and history matching.
    % 
    % REQUIRED PARAMETERS:
    %   baseCase - Instance of TestCase from the test-suite module, or a
    %              function name that generates a TestCase, that serves as
    %              a base problem for the ensemble. If will contain all 
    %              aspects and configurations that will be common for all
    %              ensemble members
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
    %                       'spmd' - Run ensemble in parallel using     
    %                           the Parallel Computing Toolbox
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
    %   `TestCase`, `BaseSamples`, `BaseQoI`
    
    properties
        qoi
        storeOutput
        figures = struct('progress', [], 'qoi', []);
    end
    
    methods
        %-----------------------------------------------------------------%
        function ensemble = MRSTEnsemble(baseCase, samples, qoi, varargin)
            % handle local options
            opt  = struct('reset',             false, ...
                          'prepareSimulation', true,  ...
                          'storeOutput',       false, ...
                          'exampleArgs',       {{}});
            [opt, other] = merge_options(opt, varargin{:});
            % Set example. This defines the base problem
            if isa(baseCase, 'TestCase')
                % Example given
                setup = baseCase;
            else
                % Example name given - set up example
                setup = TestCase(baseCase, opt.exampleArgs{:});
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
            for f = fieldnames(ensemble.figures)'
                if isgraphics(ensemble.figures.(f{1}))%isa(ensemble.figures.(f{1}), 'matlab.ui.Figure')
                    delete(ensemble.figures.(f{1}));
                end
                ensemble.figures.(f{1}) = [];
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
                    fprintf(['Simulation output for seed %d found on ', ...
                             'disk (skipping)\n'], seed               );
                end
                return
            end
            [problem, status] = ensemble.simulateEnsembleMember@BaseEnsemble(seed, varargin{:});
            if status.success
                % Compute QoI
                try
                    ensemble.qoi.getQoI(problem);
                catch me
                    % QoI computation failed, provide reason
                    warning(['Failed to compute QoI from simulation ', ...
                             'results for seed %d\n'], seed          );
                    status.success = false;
                    status.message = me;
                    ensemble.simulationStatus{seed} = status;
                end
            else
                % Simulation failed
                warning(['Could not compute QoI for seed %d due to ', ...
                         'failed simulation\n'], seed               );
            end
            % Clear output if requested
            if ~ensemble.storeOutput
                clearPackedSimulatorOutput(problem, 'prompt', false);
%                 ensemble.simulationStatus.resetData(seed);
            end
        end
        
        %-----------------------------------------------------------------%
        function flag = hasSimulationOutput(ensemble, range)
            flag = hasSimulationOutput@BaseEnsemble(ensemble, range);
            ids = ensemble.qoi.ResultHandler.getValidIds();
            flag(ismember(range, ids)) = true;
        end
        
        %-----------------------------------------------------------------%
        function flag = getSimulationStatus(ensemble, range)
            flag = getSimulationStatus@BaseEnsemble(ensemble, range);
            ids = ensemble.qoi.ResultHandler.getValidIds();
            flag(ismember(range, ids)) = 2;
        end
        
        %-----------------------------------------------------------------%
         function plotProgress(ensemble, range, plotProgress, plotIntermediateQoI, varargin)
            % Utility function for showing the progress of simulating a
            % range of ensemble members. Only available for
            % 'simulationStrategy' = 'background'.
            
            opt = struct('progressTitle', 'Simulating ensemble');
            [opt, extra] = merge_options(opt, varargin{:});
            n = 0;
            while true
                progress = ensemble.getEnsembleMemberProgress(range);
                pause(0.05);
                if plotProgress
                    ensemble.figures.progress = ...
                        plotEnsembleProgress(ensemble, progress, range, ensemble.figures.progress, ...
                                             'title', opt.progressTitle);
                    drawnow();
                    if ensemble.qoi.ResultHandler.numelData > n && plotIntermediateQoI
                        ensemble.figures.qoi = ensemble.qoi.plotEnsembleQoI(ensemble, ensemble.figures.qoi, extra{:});
                        n = ensemble.qoi.ResultHandler.numelData;
                        drawnow();
                    end
                end
                if all(isinf(progress) | isnan(progress)), break; end
            end
        end
        
 
        
        %-----------------------------------------------------------------%
        function h = plotQoI(ensemble, varargin)
            opt = struct('h', []);
            [opt, extra] = merge_options(opt, varargin{:});
            h = ensemble.qoi.plotEnsembleQoI(ensemble, opt.h, extra{:});
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
            flag = ensemble.getSimulationStatus(range);
            for i = 1:numel(range)
                if flag(i) == -1, progress(i) = nan; continue; end % Failed
                if flag(i) ==  2, progress(i) = inf; continue; end % Finished
                % Running - report fraction of completed steps
                dataDir     = fullfile(ensemble.directory(), num2str(range(i)));
                if exist(dataDir, 'dir')
                    files       = ls(dataDir);
                    progress(i) = numel(folderRegexp(files, 'state\d+\.mat', 'match'))/nsteps;
                end
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
