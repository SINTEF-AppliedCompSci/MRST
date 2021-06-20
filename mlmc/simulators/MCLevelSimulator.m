classdef MCLevelSimulator < MCSimulator
    
    properties
        levels
        levelNo
    end
    
    methods
        function mcl = MCLevelSimulator(setup, samples, qoi, levels, varargin)
            mcl = mcl@MCSimulator(setup, samples, qoi, varargin{:});
            mcl.levels = cell(numel(levels), 1);
            mcl.levelNo = levels{end}.levelNo;
            % Level QoIs will be stored in a subfolder `level-<levelNo>`
            mcl.directory = fullfile(mcl.getDataPath(), ...
                                         ['level-', num2str(mcl.levelNo)]);
            % Individual QoIs for each level will be stored in a subfolder
            % `level-<levelNo>/level-<levelNoA>-<levelNoB>/`
            dirName = @(i) fullfile(['level-', num2str(levels{end}.levelNo), ...
                                     '-', num2str(levels{i}.levelNo)     ]);
            for i = 1:mcl.numLevels
                setup    = levels{i}.constructor(mcl.setup);        % Get setup
                dataPath = fullfile(mcl.getDataPath(), dirName(i)); % Set data path
                % Set up ensemble for the level
                mcl.levels{i} = MRSTEnsemble(setup, samples, qoi      , ...
                          'solve'             , levels{i}.solver      , ...
                          'simulationStrategy', mcl.simulationStrategy, ...
                          'directory'         , dataPath              );
            end
            levelQoIs = cellfun(@(level) level.qoi, mcl.levels, 'UniformOutput', false);
            mcl.qoi = MCLevelQoI(levelQoIs);
            mcl.qoi = mcl.qoi.validateQoI(mcl.getBaseProblem());
        end
        
        %-----------------------------------------------------------------%
        function simulateEnsembleMember(mcl, seed, varargin)
            % Run simulation on MLMC level for a specific ensemble member.
            %
            % SYNOPSIS:
            %   ensemble.simulateEnsembleMember(seed);
            %
            % PARAMETERS:
            %   seed - Integer specifying which ensemble member to run.
            %
            % Get sample on finest level
            problem = mcl.levels{end}.getBaseProblem();
            sample = mcl.levels{end}.samples.getSample(seed, problem);
            % Simulate sample on all levels
            for i = 1:mcl.numLevels
                mcl.levels{i}.simulateEnsembleMember(seed, 'sample', sample);
            end
            status = mcl.levels{end}.simulationStatus{seed};
            if status.success
                try
                    mcl.qoi.getQoI(seed);
                catch me
                    % QoI computation failed, provide reason
                    warning(['Failed to compute QoI from simulation ', ...
                             'results for seed %d\n'], seed          );
                    status.success = false;
                    status.message = me;
                end
            end
            mcl.simulationStatus{seed} = status;
        end
        
        %-----------------------------------------------------------------%
        function n = numLevels(mcl)
            n = numel(mcl.levels);
        end
        
        %-----------------------------------------------------------------%
        function reset(mcl, varargin)
            reset@MCSimulator(mcl, varargin{:});
            for i = 1:mcl.numLevels
                mcl.levels{i}.reset(varargin{:});
            end
        end
        
    end
    
    methods (Access = protected)
        %-----------------------------------------------------------------%
        function progress = getEnsembleMemberProgress(mcl, range)
            % Utility function for monitoring the progression of an
            % ensemble member that is being run right now.
            progress = zeros(numel(range), mcl.numLevels);
            for i = 1:mcl.numLevels
                progress(:, i) = mcl.levels{i}.getEnsembleMemberProgress(range);
            end
            progress(progress == inf) = 1;
            progress = mean(progress, 2);
            progress(progress == 1) = inf;
        end

    end
    
end