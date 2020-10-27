classdef WellQoI < BaseQoI
    
    properties
        fldname      = {'qOs'} % Well output field names
        wellIndices            % Well indices
        wellNames              % Well names

        cumulative   = false % Cumulative production
        total        = false % Total production
        combined     = false % Combine all wells, or give per well values
        
        dt % Timestep sizes
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function qoi = WellQoI(varargin)
            qoi = qoi@BaseQoI();
            qoi = merge_options(qoi, varargin{:});
            % Check input
            assert(xor(isempty(qoi.wellIndices), isempty(qoi.wellNames)), ...
                'Please provide either well indices or well names');
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem)
            % Check that the configs that are inserted to the constructor
            % makes sense for the base problem for the ensemble, and
            % updates remaining fields.
            qoi = validateQoI@BaseQoI(qoi, problem);
            % Either, we provide the indices of the chosen wells, or their
            % names. In either case, we find the other.
            names = {problem.SimulatorSetup.schedule.control(1).W.name};
            if ~isempty(qoi.wellNames)
                % Find well indices
                ix = cellfun(@(wn) find(strcmpi(wn, names)), ...
                                    qoi.wellNames, 'UniformOutput', false);
                % Verify that all wells exist in schedule
                found = ~cellfun(@isempty, ix);
                assert(all(found), 'Did not find well %s in schedule', ...
                                                qoi.wellNames{~found});
                % Assign well indices
                qoi.wellIndices = cell2mat(ix);
            elseif ~isempty(qoi.wellIndices)
                if strcmpi(qoi.wellIndices, ':')
                    % All wells requested
                    qoi.wellIndices = 1:numel(names);
                elseif islogical(qoi.wellIndices)
                    % Logical mask - find corresponding indices
                    qoi.wellIndices = find(qoi.wellIndices);
                end
                % Verify well indices
                assert(min(qoi.wellIndices) >= 1 && max(qoi.wellIndices) <= numel(names), ...
                    'Well indices must be in the range [1, numel(W)] = [1, %d]', numel(names));
                % Find well names
                qoi.wellNames = names(qoi.wellIndices);
            else
                % We should never get here
                error('Please provide either well names or well indices');
            end
            % Get timesteps
            qoi.dt = problem.SimulatorSetup.schedule.step.val;
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem)
            % Reads the well solutions from the given problem and extract
            % the relevant data, which is then saved to file.
            % Read well solutions
            wellSols = reshape(problem.OutputHandlers.wellSols(:), [], 1);
            % Organize as matrix with dimensions:
            % (numTimesteps, numWells, numFields)
            u = getWellOutput(wellSols, qoi.fldname, qoi.wellNames);
            % Compute total or cumulative if requested
            dtProblem = getTimestepsFromProblem(problem);
            if qoi.total
                u = sum(u.*dtProblem,1);
            elseif qoi.cumulative
                u = cumsum(u.*dtProblem,1);
            end
            % Compute combined output from all wells if requested
            if qoi.combined && numel(qoi.wellNames) > 1
                u = sum(u, 2);
            end
            % If we have a time series and the problem was simulated with
            % different timesteps than the base problem, interpolate well
            % output onto base problem timesteps
            if numel(qoi.dt) ~= numel(dtProblem) && ~qoi.total
                u = qoi.interpolateWellOutput(dtProblem, u);
            end
            u = {u};
        end
        
        %-----------------------------------------------------------------%
        function wellOutput = interpolateWellOutput(qoi, dtProblem, wellOutput)
            % Get base problem and current problem timestamps
            timeProblem = [0; cumsum(dtProblem)];
            timeBase    = [0; cumsum(qoi.dt)   ];
            % Compute overlap
            t   = bsxfun(@min, timeBase(2:end)  , timeProblem(2:end)'  );
            t0  = bsxfun(@max, timeBase(1:end-1), timeProblem(1:end-1)');
            psi = max(t - t0, 0)./qoi.dt;
            % Integrate
            wellOutput = psi*wellOutput;
        end
            
        %-----------------------------------------------------------------%
        function plotQoI(qoi, ensemble, u, varargin) %#ok
            % Plot a single well QoI u in current figure.
            opt = struct('isMean'   , true, ...
                         'timescale', day , ...
                         'labels'   , true);
            [opt, extra] = merge_options(opt, varargin{:});
            color = [1,1,1]*0.8*(1-opt.isMean); % Plot mean in distinct color
            is_timeseries = true;
            if is_timeseries
                time = cumsum(qoi.dt)./opt.timescale;
                if qoi.cumulative
                    u = cumsum(u.*qoi.timesteps).*opt.timescale;
                else
                    u = u.*opt.timescale;
                end
                plot(time, u, 'lineWidth', 2, 'color', color, extra{:});
                xlim([time(1), time(end)]);
                if opt.labels
                    xlabel(sprintf('Time (%s)', formatTime(opt.timescale)));
                    ylabel(sprintf('%s', qoi.fldname{1}));
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function n = norm(qoi, u)
            if ~qoi.total
                n = cellfun(@(u) sum(u.*qoi.dt), u);
            else
                n = norm@BaseQoI(qoi, u);
            end
        end
        
    end
end

%-------------------------------------------------------------------------%
% Helpers
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function dt = getTimestepsFromProblem(problem)
    reports = reshape(problem.OutputHandlers.reports(:), [], 1);
    dt = zeros(numel(reports),1);
    for i = 1:numel(reports)
        dtstep = 0;
        stepReports = reports{i}.StepReports;
        for j = 1:numel(stepReports)
            dtstep = dtstep + stepReports{j}.Timestep;
        end
        dt(i) = dtstep;
    end     
end

%-------------------------------------------------------------------------%
function str = formatTime(timescale)
    t   = [second, minute, day, year];
    str = {'Second', 'Minute', 'Day', 'Year'};
    [d, ix] = min(abs(t - timescale));
    if d == 0
        str = str{ix};
    else
        warning('Non-standard timescale, using ''Second''');
        str = str{1};
    end
end
