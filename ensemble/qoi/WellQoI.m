classdef WellQoI < BaseQoI
    % Class for extracting well related quantity of interests from
    % (an ensemble of) simulated problems.
    %
    % DESCRIPTION:
    %   This class is used within a MRSTEnsemble to extract, store, and 
    %   work with quantities of interest related to well solutions.
    %
    % SYNOPSIS
    %   qoi = WellQoI('wellNames', {'P1', 'P2'}, ...);
    %   qoi = WellQoI('wellIndices', [3 4], ...);
    %
    % PARAMETERS
    %   'wellNames' - Names of the wells that we are interested in
    %   'wellIndices' - Indices of the wells that we are interested in
    % 
    %    Either 'wellNames' or 'wellIndices' must be provided. If both are
    %    provided, 'wellIndices' is ignored.
    %
    % OPTIONAL PARAMETERS
    %   'fldname' - cell array of the well output fields that we are
    %               interested in. Valid values are well solution field
    %               names. Default: {'qOs'}
    %   'cumulative' - Quantity of interest is cumulative production data
    %                  per well. Default: false
    %   'total' - Quantity of interest is total production data (scalar per
    %             field and well). Default: false
    %   'combined' - Quantity of interest is combined production data over
    %                the given wells. Default: false
    %
    % 
    % SEE ALSO:
    %   `BaseQoI`, `ReservoirStateQoI`, `MRSTExample`, `BaseSamples`
    
    properties
        fldname      = {'qOs'} % Well output field names
        wellIndices            % Well indices
        wellNames              % Well names

        cumulative   = false % Cumulative production
        total        = false % Total production
        combined     = false % Combine all wells, or give per well values.
        
        historyMatchDtRange
        

        dt % Timestep sizes
        % TODO: Property timestep is removed
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
            %
            % SYNOPSIS:
            %   qoi = qoi.validateQoI(problem)
            %
            % PARAMETERS:
            %   problem - MRST problem that is representative for the
            %             (ensemble member) problems that will be used with
            %             this QoI instance. Especially, the well setups
            %             should be representative, as well as the
            %             production schedule.
            %
            % NOTE:
            %   When used in an MRSTEnsemble, this function is called by
            %   the ensemble constructor.
            
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
            if isempty(qoi.historyMatchDtRange)
                qoi.historyMatchDtRange = 1:numel(qoi.dt);
            end
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem)
            % Reads the well solutions from the given problem and extract
            % the relevant data.
            %
            % SYNOPSIS:
            %   u = qoi.computeQoI(problem)
            %
            % PARAMETERS:
            %   problem - The specific problem for which to compute the QoI
            
            % Read well solutions
            wellSols = reshape(problem.OutputHandlers.wellSols(:), [], 1);
            % Organize as matrix with dimensions:
            % (numTimesteps, numWells, numFields)
            uMatrix = getWellOutput(wellSols, qoi.fldname, qoi.wellNames);
            
            % Compute total or cumulative if requested
            dtProblem = getTimestepsFromProblem(problem);
            if qoi.total
                uMatrix = sum(uMatrix.*dtProblem,1);
            elseif qoi.cumulative
                uMatrix = cumsum(uMatrix.*dtProblem,1);
            end
            % Compute combined output from all wells if requested
            if qoi.combined && numel(qoi.wellNames) > 1
                uMatrix = sum(uMatrix, 2);
            end
            
            % If we have a time series and the problem was simulated with
            % different timesteps than the base problem, interpolate well
            % output onto base problem timesteps
            if numel(qoi.dt) ~= numel(dtProblem) && ~qoi.total
                uMatrix = qoi.interpolateWellOutput(dtProblem, uMatrix);
            end
            
            % Organizing u as a cell array per well of cell array per field
            % E.g, the value for field f at well w at time t will be in 
            % u{w}{f}(t)
            numFields = numel(qoi.fldname);
            numWells = numel(qoi.wellNames);
            if qoi.combined
                numWells = 1;
            end

            u = cell(numWells,1);
            [u{:}] = deal(cell(numFields,1));
            for w = 1:numWells
                for f = 1:numFields
                    u{w}{f} = uMatrix(:, w, f);
                end
            end

        end
        
        %-----------------------------------------------------------------%
        function plotQoI(qoi, ensemble, u, varargin) %#ok
            % Plot a single well QoI u in current figure.
            opt = struct('color'         , [0,0,0], ...
                         'lineWidth'     , 2      , ...
                         'alpha'         , 0.8    , ...
                         'isMean'        , true   , ...
                         'timescale'     , day    , ...
                         'labels'        , true   , ...
                         'title'         , true   , ...
                         'cellNo'        , 1      , ...
                         'subCellNo'     , 1      , ...
                         'isObservation' , false);
            [opt, extra] = merge_options(opt, varargin{:});
            
            color = opt.color; % Plot mean in distinct color
            if ~opt.isMean
                color = color.*(1-opt.alpha) + opt.alpha;
            end
            
            
            is_timeseries = true;
            if is_timeseries
                
                time = cumsum(qoi.dt)./opt.timescale;
                if opt.isObservation
                    plot(time(qoi.historyMatchDtRange), u, 'x', 'color', [0 0 0], extra{:});
                else
                    plot(time, u, 'color'    , color, ...
                                  'lineWidth', opt.lineWidth    , ...
                                   extra{:}         );
                end
                
                xlim([time(1), time(end)]);
                box on, grid on
                if opt.title
                    if qoi.combined
                        title(sprintf('Combined produced %s', qoi.fldname{opt.subCellNo}));
                    else
                        title(sprintf('%s for well %s', qoi.fldname{opt.subCellNo}, qoi.wellNames{opt.cellNo}));
                    end
                end
                if opt.labels
                    xlabel(sprintf('Time (%s)', formatTime(opt.timescale)));
                    ylabel(sprintf('%s', qoi.fldname{opt.subCellNo}));
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function h = figure(qoi, ensemble, varargin) %#ok
            % Create figure for plotting QoI
            h = figure(varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function h = plotEnsembleWellQoI(qoi, ensemble, h, varargin)
            % Plots well properties for the ensemble and ensemble mean.
            % Creates one figure per field, and organizes results from each
            % well in a subplot.
            
            opt = struct('timescale', day);
            [opt, extra] = merge_options(opt, varargin{:});
            
            numQoIs = numel(qoi.ResultHandler.getValidIds());
            
            num_fields = numel(qoi.fldname);
            num_wells  = numel(qoi.wellNames);
            if qoi.combined
                num_wells = 1;
            end
            
            num_wells_horizontal = ceil(sqrt(num_wells));
            num_wells_vertical   = ceil(num_wells/num_wells_horizontal);
            
            [mean_qoi, qois] = qoi.computeMean();
            
            t = cumsum(qoi.dt)./opt.timescale;
            for fld = 1:num_fields
                if num_fields == 1
                    if isnan(h)
                        h = figure;
                    else
                        set(0, 'CurrentFigure', h), clf(h)
                    end
                else
                    h{fld} = figure;
                end
                for w = 1:num_wells
                    subplot(num_wells_horizontal, num_wells_vertical, w);
                    hold on
                    for i = 1:numQoIs
                        plot(t, qois{i}{w}{fld}, 'color', [1 1 1]*0.6, extra{:});
                    end
                    plot(t, mean_qoi{w}{fld}, 'color', 'red', extra{:});
                    xlabel(sprintf('Time (%s)', formatTime(opt.timescale)));
                    ylabel(sprintf('%s', qoi.fldname{fld}))
                    title(strcat(qoi.fldname{fld}, " for well ", qoi.wellNames{w}));
                end
            end
            
        end
        
        
        %-----------------------------------------------------------------%
        function n = norm(qoi, u)
            % TODO: This function doesn't really work...
            
            if ~qoi.total
                n = cellfun(@(u) sum(u.*qoi.dt), u);
            else
                n = norm@BaseQoI(qoi, u);
            end
        end
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        
        function [obs, scaling] = getObservationAndScaling(qoi, varargin) 
            
            opt = struct('vectorize', true );
            [opt, extra] = merge_options(opt, varargin{:});
            
            obs = qoi.getObservationVector('vectorize', false);
            
            for w = 1:numel(qoi.wellNames)
                for f = 1:numel(qoi.fldname)
                    scaling{w}{f} = ones(size(obs{w}{f}(:)))*max(abs(obs{w}{f}(:)));
                end
            end
            
            if opt.vectorize
                obs = qoi.qoi2vector(obs);
                scaling = qoi.qoi2vector(scaling);
            end
    
        end
        
        %-----------------------------------------------------------------%
        function obs = getObservationVector(qoi, varargin)
            
            opt = struct('vectorize', true);
            [opt, extra] = merge_options(opt, varargin{:});
            
            % Check that the observation is valid
            assert(~isempty(qoi.observationResultHandler), ...
                'qoi.observationResultHandler is missing');
            assert(numel(qoi.observationResultHandler.getValidIds) > 0, ...
                'No available data in the observationResultHandler');
            
            obs = qoi.observationResultHandler{1};
            
            if ~qoi.combined
                assert(numel(obs) == numel(qoi.wellNames), ...
                    'The observation does not match the number of wells in QoI');
            else
                assert(numel(obs) == 1, ...
                    'The observation contains several wells, but the QoI is supposed to be combined values');
            end
            
            assert(numel(obs{1}) == numel(qoi.fldname), ...
                'The observation does not match the number of fldnames in QoI');
            
            assert(numel(obs{1}{1}) <= numel(qoi.dt), ...
                'The observation has too many timesteps to match the QoI class');
            
            obs = qoi.extractHistoryMatchingTimestep(obs);
            if opt.vectorize
                obs = qoi.qoi2vector(obs);
            end
        end
        
        %-----------------------------------------------------------------%
        function u = getQoIVector(qoi, seed)
            u = qoi.ResultHandler{seed};
            u = qoi.extractHistoryMatchingTimestep(u);
            u = qoi.qoi2vector(u);
        end
            
        %-----------------------------------------------------------------%
        function R = getObservationErrorCov(qoi)
            
            assert(~isempty(qoi.observationCov), ...
                'qoi.observationCov is missing');
            
            % Check how observationCov matches the observation
            u = qoi.getObservationVector('vectorize', false);
            u_vec = qoi.qoi2vector(u);
            numObs = size(u_vec,1);
            
            if isscalar(qoi.observationCov)
                R = speye(numObs)*qoi.observationCov;
            
            elseif isvector(qoi.observationCov)
                rdiag = [];
            
                if numel(qoi.observationCov) == numel(qoi.fldnames)*numel(qoi.wellNames)
                    for w = 1:numel(qoi.wellNames)
                        for f = 1:numel(qoi.fldnames)
                            covindex = w*(numel(qoi.fldnames-1)) + f;
                            rdiag = [rdiag ; repmat(qoi.observationCov(covindex), [numel(u{w}{f}), 1])];
                        end
                    end  
                elseif numel(qoi.observationCov) == numel(qoi.fldnames)
                    % Assume that there are different covariances for each
                    % fieldName
                    for w = 1:numel(qoi.wellNames)
                        for f = 1:numel(qoi.fldnames)
                            rdiag = [rdiag ; repmat(qoi.observationCov(f), [numel(u{w}{f}), 1])];
                        end
                    end    
                elseif numel(qoi.observationCov) == numObs
                    rdiag = qoi.observationCov;
                end
                assert(numel(rdiag) == numObs, ...
                    'Wrong number of diagonal elements after mapping the vector in qoi.observationCov');
                R = sparse([1:numObs], [1:numObs], rdiag, numObs, numObs); 
                
            else % observationCov is matrix
                assert(all(size(qoi.observationCov) == [numObs, numObs]), ...
                    'The matrix in qoi.observationCov is of wrong size');
                R = qoi.observationCov;
            end
        end
        
    end % methods
    
    
    methods (Access = protected)
        
        %-----------------------------------------------------------------%
        function wellOutput = interpolateWellOutput(qoi, dtProblem, wellOutput)
            % If wellOutput is given with time intervals dtProblem, this
            % function can be used to give wellOutputs at the timesteps 
            % found in qoi.dt
            
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
        function u = extractHistoryMatchingTimestep(qoi, u)
            % u is now u{well}{field}(time)
            for w = 1:numel(qoi.wellNames)
                for f = 1:numel(qoi.fldname)
                    u{w}{f} = u{w}{f}(qoi.historyMatchDtRange);
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function u = qoi2vector(qoi, u)
            % For multiple fields and wells, this vectorization will result in
            % u = [ (well1, field1), (well1, field2), (well2, field1), (well2,
            % field2) ...]
            u_tmp = [];
            for w = 1:numel(qoi.wellNames)
                for f = 1:numel(qoi.fldname)
                    u_tmp = cat(1, u_tmp, u{w}{f});
                end
            end
            u = u_tmp;
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

%{
#COPYRIGHT#
%}