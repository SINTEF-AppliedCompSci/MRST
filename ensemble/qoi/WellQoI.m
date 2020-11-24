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
            opt = struct('color'     , [0,0,0], ...
                         'alpha'     , 0.8    , ...
                         'isMean'    , true   , ...
                         'timescale' , day    , ...
                         'labels'    , true   , ...
                         'title'     , true   , ...
                         'cellNo'    , 1      , ...
                         'subCellNo' , 1      );
            [opt, extra] = merge_options(opt, varargin{:});
            
            color = opt.color; % Plot mean in distinct color
            if ~opt.isMean
                color = color.*(1-opt.alpha) + opt.alpha;
            end
            is_timeseries = true;
            if is_timeseries
                time = cumsum(qoi.dt)./opt.timescale;
                plot(time, u, 'color'    , color, ...
                              'lineWidth', 2    , ...
                               extra{:}         );
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
            if ~qoi.total
                n = cellfun(@(u) sum(u.*qoi.dt), u);
            else
                n = norm@BaseQoI(qoi, u);
            end
        end
        
    end
    
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