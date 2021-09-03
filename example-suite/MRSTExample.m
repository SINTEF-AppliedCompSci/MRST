classdef MRSTExample
    
    properties
        % Core example properties
        name        % Example name
        description % One-line example description
        state0      % Initial state
        model       % Model
        schedule    % Simulation schedule
        options     % Options passed to example get (stored for reference)
        verbose     % Verbose flag
        % Properties for plotting
        figureProperties  % Figure propreties (defaults set by constructor)
        axisProperties    % Axis properties (defaults set by constructor)
        toolbarOptions    % Options that can be passed to plotToolbar
        visualizationGrid % Optional grid for plotting
    end
    
    methods
        %-----------------------------------------------------------------%
        function example = MRSTExample(name, varargin)
            if nargin == 0, return; end
            % Set example name
            example.name = lower(name);
            % Merge options
            opt = struct('plot', false, 'verbose', [], 'readFromDisk', true);
            [opt, varargin]  = merge_options(opt, varargin{:});
            example.verbose = opt.verbose;
            if isempty(opt.verbose), example.verbose = mrstVerbose(); end
            if example.verbose
                fprintf('Setting up %s example \n\n', example.name);
                timer = tic();
            end
            [example, extra] = merge_options(example, varargin{:});
            if opt.readFromDisk
                % Check if example is already stored on disk
                [example.description, ...
                 example.options    ] = feval(lower(name), extra{:});
                ex = example.load(false);
                if ~isempty(ex)
                    % We found it! Early return
                    ex0 = example;
                    example = ex;
                    if ex0.verbose
                        time = toc(timer);
                        fprintf(['Found a saved version of this example. '        , ...
                                'Example loaded in %s\n\n'], formatTimeRange(time));
                    end
                    return
                elseif example.verbose
                    % No luck, we need to set up from scratch
                    fprintf(['Did not find a saved version of this ', ...
                             'example, setting up from scratch \n\n']);
                end
            end
            [example.description, ... % Example description
             example.options    , ... % Options for reference
             example.state0     , ... % Initial state
             example.model      , ... % Model
             example.schedule   , ... % Simulation schedule
             plotOptions        ] ... % Plotting options
                = feval(lower(name), extra{:});
            if example.verbose
                time = toc(timer);
                fprintf('Example set up in %s\n\n', formatTimeRange(time));
            end
            % Set visualization grid if given
            [example, plotOptions] = merge_options(example, plotOptions{:});
            % Set figure properties
            [example.figureProperties, plotOptions]                ...
                = merge_options(example.defaultFigureProperties(), ...
                                plotOptions{:}                   );
            % Set axis properties
            [example.axisProperties, plotOptions]                ...
                = merge_options(example.defaultAxisProperties(), ...
                                plotOptions{:}                 );
            % Set toolbar options
            example.toolbarOptions ...
                = merge_options(example.defaultToolbarOptions(), ...
                                plotOptions{:}                 );
            names  = fieldnames(example.toolbarOptions);
            values = struct2cell(example.toolbarOptions);
            example.toolbarOptions = cell(1, 2*numel(names));
            example.toolbarOptions([1:2:end, 2:2:end]) = horzcat(names, values);
            if opt.plot
                % Plot example
                example.plot(example.model.rock, 'Name', 'rock'  , 'log10', true);
                example.plot(example.state0    , 'Name', 'state0'               );
            end
        end
                
        %-----------------------------------------------------------------%
        function props = defaultFigureProperties(example) %#ok
            % Set default figure properties
            props = struct();
            pos   = get(0, 'DefaultFigurePosition');
            props.Size = pos(3:4);
        end
        
        %-----------------------------------------------------------------%
        function h = figure(example, varargin)
            % Get example figure
            h     = figure(varargin{:});
            names = fieldnames(example.figureProperties);
            for i = 1:numel(names)
                if strcmpi(names{i}, 'Size')
                    screensize = get(groot, 'ScreenSize');
                    wsize      = example.figureProperties.(names{i});
                    wsize      = min([wsize; screensize(3:4) - [0 80]]);
                    pos        = get(h,'Position');
                    pos        = min([pos(1:2); screensize(3:4)-wsize-[0 80]]);
                    set(h, 'Position', [pos, wsize]);
                else
                    set(h, names{i}, example.figureProperties.(names{i}));
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function props = defaultAxisProperties(example)
            % Get default axis properties
            props = struct();
            % Get grid
            G = example.model.G;
            if ~isempty(example.visualizationGrid)
                G = example.visualizationGrid;
            end
            % Set XYZLim
            if isfield(G, 'nodes')
                x = G.nodes.coords;
            else
                x = G.faces.centroids;
            end
            xmin = min(x);
            xmax = max(x);
            xyz = 'XYZ';
            for i = 1:G.griddim
                props.([xyz(i), 'Lim']) = [xmin(i), xmax(i)];
            end
            if isfield(example.schedule.control(1), 'W')
                W = example.schedule.control(1).W;
                if ~isempty(W) && G.griddim == 3
                    dz = xmax(3) - xmin(3);
                    props.ZLim(1) = props.ZLim(1) - dz*0.2;
                end
            end
            % Set aspect ratio
            props.PlotBoxAspectRatio = [1,1,1];
            % Set viewpoint
            if G.griddim == 2
                props.View = [0, 90];
            elseif G.griddim == 3
                props.View = [-37.5, 30];
            end
            % Set axis projection
            if G.griddim == 2
                props.Projection = 'orthographic';
            else
                props.Projection = 'perspective';
            end
            props.Box = false;
        end
        
        %-----------------------------------------------------------------%
        function ax = setAxisProperties(example, ax)
            % Set properties to example figure axis
            names = fieldnames(example.axisProperties);
            for i = 1:numel(names)
                set(ax, names{i}, example.axisProperties.(names{i}));
            end
        end
        
        %-----------------------------------------------------------------%
        function opt = defaultToolbarOptions(example) %#ok
            % Get default toolbar options
            opt = struct('log10'        , false, ...
                         'exp'          , false, ...
                         'abs'          , false, ...
                         'filterzero'   , false, ...
                         'logical'      , false, ...
                         'outline'      , false, ...
                         'pauseTime'    , 0.150, ...
                         'lockCaxis'    , false, ...
                         'plot1d'       , false, ...
                         'plotmarkers'  , false, ...
                         'skipAugmented', false, ...
                         'field'        , 's:1', ...
                         'step_index'   , 1    , ...
                         'startplayback', false);
        end
        
        %-----------------------------------------------------------------%
        function h = plot(example, v, varargin)
            % Plot filed v on example grid using plotToolbar
            opt = struct('Name', '', 'plotWells', true, 'camlight', true);
            [opt, extra] = merge_options(opt, varargin{:});
            Name = example.name;
            if ~isempty(opt.Name)
                Name = [Name, ' ', opt.Name];
            end
            if nargin == 1
                v = example.model.rock;
            end
            h = example.figure('Name', Name);
            G = example.getVisualizationGrid();
            plotToolbar(G, v, example.toolbarOptions{:}, extra{:});
            if opt.plotWells
                example.plotWells();
            end
            example.setAxisProperties(gca);
            if opt.camlight && G.griddim == 3 ...
                    && ~all(example.axisProperties.View == [0,90])
                camlight;
            end
        end
        
        %-----------------------------------------------------------------%
        function G = getVisualizationGrid(example)
            if isempty(example.visualizationGrid)
                G = example.model.G;
            else
                G = example.visualizationGrid;
            end
        end
        
        %-----------------------------------------------------------------%
        function varargout = plotWells(example, varargin)
            if ~isfield(example.schedule.control(1), 'W')
                return;
            end
            W = example.schedule.control(1).W;
            if ~isempty(W)
                G = example.getVisualizationGrid();
                if G.griddim == 3
                    dz = example.axisProperties.ZLim(2) ...
                       - example.axisProperties.ZLim(1);
                    varargout = cell(nargout, 1);
                    [varargout{:}] = plotWell(G, W, 'color' , 'k'    , ...
                                       'height', 0.15*dz, ...
                                       varargin{:}      );
                else
                    hold on
                    x = G.cells.centroids(vertcat(W.cells), :);
                    h = plot(x(:,1), x(:,2)   , 'ok', ...
                             'markerFaceColor', 'w' , ...
                             'markerSize'     , 8   , ...
                             varargin{:}            );
                    hold off
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function hash = getExampleHash(example)
            % Get example options as string
            optnames    = fieldnames(example.options);
            optval      = struct2cell(example.options);
            keep        = cellfun(@(v) ischar(v) | isnumeric(v), optval);
            optval      = optval(keep);
            optnames    = optnames(keep);
            fix         = cellfun(@(v) ~ischar(v), optval);
            optval(fix) = cellfun(@num2str, optval(fix), 'UniformOutput', false);
            % Prepend with example name
            str = cell(1, 1 + 2*numel(optval));
            str([1, 2:2:end, 3:2:end]) = [{example.name}; optnames; optval];
            str = horzcat(str{:});
            try
                % Calculate hash value
                md = java.security.MessageDigest.getInstance('SHA-256');
                hash = sprintf('%2.2x', typecast(md.digest(uint8(str)), 'uint8')');
            catch
                % ... fallback using example string
                hash = str;
            end
        end
        
        %-----------------------------------------------------------------%
        function ok = save(example)
            % Get example hash value
            hash = example.getExampleHash();
            % Store in default data directory under 'example-suite'
            pth = fullfile(mrstDataDirectory(), 'example-suite');
            if ~exist(pth, 'dir')
                mkdir(pth)
            end
            save(fullfile(pth, hash), 'example', '-v7.3');
            ok = true;
        end
        
        %-----------------------------------------------------------------%
        function ex = load(example, throwError)
            % Attempt to load example
            if nargin < 2
                throwError = true;
            end
            % Example is stored with a filename equal to its hash value
            hash = example.getExampleHash();
            pth = fullfile(mrstDataDirectory(), 'example-suite');
            fn = [fullfile(pth, hash), '.mat'];
            ex = [];
            if isfile(fn)
                % The file exists - load it
                data = load(fn);
                ex   = data.example;
            elseif throwError
                % The file does no exists - throw an error
                error(['Could not find example '   , ...
                       'named %s with hash key %s'], example.name, hash);
            end
        end
        
        %-----------------------------------------------------------------%
        function problem = getPackedSimulationProblem(example, varargin)
            % Make packed problem with reasonable simulation setup
            opt = struct('LinearSolver'   , [], ...
                         'NonLinearSolver', []);
            [opt, extra] = merge_options(opt, varargin{:});
            has_ls  = ~isempty(opt.LinearSolver);    % Linear solver given
            has_nls = ~isempty(opt.NonLinearSolver); % Nonlinear solver given
            if has_nls
                % Check if nonlinear solver has non-default linear solver
                has_ls = has_ls || ~isa(opt.NonLinearSolver.LinearSolver, 'BackslashSolverAD');
            end
            if ~has_ls
                % Select apropriate linear solver
                rmodel = example.model;
                if isa(rmodel, 'WrapperModel')
                    rmodel = rmodel.getReservoirModel;
                end
                opt.LinearSolver = selectLinearSolverAD(rmodel);
            end
            if ~has_nls
                % Select default nonlinear solver
                opt.NonLinearSolver = NonLinearSolver('LinearSolver' , opt.LinearSolver, ...
                                                      'useRelaxation', true            );
            elseif ~has_ls
                % ... or assign linear solver
                opt.NonLinearSolver.LinearSolver = opt.LinearSolver;
            end
            % Pack problem
            problem = packSimulationProblem(                                   ...
                example.state0, example.model, example.schedule, example.name, ...
                'NonLinearSolver', opt.NonLinearSolver, extra{:}             );
        end
        
    end
    
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
