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
        function test = MRSTExample(name, varargin)
            if nargin == 0, return; end
            % Set example name
            test.name = lower(name);
            % Merge options
            opt = struct('plot', false, 'verbose', [], 'readFromDisk', true);
            [opt, varargin]  = merge_options(opt, varargin{:});
            test.verbose = opt.verbose;
            if isempty(opt.verbose), test.verbose = mrstVerbose(); end
            if test.verbose
                fprintf('Setting up %s example \n\n', test.name);
                timer = tic();
            end
            [test, extra] = merge_options(test, varargin{:});
            if opt.readFromDisk
                % Get example options and description
                setup = feval(lower(name), true, extra{:});
                test.options     = setup.options;
                test.description = setup.description;
                % Check if example exists on disk and load it
                tc = test.load(false);
                if ~isempty(tc)
                    % We found it! Early return
                    test = tc;
                    if test.verbose
                        time = toc(timer);
                        fprintf('Example loaded in %s\n\n', formatTimeRange(time));
                    end
                    return
                elseif test.verbose
                    % No luck, we need to set up from scratch
                    fprintf(['Did not find a saved version of this ', ...
                             'example, setting up from scratch \n\n']);
                end
            end
            % Set up test case and set properties
            setup = feval(lower(name), extra{:});
            test.description = setup.description;
            test.state0      = setup.state0;
            test.model       = setup.model;
            test.schedule    = setup.schedule;
            test.options     = setup.options;
            if test.verbose
                time = toc(timer);
                fprintf('Example set up in %s\n\n', formatTimeRange(time));
            end
            % Set visualization grid if given
            [test, plotOptions] = merge_options(test, setup.plotOptions{:});
            % Set figure properties
            [test.figureProperties, plotOptions]                ...
                = merge_options(test.defaultFigureProperties(), ...
                                plotOptions{:}                   );
            % Set axis properties
            [test.axisProperties, plotOptions]                ...
                = merge_options(test.defaultAxisProperties(), ...
                                plotOptions{:}                 );
            % Set toolbar options
            test.toolbarOptions ...
                = merge_options(test.defaultToolbarOptions(), ...
                                plotOptions{:}                 );
            names  = fieldnames(test.toolbarOptions);
            values = struct2cell(test.toolbarOptions);
            test.toolbarOptions = cell(1, 2*numel(names));
            test.toolbarOptions([1:2:end, 2:2:end]) = horzcat(names, values);
            if opt.plot
                % Plot example
                test.plot(test.model.rock, 'Name', 'rock'  , 'log10', true);
                test.plot(test.state0    , 'Name', 'state0'               );
            end
        end
        
        %-----------------------------------------------------------------%
        function props = defaultProperties(test, name) %#ok
        % Get default properties for given object
            name  = ['factory', name];
            props = get(groot, name);
            names = cellfun(@(pname) pname(numel(name)+1:end)        , ...
                             fieldnames(props), 'UniformOutput', false);
            props = cell2struct(struct2cell(props), names);
        end
                
        %-----------------------------------------------------------------%
        function props = defaultFigureProperties(test)
        % Set default figure properties
            % Get default properties
            props = test.defaultProperties('Figure');
            % Remove read-only properties
            props = rmfield(props, 'XDisplay');
            % Set extra property Size
            props.Size = props.Position(3:4);
        end
        
        %-----------------------------------------------------------------%
        function h = figure(test, varargin)
        % Get example figure
            h     = figure(varargin{:});
            names = fieldnames(test.figureProperties);
            for i = 1:numel(names)
                if strcmpi(names{i}, 'Size')
                    screensize = get(groot, 'ScreenSize');
                    wsize      = test.figureProperties.(names{i});
                    wsize      = min([wsize; screensize(3:4) - [0 80]]);
                    pos        = get(h,'Position');
                    pos        = min([pos(1:2); screensize(3:4)-wsize-[0 80]]);
                    set(h, 'Position', [pos, wsize]);
                else
                    set(h, names{i}, test.figureProperties.(names{i}));
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function props = defaultAxisProperties(test)
        % Get default axis properties
            props = test.defaultProperties('Axes');
            % Get grid
            G = test.model.G;
            if ~isempty(test.visualizationGrid)
                G = test.visualizationGrid;
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
                props.([xyz(i), 'LimitMethod']) = 'tight';
                props.([xyz(i), 'LimMode'])     = 'auto';
            end
            if isfield(test.schedule.control(1), 'W')
                W = test.schedule.control(1).W;
                if ~isempty(W) && G.griddim == 3
                    dz = xmax(3) - xmin(3);
                    props.ZLim(1) = props.ZLim(1) - dz*0.2;
                end
            end
            % Set aspect ratio
            props.PlotBoxAspectRatio = [1,1,0.5];
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
            % Remove read-only and other properties
            props = rmfield(props, 'Legend');
            props = rmfield(props, 'Colormap');
            pv = struct2cell(props);
            pn = fieldnames(props);
            rmtypes = {'matlab.graphics.GraphicsPlaceholder'};
            keep = true;
            for t = rmtypes
                keep = keep & cellfun(@(prop) ~isa(prop, t{1}), pv);
            end
            pv = pv(keep); pn = pn(keep);
            props = cell2struct(pv, pn);
        end
        
        %-----------------------------------------------------------------%
        function ax = setAxisProperties(test, ax)
            % Set properties to example figure axis
            names = fieldnames(test.axisProperties);
            for i = 1:numel(names)
                set(ax, names{i}, test.axisProperties.(names{i}));
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
        function hash = getTestCaseHash(test, fullHash)
            % Get hash of options struct, prepended with test case name
            hashOpt = struct2hash(test.options, test.name);
            if nargin < 2 || ~fullHash
                hash = hashOpt; return;
            else
                % Get hash of selected test case properties
                hashG        = struct2hash(test.model.G);
                hashRock     = struct2hash(test.model.rock);
                hashFluid    = struct2hash(test.model.fluid);
                hashState0   = struct2hash(test.state0);
                hashSchedule = struct2hash(test.schedule);
                % Concatenate and compute hash of combined hashes
                str = strjoin({hashOpt, hashG, hashRock, hashFluid, ...
                                           hashState0, hashSchedule}, '_');
                hash = str2hash(str);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function ok = save(example)
            % Get example hash value
            hash = example.getTestCaseHash();
            % Store in default data directory under 'example-suite'
            pth = fullfile(mrstDataDirectory(), 'example-suite');
            if ~exist(pth, 'dir')
                mkdir(pth)
            end
            save(fullfile(pth, hash), 'example', '-v7.3');
            ok = true;
        end
        
        %-----------------------------------------------------------------%
        function tc = load(test, throwError)
            % Attempt to load test case
            if nargin < 2
                throwError = true;
            end
            % Test case is stored with a filename equal to its hash value
            hash = test.getTestCaseHash();
            pth  = fullfile(mrstDataDirectory(), 'test-suite', test.name);
            fn   = [fullfile(pth, hash), '.mat'];
            tc   = [];
            if isfile(fn)
                % The file exists - load it
                if test.verbose
                    fprintf('Found a saved version of this example. Loading, ...');
                end
                data = load(fn);
                tc   = data.example;
            elseif throwError
                % The file does no exists - throw an error
                error(['Could not find example '   , ...
                       'named %s with hash key %s'], test.name, hash);
            end
        end
        
        %-----------------------------------------------------------------%
        function problem = getPackedSimulationProblem(test, varargin)
        % Make packed problem with reasonable simulation setup
            opt = struct('Name'           , [], ...
                         'LinearSolver'   , [], ...
                         'NonLinearSolver', []);
            [opt, extra] = merge_options(opt, varargin{:});
            if isempty(opt.Name)
                opt.Name = test.getTestCaseHash(true);
            end
            has_ls  = ~isempty(opt.LinearSolver);    % Linear solver given
            has_nls = ~isempty(opt.NonLinearSolver); % Nonlinear solver given
            if has_nls
                % Check if nonlinear solver has non-default linear solver
                has_ls = has_ls || ...
                    ~isa(opt.NonLinearSolver.LinearSolver, 'BackslashSolverAD');
            end
            if ~has_ls
                % Select apropriate linear solver
                rmodel = test.model;
                if isa(rmodel, 'WrapperModel')
                    rmodel = rmodel.getReservoirModel;
                end
                opt.LinearSolver = selectLinearSolverAD(rmodel);
            end
            if ~has_nls
                % Select default nonlinear solver
                opt.NonLinearSolver = NonLinearSolver( ...
                                    'LinearSolver' , opt.LinearSolver, ...
                                    'useRelaxation', true            );
            elseif ~has_ls
                % ... or assign linear solver
                opt.NonLinearSolver.LinearSolver = opt.LinearSolver;
            end
            % Pack problem
            desc = test.description(1:min(numel(test.description), 113));
            problem = packSimulationProblem(                          ...
                test.state0, test.model, test.schedule, test.name,    ...
                    'Name'           , opt.Name                     , ...
                    'Description'    , desc                         , ...
                    'NonLinearSolver', opt.NonLinearSolver, extra{:});
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
