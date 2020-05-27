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
        figureProperties % Figure propreties (defaults set by constructor)
        axisProperties   % Axis properties (defaults set by constructor)
        toolbarOptions   % Options that can be passed to plotToolbar
    end
    
    methods
        %-----------------------------------------------------------------%
        function example = MRSTExample(name, varargin)
            % Set example name
            example.name = lower(name);
            % Merge options
            opt = struct('plot', false, 'verbose', []);
            [opt, varargin]  = merge_options(opt, varargin{:});
            example.verbose = opt.verbose;
            if isempty(opt.verbose), example.verbose = mrstVerbose(); end
            if example.verbose
                fprintf('Setting up %s example \n\n', example.name);
                timer = tic();
            end
            [example, extra] = merge_options(example, varargin{:});
            [example.description, ... % Example description
             example.state0     , ... % Initial state
             example.model      , ... % Model
             example.schedule   , ... % Simulation schedule
             example.options    , ... % Options for reference
             plotOptions        ] ... % Plotting options
                = feval(lower(name), extra{:});
            if example.verbose
                time = toc(timer);
                fprintf('Example set up in %s\n\n', formatTimeRange(time));
            end
            % Set figure properties
            [example.figureProperties, plotOptions]                ...
                = merge_options(example.defaultFigureProperties(), ...
                                plotOptions{:}                   );
            % Set axis properties
            [example.axisProperties, plotOptions]                ...
                = merge_options(example.defaultAxisProperties(), ...
                                plotOptions{:}                 );
            % Set toolbar options
            example.toolbarOptions = merge_options(example.defaultToolbarOptions, plotOptions{:});
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
                    position = [h.Position(1:2), example.figureProperties.(names{i})];
                    set(h, 'Position', position);
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
            if isfield(example.options, 'Gviz')
                G = example.options.Gviz;
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
            W = example.schedule.control(1).W;
            if ~isempty(W) && G.griddim == 3
                dz = xmax(3) - xmin(3);
                props.ZLim(1) = props.ZLim(1) - dz*0.2;
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
        function plot(example, v, varargin)
            % Plot filed v on example grid using plotToolbar
            opt = struct('Name', '');
            [opt, extra] = merge_options(opt, varargin{:});
            Name = example.name;
            if ~isempty(opt.Name)
                Name = [Name, ' ', opt.Name];
            end
            if nargin == 1
                v = example.model.rock;
            end
            G = example.model.G;
            if isfield(example.options, 'Gviz')
                G = example.options.Gviz;
            end
            example.figure('Name', Name);
            plotToolbar(G, v, example.toolbarOptions{:}, extra{:});
            example.setAxisProperties(gca);
            W = example.schedule.control(1).W;
            if ~isempty(W)
                if G.griddim == 3
                    plotWell(G, W, 'color', 'k');
                else
                    hold on
                    x = G.cells.centroids(vertcat(W.cells), :);
                    plot3(x(:,1), x(:,2), ones(size(x,1),1), 'ok', ...
                         'markerFaceColor', 'w', ...
                         'markerSize'     , 8  );
                    hold off
                end
            end
            if G.griddim == 3 && ~all(example.axisProperties.View == [0,90])
                camlight;
            end
        end
        
        %-----------------------------------------------------------------%
        function problem = getPackedSimulationProblem(example, varargin)
            % Make packed problem with reasonable simulation setup
            opt = struct('LinearSolver'   , [], ...
                         'NonlinearSolver', []);
            [opt, varargin] = merge_options(opt, varargin{:});
            ix = strcmpi('ExtraArguments', varargin);
            extra = {};
            if ~isempty(ix)
                extra = varargin{ix+1};
                opt   = merge_options(opt, extra{:});
                varargin(ix:ix+1) = [];
            end
            if isempty(opt.LinearSolver)
                rmodel = example.model;
                if isa(rmodel, 'WrapperModel')
                    rmodel = rmodel.getReservoirModel;
                end
                lsolver = selectLinearSolverAD(rmodel);
                opt.LinearSolver = lsolver;
            end
            extra = horzcat(extra, {'LinearSolver'   , opt.LinearSolver   , ...
                                    'NonlinearSolver', opt.NonlinearSolver});
            problem = packSimulationProblem(                                   ...
                example.state0, example.model, example.schedule, example.name, ...
                varargin{:}, 'ExtraArguments', extra);
        end
        
    end
    
end