classdef TestCase
    
    properties
        % Core test case properties
        name        % Example name
        description % One-line test case description
        state0      % Initial state
        model       % Model
        schedule    % Simulation schedule
        options     % Options passed to setup function
        verbose     % Verbose flag
        extra       % Anything that does not fit into the other categories 
                    % above (e.g., reference data)
        % Properties for plotting
        figureProperties  % Figure propreties (defaults set by constructor)
        axisProperties    % Axis properties (defaults set by constructor)
        toolbarOptions    % Options that can be passed to plotToolbar
        visualizationGrid % Optional grid for plotting
    end

    properties (Access = private)
        baseName  % Name of the setup function called by the constructor
    end
    
    methods
        %-----------------------------------------------------------------%
        function test = TestCase(name, varargin)
            
            if nargin == 0, return; end
            % Set test name
            test.baseName = lower(name);
            % Merge options
            opt = struct('plot', false, 'readFromDisk', true);
            [opt , varargin] = merge_options(opt, varargin{:});
            [test, extra   ] = merge_options(test, varargin{:});
            if isempty(test.verbose), test.verbose = mrstVerbose(); end
            if test.verbose
                fprintf('Setting up %s test \n\n', test.baseName);
                timer = tic();
            end
            if opt.readFromDisk
                % Get test options and description
                setup = feval(lower(name), false, extra{:});
                test.options     = setup.options;
                test.description = setup.description;
                % Check if test exists on disk and load it
                tc = test.load(false);
                if ~isempty(tc)
                    % We found it! Early return
                    test = tc;
                    if test.verbose
                        time = toc(timer);
                        fprintf('Test case loaded in %s\n\n', formatTimeRange(time));
                    end
                    return
                elseif test.verbose
                    % No luck, we need to set up from scratch
                    fprintf(['Did not find a saved version of this ', ...
                             'test case, setting up from scratch \n\n']);
                end
            end
            % Set up test case
            setup = feval(lower(name), extra{:});
            % Set properties to TestCase object
            if isempty(test.name), test.name = setup.name; end
            test.description = setup.description;
            test.state0      = setup.state0;
            test.model       = setup.model;
            test.schedule    = setup.schedule;
            test.options     = setup.options;
            test.extra       = setup.extra;
            % Print progress
            if test.verbose
                time = toc(timer);
                fprintf('Test case set up in %s\n\n', formatTimeRange(time));
            end
            % Set visualization grid if given
            [test, plotOptions] = merge_options(test, setup.plotOptions{:});
            % Set figure properties
            [test.figureProperties, plotOptions] ...
                = test.processFigureProps(plotOptions{:});
            % Set axis properties
            [test.axisProperties, plotOptions] ...
                = test.processAxisProps(plotOptions{:});
            % Set toolbar options
            test.toolbarOptions ...
                = merge_options(test.defaultToolbarOptions(), ...
                                plotOptions{:}              );
            names  = fieldnames(test.toolbarOptions);
            values = struct2cell(test.toolbarOptions);
            test.toolbarOptions = cell(1, 2*numel(names));
            test.toolbarOptions([1:2:end, 2:2:end]) = horzcat(names, values);
            if opt.plot
                % Plot test case
                test.plot(test.model.rock, 'Name', 'rock'  , 'log10', true);
                test.plot(test.state0    , 'Name', 'state0'               );
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [props, other] = validProperties(test, type, args, extra, readOnly) %#ok
        % Get default properties for given object

            type  = ['factory', type];
            props = get(groot, type);
            validNames = cellfun(@(pname) pname(numel(type)+1:end)        , ...
                             fieldnames(props), 'UniformOutput', false);
            if nargin > 2, validNames = [validNames; extra]; end
            if nargin > 3, validNames = setdiff(validNames, readOnly); end

            names = args(1:2:end);
            keep  = reshape(ismember(names, validNames), 1, []);
            keep  = reshape(repmat(keep,2,1), [], 1);
            props = args(keep);
            other = args(~keep);
            if ~any(keep), props = struct(); return; end
            props = cell2struct(props(2:2:end), props(1:2:end), 2);
            
        end
                
        %-----------------------------------------------------------------%
        function [props, other] = processFigureProps(test, varargin)
        % Set default figure properties
        
            % Get default properties
            [props, other] = test.validProperties('Figure', varargin, 'Size', 'XDisplay');

        end
        
        %-----------------------------------------------------------------%
        function h = figure(test, varargin)
        % Get test case figure
        
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
        function [props, other] = processAxisProps(test, varargin)
        % Get default axis properties
        
            [props, other] = test.validProperties('Axes', varargin, [], ...
                {'Legend', 'Colormap', 'ZDir', 'PlotBoxAspectRatioMode'});
            % Get grid
            G = test.model.G;
            if ~isempty(test.visualizationGrid)
                G = test.visualizationGrid;
            end
            xyz = 'XYZ';
            for i = 1:G.griddim
                % Set axis tight
                nm = [xyz(i), 'LimSpec'];
                if ~isfield(props, nm), props.(nm) = 'tight'; end
                nm = [xyz(i), 'LimitMethod'];
                if ~isfield(props, nm), props.(nm) = 'tight'; end
                nm = [xyz(i), 'LimMode'];
                if ~isfield(props, nm), props.(nm) = 'auto' ; end
            end
            % Set viewpoint
            if ~isfield(props, 'View')
                if G.griddim == 2
                    props.View = [0, 90];
                elseif G.griddim == 3
                    props.View = [-37.5, 30];
                end
            end
            % Set axis projection
            if ~isfield(props, 'Projection')
                if G.griddim == 2
                    props.Projection = 'orthographic';
                else
                    props.Projection = 'perspective';
                end
            end
            
        end
        
        %-----------------------------------------------------------------%
        function ax = setAxisProperties(test, ax)
            % Set properties to test case figure axis
            
            names = fieldnames(test.axisProperties);
            for i = 1:numel(names)
                set(ax, names{i}, test.axisProperties.(names{i}));
            end
            
        end
        
        %-----------------------------------------------------------------%
        function opt = defaultToolbarOptions(test) %#ok
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
                         'surf'         , false, ...
                         'plotmarkers'  , false, ...
                         'skipAugmented', false, ...
                         'field'        , 's:1', ...
                         'step_index'   , 1    , ...
                         'startplayback', false);
        
        end
        
        %-----------------------------------------------------------------%
        function h = plot(test, varargin)
        % Plot field v on test case grid using plotToolbar
            
            if nargin == 1 || ischar(varargin{1})
                v = test.model.rock;
            else
                v = varargin{1};
                varargin = varargin(2:end);
            end
            opt = struct('Name'     , ''  , ...
                         'plotWells', true, ...
                         'wellOpts' , {{}}, ...
                         'camlight' , true);
            [opt, varargin] = merge_options(opt, varargin{:});
            Name = test.name;
            if ~isempty(opt.Name)
                Name = [Name, ' ', opt.Name];
            end
            h = test.figure('Name', Name);
            G = test.getVisualizationGrid();
            plotToolbar(G, v, test.toolbarOptions{:}, varargin{:});
            if opt.plotWells
                test.plotWells(opt.wellOpts{:});
            end
            test.setAxisProperties(gca);
            if opt.camlight && G.griddim == 3 ...
                    && ~all(test.axisProperties.View == [0,90])
                camlight;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function G = getVisualizationGrid(test)
            
            if isempty(test.visualizationGrid)
                G = test.model.G;
            else
                G = test.visualizationGrid;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function varargout = plotWells(test, varargin)
            
            if ~isfield(test.schedule.control(1), 'W')
                return;
            end
            W = test.schedule.control(1).W;
            if ~isempty(W)
                G = test.getVisualizationGrid();
                if G.griddim == 3
                    ax = gca;
                    dz = ax.ZLim(2) ...
                       - ax.ZLim(1);
                    varargout = cell(nargout, 1);
                    [varargout{:}] = plotWell(G, W, 'color' , 'k'    , ...
                                       'height', 0.15*dz, ...
                                       varargin{:}      );
                else
                    opt = struct('fontsize', 12 , ...
                                 'color'   , 'k', ...
                                 'color2'  , 'k');
                    [opt, varargin] = merge_options(opt, varargin{:});
                    nw = numel(W);
                    [hw, htext] = deal(zeros(nw,1));
                    gray = [1,1,1]*0.8;
                    hold on
                    for w = 1:numel(W)
                        color = opt.color;
                        if W(w).sign < 0, color = opt.color2; end
                        x  = G.cells.centroids(vertcat(W(w).cells), :);
                        hw(w) = plot(x(1), x(2)      , 'o'   , ...
                                     'Color'          , color, ...
                                     'markerFaceColor', gray , ...
                                     'markerSize'     , 8    , ...
                                     varargin{:}             );
                        if opt.fontsize > 0
                            htext(w) = text(x(1), x(2), W(w).name, ...
                                'FontSize'         , opt.fontsize, ...
                                'VerticalAlignment', 'bottom'    , ...
                                'interp'           , 'latex'     , ...
                                'Color'            , color       );
                        end
                    end
                    hold off
                end
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [hash, s] = getTestCaseHash(test, fullHash)
        % Get hash of options struct, prepended with test case name
        
            skip = {'figureProperties' , ...
                    'axisProperties'   , ...
                    'toolbarOptions'   , ...
                    'visualizationGrid', ...
                    'verbose'          };
            if nargin > 1 && ~fullHash
                skip = [skip, {'model', 'state0', 'schedule', 'extra'}];
            end
            [hash, s] = obj2hash(test, 'skip', skip);
            
        end
        
        %-----------------------------------------------------------------%
        function ok = save(test)
        % Save test case to disk
        
            % Get test case hash value
            hash = test.getTestCaseHash(false);
            % Store in default data directory under 'test-suite'
            pth = fullfile(mrstDataDirectory(), 'test-suite', ...
                           test.baseName, test.name);
            if ~exist(pth, 'dir')
                mkdir(pth)
            end
            save(fullfile(pth, hash), 'test', '-v7.3');
            ok = true;
            
        end
        
        %-----------------------------------------------------------------%
        function tc = load(test, throwError)
        % Load test case from disk
        
            % Attempt to load test case
            if nargin < 2
                throwError = true;
            end
            % Test case is stored with a filename equal to its hash value
            hash = test.getTestCaseHash(false);
            pth  = fullfile(mrstDataDirectory(), 'test-suite', ...
                            test.baseName, test.name);
            fn   = [fullfile(pth, hash), '.mat'];
            tc   = [];
            if isfile(fn)
                % The file exists - load it
                if test.verbose
                    fprintf('Found a saved version of this test case. Loading, ...');
                end
                data = load(fn);
                tc   = data.test;
            elseif throwError
                % The file does no exists - throw an error
                error(['Could not find test case '   , ...
                       'named %s with hash key %s'], test.name, hash);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function problem = getPackedSimulationProblem(test, varargin)
        % Make packed problem with reasonable simulation setup
        
            opt = struct('Name'           , test.name, ...
                         'useHash'        , true     , ...
                         'LinearSolver'   , []       , ...
                         'NonLinearSolver', []       );
            [opt, varargin] = merge_options(opt, varargin{:});
            if opt.useHash
                if ~isempty(opt.Name)
                    opt.Name = [opt.Name, '_'];
                end
                opt.Name = [opt.Name, test.getTestCaseHash()];
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
                try
                    opt.LinearSolver = selectLinearSolverAD(rmodel);
                catch ex
                    warning(ex.identifier, ...
                            ['Unable to select linear solver from ' , ...
                             'model. Reason: %s \n'                 , ...
                             'Proceeding with BackslashSolverAD'   ], ...
                             ex.message                             );
                    opt.linearSolver = BackslashSolverAD();
                end
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
            if numel(desc) < numel(test.description)
                desc = [desc, ' ...'];
            end
            problem = packSimulationProblem(test.state0, test.model, ...
                    test.schedule, test.baseName,  ...
                    'Name'           , opt.Name,   ...
                    'Description'    , desc,       ...
                    'NonLinearSolver', opt.NonLinearSolver, varargin{:});
        
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
