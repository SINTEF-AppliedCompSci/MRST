function interactiveDiagnostics(G, rock, W, varargin)
%Launch an interactive diagnostics session
%
% SYNOPSIS:
%   interactiveDiagnostics(G, rock, W);
%   interactiveDiagnostics(G, rock, W, 'state', state)
%   interactiveDiagnostics(G, rock, W, dataset, 'state', state)
%
% DESCRIPTION:
%   This function launches an interactive session for doing flow
%   diagnostics. The functionality differs slightly based on the input
%   arguments given:
%     - If a dataset is given, this set of cellwise data will be available
%       for visualization.
%     - If a state is given, this state will allow for visualization of
%       the component ratios inside a drainage volume. The flux from this
%       state can also be used to calculate time of flight if "computeFlux"
%       is disabled, for instance if the user has some external means of
%       computing fluxes.
%
%   Once the initialization is complete, two windows will be produced:
%     - A plotting window, showing the reservoir along with the wells and
%       visualized quantities. In the plotting window, it is possible to
%       click wells to get additional information, such as the allocation
%       factors per perforation in the well, the phase distribution grouped
%       by time of flight and the corresponding pore volumes swept/drained.
%
%     - A controller window which is used to alter the state of the
%       plotting window:
%
%       1) Wells can be selected for display (if a well is selected, the
%       corresponding drainage (producer) or flooding (injector) volumes
%       will be visualized. A simple playback function can be used to show
%       propagation of time of flight. Different quantities can be
%       visualized to get a better understanding of the system.
%
%       2) A set of buttons provide (experimental) well editing, access to
%       visualization of the well pair connections, Phi/F diagram with
%       Lorenz' coefficient etc.
%
%       3) A player controller (experimental), if a time-series of wells
%       (and states) are supplied.
%
%
% REQUIRED PARAMETERS:
%
%   G    - Valid grid structure.
%
%   rock - Rock with valid permeability and porosity fields.
%
%   W    - A set of wells which are compatible with the incompTPFA solver.
%
%
% OPTIONAL PARAMETERS:
%
%  'state' - Reservoir state containing fluid saturations and optionally
%       flux and pressure (if 'computeFlux' is false)
%
%  'computeGrid' - Use G for plotting, but computeGrid for computing
%       time-of-flight, tracer partitions, etc. 
%
%  'tracerfluid' - Fluid used for tracer computation. This defaults to a
%       trivial single-phase fluid.
%
%  'fluid' - Fluid used for mobility calculations if 'useMobilityArrival'
%       is enabled.
%
%  'LinSolve' - Linear solver for pressure systems. Defaults: mldivide
%
%  'computeFlux' - If set to false, fluxes are extracted from the provided
%       state keyword argument. This requires a state to be provided and
%       can be used if the fluxes are computed externally (for instance
%       from a expensive full-physics simulation)
%
%  'useMobilityArrival' - If the well plot showing nearby saturations
%       should plot mobility instead of saturations. This may be
%       interesting in some cases because the mobile fluids are more likely
%       to be extracted. However, this plot is often dominated by very
%       mobile gas regions.
%
%  'daspect' - Data aspect ratio, in a format understood by daspect()
%
%  'name' - Name to use for windows
%
%  'leaveOpenOnClose' - Default false. Leaves all figures open when closing
%       the controller.
%
%  'fastRotate' - Enables fast rotate if true, disables if false. Default
%       is enabled for models with > 50000 cells, disabled otherwise.
%
% RETURNS:
%
%   Nothing. Creates two figures.
%
% EXAMPLE:
%     G = computeGeometry(cartGrid([10, 10, 2]));
%     rock = struct('poro', ones(G.cells.num, 1), 'perm', ones(G.cells.num, 1)*darcy)
%
%     W = verticalWell([], G, rock, 1,  1, [], 'Val', 1)
%     W = verticalWell(W, G, rock,  1,  10, [], 'Val', 1)
%     W = verticalWell(W, G, rock,  10, 1, [], 'Val', 1)
%     W = verticalWell(W, G, rock,  10, 10, [], 'Val', 1)
%     W = verticalWell(W, G, rock,  5,  5, [], 'Val', 0)
%     interactiveDiagnostics(G, rock, W);
%
% SEE ALSO:
%   `Diagnostics` `examples`

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
    try
       require mrst-gui incomp
    catch %#ok<CTCH>
       mrstModule add mrst-gui incomp
    end

    if nargin > 3 && mod(numel(varargin), 2) == 1
        dsname = inputname(4);
        datasets = varargin{1};
        varargin = varargin(2 : end);
    else
        dsname = '';
        datasets = [];
    end
    
    water    = initSingleFluid('mu' , 1*centi*poise, 'rho', 1000*kilogram/meter^3);
    oilwater = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                               'rho', [1000, 859]*kilogram/meter^3, ...
                               'n'  , [   2,   2]);
    oilwater.names = {'Water', 'Oil', 'Gas'};
    
    opt = struct('state',               [],...
                 'computeGrid',         [],...
                 'tracerfluid',         water, ...
                 'LinSolve',            @mldivide, ...
                 'computeFlux',         true, ...
                 'fluid',               oilwater, ...
                 'wellPlotFn',          [], ...
                 'name',                [], ...
                 'daspect',             [], ...
                 'leaveOpenOnClose',    false, ...
                 'lineWells',           false, ...
                 'maxTOF',              [], ...
                 'useLight',            false, ...
                 'fastRotate',          [], ...
                 'computeWellTOFs',     false, ...
                 'showGrid',            false, ...
                 'ignoreWellRate',      .01 ...
    );

    opt = merge_options(opt, varargin{:});
    
    if (isempty(opt.name))
        opt.name = {'Interactive Diagnostics'};
    elseif (isa(opt.name, 'char'))
        opt.name = {opt.name};
    end

    assert(opt.computeFlux || ~isempty(opt.state),...
        'If computeFlux is off a state must be provided!')
    
    if (isstruct(W))
        W = {W};
    end
    
    %Create anonymous function for remapping
    %values from the compute grid to fit the plotting grid G
    cdataToPlotGrid = @(f) f;
    wellsToPlotGrid = @(f) f;
    computeGrid = G;
    if (~isempty(opt.computeGrid))
        assert(isfield(G.cells, 'eMap'), ...
            'G must have an eMap field when using a separate computeGrid');
        cdataToPlotGrid = @(cdata) cdata(G.cells.eMap,:);
        wellsToPlotGrid = @(wells) remapWells(G, wells);
        computeGrid = opt.computeGrid;
    end
        
    state = {[]};
    state_idx = 1;
    if (~isempty(opt.state))        
        if (numel(opt.state) == 1)
            state = {opt.state};
        else
            state = opt.state;
        end
        
        assert(numel(W) == numel(state), ...
            'W and state must have equal number of elements');
    end
    

    if (~isempty(datasets))
        if (numel(datasets) == 1)
            datasets = {datasets};
        end
        
        assert(numel(W) == numel(datasets), ...
            'W and datasets must have equal number of elements');
    elseif(isempty(datasets) && ~isempty(state))
        dsname = 'State';
        datasets = state;
    else
        dsname = '';
        datasets = cell(1,1);
    end
    


    % Main scope variables for tracers and tof
    [D, D_int, WP] = deal([]);

    % Main figures
    [fig_main,... 
     fig_ctrl,...
     fig_well,...
     fig_phi,...
     fig_well_alloc,...
    ] = deal(NaN);

    %Gui elements
    [ctrl_drain_vols,...
     ctrl_flood_vols,...
     ctrl_well_conn,...
     ctrl_use_average,...
     ctrl_show_grid,...
     ctrl_plot_all_wells,...
     speedsh,...
     tofsh,...
     nwtofsh,...
     mtofsh,...
     mtofeh,...
     Mtofsh,...
     Mtofeh,...
     alfash,...
     tofext,...
     extents,...
     hdataset,...
     hset_op,...
     mrst_ds,...
     ds_panel,...
     lights] = deal(NaN);
    
    %Handles for well plot in main figure
    fig_main_wells = {};
    fig_main_wells.hwells = [];
    fig_main_wells.dirty = true;
    
    %Handles for cell data in main figure
    fig_main_cell_data = []; % "Grid with cell values"
    fig_main_grid = []; %Grid outline
    
    %Well connection pairs in main figure
    fig_main_well_conns = []; 
    
    %Handles that should be enabled / disabled
    %according to if we're using averaged or single 
    %timestep values
    average_gui_handles = [];
    non_average_gui_handles = [];
    
    %"GUI State" variables
    window_name = opt.name{state_idx};
    compute_average = false;
    playback = false; %Play back TOF?
    time_varying = (~isempty(state) && numel(state) > 1) || (~isempty(W) && numel(W) > 1);
    selection = []; % Selection of cells
    
    % Storage for diagnostics and wellpair computations
    [WellPairs, Diagnostics] = deal(cell(numel(state), 1));
    
    
    %Create main figure
    fig_main = figure('Name', window_name);    

    %Set scaling
    max_G = max(G.faces.centroids);
    min_G = min(G.faces.centroids);
    extents = zeros(1,6);
    extents(1:2:end) = min_G - 0.05*(max_G-min_G);
    extents(2:2:end) = max_G + 0.05*(max_G-min_G);
    
    if (~isempty(opt.daspect))
        daspect(opt.daspect);
    else
        da = max_G - min_G;
        da(da == 0) = 1;
        daspect(da);
    end
    % Create control panel
    % Precompute TOF etc. for creating main control
	pv = poreVolume(computeGrid, rock);
    [D, WP] = getDiagnostics(state_idx);
    tofext = getTOFRange(D, opt);
    createMainControl();
        
    % Trigger plot initial setup
    plotMain();
    axis off vis3d;
    axis(extents);
    %view(3);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Only function definitions from here on %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Clamp real: Clamps data to real valued numbers (non-inf)
    function val = clamp_real(in_val)
        val = max(-realmax, min(realmax, in_val));
    end

    %
    function tofext = getTOFRange(D, opt)
        if (D.isvalid)
            tof = D.tof(:);
            if ~isempty(opt.maxTOF)
               tofext = convertTo([min(tof), opt.maxTOF], year);
            else
               isNeg = (tof <= 0);
               tof(isNeg) = min(tof(~isNeg));

               tofext = convertTo(([min(tof), 5*10^(mean(log10(tof)))]), year);
               tofext(2) = clamp_real(max(tofext(2), 15));
            end
        else
            tofext = [0 1];
        end
    end
        
    function setTOFSliderExtents()
        %Update min/max and current value
        min_val = tofext(1);
        max_val = tofext(2);
        cur_m_val = str2double(get(mtofeh, 'String'));
        cur_m_val = max(min(cur_m_val, max_val), min_val);

        cur_M_val = str2double(get(Mtofeh, 'String'));
        cur_M_val = max(min(cur_M_val, max_val), min_val);

        set(mtofsh, 'Min', tofext(1));
        set(mtofsh, 'Max', tofext(2));
        set(mtofsh, 'Value', cur_m_val);
        %set(mtofeh, 'String', num2str(cur_m_val));
        mtofs_callback([], []);

        set(Mtofsh, 'Min', tofext(1));
        set(Mtofsh, 'Max', tofext(2));
        set(Mtofsh, 'Value', cur_M_val);
        %set(Mtofeh, 'String', num2str(cur_M_val));
        Mtofs_callback([], []);
    end

    % Set special functions for the min/max time of flight slider
    % handle
    function mtofs_callback(src, event)
        min_val = get(mtofsh, 'Min');
        max_val = get(mtofsh, 'Max');
        cur_val = get(mtofsh, 'Value');
        if (cur_val <= min_val)
            set(mtofsh, 'Value', min_val);
            set(mtofeh, 'String', sprintf('%.1f', 0));
        elseif (cur_val >= max_val)
            set(mtofsh, 'Value', max_val);
            set(mtofeh, 'String', sprintf('%.1f', Inf));
        else
            set(mtofeh, 'String', sprintf('%.1f', cur_val));
        end

        %Only plot if src is nonempty (true callback)
        if (~isempty(src))
            if (compute_average)
                D_int = [];
                computeAverage();
            end
            plotMain(src, event);
        end
    end

    function Mtofs_callback(src, event)
        min_val = get(Mtofsh, 'Min');
        max_val = get(Mtofsh, 'Max');
        cur_val = get(Mtofsh, 'Value');
        if (cur_val <= min_val)
            set(Mtofsh, 'Value', min_val);
            set(Mtofeh, 'String', sprintf('%.1f', -Inf));
        elseif (cur_val >= max_val)
            set(Mtofsh, 'Value', max_val);
            set(Mtofeh, 'String', sprintf('%.1f', Inf));
        else
            set(Mtofeh, 'String', sprintf('%.1f', cur_val));
        end

        %Only plot if src is nonempty (true callback)
        if (~isempty(src))
            if (compute_average)
                D_int = [];
                computeAverage();
            end
            plotMain(src, event);
        end
    end


    function addVisualizationControls(parent, pos)
        vis_N = 1;
        vis_H = 1 / vis_N;
        vis_h = vis_H*0.9;
        bwidth = .2;
        
        vis_control = uipanel('Parent', parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'None');
        
        confp = uipanel(vis_control, 'Units', 'normalized', 'Position', [0 0*vis_H 1 vis_h], 'BorderType', 'None');
        uicontrol(confp, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 0 bwidth 1],...
                   'String', 'Phi/F diagram',...
                   'Callback', @(src, event) plotPhi(src, event)...
                   );
               
        uicontrol(confp, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [1*bwidth 0 bwidth 1],...
                   'String', 'Well allocations',...
                   'Callback', @plotWellAllocations ...
                   ); 

        if opt.computeFlux
            uicontrol(confp, 'Style', 'pushbutton',...
                       'Units', 'normalized',...
                       'Position', [2*bwidth 0 bwidth 1],...
                       'String', 'Edit wells',...
                       'Callback', @(src, event) changeWells()...
                       );
        end
    end


    function addTOFControls(parent, pos)
        %Function that forces update of time integral TOF
        function dirtyTOFPlotMain(src, event)
            if (compute_average)
                D_int = [];
                computeAverage();
            end
            plotMain(src, event);
        end
        
        % Panel for Time of flight selection parameters
        tof_panel = uipanel('Parent', parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'none');
        tof_N = 7;
        tof_H = 1/tof_N;
        tof_h = tof_H*0.8;
                               
        tof_N = tof_N - 1;
        uicontrol(tof_panel, 'Style', 'text', ...
                             'string', 'Selection:', ...
                             'HorizontalAlignment', 'left', ...
                             'Units', 'normalized',...
                             'Position', [.0 tof_N*tof_H .25 tof_h]); 
        hset_op = uicontrol(tof_panel, 'Style', 'popupmenu',...
                                   'Units', 'normalized',...
                                   'Position', [.4 tof_N*tof_H .6 tof_h],...
                                   'Callback', @plotMain,...
                                   'String', {'Union {Flood, Drain} volumes',...
                                              'Intersection {Flood, Drain} volumes',...
                                              'Flood volumes',...
                                              'Drain volumes'}...
                                   );
        
        tof_N = tof_N - 1;
        cellfields = getStructFields(G, datasets{state_idx}, dsname);
        perm = strcat({'X', 'Y', 'Z'}, ' permeability');
        uicontrol(tof_panel, 'Style', 'text', ...
                             'string', 'Display:', ...
                             'HorizontalAlignment', 'left', ...
                             'Units', 'normalized',...
                             'Position', [.0 tof_N*tof_H .25 tof_h]);
        % include well-tofs if computed
        if all(isfield(D, {'itof', 'ptof'}))
            wtofNms = {'TOF selected injectors', 'TOF selected producers'};
        else
            wtofNms = {};
        end 
        hdataset = uicontrol(tof_panel, 'Style', 'popupmenu',...
                                    'Units', 'normalized',...
                                    'Position', [.4 tof_N*tof_H .6 tof_h],...
                                    'Callback', @plotMain,...
                                    'String', {'Forward TOF',...
                                               'Backward TOF',...
                                               'Sum of TOFs',...
                                               'Drainage region',...
                                               'Flooding region',...
                                               'Drainage blend',...
                                               'Flooding blend',...
                                               'Tracer selected injectors',...
                                               'Tracer selected producers',...
                                               wtofNms{:}, ...
                                               'Porosity', ...
                                               perm{1:size(rock.perm, 2)},...
                                               cellfields{:}}...
                                   ); %#ok<CCAT>
        
        tof_N = tof_N - 1;
        [nwtofsh, ~, ~]   = linkedSlider(tof_panel, [0 tof_N*tof_H 1 tof_h], .15, [0 1], 0.0, 'Near well', @dirtyTOFPlotMain);
        
        tof_N = tof_N - 1;
        [Mtofsh, Mtofeh, ~]   = linkedSlider(tof_panel, [0 tof_N*tof_H 1 tof_h], .15, tofext, tofext(2), 'Max TOF', @dirtyTOFPlotMain);
        
        tof_N = tof_N - 1;
        [mtofsh, mtofeh, ~]   = linkedSlider(tof_panel, [0 tof_N*tof_H 1 tof_h], .15, tofext, tofext(1), 'Min TOF', @dirtyTOFPlotMain);
        set(mtofsh, 'Callback', @mtofs_callback);
        set(mtofeh, 'String', num2str(-Inf))
        set(Mtofsh, 'Callback', @Mtofs_callback);
        set(Mtofeh, 'String', num2str(+Inf))
        
        if (time_varying)
            tof_N = tof_N - 1;
            [tofsh, tofeh, ~]   = linkedSlider(tof_panel, [0 tof_N*tof_H 1 tof_h], .15, [0 1], 0.25, 'TOF Freq.', @plotMain);
            average_gui_handles = [average_gui_handles(:); tofsh; tofeh];
        end
        
        
        tof_N = tof_N - 1;
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 tof_N*tof_H .2 tof_h],...
                   'String', 'Play TOF',...
                   'Callback', @playBackTof...
                   );
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [.2 tof_N*tof_H .2 tof_h],...
                   'String', 'Stop TOF',...
                   'Callback', @stopPlayBackTof...
                   );
        [speedsh, ~, ~] = linkedSlider(tof_panel, [.45 tof_N*tof_H .55 tof_h], .15, [1 1000], 50, 'Steps', []);
        uicontrol(tof_panel,'Style', 'pushbutton', ...
                   'Units', 'normalized', ...
                   'Position', [.8 0 .2 tof_h], ...
                   'String', 'Exit', ...
                   'Callback', @stopApplication);
    end

        
    function addAdvancedControls(parent, pos)
        adv_panel = uipanel('Parent', parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'none');
        adv_N = 4;
        adv_H = 1/adv_N;
        adv_h = adv_H*0.8;
        bwidth = .2;
        
        [alfash, ~, ~]   = linkedSlider(adv_panel, [0 3*adv_H 1 adv_h], .15, [0 1], 1, 'Alpha', @plotMain);
        
        
               
        ctrl_show_grid = uicontrol(adv_panel, 'Style', 'checkbox',...
                   'Units', 'normalized',...
                   'Position', [0 2*adv_H bwidth adv_h],...
                   'String', 'Show grid',...
                   'Callback', @plotMain, ...
                   'Value', opt.showGrid ...
                   );
               
        ctrl_well_conn = uicontrol(adv_panel, 'Style', 'checkbox',...
                   'Units', 'normalized',...
                   'Position', [bwidth 2*adv_H bwidth adv_h],...
                   'String', 'Well pairs',...
                   'Callback', @plotWellConnections ...
                   );
               
        function allWellsChange(src, event)
            fig_main_wells.dirty = true;
            plotMain(src, event)
        end
        ctrl_plot_all_wells = uicontrol(adv_panel, 'Style', 'checkbox',...
                   'Units', 'normalized',...
                   'Position', [2*bwidth 2*adv_H bwidth adv_h],...
                   'String', 'Plot all wells',...
                   'Callback', @allWellsChange, ...
                   'Value', 1 ...
                   );
               
               
        % lighting (set both above and beneath..)
        function setLights(src, event)
            if (get(src, 'Value'))
                set(0, 'CurrentFigure', fig_main);
                ax = axis;
                [p1, p2] = deal(ax([2,3,5]), ax([1,4,6]));
                l1 = light('Position', p1 + 3*(p2-p1));
                [p1, p2] = deal(ax([2,3,6]), ax([1,4,5]));
                l2 = light('Position', p1 + 3*(p2-p1));
                lights = [l1, l2];
            elseif (numel(lights) > 0)
                delete(lights);
            end
        end
        uicontrol(adv_panel, 'Style', 'checkbox',...
                   'Units', 'normalized',...
                   'Position', [3*bwidth 2*adv_H bwidth adv_h],...
                   'String', 'Enable lighting',...
                   'Callback', @setLights, ...
                   'Value', opt.useLight ...
                   );
        
        % Function which saves variables to the workspace on demand
        function saveDataToWorkSpace(src, event) %#ok<*INUSD>
            indices = find(selection);
            
            labels = {'Sector model', 'Sector model (index map)'};
            varnames = labels;
            vars = {indices, computeGrid.cells.indexMap(indices)};
            
            %Clean the varnames of spaces etc
            clean = @(x) lower(regexprep(x,'[^a-zA-Z0-9]','_'));
            varnames = cellfun(@(x) clean(x), varnames, 'UniformOutput', false);

            export2wsdlg(labels, varnames, vars)
        end
        
        function saveDataToFile(src, event)
            title = 'Export selection to FLUXNUM';
            num_lines = 1;
            prompt = {'Enter filename:'};
            default = {'FLUXNUM.PROP'};
            
            while (true)
                answer = inputdlg(prompt, title, num_lines, default);
                if (numel(answer) == 0)
                    return;
                end

                filename = answer{1};
                if (exist(filename, 'file'))
                    choice = questdlg('File already exists. Do you want to overwrite?', ...
                        'File exists', ...
                        'Yes','No', 'No');
                    if (strcmp(choice, 'Yes'))
                        break;
                    end
                else
                    break;
                end
            end
            
            [fid, msg] = fopen(filename, 'wt');
            if fid < 0
                errordlg(['Failed to open output file ', msg, '. Aborting.']);
                return
            end
            
            fprintf(fid, 'FLUXNUM\n');
            var = zeros(prod(computeGrid.cartDims), 1);
            var( computeGrid.cells.indexMap(selection) ) = 1;
            out_format = [repmat('%d\t', 1, 20), '\n'];
            fprintf(fid, out_format, var);
            fprintf(fid, '/\n');
            fclose(fid);
        end
        
        save_panel = uipanel('Parent', adv_panel, 'Units', 'normalized', 'Position', [0 0*adv_H 1 2*adv_h], 'Title', 'Export selection');
        uicontrol(save_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 0.2 0.3 .8],...
                   'String', 'Save to workspace', ...
                   'Callback', @saveDataToWorkSpace...
                   );
        uicontrol(save_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0.4 0.2 0.3 0.8],...
                   'String', 'Save to file', ...
                   'Callback', @saveDataToFile...
                   );
    end


    function addTimestepControls(parent, pos)
        ds_panel = uipanel('Parent', parent,...
                'Units', 'normalized',...
                'Position', pos,...
                'BorderType', 'none');
            
        uicontrol(ds_panel,'Style','text',...
            'Units', 'normalized',...
            'HorizontalAlignment', 'left',...
            'String','Simulation timestep',...
            'Position', [0 0.5 0.4 .5]);

        if (numel(state) > 1)
            selector_datasets = state;
        elseif(numel(W) > 1)
            selector_datasets = W;
        else
            error('Programmer error');
        end
        
        function selectTimeStep(datas, index, fullname)                    %#ok<INUSL>
            %Get ID of currently selected wells
            drain_wells = get(ctrl_drain_vols, 'Value');
            flood_wells = get(ctrl_flood_vols, 'Value');
            drain_selection = D.prod(drain_wells);
            flood_selection = D.inj(flood_wells);
            well_selection = [drain_selection, flood_selection];
            fig_main_wells.dirty = true;
            
            state_idx = index;
            [D, WP] = getDiagnostics(state_idx);
            
            %Update main plot
            plotMain();
            
            %Update producers / injectors selection box text
            set(ctrl_drain_vols, 'String', {W{state_idx}(D.prod).name});
            set(ctrl_flood_vols, 'String', {W{state_idx}(D.inj).name});
            
            %Find out which items we should be selecting
            drain_wells = find(ismember(D.prod, well_selection));
            flood_wells = find(ismember(D.inj, well_selection));
            
            %Set selections
            set(ctrl_drain_vols, 'Value', drain_wells);
            set(ctrl_flood_vols, 'Value', flood_wells);
            
            %Update aux plots
            if (ishandle(fig_phi))
                plotPhi([], []);                
            end
            if (ishandle(fig_well_alloc))
                plotWellAllocations([], [])
            end
        end

        datasetSelector(G, selector_datasets, 'Parent', ds_panel, 'Location', ...
            [0.4, 0.5, 0.6, .5], 'Callback', @selectTimeStep, ...
            'active', 1, 'Nofields', true,...
            'Tag', 'datasetselector');
        mrst_ds = findobj('Tag', 'datasetselector', '-and', 'Parent', ds_panel);
        set(mrst_ds, 'BorderType', 'None');
        mrst_ds_children = allchild(mrst_ds);
        non_average_gui_handles = [non_average_gui_handles(:); mrst_ds_children];

        %Control for showing / computing average over time of e.g., TOF
        ctrl_use_average = uicontrol(ds_panel, 'Style', 'checkbox',...
                                   'Units', 'normalized',...
                                   'Position', [0 0 1 0.5],...
                                   'Callback', @handleIntegralChange,...
                                   'HorizontalAlignment', 'right',...
                                   'String', 'Average over timesteps'...
                                   );                               
                               
        function handleIntegralChange(src, event)
            %Options specific when using the average over timesteps of the TOF
            time_integral_options = {
                'Forward TOF frequency',...
                'Backward TOF frequency',...
                'Sum TOF frequency'};

            %Check if we are to compute the average over time
            %or for an instantaneous snapshot.
            if (ishandle(ctrl_use_average))
                compute_average = logical(get(ctrl_use_average, 'Value'));
                if (compute_average)
                    set(average_gui_handles, 'Enable', 'on');
                    set(non_average_gui_handles, 'Enable', 'off');
                    
                    %Add new integral options
                    old_options = get(hdataset, 'String');
                    new_options = {old_options{:}, time_integral_options{:}}; %#ok<CCAT>
                    set(hdataset, 'String', new_options);
                    
                    computeAverage();
                else
                    set(average_gui_handles, 'Enable', 'off');
                    set(non_average_gui_handles, 'Enable', 'on');

                    %Remove integral options
                    options = get(hdataset, 'String');
                    for idx=1:numel(time_integral_options)
                        mask = ~ismember(options, time_integral_options{idx});
                        options = options(mask);
                    end
                    set(hdataset, 'String', options);
                    curr_val = get(hdataset, 'Value');
                    if (curr_val > numel(options))
                        set(hdataset, 'Value', 1);
                    else
                        set(hdataset, 'Value', curr_val);
                    end
                    
                    [D, WP] = getDiagnostics(state_idx);
                end
                
            end
            
            %Only plot if src is nonempty (true callback)
            if (~isempty(src))
                plotMain(src, event);
            end
        end
    end



    function createMainControl()
        % Set up figure handles
        if ~ishandle(fig_ctrl)
            pos = get(fig_main, 'OuterPosition');
            size_xy = [475 550];
            pos_xy = pos(1:2) + pos(3:4) - size_xy - [400 0];
            fig_ctrl = figure('Position',[pos_xy,  size_xy], ...
                'Toolbar','none',...
                'MenuBar', 'none',...
                'CloseRequestFcn', @close_all_func);
            set(fig_ctrl, 'Name', ['Controller ', window_name]);
        else
            clf(fig_ctrl)
        end

        % Drainage / flooding controls
        np = numel(D.prod);
        ni = numel(D.inj);

        ctrl_drain_vols = zeros(np, 1);
        ctrl_flood_vols = zeros(ni, 1);
        
        gpy = 0.45;
        gph = .55;
        if (ishandle(mrst_ds))
            gpy = gpy + 0.15;
            gph = gph - 0.15;
        end

        drainp = uipanel(fig_ctrl, 'Title', 'Drainage volumes',...
                                   'Units', 'normalized',...
                                   'Position', [0 gpy .5 gph]);
        floodp = uipanel(fig_ctrl, 'Title', 'Flooding volumes',...
                                   'Units', 'normalized',...
                                   'Position', [.5 gpy .5 gph]);

        function dirtyWellsPlotMain(src, event)
            fig_main_wells.dirty = true;
            plotMain(src, event);
        end
        
        ctrl_drain_vols = uicontrol(drainp, 'Style', 'listbox', ...
                                          'Units', 'normalized',...
                                          'Max', 2, ...
                                          'Min', 0, ...
                                          'String', {W{state_idx}(D.prod).name},...
                                          'Value', [], ...
                                          'Callback', @dirtyWellsPlotMain, ...
                                          'Position', [0 0 1 1]);

        ctrl_flood_vols = uicontrol(floodp, 'Style', 'listbox', ...
                                          'Units', 'normalized',...
                                          'Max', 2, ...
                                          'Min', 0, ...
                                          'String', {W{state_idx}(D.inj).name},...
                                          'Value', [], ...
                                          'Callback', @dirtyWellsPlotMain, ...
                                          'Position', [0 0 1 1]);
        
        warning('off', 'MATLAB:uitabgroup:OldVersion');
        tabgp = uitabgroup(fig_ctrl,'Position', [0 0 1 .4]);
        
        tof_controls_tab = uitab(tabgp, 'Title', 'Region selection');
        addTOFControls(tof_controls_tab, [0 0 1 1]); 
        
        %Add time-varying dataset slider
        if (time_varying)
            data_controls_tab = uitab(tabgp, 'Title', 'Dataset selection');
            addTimestepControls(data_controls_tab, [0 0.7 1 0.3]);
        end
        
        vis_controls_tab = uitab(tabgp, 'Title', 'Plots');
        addVisualizationControls(vis_controls_tab, [0 .8 1 .2]); 
        
        advanced_controls_tab = uitab(tabgp, 'Title', 'Advanced');
        addAdvancedControls(advanced_controls_tab, [0 .4 1 .6]); 
        
        %Enable/disable the correct controls by firing callback on
        %time_integral button
        if (ishandle(ctrl_use_average))
            f = get(ctrl_use_average, 'Callback');
            f([], []);
        end
    end

    function computeAverage()
        if (~isempty(D_int))
            D = D_int;
        else
            extra_args = {};
            if (isfield(state{1}, 'time'))
                times = cellfun(@(x) x.time, state);
                times = [0; times];
                dt = (times(2:end) - times(1:end-1)) / times(end);
                extra_args = {extra_args{:}, 'dt', dt}; %#ok<CCAT>
            end
            
            min_tof = convertFrom(str2double(get(mtofeh, 'String')), year);
            max_tof = convertFrom(str2double(get(Mtofeh, 'String')), year);

            % Precompute all diagnostics states
            computeAllSteps();
            
            D_int = computeTOFandTracerAverage(state, computeGrid, rock, ...
                'wells', W, ...
                'max_tof', max_tof,...
                'min_tof', min_tof,...
                'diagnostics', Diagnostics, ...
                extra_args{:});
                
            D = D_int;
        end
    end

    function selection = getTOFSelection()
        %Get options for TOF region selection
        drain_wells = get(ctrl_drain_vols, 'Value');
        flood_wells = get(ctrl_flood_vols, 'Value');
        near_well_max_tof = convertFrom(get(nwtofsh, 'Value'), year);
        min_tof = convertFrom(str2double(get(mtofeh, 'String')), year);
        max_tof = convertFrom(str2double(get(Mtofeh, 'String')), year);
        switch(get(hset_op, 'Value'))
            case 1
                set_op = 'union';
            case 2
                set_op = 'intersection';
            case 3
                set_op = 'flood';
            case 4
            	set_op = 'drain';
        end
        
        %Get TOF region
        if (compute_average)
            % Get current TOF frequency threshold
            tof_freq = get(tofsh, 'Value');
            
            isubset = D.isubset;
            psubset = D.psubset;
            
            frequency = (psubset+isubset) / 2.0;
            psubset = psubset & (frequency >= tof_freq);
            isubset = isubset & (frequency >= tof_freq);
            
            selection = selectTOFRegion(D, max_tof, min_tof, ...
                                        'drain_wells', drain_wells, ...
                                        'flood_wells', flood_wells, ...
                                        'set_op', set_op, ...
                                        'near_well_max_tof', near_well_max_tof, ...
                                        'psubset', psubset, ...
                                        'isubset', isubset);
        else
            selection = selectTOFRegion(D, max_tof, min_tof, ...
                                        'drain_wells', drain_wells, ...
                                        'flood_wells', flood_wells, ...
                                        'set_op', set_op, ...
                                        'near_well_max_tof', near_well_max_tof);
        end
    end

    function plotMain(src, event)
        if isa(gcf, 'double') %pre 2014b
            grid_plot_opts = {'HitTest', 'off'};
        else
            grid_plot_opts = {'HitTest', 'off', 'PickableParts', 'none'};
        end
        if (compute_average)
            window_name = 'Average over timesteps';
        else
            if (numel(opt.name) == 1)
                window_name = [opt.name{1}, ' timestep ', num2str(state_idx)];
            else
                window_name = opt.name{state_idx};
            end
        end
        
        if(ishandle(fig_ctrl))
            set(fig_ctrl, 'Name', ['Controller ', window_name]);
        end            
            
        if ishandle(fig_main)
            set(0, 'CurrentFigure', fig_main);
        else
            fig_main = figure();
        end
        set(fig_main, 'Name', window_name);

        plot_grid = get(ctrl_show_grid, 'Value');
        if (plot_grid)
            if ~any(ishandle(fig_main_grid))
                fig_main_grid = plotGrid(G, 'facec', 'none', 'edgea', .05, 'edgec', 'black');
                set(fig_main_grid, 'UserData', 'gridoutline', grid_plot_opts{:});
            end
        else
            if (any(ishandle(fig_main_grid)))
                delete(fig_main_grid);
            end
        end
        
        if (~D.isvalid)
            warning('Time of flight is invalid.');
            clf();
            return;
        end
        
        tofext = getTOFRange(D, opt);
        setTOFSliderExtents();
        selection = getTOFSelection();
        
        dataind = get(hdataset, 'Value');
        [cdata, clim, ~] = selectDataset(dataind);

        alpha   = get(alfash, 'Value');

        % Plot current selection
        if any(ishandle(fig_main_cell_data))
            delete(fig_main_cell_data);
        end
        [az,el] = view;
        if any(selection)
            if (dataind == 6 || dataind == 7)
                fig_main_cell_data = plotTracerBlend(G, ...
                    cdataToPlotGrid(cdata{1}), ...
                    cdataToPlotGrid(cdata{2}), ...
                    'cells', cdataToPlotGrid(selection), ...
                    'FaceAlpha', alpha, grid_plot_opts{:});
            else
                fig_main_cell_data = plotCellData(G, ...
                    cdataToPlotGrid(cdata), cdataToPlotGrid(selection), ...
                    'EdgeColor', 'none', 'FaceAlpha', alpha, grid_plot_opts{:});
            end
        end

        if isempty(opt.fastRotate)
            if (G.cells.num > 50000)
                fastRotateButton();
            end
        elseif opt.fastRotate
            fastRotateButton();
        end

        plotWells();
        
        if (clim(2) > clim(1))
            caxis(clim);
        end
        
        %Update plot of wells (if enabled)
        plotWellConnections();
        
        %Try forcing manual axis when plotting to avoid rescaling
        axis(extents);
        view(az,el);
    end

    function stopApplication(src, event)
     if ishghandle(fig_well_alloc), close(fig_well_alloc); end
     if ishghandle(fig_phi),  close(fig_phi);  end
     if ishghandle(fig_well), close(fig_well); end
     if ishghandle(fig_main), close(fig_main); end
     if ishghandle(fig_ctrl), close(fig_ctrl); end
     return
    end

    function playBackTof(src, event)
        if playback
            return
        end
        N = round(get(speedsh, 'Value'));
        
        if (compute_average)
            t_0 = 1;
            t_end = get(tofsh, 'Value');
            callback = get(tofsh, 'Callback');
        else
            t_0 = get(mtofsh, 'Value');
            t_end = get(Mtofsh, 'Value');
            callback = get(Mtofsh, 'Callback');
            if t_0 >= t_end
                disp 'Minimum TOF is larger than maximum TOF...'
                return
            end
        end
        playback = true;
        % Pop figure to front explicitly before loop
        figure(fig_main);
        
        for i = 0:N
            if ~playback
                return
            end
            
            if (compute_average)
                set(tofsh, 'Value', t_0 + i*(t_end - t_0)/N);
                callback(tofsh, []);
            else
                set(Mtofsh, 'Value', min(t_0 + i*(t_end - t_0)/N, tofext(2)));
                callback(Mtofsh, []);
            end
            
            timer = tic();
            plotMain();
            drawnow
            % The whole event should take a minimum of 10 seconds
            pause(min(0, 10/N - toc(timer)));
        end
        
        playback = false;
    end

    function stopPlayBackTof(src, event)
        playback = false;
    end

    function onClickWell(src, event, wk)
        clickType = get(gcf, 'SelectionType');
        if ~isempty(opt.wellPlotFn) && strcmpi(clickType, 'alt')
            t = state_idx;
            if isfield(state{state_idx}, 'time')
                t = state{state_idx}.time;
            end
            opt.wellPlotFn(src, event, W{state_idx}(wk).name, t);
            return
        end
        
        %np = numel(D.prod);
        %ni = numel(D.inj);
            
        %winj = repmat(D.inj,[1,np]);
        %wpro = rldecode(D.prod,ni,2);
        
        wp = [];

        isInj = ismember(wk, D.inj);
        if isInj
            ik = find(D.inj == wk);
            %sub = winj == wk;

            % well pair stuff
            if (numel(WP) > 0)
                wp = @(x) WP.inj(x);
            end
            otherNames = {W{state_idx}(D.prod).name, 'reservoir'};

            % set plots to match piecharts
            v = find(strcmpi(get(hdataset, 'String'), 'drainage region'));
            set(hdataset, 'Value', v(1))
            set(ctrl_flood_vols, 'Value', ik)
            set(ctrl_drain_vols, 'Value', []);
        else
            ik = find(D.prod == wk);
            %sub = wpro == wk;

            % well pair stuff
            if (numel(WP) > 0)
                wp = @(x) WP.prod(x);
            end
            otherNames = {W{state_idx}(D.inj).name, 'reservoir'};

            v = find(strcmpi(get(hdataset, 'String'), 'flooding region'));
            set(hdataset, 'Value', v(1))
            set(ctrl_drain_vols, 'Value', ik)
            set(ctrl_flood_vols, 'Value', []);
        end

        if isempty(fig_well) || ~ishandle(fig_well)
            fig_well = figure();
        else
            set(0, 'CurrentFigure', fig_well); clf
        end
        set(fig_well, 'name', ['Well ', W{state_idx}(wk).name]);

        plotArrival = ~isempty(state{state_idx}) && ~isInj;

        if plotArrival;
            subplot(2, 2, 1);
        else
            subplot(2, 2,[1 3])
        end

        if (numel(WP) > 0)
            tmp = wp(ik);
            %salloc = abs(sum(tmp.alloc,1));
            alloc = [tmp.alloc, tmp.ralloc];
            salloc = abs(sum(alloc,1));
            [sa, tmp_ix] = sort(salloc, 'descend');
            ssa = sum(sa);
            ix = find(cumsum(sa) >= (1-opt.ignoreWellRate)*ssa, 1, 'first');
            [sa, tmp_ix] = deal(sa(ix:-1:1), tmp_ix(ix:-1:1));
            if ~isempty(sa) && ssa~=0
                pie(sa/ssa,otherNames(tmp_ix));
                %pie(max(WP.vols(sub), eps), ones(size(WP.vols(sub))))
                title('Well allocation factors')
                if isInj
                    set(ctrl_drain_vols, 'Value', tmp_ix(tmp_ix<=numel(D.prod)));
                else
                    set(ctrl_flood_vols, 'Value', tmp_ix(tmp_ix<=numel(D.inj)));
                end
            end
            fig_main_wells.dirty = true;
        end
        if plotArrival && ~isempty(tmp_ix)
            subplot(2,2,3);  cla;
            plotTOFArrival(state{state_idx}, W{state_idx}, pv, opt.fluid, find(D.prod == wk), D, 'inj_ix', tmp_ix, 'maxTOF', opt.maxTOF);
        end

        
        if (numel(WP) > 0)
            if numel(tmp.z) > 1
                subplot(2,2,[2 4])
                % Allocation factors by connection does not make sense for only one
                % perforation!
                %[z, zind] = sort(tmp.z, 'descend');
                %alloc  = abs(tmp.alloc(:, tmp_ix));
                alloc  = abs(alloc(:, tmp_ix));
                na = size(alloc,1);
                % reverse cumsum
                calloc = cumsum(flipud(alloc), 1);
                calloc = flipud(calloc);
                % Use standard units
                calloc = convertTo(calloc, 1/day);
                if ~isempty(calloc)
                    area((1:na)', calloc); axis tight;
                    hold on
                    plot((1:na)', zeros(na,1), '>', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'MarkerSize', 5);
                    % Flip it around, xlabel is really ylabel
                    view(90, -90);
                    set(gca, 'XDir', 'reverse')
                    xlabel('Connection #')
                    ylabel('Accumulated flux [m^3/day]')
                    title('Allocation by connection');
                end
            end
        end

        plotMain();
    end

    function plotPhi(src, event)
        if (~D.isvalid)
            warndlg('It appears that the time of flight is invalid. Aborting.');
            return
        end
        
        if isempty(fig_phi) || ~ishandle(fig_phi)
            fig_phi = figure();
        else
            set(0, 'CurrentFigure', fig_phi);
        end
        set(fig_phi, 'Name', ['Phi/F diagram (', window_name, ')']);
        [F,Phi] = computeFandPhi(pv,D.tof);
        plot(Phi,F,'.');
        title(sprintf('Lorenz coefficient: %f\n', computeLorenz(F,Phi)));
    end

    function plotWells()
        % Plot wells
        if (fig_main_wells.dirty)
            if isempty(fig_main_wells.hwells) || ~ishandle(fig_main_wells.hwells(1).label)
                fig_main_wells.hwells = plotWellData(G, wellsToPlotGrid(W{state_idx}), ...
                    'color', [0 0 0], 'linePlot', opt.lineWells);
                for j = 1:numel(W{state_idx})
                    color = colorizeWell('global', j, D);
                    set([fig_main_wells.hwells(j).label, ...
                         fig_main_wells.hwells(j).connector,...
                         fig_main_wells.hwells(j).body], 'ButtonDownFcn', @(src, event) onClickWell(src, event, j));
                    set([fig_main_wells.hwells(j).label, fig_main_wells.hwells(j).connector], 'Color', color);
                    set(fig_main_wells.hwells(j).label, 'FontWeight', 'bold', 'Interpreter', 'none')
                    for k = 1:numel(fig_main_wells.hwells(j).body)
                        if any(strcmp(get(fig_main_wells.hwells(2).body(k), 'type'), {'patch', 'surface'}))
                            set(fig_main_wells.hwells(j).body(k), 'FaceColor', color, 'EdgeColor', 'none');
                        else
                            set(fig_main_wells.hwells(j).body(k), 'Color', color);
                        end
                    end
                end
            end
    
            %Plot all wells
            if (get(ctrl_plot_all_wells, 'Value'))
                sel_wells = [D.prod, D.inj];
            else %Or the selected wells only
                drain_wells = D.prod(get(ctrl_drain_vols, 'Value'));
                flood_wells = D.inj(get(ctrl_flood_vols, 'Value'));
                sel_wells = [drain_wells, flood_wells];
            end
            
            vis = repmat({'off'}, [numel(W{state_idx}), 1]);
            [vis{sel_wells}] = deal('on');
            for j = 1: numel(W{state_idx})
                set([fig_main_wells.hwells(j).label, ...
                     fig_main_wells.hwells(j).connector, ...
                     fig_main_wells.hwells(j).body], 'Visible', vis{j});
            end          
            fig_main_wells.dirty = false;
        end
    end

    function plotWellAllocations(~, ~)
        if (~D.isvalid)
            warndlg('It appears that the time of flight is invalid. Aborting');
            return
        end
        
        if isempty(fig_well_alloc) || ~ishandle(fig_well_alloc)
            fig_well_alloc = figure();
        else
            set(0, 'CurrentFigure', fig_well_alloc);
        end
        set(fig_well_alloc, 'Name', ['Well allocations (', window_name, ')']);
        
        %Plot all wells
        if (get(ctrl_plot_all_wells, 'Value'))
            plotWellAllocationPanel(D, WP);
        else %Or the selected wells only
            drain_wells = D.prod(get(ctrl_drain_vols, 'Value'));
            flood_wells = D.inj(get(ctrl_flood_vols, 'Value'));
            sel_wells = [drain_wells, flood_wells];
            plotWellAllocationPanel(D, WP, 'plotOnly', sel_wells);
        end
    end

    function plotWellConnections(~, ~)
        if get(ctrl_well_conn, 'Value')
            if ishandle(fig_main)
                set(0, 'CurrentFigure', fig_main);
            else
                fig_main = figure;
                axis tight off
                plotMain();
            end
            
            if any(ishandle(fig_main_well_conns))
                delete(fig_main_well_conns)
                fig_main_well_conns = [];
            end
            
            if (D.isvalid)
                fig_main_well_conns = plotWellPairConnections(computeGrid, WP, D, W{state_idx}, pv);
            else
                return
            end
        else
            delete(fig_main_well_conns)
            fig_main_well_conns = [];
        end
    end

    function [cdata, clim, cmap] = selectDataset(dataind)
        datanames = get(hdataset, 'String');
        clim = [];
        cmap = 'jet';
        
        switch lower(datanames{dataind})
            case 'forward tof'
                cdata = log10(D.tof(:,1));
                clim = log10(clamp_real(convertFrom(tofext, year)));
            case 'backward tof'
                cdata = log10(D.tof(:,2));
                clim = log10(clamp_real(convertFrom(tofext, year)));
            case 'sum of tofs'
                data = sum(D.tof(:,1:2), 2);
                cdata = log10(data);
                clim = log10(clamp_real(convertFrom([min(data), tofext(2)], year)));
            case 'drainage region'
                cmap = 'gray';
                cdata = D.ppart;
                clim = clamp_real([1, max(D.ppart)]);
            case 'flooding region'
                cdata = D.ipart;
                clim = clamp_real([1, max(D.ipart)]);
            case 'drainage blend'
                cdata = {D.ppart, max(D.ptracer, [], 2)};
                clim = [1, max(D.ipart)];
            case 'flooding blend'
                cdata = {D.ipart, max(D.itracer, [], 2)};
                clim = [1, max(D.ipart)];
            case 'tracer selected injectors'
                cdata  = sum(D.itracer(:, get(ctrl_flood_vols, 'Value')), 2);
            case 'tracer selected producers'
                cdata  = sum(D.ptracer(:, get(ctrl_drain_vols, 'Value')), 2);
            case 'tof selected injectors'
                cdata = log10(tracerAveragedTOF(1));
                clim =  log10(clamp_real(convertFrom(tofext, year)));
            case 'tof selected producers'
                cdata = log10(tracerAveragedTOF(-1));
                clim =  log10(clamp_real(convertFrom(tofext, year)));
            case 'porosity'
                cdata = rock.poro;
            case 'x permeability'
                cdata = log10(rock.perm(:,1));
            case 'y permeability'
                cdata = log10(rock.perm(:,2));
            case 'z permeability'
                cdata = log10(rock.perm(:,3));
            case 'forward tof frequency'
                cdata = D.psubset;
            case 'backward tof frequency'
                cdata = D.isubset;
            case 'sum tof frequency'
                cdata = D.isubset + D.psubset;
            case 'state.s'
                if size(state{state_idx}.s, 2) == 3
                    cdata = state{state_idx}.s(:, [2 3 1]);
                else
                    cdata = state{state_idx}.s;
                end
            case 'state.mob'
                if size(state{state_idx}.s, 2) == 3
                    cdata = state{state_idx}.mob(:, [2 3 1]);
                else
                    cdata = state{state_idx}.mob;
                end                 
            otherwise
                cdata = readStructField(datasets{state_idx}, datanames{dataind});
        end
        
        if isempty(clim)
            m = min(cdata(:));
            M = max(cdata(:));

            clim = [m - eps(m), M + eps(m)];
        end
    end

    function tof = tracerAveragedTOF(sgn)
        if sgn > 0
            tr  = D.itracer(:, get(ctrl_flood_vols, 'Value'));
            tof = sum(D.itof(:, get(ctrl_flood_vols, 'Value')).*tr, 2);
        else
            tr  = D.ptracer(:, get(ctrl_drain_vols, 'Value'));
            tof = sum(D.ptof(:, get(ctrl_drain_vols, 'Value')).*tr, 2);
        end
        str = sum(tr, 2); 
        ok  = and(tof < str*opt.maxTOF, str > sqrt(eps));
        tof(ok)  = tof(ok)./str(ok);
        tof(~ok) = opt.maxTOF;
    end
        
    function [D, WP] = computeTOFAndTracerAndWellPairs(W_step, state_step)
        if opt.computeFlux
            rS = initState(G, W_step, 0);
            T  = computeTrans(G, rock);
            rS = incompTPFA(rS, G, T, opt.tracerfluid, 'wells', W_step, 'LinSolve', opt.LinSolve);
        else
            rS = initState(G, W_step, 0);
            if isfield(state_step, 'wellSol')
                rS.wellSol = state_step.wellSol;
            end
            rS.flux     = state_step.flux;
            rS.pressure = state_step.pressure;
        end
        
        D = computeTOFandTracer(rS, computeGrid, rock, 'wells', W_step, 'processCycles', true, 'maxTOF', opt.maxTOF, ...
                                'computeWellTOFs', opt.computeWellTOFs);
        D.itracer(isnan(D.itracer)) = 0;
        D.ptracer(isnan(D.ptracer)) = 0;
        
        % Cap tof to maximum tof for unreachable areas for the time being
        tf = D.tof(:);
        if (all(isinf(tf)))
            warndlg('TOF invalid: all values are inf.');
            D.isvalid = false;
            WP = [];
        elseif (numel(D.inj) * numel(D.prod) == 0)
            warndlg('Missing at least one injecting or producing well.');
            D.isvalid = false;
            WP = [];
        else
            D.tof(isinf(D.tof)) = max(tf(isfinite(tf)));
            D.isvalid = true;
            WP = computeWellPairs(rS, computeGrid, rock, W_step, D);
        end
    end

    function computeValues(stepNo, recompute)
        if isempty(Diagnostics{stepNo}) || recompute
            disp('New state encountered, computing diagnostics...');
            [d, wp] = ...
                computeTOFAndTracerAndWellPairs(W{stepNo}, state{stepNo});
            Diagnostics{stepNo} = d;
            WellPairs{stepNo}   = wp;
            if (~d.isvalid || isempty(wp))
                warning('Time of flight returned inf. Are there both active injectors and producers present?')
            end
        end
    end

    function computeAllSteps()
        h = waitbar(0, 'Computing diagnostics...');
        for ix = 1:numel(state)
            if isempty(Diagnostics{ix})
                waitbar(ix/numel(state), h, 'Computing diagnostics...');
            end
            computeValues(ix, false);
        end
        if ishandle(h)
            close(h);
        end
    end

    function [D, WP] = getDiagnostics(stepNo, recompute)
        if nargin == 1
            recompute = false;
        end
        computeValues(stepNo, recompute);
        D = Diagnostics{stepNo};
        WP = WellPairs{stepNo};
    end

    function changeWells()
        W{state_idx} = editWells(G, W{state_idx}, rock);
        [D, WP] = getDiagnostics(state_idx, true);
        createMainControl();
        
        %If we changed the number of wells, reset the figure handles
        if (numel(fig_main_wells.hwells) ~= numel(W{state_idx}))
            for i=1:numel(fig_main_wells.hwells)
                delete(fig_main_wells.hwells(i).body);
                delete(fig_main_wells.hwells(i).connector);
                delete(fig_main_wells.hwells(i).label);
            end
            fig_main_wells = {};
            fig_main_wells.hwells = [];
        end
        fig_main_wells.dirty = true;
        
        plotMain();
    end



    function [sliderhandle, edithandle, grouphandle] = linkedSlider(parent, pos, fieldsize, ext, defaultval, title, callback)
        minval = abs(ext(1));
        maxval = abs(ext(2));
        defaultval = abs(defaultval);

        grouphandle = uipanel(parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'None');
        uicontrol(grouphandle, 'Style', 'text', 'Units', 'normalized', 'Position', [0, 0, fieldsize 1], 'string', title)

        cap = @(x) max(minval, min(x, maxval));

        edithandle = uicontrol(grouphandle, 'Style', 'edit', 'Units', 'normalized', 'Value', defaultval, 'Position', [0 + (1-fieldsize), 0, fieldsize 1], 'String', sprintf('%.1f', defaultval));
        sliderhandle = uicontrol(grouphandle, 'Style', 'slider', 'Units', 'normalized', 'Position', [fieldsize, 0, (1-2*fieldsize) 1], 'Min', minval, 'Max', maxval, 'Value', defaultval);

        slidercallback2 = @(src, event) set(edithandle, 'String', sprintf('%.1f', get(src, 'Value')));
        editcallback2 = @(src, event) set(sliderhandle, 'Value', cap(sscanf(get(src, 'String'), '%f')));

        function slidercallback3(src, event)
            callback(src, event);
            slidercallback2(src, event);
        end
        function editcallback3(src, event)
            callback(src, event);
            editcallback2(src, event);
        end

        if (~isempty(callback))
            set(sliderhandle, 'Callback', @slidercallback3);
        else 
            set(sliderhandle, 'Callback', slidercallback2);
        end

        if (~isempty(callback))
            set(edithandle, 'Callback', @editcallback3);
        else
            set(edithandle, 'Callback', editcallback2);
        end
    end


    % Private helpers
    function close_all_func(src, event)
        if (~opt.leaveOpenOnClose)
            if any(ishandle(fig_ctrl))
                delete(fig_ctrl);
            end
            if any(ishandle(fig_well_alloc))
                delete(fig_well_alloc);
            end
            if any(ishandle(fig_phi))
                delete(fig_phi);
            end
            if any(ishandle(fig_main))
                delete(fig_main);
            end
            if any(ishandle(fig_well))
                delete(fig_well);
            end
        end
    end

end
