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
%     for visualization.
%     - If a state is given, this state will allow for visualization of
%     the component ratios inside a drainage volume. The flux from this
%     state can also be used to calculate time of flight if "computeFlux"
%     is disabled, for instance if the user has some external means of
%     computing fluxes.
%
%   Once the initialization is complete, two windows will be produced:
%     - A plotting window, showing the reservoir along with the wells and
%     visualized quantitites.
%         In the plotting window, it is possible to click wells to get
%         additional information, such as the allocation factors per
%         perforation in the well, the phase distribution grouped by time
%         of flight and the corresonding pore volumes swept/drained.
%
%     - A controller window which is used to alter the state of the
%     plotting window:
%     Wells can be selected for display (if a well is selected, the
%     corresponding drainage (producer) or flooding (injector) volumes will
%     be visualized. A simple playback function can be used to show
%     propagation of time of flight. Different quantitites can be
%     visualized to get a better understanding of the system.
%
%     A set of buttons provide (experimental) well editing, access to
%     visualization of the well pair connections, Phi/F diagram with
%     Lorenz' coefficient etc.
%
%     A player controller (experimental), if a time-series of wells (and
%     states) are supplied.
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
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%  'state' - Reservoir state containing fluid saturations and optionally
%            flux and pressure (if 'computeFlux' is false)
%
%  'computeGrid' - Use G for plotting, but computeGrid for computing time
%                  of flight etc.
%
%  'tracerfluid' - Fluid used for tracer computation. This defaults to a
%            trivial fluid.
%
%  'fluid'  - Fluid used for mobility calculations if 'useMobilityArrival'
%             is enabled.
%
%  'LinSolve' - Linear solver for pressure systems. Defaults to mldivide
%  (backslash)
%
%  'computeFlux' - If set to false, fluxes are extracted from the provided
%                  state keyword argument. This requires a state to be
%                  provided and can be used if the fluxes are computed
%                  externally (for instance from a expensive full-physics
%                  simulation)
%
%  'useMobilityArrival' - If the well plot showing nearby saturations
%                         should plot mobility instead of saturations. This
%                         may be interesting in some cases because the
%                         mobile fluids are more likely to be extracted.
%                         However, this plot is often dominated by very
%                         mobile gas regions.
%
%  'daspect' - Data aspect ratio, in a format understood by daspect()
%
%  'name' - Name to use for windows
%
%  'leaveOpenOnClose' - Default false. Leaves all figures open when closing
%                       the controller.
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
%   Diagnostics examples

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
                 'useMobilityArrival',  false,...
                 'fluid',               oilwater, ...
                 'name',                [], ...
                 'daspect',             [], ...
                 'leaveOpenOnClose',    false ...
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
        cdataToPlotGrid = @(cdata) cdata(G.cells.eMap);
        wellsToPlotGrid = @(wells) remapWells(G, wells);
        computeGrid = opt.computeGrid;
    end
        
    state = [];
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
    

    % Currently playing back
    playback = false;

    % Main scope variables for tracers and tof
    [D, WP] = deal([]);

    % Main scope variables for control panel
    [fig_ctrl,...
     wellPlot,...
     draintoggle,...
     floodtoggle,...
     wctoggle,...
     speedsh,...
     tofsh,...
     nwtofsh,...
     mtofsh,...
     mtofeh,...
     Mtofsh,...
     Mtofeh,...
     alfash,...
     tofext,...
     hdataset,...
     hset_op,...
     time_integral,...
     mrst_ds,...
     ds_panel,...
     ni, np] = deal(NaN);
    
    % Setup "persistent" variables
    pm_htop = [];
    pm_htext = [];
    pm_hs = [];
    pm_hline = [];
    pm_mainph = [];
    pm_outlineph = [];
    phiPlot = []; % Used in plotPhi
    pwc_wch = []; % Used in plotWellConnections
    wah_fig = []; % plotWellAllocations
    selection = []; % Selection of cells
    
    time_integral_handles = [];
    static_handles = [];
    
    %State variables
    compute_time_integral = false;
    D_int = [];
    window_name = opt.name{1};
    
    %Create main figure
    fig_main = figure('Name', window_name);
    
    axis tight off
    if (~isempty(opt.daspect))
        daspect(opt.daspect);
    end
    view(3);

    % Create control panel
    % Precompute TOF etc. for creating main control
	pv = poreVolume(computeGrid, rock);
    [D, WP] = computeTOFAndTracerAndWellPairs(W{state_idx}, state{state_idx});
    tofext = getTOFRange(D);
    createMainControl();
    [D.isubset, D.psubset] = getSubset(D);
    
    % Trigger plot initial setup
    plotMain();

    % Private helpers
    function ctrl_close_func(src, event)
        if (~opt.leaveOpenOnClose)
            if any(ishandle(pm_mainph))
                delete(pm_mainph);
            end
            if any(ishandle(fig_ctrl))
                delete(fig_ctrl);
            end
            if any(ishandle(pwc_wch))
                delete(pwc_wch);
            end
            if any(ishandle(wah_fig))
                delete(wah_fig);
            end
            if any(ishandle(phiPlot))
                delete(phiPlot);
            end
            if any(ishandle(fig_main))
                delete(fig_main);
            end
            if any(ishandle(wellPlot))
                delete(wellPlot);
            end
        end
    end

    function val = clamp_real(in_val)
        val = max(-realmax, min(realmax, in_val));
    end

    function tofext = getTOFRange(D)
        if (D.isvalid)
            tof = D.tof(:);
            isNeg = (tof <= 0);
            tof(isNeg) = min(tof(~isNeg));

            tofext = convertTo(([min(tof), 5*10^(mean(log10(tof)))]), year);
            tofext(2) = clamp_real(max(tofext(2), 15));
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
            if (compute_time_integral)
                D_int = [];
                computeTimeIntegral();
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
            if (compute_time_integral)
                D_int = [];
                computeTimeIntegral();
            end
            plotMain(src, event);
        end
    end


    function addVisualizationControls(parent, pos)
        vis_N = 3;
        vis_H = 1 / vis_N;
        vis_h = vis_H*0.9;
        bwidth = .2;
        
        vis_control = uipanel('Parent', parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'None');
        
        [alfash, ~, ~]   = linkedSlider(vis_control, [0 2*vis_H 1 vis_h], .15, [0 1], 1, 'Alpha', @plotMain);
        
        
        wctoggle = uicontrol(vis_control, 'Style', 'checkbox',...
                   'Units', 'normalized',...
                   'Position', [0 1*vis_H bwidth vis_h],...
                   'String', 'Well pairs',...
                   'Callback', @plotWellConnections ...
                   );
               
               
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
        function dirtyPlotMain(src, event)
            if (compute_time_integral)
                D_int = [];
                computeTimeIntegral();
            end
            plotMain(src, event);
        end
        
        
        % Panel for Time of flight selection parameters
        tof_panel = uipanel('Parent', parent, 'Units', 'normalized', 'Position', pos, 'BorderType', 'none');
        tof_N = 7;
        tof_H = 1/tof_N;
        tof_h = tof_H*0.9;
        
        cellfields = getStructFields(G, datasets{state_idx}, dsname);
        perm = strcat({'X', 'Y', 'Z'}, ' permeability');
        uicontrol(tof_panel, 'Style', 'text', ...
                             'string', 'Display:', ...
                             'HorizontalAlignment', 'left', ...
                             'Units', 'normalized',...
                             'Position', [.0 6*tof_H .25 tof_h]);
                
        hdataset = uicontrol(tof_panel, 'Style', 'popupmenu',...
                                    'Units', 'normalized',...
                                    'Position', [.2 6*tof_H .33 tof_h],...
                                    'Callback', @plotMain,...
                                    'String', {'Forward TOF',...
                                               'Backward TOF',...
                                               'Sum of TOFs',...
                                               'Drainage region',...
                                               'Flooding region',...
                                               'Drainage blend',...
                                               'Flooding blend',...
                                               'Porosity', ...
                                               perm{1:size(rock.perm, 2)},...
                                               cellfields{:}}...
                                   );
        hset_op = uicontrol(tof_panel, 'Style', 'popupmenu',...
                                   'Units', 'normalized',...
                                   'Position', [.66 6*tof_H .33 tof_h],...
                                   'Callback', @plotMain,...
                                   'String', {'Union {Flood, Drain} volumes',...
                                              'Intersection {Flood, Drain} volumes',...
                                              'Flood volumes',...
                                              'Drain volumes'}...
                                   );
        
        [nwtofsh, ~, ~]   = linkedSlider(tof_panel, [0 5*tof_H 1 tof_h], .15, [0 1], 0.2, 'Near well TOF', @dirtyPlotMain);
        [Mtofsh, Mtofeh, ~]   = linkedSlider(tof_panel, [0 4*tof_H 1 tof_h], .15, tofext, tofext(2), 'Max TOF', @dirtyPlotMain);
        [mtofsh, mtofeh, ~]   = linkedSlider(tof_panel, [0 3*tof_H 1 tof_h], .15, tofext, tofext(1), 'Min TOF', @dirtyPlotMain);
        set(mtofsh, 'Callback', @mtofs_callback);
        set(mtofeh, 'String', num2str(-Inf))
        set(Mtofsh, 'Callback', @Mtofs_callback);
        set(Mtofeh, 'String', num2str(+Inf))
        
        [tofsh, tofeh, ~]   = linkedSlider(tof_panel, [0 2*tof_H 1 tof_h], .15, [0 1], 0.25, 'TOF Freq.', @plotMain);
        time_integral_handles = [time_integral_handles(:); tofsh; tofeh];
        
        
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 1*tof_H .2 tof_h],...
                   'String', 'Play TOF',...
                   'Callback', @playBackTof...
                   );
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [.2 1*tof_H .2 tof_h],...
                   'String', 'Stop TOF',...
                   'Callback', @stopPlayBackTof...
                   );
        [speedsh, ~, ~] = linkedSlider(tof_panel, [.45 1*tof_H .55 tof_h], .15, [1 1000], 50, 'Steps', []);  
        
        
        
        % Function which saves variables to the workspace on demand
        function saveDataToWorkSpace(src, event)
            indices = find(selection);
            
            labels = {'Sector model', 'Sector model (index map)'};
            varnames = labels;
            vars = {indices, computeGrid.cells.indexMap(indices)};
            
            %Clean the varnames of spaces etc
            clean = @(x) lower(regexprep(x,'[^a-zA-Z0-9]','_'));
            varnames = cellfun(@(x) clean(x), varnames, 'UniformOutput', false);

            export2wsdlg(labels, varnames, vars)
        end
        
        uicontrol(tof_panel, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'Position', [0 0*tof_H 0.3 tof_h],...
                   'String', 'Export selection', ...
                   'Callback', @saveDataToWorkSpace...
                   );
    end


    function addTimestepSelector(parent, pos)
        ds_panel = uipanel('Parent', parent,...
                'Units', 'normalized',...
                'Position', pos,...
                'BorderType', 'none');
            
        uicontrol(ds_panel,'Style','text',...
            'Units', 'normalized',...
            'HorizontalAlignment', 'left',...
            'String','Simulation timestep',...
            'Position', [0 0.5 0.4 .5]);

        selector_datasets = [];
        if (numel(state) > 1)
            selector_datasets = state;
        elseif(numel(W) > 1)
            selector_datasets = W;
        else
            error('Programmer error');
        end
        
        function selecttimestep(datas, index, fullname)
            state_idx = index;
            computeValues(state_idx);
            
            %Update main plot
            plotMain();
            
            %Update aux plots
            if (ishandle(phiPlot))
                plotPhi([], []);                
            end
            if (ishandle(wah_fig))
                plotWellAllocations([], [])
            end
        end

        datasetSelector(G, selector_datasets, 'Parent', ds_panel, 'Location', ...
            [0.4, 0.5, 0.6, .5], 'Callback', @selecttimestep, ...
            'active', 1, 'Nofields', true);
        mrst_ds = findobj('Tag', 'mrst-datasetselector');
        set(mrst_ds, 'BorderType', 'None');
        mrst_ds_children = allchild(mrst_ds);
        static_handles = [static_handles(:); mrst_ds_children];

        %Control for showing / computing average over time of e.g., TOF
        time_integral = uicontrol(ds_panel, 'Style', 'checkbox',...
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
            if (ishandle(time_integral))
                compute_time_integral = logical(get(time_integral, 'Value'));
                if (compute_time_integral)
                    set(time_integral_handles, 'Enable', 'on');
                    set(static_handles, 'Enable', 'off');
                    
                    %Add new integral options
                    old_options = get(hdataset, 'String');
                    new_options = {old_options{:}, time_integral_options{:}};
                    set(hdataset, 'String', new_options);
                    
                    computeTimeIntegral();
                else
                    set(time_integral_handles, 'Enable', 'off');
                    set(static_handles, 'Enable', 'on');

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
                    
                    computeValues(state_idx);
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
            fig_ctrl = figure('Position',[pos_xy,  size_xy], 'Toolbar','none', 'MenuBar', 'none', 'CloseRequestFcn', @ctrl_close_func);
            set(fig_ctrl, 'Name', ['Controller ', window_name]);
        else
            clf(fig_ctrl)
        end

        % Drainage / flooding controls
        pp = unique(D.ppart);
        ip = unique(D.ipart);

        np = numel(pp);
        ni = numel(ip);

        draintoggle = zeros(np, 1);
        floodtoggle = zeros(ni, 1);
        
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

        draintoggle = uicontrol(drainp, 'Style', 'listbox', ...
                                          'Units', 'normalized',...
                                          'Max', 2, ...
                                          'Min', 0, ...
                                          'String', {W{state_idx}(D.prod).name},...
                                          'Value', [], ...
                                          'Callback', @plotMain, ...
                                          'Position', [0 0 1 1]);

        floodtoggle = uicontrol(floodp, 'Style', 'listbox', ...
                                          'Units', 'normalized',...
                                          'Max', 2, ...
                                          'Min', 0, ...
                                          'String', {W{state_idx}(D.inj).name},...
                                          'Value', 1:numel(D.inj), ...
                                          'Callback', @plotMain, ...
                                          'Position', [0 0 1 1]);
        
                                      
        tabgp = uitabgroup(fig_ctrl,'Position', [0 0 1 .4]);
        
        tof_controls_tab = uitab(tabgp, 'Title', 'Time of flight controls');
        addTOFControls(tof_controls_tab, [0 .1 1 .9]); 
        
        vis_controls_tab = uitab(tabgp, 'Title', 'Visualization config');
        addVisualizationControls(vis_controls_tab, [0 .6 1 .4]); 
        
        %Add time-varying dataset slider
        if ((~isempty(state) && numel(state) > 1) || (~isempty(W) && numel(W) > 1))
            data_controls_tab = uitab(tabgp, 'Title', 'Dataset selection');
            addTimestepSelector(data_controls_tab, [0 0.7 1 0.3]);
        end
        
        %Enable/disable the correct controls by firing callback on
        %time_integral button
        f = get(time_integral, 'Callback');
        f([], []);
    end

    function [isubset, psubset] = getSubset(D)
        % Get current TOF settings
        min_tof = convertFrom(str2double(get(mtofeh, 'String')), year);
        max_tof = convertFrom(str2double(get(Mtofeh, 'String')), year);
        
        % Drainage volumes
        data = D.tof(:,2);
        psubset = data >= min_tof & data <= max_tof;

        % Flooding volumes
        data = D.tof(:,1);
        isubset = data >= min_tof & data <= max_tof;
    end

    function computeTimeIntegral()
        if (~isempty(D_int))
            D = D_int;
        else
            extra_args = {};
            if (isfield(state{1}, 'time'))
                times = cellfun(@(x) x.time, state);
                times = [0; times];
                dt = (times(2:end) - times(1:end-1)) / times(end);
                extra_args = {extra_args{:}, 'dt', dt};
            end
            
            min_tof = convertFrom(str2double(get(mtofeh, 'String')), year);
            max_tof = convertFrom(str2double(get(Mtofeh, 'String')), year);
            D_int = computeTOFandTracerAverage(state, computeGrid, rock, ...
                'wells', W, ...
                'max_tof', max_tof,...
                'min_tof', min_tof,...
                extra_args{:});
                
            D = D_int;
        end
    end


    function plotMain(src, event)
        if (compute_time_integral)
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

        if any(ishandle(pm_mainph))
            delete(pm_mainph);
        end

        if ~any(ishandle(pm_outlineph))
            pm_outlineph = plotGrid(G, 'facec', 'none', 'edgea', .05, 'edgec', 'black');
            set(pm_outlineph, 'UserData', 'gridoutline');
        end
        
        if (compute_time_integral)
            tofext = getTOFRange(D);
            isubset = D.isubset;
            psubset = D.psubset;
            
            % Get current TOF frequency threshold
            frequency = (psubset+isubset) / 2.0;
            tof_freq = get(tofsh, 'Value');
            psubset = psubset & (frequency >= tof_freq);
            isubset = isubset & (frequency >= tof_freq);
            
            setTOFSliderExtents();
        else        
            [isubset, psubset] = getSubset(D);
            tofext = getTOFRange(D);
            setTOFSliderExtents();
        end
        
        %Get setting for drain/flood wells
        drain_wells = get(draintoggle, 'Value');
        flood_wells = get(floodtoggle, 'Value');
        
        psubs = ismember(D.ppart, drain_wells);
        isubs = ismember(D.ipart, flood_wells);
        
        % Use logical operation on the selection
        switch(get(hset_op, 'Value'))
            case 1
            selection = (isubset & isubs) | (psubset & psubs);
            case 2
            selection = (isubset & isubs) & (psubset & psubs);
            case 3
            selection = (isubset & isubs);
            case 4
            selection = (psubset & psubs);
        end
        
        %Grow selection to include near-well regions
        near_well_tof = convertFrom(get(nwtofsh, 'Value'), year);
        isubset = D.tof(:,1) <= near_well_tof;
        psubset = D.tof(:,2) <= near_well_tof;
        near_well_selection = (isubset & isubs) | (psubset & psubs);
        selection  = near_well_selection | selection;
        
        
        dataind = get(hdataset, 'Value');
        [cdata, clim, ~] = selectDataset(dataind);

        alpha   = get(alfash, 'Value');

        % Plot current selection
        if any(selection)
            if (dataind == 6 || dataind == 7)
                pm_mainph = plotTracerBlend(G, ...
                    cdataToPlotGrid(cdata{1}), ...
                    cdataToPlotGrid(cdata{2}), ...
                    'cells', cdataToPlotGrid(selection), ...
                    'FaceAlpha', alpha);
            else
                pm_mainph = plotCellData(G, ...
                    cdataToPlotGrid(cdata), cdataToPlotGrid(selection), ...
                    'EdgeColor', 'none', 'FaceAlpha', alpha);
            end
        end

        fastRotateButton();

        % Plot wells
        if ~all(ishandle([pm_htop, pm_htext, pm_hs])) || isempty([pm_htop, pm_htext, pm_hs])

            [pm_htop, pm_htext, pm_hs, pm_hline] = ...
                plotWell(G, wellsToPlotGrid(W{state_idx}), ...
                    'color', 'red', 'height',  0);
            for i = 1:numel(W{state_idx})
                color = colorizeWell('global', i, D);
                set([pm_htop(i) pm_htext(i) pm_hs(i)], 'ButtonDownFcn', @(src, event) onClickWell(src, event, i));
                set([pm_htop(i) pm_hs(i)], 'FaceColor', color, 'EdgeColor', color)
                set([pm_htext(i) pm_hline(i)], 'Color', color)
                set(pm_htext(i), 'FontWeight', 'bold', 'Interpreter', 'none')
            end

        end
        axis tight
        if (clim(2) > clim(1))
            caxis(clim);
        end
        
        %Update plot of wells (if enabled)
        plotWellConnections();
    end

    function playBackTof(src, event)
        if playback
            return
        end
        N = round(get(speedsh, 'Value'));
        
        if (get(time_integral, 'Value'))
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
            
            if (get(time_integral, 'Value'))
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

        winj = repmat(D.inj,[1,np]);
        wpro = rldecode(D.prod,ni,2);
        
        wp = [];

        isInj = ismember(wk, D.inj);
        if isInj
            ik = find(D.inj == wk);
            sub = winj == wk;

            % well pair stuff
            if (numel(WP) > 0)
                wp = @(x) WP.inj(x);
            end
            otherNames = {W{state_idx}(D.prod).name};

            % set plots to match piecharts
            v = find(strcmpi(get(hdataset, 'String'), 'drainage region'));
            set(hdataset, 'Value', v(1))
            set(floodtoggle, 'Value', ik)
            set(draintoggle, 'Value', []);
        else
            ik = find(D.prod == wk);
            sub = wpro == wk;

            % well pair stuff
            if (numel(WP) > 0)
                wp = @(x) WP.prod(x);
            end
            otherNames = {W{state_idx}(D.inj).name};

            v = find(strcmpi(get(hdataset, 'String'), 'flooding region'));
            set(hdataset, 'Value', v(1))
            set(draintoggle, 'Value', ik)
            set(floodtoggle, 'Value', []);
        end

        if isempty(wellPlot) || ~ishandle(wellPlot)
            wellPlot = figure();
        else
            set(0, 'CurrentFigure', wellPlot); clf
        end
        set(wellPlot, 'name', ['Well ', W{state_idx}(wk).name]);

        plotArrival = ~isempty(state) && ~isInj;

        if plotArrival;
            subplot(2, 2, 1);
        else
            subplot(2, 2,[1 3])
        end

        if (numel(WP) > 0) 
            pie(max(WP.vols(sub), eps), ones(size(WP.vols(sub))))
            title('Pore volumes')
        end
        if plotArrival
            subplot(2,2,3);  cla;
            plotTOFArrival(state{state_idx}, W{state_idx}, pv, opt.fluid, find(D.prod == wk), D, opt.useMobilityArrival)
        end

        subplot(2,2,[2 4])
        if (numel(wp) > 0)
            tmp = wp(ik);
            if numel(tmp.z) > 1
                % Allocation factors by depth does not make sense for only one
                % perforation!
                [z, zind] = sort(tmp.z, 'descend');
                alloc = tmp.alloc(zind, :);
                % Always some allocation - avoid division by zero
                alloc(alloc == 0) = eps;
                walloc = bsxfun(@rdivide, cumsum(alloc,1), sum(alloc(:)));
                area(z, walloc, eps); axis tight;
                hold on
                plot(z, zeros(numel(z, 1)), '>', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'MarkerSize', 5);
                % Flip it around, xlabel is really ylabel
                view(90, -90);
                set(gca, 'XDir', 'reverse')
                legend(otherNames, 'Location', 'EastOutside');
                xlabel('Depth')
                title(['Allocation factors by depth for ', W{state_idx}(wk).name]);
            end
        end

        plotMain();
    end

    function plotPhi(src, event)
        if (~D.isvalid)
            warning('It appears that the time of flight is invalid');
            return
        end
        
        if isempty(phiPlot) || ~ishandle(phiPlot)
            phiPlot = figure();
        else
            set(0, 'CurrentFigure', phiPlot);
        end
        set(phiPlot, 'Name', ['Phi/F diagram (', window_name, ')']);
        [F,Phi] = computeFandPhi(pv,D.tof);
        plot(Phi,F,'.');
        title(sprintf('Lorenz coefficient: %f\n', computeLorenz(F,Phi)));
    end

    function plotWellAllocations(src, event)
        if (~D.isvalid)
            warning('It appears that the time of flight is invalid');
            return
        end
        
        if isempty(wah_fig) || ~ishandle(wah_fig)
            wah_fig = figure();
        else
            set(0, 'CurrentFigure', wah_fig);
        end
        set(wah_fig, 'Name', ['Well allocations (', window_name, ')']);
        plotWellAllocationComparison(D, WP, [], []);
    end

    function plotWellConnections(src, event)        
        if get(wctoggle, 'Value')            
            if ishandle(fig_main)
                set(0, 'CurrentFigure', fig_main);
            else
                fig_main = figure;
                axis tight off
                plotMain();
            end
            
            if any(ishandle(pwc_wch))
                delete(pwc_wch)
                pwc_wch = [];
            end
            
            if (D.isvalid)
                pwc_wch = plotWellPairConnections(computeGrid, WP, D, W{state_idx}, pv);
            else
                return
            end
        else
            delete(pwc_wch)
            pwc_wch = [];
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
            otherwise
                assert(isfield(datasets{state_idx}, datanames{dataind}), 'Trying to access non-existent field');
                cdata = readStructField(datasets{state_idx}, datanames{dataind});
        end
        
        if isempty(clim)
            m = min(cdata(:));
            M = max(cdata(:));

            clim = [m - eps(m), M + eps(m)];
        end
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
        
        D = computeTOFandTracer(rS, computeGrid, rock, 'wells', W_step);
        D.itracer(isnan(D.itracer)) = 0;
        D.ptracer(isnan(D.ptracer)) = 0;
        
        % Cap tof to maximum tof for unreachable areas for the time being
        tf = D.tof(:);
        if (all(isinf(tf)))
            D.isvalid = false;
            WP = [];
        else
            D.tof(isinf(D.tof)) = max(tf(isfinite(tf)));
            D.isvalid = true;
            WP = computeWellPairs(rS, computeGrid, rock, W_step, D);
        end
    end

    function computeValues(idx)
        [D, WP] = computeTOFAndTracerAndWellPairs(W{idx}, state{idx});

        if (~D.isvalid || isempty(WP))
            warning('Time of flight returned inf. Are there both active injectors and producers present?')
        end
    end

    function changeWells()
        W{state_idx} = editWells(G, W{state_idx}, rock);
        recomputeValues();
        createMainControl();
    end

    function recomputeValues()
        computeValues(state_idx)
        if ishandle(fig_main)
            close(fig_main)
        end
        plotMain();
    end
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

