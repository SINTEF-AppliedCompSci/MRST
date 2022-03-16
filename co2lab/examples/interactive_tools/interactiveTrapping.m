function varargout = interactiveTrapping(inp, varargin)
%Create an interactive figure showing trapping structure for a top surface grid
%
% For a detailed list of functionality and controls, please see the
% showTrapsInteractively example.
%
% SYNOPSIS:
%       interactiveTrapping('Johansenfm')
%       interactiveTrapping(G_top)
%
% PARAMETERS:
%   inp     - Either a valid top surface grid as defined by
%             topSurfaceGrid(G) or a string which is valid input for
%             getAtlasGrid.
%
%   'pn'/pv - List of optional property names/property values:
%                   
%               - coarsening: If inp is the name of a Atlas grid, this
%                 argument will be passed onto getAtlasGrid to coarsen the
%                 surface. Default: 1 for the full grid.
%
%               - light: toggle the lighting of the surface grid.
%                 Default: false
%
%               - spillregions: toggle the display of spill regions on the
%                 surface grid
%                 Default: false
%
%               - traps: toggle the display of traps on the surface grid
%                 Default: true
%
%               - colorpath: toggle red/gray color scheme to mark the traps
%                 encountered along a spill path; same color is applied on
%                 the histogram of trapping volumes.
%                 Default: true
%
%               - method: Choose cell or node based algorithm. Valid
%                 inputs: 'node', 'cell'. Default: 'cell'.
%
%               - injpt: Choose this cell number as injection point at
%                 startup. If zero, no injection point is selected.
%                 Default: zero
%
%               - contour: toggle drawing of contour lines of height data
%                 if these are available as part of the data set
%
% RETURNS:
%   h  - Handle to resulting figure object.
%
%
% SEE ALSO:
%   `trapAnalysis`, `showTrapsInteractively`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   opt = struct('coarsening',    1,     ...
                'light',        false,  ...
                'spillregions', false,  ...
                'traps',        true,   ...
                'method',       'cell', ...
                'colorpath',    true,   ...
                'injpt',        0,      ...
                'contour',      false  ...
                );
   
   global selectedCell;
   
   if nargout > 0
       varargout{1} = NaN;
   end
   
   if strcmpi(opt.method, 'cell')
       require('coarsegrid');
   end
   
   opt = merge_options(opt, varargin{:});
   selectedCell = opt.injpt;

    if isfield(inp, 'cells')
        % Assume input is a top surface grid
        Gt = inp;
        assert(isfield(Gt, 'columns'), 'Please provide a top surface grid');
        top = [];

    elseif ischar(inp)
        % Assume input is the name of a CO2 Storage Atlas grid

        [gr, data, petroinfo] = getAtlasGrid(inp, 'coarsening', opt.coarsening, 'nz', 1);

        % Grab topgrid
        top = data(cellfun(@(x) strcmpi(x.variant, 'top'), data));
        if numel(top)>0
            top = top{1};
        else
            top = [];
        end
        % Process the grid
        G = processGRDECL(gr{1});

        % The coarsening may disconnect small components, take only the largest and
        % first grid produced by processGRDECL and add geometry data to it.
        try
            mlist = mrstModule();
            mrstModule add libgeometry
            G = mcomputeGeometry(G(1));
            mrstModule('reset', mlist{:})
        catch %#ok
            G = computeGeometry(G(1));
        end

        Gt = topSurfaceGrid(G);
        % Pass it silently along with the grid
        Gt.petroinfo = petroinfo{1};
    end

    res = trapAnalysis(Gt, strcmpi(opt.method, true));

    if isempty(res.trap_adj) || all(res.trap_adj(:)==0)
        disp('No trap adjacency graph produced. Does the grid have any traps?')
        return
    end
    
    res.Gflat = flattenTraps(Gt, res);
    % Store a list of traps to avoid slicing sparse matrices
    % Commented out because it wasn't needed. May be used for very large
    % grids.
%     [i, j] = find(res.trap_adj);
%     T = sortrows([i,j]);
%     i = T(:,1);
%     j = T(:,2);
%     p = cumsum([1; accumarray(i, 1)]);
%     res.p = p;
% 
%     res.i = i;
%     res.j = j;


    defaultpos = get(0, 'DefaultFigurePosition');
    % Try to be sort of not stupid when creating a large figure by going via
    % the default setting.

    h = figure('Position', [defaultpos(1:2) - [0 100], 1000 500]);
    bf = boundaryFaces(Gt);

    % Create a toolbar for toggling different displays
    ht = uitoolbar(h);
    
    uitoggletool(ht, 'TooltipString', 'Plot traps',...
                     'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top),...
                     'cdata', geticon('traps'),...
                     'Tag', 'toggleTraps', ...
                     'State', boolToToggle(opt.traps));

    uitoggletool(ht, 'TooltipString', 'Plot spill regions',...
                     'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top),...
                     'cdata', geticon('spill'),...
                     'Tag', 'toggleSpillRegions',...
                     'State', boolToToggle(opt.spillregions));
    
                 
    uitoggletool(ht, 'TooltipString', 'Lighting',...
                     'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top), ...
                     'cdata', geticon('light'),...
                     'Tag', 'toggleLight', 'State', boolToToggle(opt.light));

    uitoggletool(ht, 'TooltipString', 'Toggle contourlines',...
                     'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top),...
                     'cdata', geticon('contour'),...
                     'Tag', 'toggleContour', 'State', boolToToggle(opt.contour));
                  
    uitoggletool(ht, 'TooltipString', 'Red/gray color along spill paths', ...
                      'ClickedCallback', @(src, event) plotMain(Gt, res, bf, top), ...
                      'cdata', geticon('colorpath'), ...
                      'Tag', 'toggleColorpath', 'State', boolToToggle(opt.colorpath));

    uipushtool(ht, 'TooltipString', 'Reset view',...
                   'cdata', geticon('reset'),...
                   'ClickedCallback', @(src, event) view(0, 90));

    uipushtool(ht, 'TooltipString', 'Setup VE simulation',...
                   'cdata', geticon('simulate'),...
                   'ClickedCallback', @(src, event) runSimulation(Gt, res, src, event));
                
    plotMain(Gt, res, bf, top);
    
    if nargout > 0
       % Return a handle to the figure if it was asked for.
       varargout{1} = h;
    end
end

function clickHandler(Gt, res, bf, data, src, event)                       %#ok
    global selectedCell;

    % Get data
    pts = get(gca, 'CurrentPoint');
    [c f] = nearestCellLine(Gt, 1:Gt.faces.num, pts);                      %#ok
   
    if strcmpi(get(gcf, 'SelectionType'), 'alt')
        figure; clf;
        
        trap = res.trap_regions(c(1));
        currTrap = res.traps == trap;
        
        volume = volumesOfTraps(Gt, res, trap);
        title(sprintf('Trap %d, containing %1.3g m^3 volume', trap, volume));
        gc = Gt.cells.centroids(currTrap, :);
        
        m = min(gc);
        M = max(gc);

        % Find the general area around the trap
        currRegion = findRegion(Gt, m, M, .3);

        % Actual plot
        
        % Plot bottom of trap
        plotGrid(res.Gflat, currTrap, 'facec', 'red')
        % Plot ceiling of trap
        plotGrid(Gt, currTrap & currRegion, 'facec', 'blue', 'facea', .8, 'edgec', 'none')
        plotCellData(Gt, Gt.cells.z, ~currTrap & currRegion, 'facea', .3, 'edgec', 'none');
        
        axis off tight
        view(-20, 45)
        legend({'Trap floor', 'Trap ceiling'})
        cbar = colorbar('South');
        set(get(cbar, 'YLabel'), 'String', 'Depth')
    else
       selectedCell  = c(1);
       
       reptext = findobj('Tag', 'gridReport');
       if ~isempty(reptext)
          set(reptext(1), 'string', getReport(Gt, selectedCell));
       end
       plotMain(Gt, res, bf, data);
    end
end

function plotMain(Gt, res, bf, atlasdata)
% Main plotting of grid

    global selectedCell;

    
    isOn = @(tag) strcmpi(get(findobj(gcf, 'tag', tag), 'State'), 'on');
    
    showTraps = isOn('toggleTraps');
    showSpill = isOn('toggleSpillRegions');
    showLight = isOn('toggleLight');  
    showContour = isOn('toggleContour');  
    
    fn = @(src, event) clickHandler(Gt, res, bf, atlasdata, src, event);
    subplot('position', [.025 .025 .65 .95]);

    [az, el] = view();
    cla; ax = gca;
    hold on

    zdata = Gt.cells.z;
    zrange = (zdata - min(zdata))./(max(zdata) - min(zdata));
    
    if isfield(res, 'trap_regions') && showSpill
        % Plot accumulation regions
        if max(res.traps)<128
           cmap = lines(max(res.traps) + 1);
        else
           cmap = jet(max(res.traps) + 1);
        end
        colormap(cmap);
        map = greyMask(cmap);
        map(1,:) = get(gcf, 'Color');
        tmp = res.trap_regions;
        tmp(tmp>max(res.traps)) = max(res.traps);
        
        h = plotCellData(Gt, ones(Gt.cells.num, 1), 'ButtonDownFcn', fn, 'EdgeColor', 'none');
        set(h, 'FaceVertexCData', map(tmp + 1, :))
    else
        % Plot heightdata
        h = plotCellData(Gt, zdata, 'ButtonDownFcn', fn, 'EdgeColor','none');
        map = jet(64);
        v = interp1(linspace(0,1,64), map, zrange );
        set(h, 'FaceVertexCData', v);
    end
    
    if showTraps
        plotCellData(Gt, res.traps, res.traps ~= 0, ...
           'ButtonDownFcn', fn, 'edgecolor', [.3 .3 .3], 'edgea', .05)
        plotPartitionOutlineTopsurface(Gt, res.traps, 'LineWidth', 1.5)
        caxis([0 max(res.traps)]);
    end
    
    % Determine whether we have forward or backward mode
    if strcmpi(get(gcf, 'SelectionType'), 'extend')
        outline = 'green';
        A = res.trap_adj';
    else
        outline = 'blue';
        A = res.trap_adj;
    end
    
    % Do we have any cell selected?
    if selectedCell>0
       trap = res.trap_regions(selectedCell);
       wpos = Gt.cells.centroids(selectedCell,:);
    else
       trap = 0;
       wpos = [];
    end
    zm = min(Gt.cells.z);
    zM = max(Gt.cells.z);
    
    if trap==0,
       disp(['Current position is not downstream from any trap, '...
              'or spillpoint data was not produced.'])
       if ~isempty(wpos)
          plotGrid(Gt, selectedCell, 'facec', 'yellow', 'ButtonDownFcn', fn)
          plot3(wpos([1 1],1), wpos([1 1],2), [zm zM],'r','LineWidth',2);
       end
    else
    
       plotGrid(Gt, res.traps == trap, 'edgec', 'none', 'facec', 'red', 'ButtonDownFcn', fn);
       subt = getMigrationTree(Gt, A, trap, 0, fn);
    
       t = [trap subt];
       colorizePath = strcmpi(get(findobj(gcf, 'tag', 'toggleColorpath'), 'State'), 'on');
       for i = 1:numel(t)
          ci = res.traps == t(i);
        
          plotGrid(Gt, ci, 'facec', colorizeDepth(i, numel(t), colorizePath),...
             'edgec', 'none', 'ButtonDownFcn', fn);
       end
       partition = res.traps;
       partition(~ismember(partition, t)) = max(t)+1;
    
       Gtmp = Gt;
       Gtmp.nodes.z = Gt.nodes.z - (zM- zm)/1000;
       plotPartitionOutlineTopsurface(Gtmp, partition, 'LineWidth', 1.5, 'EdgeColor', outline)
       
       if isfield(res, 'cell_lines')
          [x, y, z] = getPlotRivers(Gt, [res.cell_lines{[unique(subt) trap]}], res);
       end
       plot3(x, y, z, 'k', 'LineWidth', 2);
    
       % Highlight current selection
       plotGrid(Gt, selectedCell, 'facec', 'yellow', 'ButtonDownFcn', fn)
       plot3(wpos([1 1],1), wpos([1 1],2), [zm zM],'r','LineWidth',2);
    end
    hold off
    %colormap jet
    axis tight off
    view(az,el);

    if ~isempty(atlasdata) && showContour
        contourAtlas(atlasdata, 10, 1)
    end
    
    if showLight
        % Account for lighting bug in Matlab due to reversed Z-axis
        if exist('verLessThan', 'file') && ~verLessThan('matlab', '9.0')
            sgn = -1;
        else
            sgn = 1;
        end
        light('position', [max(Gt.cells.centroids) sgn*4*max(Gt.cells.z)],'Style','local');
        lighting phong
        material dull
    end
    if trap==0, return, end;
    
    %%% Plot piechart
    subplot('position', [.75 .525 .225 .425]); cla reset; hold off;
    subtraps = unique(subt);
    subtraps(subtraps == trap) = [];
    
    % Cap above zero since pie refuses to plot zero values
 
    vprimary   = volumesOfTraps(Gt, res, trap);
    vsecondary = volumesOfTraps(Gt, res, subtraps);
    if isempty(subtraps) || sum(vsecondary) < 1;
        % To get nice plot, set zero values to 1 m^3 volume. This is small
        % at the scales we are operating on, and will give a better plot.
        vsecondary = 1;
    end
    
    trapped = vprimary + sum(vsecondary);
    
    total = sum(volumesOfTraps(Gt, res, 1:max(res.traps))) - trapped;
    total = max(total,eps);  % ensure a small positive number
    
    hp = pie([vprimary, sum(vsecondary), total], [1,1,1]);
    pcm = [0 0 .7; .1 .8 .1; .6 0 0];
    for k=1:3, set(hp(2*k-1), 'FaceColor', pcm(k,:)); end
    hl = legend('Primary', 'Migration', 'Not Filled', ...
       'Location','SouthOutside','Orientation','horizontal');
    set(hl, 'FontSize', 8, 'Color', min(get(gcf,'Color')+.1, 1));

    %%%  Plot barchart
    subplot('position', [.775 .025 .2 .45]); cla
    % While setting the axis to logarithmic would be the best idea, this
    % unfortunately breaks opengl support for ALL subplots. We instead
    % manually create a logarithmic axes
    data = [vprimary vsecondary(end:-1:1)];
    ldata = log10(data);
    ldata(ldata < 0) = 0;
    % Transform and plot the data to the unit axis
    % trick: By using hist we get tighter spacing *and* the bars are for
    % some reason drawn as patches. Moreover, if N>150, we need to set
    % EdgeColor to 'none' to avoid causing problems when 'bar' sets
    % EdgeColor equal FaceColor.
    hbar = bar(ldata./max(ldata), 'hist');
    N = numel(subtraps) + 1;
    if N>150
       set(hbar, 'FaceVertexCData', colorizeDepth(1:N, N, colorizePath), 'EdgeColor', 'none');
    else
       set(hbar, 'FaceVertexCData', colorizeDepth(1:N, N, colorizePath));
    end
    yt = get(gca, 'YTick');
    ytl = arrayfun(@(y) sprintf('%2.3g', 10.^(y*max(ldata))), yt, 'UniformOutput', false);
    set(gca, 'YTickLabel', ytl)
    set(gca, 'YTickMode', 'manual')
    axis tight
    xlabel('Trap graph distance')
    ylabel('Trapped volume (m^3)')
    
    % Set axis to the main plot to ensure that immediate changes to view
    % etc. go to the actual 3D axis.
    axes(ax); %#ok<MAXES>
end

function subtraps = getMigrationTree(G, A, trap, depth, fn)
% Recursively traverse and find the full migration tree
    subtraps = find(A(trap, :));
    tmp = [];
    for i = 1:numel(subtraps);
        trp = getMigrationTree(G, A, subtraps(i), depth + 1, fn);
        tmp = [tmp trp]; %#ok
    end
    subtraps = [reshape(subtraps, 1, []) reshape(tmp, 1, [])];
end

function [c, f] = nearestCellLine(G, bf, pts)
% Find a point on the surface based on a click.
    N = numel(bf);
    x1 = repmat(pts(1,:), N, 1);
    x2 = repmat(pts(2,:), N, 1);
    x0 = G.faces.centroids(bf,:);
    if size(x0, 2) == 2
        x0 = [x0 G.faces.z(bf)];
    end
    % Find distance projected on line with a fudge factor for strange grids
    linedist = sum(cross(x0-x1, x0-x2 ,2).^2, 2);
%     ptdist = sum((x0-x1).^2, 2);
%     dist = sqrt(ptdist - 10*linedist);
    
    dist = linedist;
    
    [val ind] = min(dist); %#ok backwards compatability
    f = bf(ind);
    c = G.faces.neighbors(f, :);
    c = c(c~=0);
end


function c = findRegion(G, m, M, radius)
% Find a square region in some radius given minimum and maximum values for
% x and y.
    nc = G.cells.num;
    r = radius*(M-m);
    
    c = true(nc, 1);
    for i = 1:2
        c = c & G.cells.centroids(:,i) < M(i) + r(i);
        c = c & G.cells.centroids(:,i) > m(i) - r(i);
    end
end

function [x, y, z] = getPlotRivers(Gt, cell_lines, res)
% Plot rivers
    cc=[Gt.cells.centroids, Gt.cells.z];
    [x,y,z] = deal([]);
    
    for i=1:numel(cell_lines);
       cl = unique(cell_lines{i});
       cl = cl(res.traps(cl) == 0);
       x = [x; cc(double(cl),1); NaN]; %#ok
       y = [y; cc(double(cl),2); NaN]; %#ok
       z = [z; cc(double(cl),3); NaN]; %#ok
    end

end

function c = colorizeDepth(depth, N, colorize)
% Colorize from bright red to nearly black
    depth = reshape(depth, [], 1);
    c= max(1-eps - (depth-1)/(N+1), 0);
    if colorize
       c = [.25+.75*c .15*ones(length(c),2)];
    else
       c = .25 + .5*repmat(c, 1, 3);
    end
end

function cdata = geticon(name)
    n = mfilename('fullpath');
    icon = fullfile(n(1:end-length(mfilename)), 'icons', [name '.gif']);
    [cdata, map] = imread(icon);
    map(find(map(:,1)+map(:,2)+map(:,3)==3)) = NaN;     %#ok trivial lookup
    cdata = ind2rgb(cdata,map);
end

function c = greyMask(c)
    c = (c + repmat(get(gcf, 'Color'), size(c, 1), 1))./2;
end

function plotPartitionOutlineTopsurface(Gt, p, varargin)
    cg = generateCoarseGrid(Gt, p);
    flag = true(cg.faces.num,1);
    flag(boundaryFaces(cg))=false;
    if sum(flag)==0
       plotFaces(cg, varargin{:});
    else
       plotFaces(cg, flag, varargin{:});
    end
end

function runSimulation(Gt, res, src, event) %#ok<INUSD>
    global selectedCell veSimAborted
    persistent fi
    
    if selectedCell==0
       disp('Please select an injection point first')
       return;
    end
    pos = get(gcf, 'OuterPosition');
    if isempty(fi) || ~ishandle(fi) 
        fi = figure('Position',[pos(1:2)-[300 0],  [510 400]], 'Toolbar','none', 'MenuBar', 'none');
    else
        figure(fi);
        clf;
    end
    
    set(fi, 'Name', 'Set up simulation');
    
    bottom = 50;

    ph = uipanel(fi, 'Position', [0 0 1 1], 'Title', 'Set up VE simulation');
    
    % Petrophysical stuff
    poro = .3;
    perm = 300;
    if isfield(Gt, 'petroinfo')
        if ~isnan(Gt.petroinfo.avgperm)
            perm = convertTo(Gt.petroinfo.avgperm, milli*darcy);
        end
        if ~isnan(Gt.petroinfo.avgporo)
            poro = Gt.petroinfo.avgporo;
        end
    end
    sporo = linkedSlider(ph, [10 bottom], 0, 1, poro, 'Porosity');
    
    sperm = linkedSlider(ph, [10 bottom + 1*25], min(.1, perm), max(1000, perm), perm, 'Permeability (mD)');
    
    % Timesteps etc
    sinj  = linkedSlider(ph, [10 bottom + 3*25 + 10], 0, 200, 50, 'Injection time (years)');
    stni  = linkedSlider(ph, [10 bottom + 2*25 + 10], 1, 400, 25, 'Time steps');
    smig  = linkedSlider(ph, [10 bottom + 5*25 + 10], 0, 10000, 2450, 'Migration time (years)');
    stnm  = linkedSlider(ph, [10 bottom + 4*25 + 10], 1, 1000, 245, 'Time steps');
    spres = linkedSlider(ph, [10 bottom + 6*25 + 20], 50, 1000, 200, 'Pressure (bar)');

    % Rough CO2 volume estimate for the reservoir to ensure that the
    % scrollbar shows something (numerically) sensibile
    estco2vol = 500*.3*sum(Gt.cells.volumes.*Gt.cells.H)/1e9;
    maxrate = min(20,estco2vol/50);
    
    sco2 = linkedSlider(ph, [10 bottom + 7*25 + 20], 0, maxrate, maxrate*0.05, 'CO2 rate (Mt/year)');
    
    uicontrol(ph, 'Style', 'pushbutton', 'Position', [100, 5, 100, 40], 'string', 'Simulate!', 'callback', @(src, event) uiresume(fi))
    uicontrol(ph, 'Style', 'pushbutton', 'Position', [240, 5, 100, 40], 'string', 'Abort', 'callback', @(src, event) abortSim()) 

    advancedph = uicontrol(ph, 'Style', 'checkbox', 'Position', [350, 5, 120, 40], 'string', 'Advanced plot', 'Value', true);
    
    uicontrol(ph, 'Style', 'text', 'Tag', 'gridReport', 'Position', [10 bottom + 7*25 + 50, 455, 100], 'string', getReport(Gt, selectedCell))
    
    while ishandle(fi)
        uiwait(fi);

        if ~ishandle(fi)
            % User closed window
            return
        end
        
        petrodata.avgperm = get(sperm, 'Value')*milli*darcy;
        petrodata.avgporo = get(sporo, 'Value');
        
        veSimAborted = false;
        mrstModule add incomp
        migrateInjection(Gt, res, petrodata, selectedCell,...
                         'amount',      get(sco2, 'Value'),...
                         'T_injection', get(sinj, 'Value')*year,...
                         'T_migration', get(smig, 'Value')*year,...
                         'topPressure', get(spres, 'Value')*barsa,...
                         'Ni',          round(get(stni, 'Value')),...
                         'Nm',          round(get(stnm, 'Value')), ...-
                         'plotPanel',   get(advancedph, 'Value') == 1);
    end
end

function [s, e] = linkedSlider(fi, pos, md, Md, val, title)
    uicontrol(fi, 'Style', 'text', 'Position', [pos(1:2) 150, 25], 'string', title)
    size1 = [pos(1)+150, pos(2), 250, 25];
    sizee1 = [size1(1:2) + [size1(3) + 10, 0], [50 25]];
    cap = @(x) max(md, min(x, Md));
    
    e = uicontrol(fi, 'Style', 'edit', 'Value', val, 'Position', sizee1, 'String', sprintf('%.2f', val));
    fun = @(src, event) set(e, 'String', sprintf('%.2f', get(src, 'Value')));
    s = uicontrol(fi, 'Style', 'slider', 'Position', size1, 'Min', md, 'Max', Md, 'Value', val, 'Callback', fun);
    fun2 = @(src, event) set(s, 'Value', cap(sscanf(get(src, 'String'), '%f')));
    set(e, 'Callback', fun2);
end

function abortSim()
    global veSimAborted
    veSimAborted = true;
end

function v = boolToToggle(boolvar)
    if boolvar
        v = 'on';
    else
        v = 'off';
    end
end

function s = getReport(Gt, c)
    [i, j] = ind2sub(Gt.cartDims, c);
    s = '';
    s = [s sprintf('Reservoir with geometric volume of %2.4g m^3 and %d cells\n', sum(Gt.cells.volumes.*Gt.cells.H), Gt.cells.num)];
    s = [s sprintf('Injection site located at cell %d at \n\tx: %1.2f   y: \t%1.2f      i: %d   j: %d', c, Gt.cells.centroids(c, 1:2), i, j)];
end
