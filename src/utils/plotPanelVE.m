function varargout = plotPanelVE(G, Gt, W, sol, t, vol, varargin)
%Make a panel plot used in the examples of the VE module
%
% SYNOPSIS
%       plotPanelVE(G, Gt, sol, t, vol)
%       plotPanelVE(G, Gt, sol, t, vol, 'pn1', pv1, ...)
%   h = plotPanelVE(...)
%
% PARAMETERS:
%   G       - Grid data structure for 3D grid
%   Gt      - Grid data structure for topsurface grid
%   W       - Well
%   sol     - Solution data structure
%   t       - Time
%   vol     - Vector of (free, trapped, total) volumes
%   opt     - Various options
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.
%
% DESCRIPTION:
%   The function will make a composit plot that consists of several parts
%    - a 3D plot of the plume
%    - a pie chart of trapped versus free volume
%    - a plane view of the plume from above
%    - two cross-sections in the x/y directions through the well
%
%   Supported options for the 3D plot of saturations
%      'maxH'   -- Double, maximal height value used for scaling the
%                  colorbar.  Default value: 10
%      'view'   -- Vector, view angle. Default value: [30,60]
%      'wireS'  -- Bool, if true, plot only positive saturation values
%                  on top of a wireframe grid. Default value: false
%      'Wadd'   -- Double, height added to the well. Default value: 100
%      'Saxis'  -- Vector, scaling of caxis. Default value: [0 1]
%
%   Supported options for the plot of CO2 height:
%      'wireH'  -- Bool, if true, plot only positive saturation values
%                  on top of a wireframe grid. Otherwise, make a color plot
%                  with height shown as contours. Default value: false
%      'ptol'   -- Double, smallest number that is considered as positive
%                  in the wireframe plot. Default value: 1e-3
%
%   Supported options for the slice plots:
%      'slice'  -- Vector, index along which to pick the two 2D slices
%                  Default value: slice=[10 10]
%      'plotPlume' -- Boolean indicating if the plume should be plotted as
%                     a sharp interface.
%
%   In addition, the function may plot the CO2 inventory as a function of
%   time if the user clicks the pie chart or if the option 'plotHist' is
%   set to true.
%
% SEE ALSO:
%   runIGEMS

persistent G_x G_y;      % Grids for 2D slices
persistent h h1 h2 h3;   % Figure and graphics handles
persistent x y;          % Node coordinates
persistent cells_sub_x;  % Cell numbers in x-slice
persistent cells_sub_y;  % Cell numbers in y-slice
persistent volHistory;   % Pie chart historic values
persistent thistory;     % Previous timesteps

prm = struct('maxH',      10, ...
             'view',      [30, 60],...
             'wireS',     false, ...
             'wireH',     false, ...
             'Wadd',      100, ...
             'Saxis',     [0 1], ...
             'ptol',      1e-3, ...
             'slice',     [10 10],...
             'plotPlume', false,  ...
             'plotHist',  false);

prm = merge_options(prm, varargin{:});

%% Prepare for subsequent plotting
if (t==0)
  clear ijk
  [ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap);
  ijk = [ijk{:}];

  cells_sub_x = find(ijk(:,2) == prm.slice(2));
  G_x = extractSubgrid(G, cells_sub_x);
  G_x = computeGeometry(G_x);
  G_x.cartDims = G.cartDims;

  cells_sub_y = find(ijk(:,1) == prm.slice(1));
  G_y = extractSubgrid(G, cells_sub_y);
  G_y = computeGeometry(G_y);
  G_y.cartDims = G.cartDims;

  scrsz  = get(0,'ScreenSize');
  fsize = min([scrsz(3:4); 1024 700]);
  h = figure('Position', [scrsz(3)-fsize(1), scrsz(4)-fsize(2)-75, fsize]);
  if nargout > 0; varargout{1} = h; end;
  
  set(h, 'color', 'white');

  if ~prm.wireH
     x = reshape(Gt.nodes.coords(:,1),Gt.cartDims+1);
     y = reshape(Gt.nodes.coords(:,2),Gt.cartDims+1);
  else
     h2 = [];
  end
  volHistory = vol;
  thistory   = 0;
  if prm.plotHist
     h3 = figure;
  else
     h3 = [];
  end
end


%% Height of the CO2 column
set(0,'CurrentFigure', h);
subplot(2,3,1:2);
if t==0
   if prm.wireS
      plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
   else
      h1 = plotCellData(G, sol.s, 'EdgeColor','w','EdgeAlpha', .15);
   end
   plotGrid(G_x, 'FaceColor', 'none', 'EdgeColor', 'c');
   plotGrid(G_y, 'FaceColor', 'none', 'EdgeColor', 'm');
   plotWell(G, W, 'height', prm.Wadd, 'color', 'r');
   axis tight off, view(prm.view);
else
   if ishandle(h1), delete(h1); end
   if prm.wireS
      h1 = plotCellData(G, sol.s, sol.s>prm.ptol, 'EdgeAlpha', 0.1);
   else
      h1 = plotCellData(G, sol.s, 'EdgeColor','w','EdgeAlpha', .15);
   end
end
title(['CO2 saturation at ', num2str(convertTo(t,year)), ' years']);
caxis(prm.Saxis);

%% Pie chart of trapped and free CO2
subplot(2,3,3), cla
if numel(vol)==3
   str1 = ['trapped ' num2str(round(vol(1)/vol(3)*100)) ' %'];
   str2 = ['free ' num2str(round(vol(2)/vol(3)*100)) ' %'];
   ph = pie([max(vol(1),eps) max(vol(2),eps)],{str1, str2});
   title(['Total volume: ', ...
      num2str(round(convertTo(vol(3),mega*meter^3))),' M m^3']);
elseif numel(vol)==6
   vplot = max([vol(1:5) vol(6)-sum(vol(1:5))], eps);
   ph = pie(vplot);
   % trick auto-shrink by first placing legend below and then moving it
   lnames = {'Struct. residual', 'Residual', 'Free residual', 'Struct. movable', ...
      'Free movable', 'Leaked'};
   hl = legend(lnames, 'Location','SouthOutside', 'orientation','horizontal');
   set(hl,'Location','SouthEastOutside','orientation','vertical');
   pl = get(hl, 'Position'); set(hl, 'Position', [.85 pl(2)-.035 pl(3:4)]);
   title(['Total injected volume: ', ...
      num2str(round(convertTo(vol(6),mega*meter^3))),' M m^3']);

   if t ~= 0
        volHistory = [volHistory; vplot];
        thistory = [thistory; t];
        set(ph, 'ButtonDownFcn', @(src, event) makePlotHistory);
        plotHistory()
        set(0,'CurrentFigure', h);
    end
else
   disp('Incompatible volume vector. No pie chart created');
end



%% Height of the CO2 plume plotted on the 2D top-surface grid
subplot(2,3,4)
if prm.wireH
   if t==0
      plotGrid(Gt, 'FaceColor', 'none', 'EdgeAlpha', 0.1);
      axis tight off, set(gca,'zdir','reverse'); view([-90, 90])
      colorbar, caxis([0 prm.maxH]);
   else
      if ishandle(h2), delete(h2); end
      h2 = plotCellData(Gt, sol.h, sol.h(:,1)>prm.ptol,'EdgeAlpha', 0.1);
   end
else
   cla
   hplot = reshape(sol.h(:,1),Gt.cartDims);
   hplot = hplot([1:end end],:);
   hplot = hplot(:,[1:end end]);
   pcolor(x,-y,hplot); shading flat; caxis([0 prm.maxH])
   hold on;
   colorbar
   contour(x,-y,reshape(Gt.nodes.z,Gt.cartDims+1),5,'w');
   hold off, axis tight off
end
title('Height of CO2-column');

%% Saturation along the two slices
subplot(4,3,8:9); cla % x-slice
if prm.plotPlume
    plotPlume(G, Gt, sol.h, cells_sub_x,'EdgeColor','w','EdgeAlpha',.1);
    plotGrid(G, cells_sub_x,'EdgeColor', 'k','EdgeAlpha',.1, 'FaceColor', 'none');
    caxis([0 prm.maxH]);
else
    plotCellData(G_x, sol.s(cells_sub_x),'EdgeColor','w','EdgeAlpha',.1);
    caxis(prm.Saxis);
end
axis tight off, view([0 0] + prm.plotPlume*4), title('CO2-saturation, x-slice (cyan)');


subplot(4,3,11:12); cla % yslice
if prm.plotPlume
    plotPlume(G, Gt, sol.h, cells_sub_y,'EdgeColor','w','EdgeAlpha',.1);
    plotGrid(G, cells_sub_y,'EdgeColor','k','EdgeAlpha',.1, 'FaceColor', 'none');
    caxis([0 prm.maxH]);
else
    plotCellData(G_y, sol.s(cells_sub_y),'EdgeColor','w','EdgeAlpha',.1);
    caxis(prm.Saxis);
end
axis tight off, view([90 0] + prm.plotPlume*4), title('CO2-saturation, y-slice (magenta)');
drawnow

function makePlotHistory()
    if isempty(h3)
        h3 = NaN;
    end
    
    if ~ishandle(h3)
        h3 = figure;
    end
    plotHistory();
end

function plotHistory()
    if ishandle(h3)
        set(0, 'CurrentFigure', h3);
        area(convertTo(thistory, year), convertTo(volHistory, mega*meter^3));
        xlabel('Years since simulation start');
        ylabel('Volume (M m^3)')
        axis tight
        legend(lnames, 'location', 'eastoutside', 'orientation', 'vertical');
    end
end

end
