function plotWellAllocationComparison(D1, WP1, D2, WP2, varargin)
%Plot a panel comparing well-allocation from models with different resolution
%
% SYNOPSIS
%   plotWellAllocationComparision(D1, WP1, D2, WP2)
%
% PARAMETERS:
%   D1, D2   - data structure with basic data for flow diagnostics computed
%              by a call to 'computeTOFandTracer' for model 1 and model 2
%
%   WP1, WP2 - data structure containing information about well pairs,
%              computed by a call to 'computeWellPairs'
%
% DESCRIPTION:
%   The routine makes a bar plot for each well segment that is represented
%   in the input data D/WP. (There is no check that the well segments in
%   D1 and D2 are the same). For injection wells, the bars represent
%   the cumulative outfluxes, from the bottom to the top of the segment,
%   that have been attributed to the different producers. The bars are
%   shown in color for model 1, with a unique color representing each of
%   the segments of the producers,  and in solid black lines for model 2.
%   For the production wells, the bars represent influxes that can be
%   attributed to different injector segments. If the two models predict
%   the same flux allocation, the color bars and the solid lines should be
%   matching.
%
%   If D2 and WP2 are empty, a graph is produced that shows the well
%   allocation for the wells in WP1.
%
% SEE ALSO:
%   `computeTOFandTracer`, `computeWellParis`, `expandWellCompletions`

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
opt = struct('useZ', true,...
             'bwidth', .98, ...
             'plotrow', false);
opt = merge_options(opt, varargin{:});

% -------------------------------------------------------------------------
% check what we should plot
if (isempty(D2) || isempty(WP2))
    plotFlag = false;
else
    plotFlag = true;
end

% -------------------------------------------------------------------------
% Parameterization along well path
% Create fake z values based on the cell numbering in the well. Useful
% whenever the well isn't vertical.
WP1 = addDummyData(WP1);
if ~opt.useZ
    for i = 1:numel(WP1.inj)
        n = numel(WP1.inj(i).z);
        WP1.inj(i).z = (1:n)./n;
    end
    for i = 1:numel(WP1.prod)
        n = numel(WP1.prod(i).z);
        WP1.prod(i).z = (1:n)./n;
    end
end
if plotFlag
    WP2 = addDummyData(WP2);
    if ~opt.useZ
        for i = 1:numel(WP1.inj)
            n = numel(WP2.inj(i).z);
            WP2.inj(i).z = (1:n)./n;
        end
        for i = 1:numel(WP1.prod)
            n = numel(WP2.prod(i).z);
            WP2.prod(i).z = (1:n)./n;
        end
    end
end

% -------------------------------------------------------------------------
% Extract and format the well-allocation factors for model 1
nit = numel(D1.inj);
npt = numel(D1.prod);
for i=1:numel(WP1.inj)
   nc   = numel(WP1.inj(i).z);
   atmp = zeros(nc, nit+npt);
   atmp(:,D1.prod) = WP1.inj(i).alloc;
   WP1.inj(i).alloc = atmp;
end
for i=1:numel(WP1.prod)
   nc = numel(WP1.prod(i).z);
   atmp = zeros(nc, nit+npt);
   atmp(:,D1.inj) = WP1.prod(i).alloc;
   WP1.prod(i).alloc = atmp;
end
wp1(D1.inj)  = WP1.inj(:);
wp1(D1.prod) = WP1.prod(:);

% -------------------------------------------------------------------------
% Extract and format the well-allocation factors for model 2
if plotFlag
    nit = numel(D2.inj);
    npt = numel(D2.prod);
    for i=1:numel(WP2.inj)
        nc   = numel(WP2.inj(i).z);
        atmp = zeros(nc,nit+npt);
        atmp(:,D2.prod) = WP2.inj(i).alloc;
        WP2.inj(i).alloc  = atmp;
    end
    for i=1:numel(WP2.prod)
        nc   = numel(WP2.prod(i).z);
        atmp = zeros(nc,nit+npt);
        atmp(:,D2.inj) = WP2.prod(i).alloc;
        WP2.prod(i).alloc  = atmp;
    end
    wp2(D2.inj)  = WP2.inj(:);
    wp2(D2.prod) = WP2.prod(:);
end

% -------------------------------------------------------------------------
% Calculate number of plots in 2D grid
num_plots = numel(wp1);
ny = ceil(sqrt(num_plots));
nx = ny-1;
if (ny*nx < num_plots)
    nx = ny;
end
if opt.plotrow, [nx,ny] = deal(num_plots,1); end

current_fig = gcf();
wb = waitbar(0, 'Well 1');
set(0, 'CurrentFigure', current_fig);

% -------------------------------------------------------------------------
% Plot the well-allocation factors for fine/coarse scale models
for i=1:num_plots
   subplot(ny,nx,i);
   
   % --------------------------------
   % Plot data for model 1
   [~,ix]   = sort(wp1(i).z);
   wp1(i).z = wp1(i).z(ix);
   wp1(i).alloc = wp1(i).alloc(ix,:);
   nseg = numel(wp1(i).z);
   if nseg == 1  % Need to trick Matlab
      z = [wp1(i).z; wp1(i).z+1];
      a = [cumsum(wp1(i).alloc,1); zeros(1,numel(wp1(i).alloc))];
      bwidth = 2*numel(wp2(i).z)+3;
      h = barh(z, a,'stacked','BarWidth',bwidth, 'EdgeColor','none');
   else
      z = flipud(wp1(i).z);
      a = cumsum(flipud(wp1(i).alloc),1);
      if nseg<21
          aargs = {'YDir', 'YLim'};
          h = barh(z, a,'stacked','BarWidth',opt.bwidth, 'EdgeColor','none');
      else
          h = area(z, a, eps,'EdgeColor','none'); axis tight;
          view(90,-90);
          aargs = {'XDir', 'XLim'};
      end
   end
   zm = min(wp1(i).z);
   zM = max(wp1(i).z);
   
   % --------------------------------
   % Plot data for model 2
   if plotFlag
       hold on
       nseg = numel(wp2(i).z);
       if nseg == 1  % Need to trick Matlab
           z = [wp2(i).z; wp2(i).z+1];
           a = [cumsum(wp2(i).alloc,1); zeros(1,numel(wp2(i).alloc))];
           barh(z,a,'stacked','FaceColor','none','BarWidth',opt.bwidth);
       else
           z = flipud(wp2(i).z);
           a = cumsum(flipud(wp2(i).alloc),1);
           if nseg<21
               aargs = {'YDir', 'YLim'};
               barh(z,a,'stacked','FaceColor','none','BarWidth',opt.bwidth);
           else
               h = area(z, a, eps,'FaceColor','none'); axis tight;
               view(90,-90);
               aargs = {'XDir', 'XLim'};
           end
       end
       hold off
       zm = min(zm,min(wp2(i).z));
       zM = max(zM,max(wp2(i).z));
   end
   
   % ---------------------------------
   % Set axis and add legend
   set(gca,aargs{1},'reverse', aargs{2}, [zm zM + sqrt(eps)] + [-.1 .1]*(zM-zm));

   if max(wp1(i).alloc(:))>0
      hl=legend(h(i),wp1(i).name,'Location','SouthEast');
      set(hl,'FontSize',14,'FontWeight','demi'); legend boxoff
   else
      hl=legend(h(i),wp1(i).name,'Location','SouthWest');
      set(hl,'FontSize',14,'FontWeight','demi'); legend boxoff
   end
   waitbar(i/num_plots,wb,['Well ', num2str(i+1)])
end

% -------------------------------------------------------------------------
% Set colormap
cmap = jet(nit+npt);
c = 0.9*cmap + .1*ones(size(cmap)); colormap(c);
drawnow;
delete(wb);

end


% =========================================================================

function WP = addDummyData(WP)
    ni = numel(WP.inj);
    np = numel(WP.prod);
    for i = 1:ni
        if isempty(WP.inj(i).alloc)
            WP.inj(i).alloc = zeros(1, np);
            WP.inj(i).z = 1;
        end
    end
    for i = 1:np
        if isempty(WP.prod(i).alloc)
            WP.prod(i).alloc = zeros(1, ni);
            WP.prod(i).z = 1;
        end
    end
end
