function plotWellAllocationPanel(D, WP, varargin)
%Plot a panel comparing well-allocation from models with different resolution
%
% SYNOPSIS
%   plotWellAllocationPanel(D, WP)
%
% PARAMETERS:
%   D  - data structure with basic data for flow diagnostics computed
%        by a call to 'computeTOFandTracer'
%
%   WP - data structure containing information about well pairs,
%        computed by a call to 'computeWellPairs'
%
% DESCRIPTION:
%   The routine makes a bar plot for each well segment that is represented
%   in the input data D/WP.  For injection wells, the bars represent
%   the cumulative outfluxes, from the bottom to the top of the segment,
%   that have been attributed to the different producers. The bars are
%   shown in color with a unique color representing each of
%   the segments of the producers.
%   For the production wells, the bars represent influxes that can be
%   attributed to different injector segments.
%
% SEE ALSO:
%   `computeTOFandTracer`, `computeWellPairs`, `expandWellCompletions`

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
             'plotOnly', NaN, ...
             'bwidth', .98);
opt = merge_options(opt, varargin{:});

% -------------------------------------------------------------------------
% check what we should plot
if (isnan(opt.plotOnly))
    num_wells = numel(WP.inj) + numel(WP.prod);
    opt.plotOnly = 1:num_wells;
else
    if (numel(opt.plotOnly) == 0)
        warndlg('plotWellAllocationComparison asked to plot zero wells... Aborting.');
        return;
    end
end

% -------------------------------------------------------------------------
% Parameterization along well path
% Create fake z values based on the cell numbering in the well. Useful
% whenever the well isn't vertical.
WP = addDummyData(WP);
if ~opt.useZ
    for k = 1:numel(WP.inj)
        n = numel(WP.inj(k).z);
        WP.inj(k).z = (1:n)./n;
    end
    for k = 1:numel(WP.prod)
        n = numel(WP.prod(k).z);
        WP.prod(k).z = (1:n)./n;
    end
end

% -------------------------------------------------------------------------
% Extract and format the well-allocation factors
nit = numel(D.inj);
npt = numel(D.prod);
for k=1:numel(WP.inj)
   nc   = numel(WP.inj(k).z);
   atmp = zeros(nc, nit+npt);
   atmp(:,D.prod) = WP.inj(k).alloc;
   WP.inj(k).alloc = atmp;
end
for k=1:numel(WP.prod)
   nc = numel(WP.prod(k).z);
   atmp = zeros(nc, nit+npt);
   atmp(:,D.inj) = WP.prod(k).alloc;
   WP.prod(k).alloc = atmp;
end
wp(D.inj)  = WP.inj(:);
wp(D.prod) = WP.prod(:);

% -------------------------------------------------------------------------
% Calculate number of plots in 2D grid
num_plots = numel(opt.plotOnly);
clf();
nx = ceil(sqrt(num_plots));
ny = nx-1;
if (nx*ny < num_plots)
    ny = nx;
end

current_fig = gcf();
wb = waitbar(0, 'Well 1');
set(0, 'CurrentFigure', current_fig);
[xs,ys]=meshgrid(linspace(0,.85,nx+1),linspace(0,1,ny+1));
dx=diff(xs,[],2); mx=.05*max(dx(:));
dy=diff(ys,[],1); my=.05*max(dy(:));

% -------------------------------------------------------------------------
% Add color legend
axes('position',[.87 .025 .11 .95]);
nw = numel(wp);
nc = 2;
nr = ceil(nw/nc);
[x,y]=meshgrid(1:nc,1:nr);
im=nan(1,nr*nc); im(1:nw)=1:nw; im = reshape(im,nc,nr)';
pcolor(im([end:-1:1 1],[1:end end]));
txt = cell(nc,nr); txt(1:nw) = {wp.name}; txt=flipud(txt'); txt=txt(:);
ht=text(x(:)+.5,y(:)+.5,txt);
set(ht,'Color','w','HorizontalAlignment','center','FontWeight','bold','FontSize',8);
%ht2=text(x(:)+.5,y(:)+.5,txt); 
%set(ht2,'Color','w','HorizontalAlignment','center','FontSize',8);
axis off


% -------------------------------------------------------------------------
% Plot the well-allocation factors for fine/coarse scale models
for n=1:numel(opt.plotOnly)
   set(0, 'CurrentFigure', current_fig);
   i = rem(n-1,nx)+1;
   j=ny-fix((n-1)/nx);
   axes('position', ...
       [xs(j,i)+mx ys(j,i)+my dx(j,i)-2*mx  dy(j,i)-2*my]);
   %subplot(nx,ny+1,n);
   k = opt.plotOnly(n);
   
   % --------------------------------
   % Plot data
   [~,ix]   = sort(wp(k).z);
   wp(k).z = wp(k).z(ix);
   wp(k).alloc = wp(k).alloc(ix,:);
   nseg = numel(wp(k).z);
   if nseg == 1  % Need to trick Matlab
      z = [wp(k).z; wp(k).z+1];
      a = [cumsum(wp(k).alloc,1); zeros(1,numel(wp(k).alloc))];
      bwidth = 2*numel(wp(k).z)+3;
      h = barh(z, a,'stacked','BarWidth',bwidth, 'EdgeColor','none');
      args = {'YDir', 'YLim'};
   else
      z = flipud(wp(k).z);
      a = cumsum(flipud(wp(k).alloc),1);
      if nseg<21
          args = {'YDir', 'YLim'};
          h = barh(z, a,'stacked','BarWidth',opt.bwidth, 'EdgeColor','none');
      else
          h = area(z, a, eps,'EdgeColor','none'); axis tight;
          view(90,-90);
          args = {'XDir', 'XLim'};
      end
   end
   zm = min(wp(k).z);
   zM = max(wp(k).z);
   
   % ---------------------------------
   % Set axis and add legend
   set(gca,args{1},'reverse', args{2}, [zm zM + sqrt(eps)] + [-.1 .1]*(zM-zm));

   if max(wp(k).alloc(:))>0
      hl=legend(h(k),wp(k).name,'Location','SouthEast');
      set(hl,'FontSize',8); legend boxoff
   else
      hl=legend(h(k),wp(k).name,'Location','SouthWest');
      set(hl,'FontSize',8); legend boxoff
   end
   set(gca,'XTickLabel',[],'YTickLabel',[]);
   waitbar(n/numel(opt.plotOnly),wb,['Well ', num2str(k+1)])
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
