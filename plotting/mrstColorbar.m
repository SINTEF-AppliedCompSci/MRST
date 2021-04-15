function [hcb, hax] = mrstColorbar(varargin)
%Append a colorbar with an accompanying histogram to the current axis
%
% SYNOPSIS:
%   mrstColorbar
%   mrstColorbar(ax)
%   mrstColorbar(..,values)
%   [hc,hh] = mrstColorbar(..,location)
%   [hc,hh] = mrstColorbar(..,location, logscale)
%   [hc,hh] = mrstColorbar(..,location, logscale, limits)
%
% PARAMETERS:
%    ax       - add colorbar to axes AX instead of current axis 
%    values   - create the accompanying histogram using the VALUES vector.
%               If not specified, the routine will pick values from CData
%               of the current axes. Notice that many 3D plots in MRST only
%               show the outer surface of a grid and hence the histogram
%               will not represent the full 3D dataset unless this is
%               explicitly specified.
%    location - location of colorbar, same as for MATLAB's colorbar except
%               that the colorbar is always placed outside of the plot 
%    logscale - flag indicating that the data displayed are shown on a
%               logarithmic scale. This will manipulate the the tick marks
%               on the colorbar so that they are given in scientific
%               notation
%    limits   - lower and upper limits for the histogram bins. NB! Setting
%               this parameter will also reset the caxis accordingly.
%
% RETURNS:
%    hc  - graphics handle to colorbar
%    hh  - graphics handle to histogram
%
% SEE ALSO:
%   `hist`

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
if ~exist('verLessThan','file') || verLessThan('matlab','8.4')
    error('Not implemented for your MATLAB version. Try colorbarHist instead');
end
nargs = nargin;

remaxes = gca;

% Is first input argument a graphics handle? If not, assume the colorbar is
% to be added to the current axes
if (nargs>0) && (numel(varargin{1})==1) && ishandle(varargin{1})
    ha       = varargin{1};
    varargin = varargin(2:end);
    nargs    = nargs-1;
else
    ha = gca;
end

% Is the first (remaining) argument an array? If not, assume that the data
% to be used in the histogram should be picked from the children of the
% current axes.
dataGiven = false;
if (nargs>0) && isa(varargin{1},'numeric')
    CData    = varargin{1};
    varargin = varargin(2:end);
    nargs    = nargs-1;
    dataGiven = true;
else
    h=get(ha,'Children');
    CData = [];
    for i=1:numel(h)
        if isfield(h(i),'CData') && ~isempty(h(i).CData)
            CData=[CData; h(i).CData]; %#ok<AGROW>
        end
    end
end

%opos = get(gca,'OuterPosition');
oapos = get(ha,'Position');

% Is the first argument a location specifier? If not, use the default
% value, which is to the south
if nargs>0
    location = lower(varargin{1});
    if isempty(strfind(location,'outside')) %#ok<STREMP>
       location=[location 'outside']; 
    end
else
    location = 'southoutside';
end

% Is the second argument true? If so, the data are on logaritmic scale
if nargs>1
    logflag = varargin{2};
    if dataGiven && logflag
        CData = log10(CData);
    end
else
    logflag = false;
end

% If there it is a 
if nargs>2
    limits = varargin{3};
    assert(numel(limits)==2,'Need two values to specify color axis');
    if logflag
        limits = log10(limits);
    end
else
    limits = caxis(ha);
end

% Specify where to put the histogram and its orientation
switch location
    case 'southoutside'
        [az,el] = deal(  0, 90);
        histogramPos = @(x) [x(1)+.05*x(3) x(2) .9*x(3) x(4)];
    case 'northoutside'
        [az,el] = deal(  0,-90);
        histogramPos = @(x) [x(1)+.05*x(3) x(2) .9*x(3) x(4)];  
    case 'eastoutside'
        [az,el] = deal(-90, 90);
        histogramPos = @(x) [x(1) x(2)+.05*x(4) x(3) .9*x(4)];  
    case 'westoutside'
        [az,el] = deal( 90,-90);
        histogramPos = @(x) [x(1) x(2)+.05*x(4) x(3) .9*x(4)];
    otherwise
        error('Unknown location: %s', location);
end

% Set up a colorbar, insert the histogram in its place, and move the
% colorbar slightly below/above/left/right of the histogram
hcb   = colorbar(ha, location); drawnow
apos  = get(ha, 'Position');
cpos  = get(hcb,'Position');
hhpos = histogramPos(cpos);
cbpos = colorbarPos(hhpos,oapos,location);
set(hcb,'Position',cbpos);
if logflag
    ticks = get(hcb, 'XTick');
    ticks = ceil(min(ticks)):floor(max(ticks));
    newticks = arrayfun(@(x) ['1e', num2str(x)], ticks, 'UniformOutput', false);
    set(hcb, 'XTick', ticks, 'XTickLabel', newticks)
end
set(ha,'Position',apos)
hax = axes('Position', hhpos);
hh = histogram(CData); view(az,el);
hh.BinLimits = limits;
hh.NumBins   = 50;
set(hax,'XLim', limits);
axis off;
axes(remaxes);
% if nargout<2,  clear hax, end
% if nargout~=1, clear hcb, end
end

function newpos=colorbarPos(c, a, location)
  switch location
      case 'southoutside'
          dy = c(2)-a(2);
          newpos = [c(1) a(2) c(3) .85*dy];
      case 'northoutside'
          dy = a(2)+a(4)-c(2)-c(4);
          newpos = [c(1) c(2)+c(4)+.15*dy c(3) .85*dy];
      case 'westoutside'
          dx = c(1)-a(1);
          newpos = [a(1) c(2) .85*dx c(4)];
      case 'eastoutside'
          dx = a(1)+a(3)-c(1)-c(3);
          newpos = [c(1)+c(3)+.15*dx c(2) .85*dx c(4)];
  end
end