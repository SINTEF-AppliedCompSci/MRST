function [hcb, hax] = colorbarHist(q, lim, varargin)
%COLORBARHIST Make colorbar with histogram on top
%
% SYNOPSIS:
%   colorbarHist(q, lim)
%   [hc,hh] = colorbarHist(q, lim, loc)
%   [hc,hh] = colorbarHist(q, lim, loc, n)
%
% PARAMETERS:
%    q   - vector for which histogram is to be defined
%    lim - defines the range of values used to set limits of colorbar and
%          axis of histogram for q
%    loc - location of colorbar: 'East','West', or 'South' (default)
%    n   - number of bins in histogram (default: 50). See the documentation
%          of hist for the interpretation of this parameter.
%
% RETURNS:
%    hc  - graphics handle to colorbar
%    hh  - graphics handle to histogram
%

assert(any(nargin==[2,3,4]), 'Incorrect number of input arguments');

if nargin==4, 
   nbins = varargin{2}+1;
else
   nbins = 51;
end
if numel(nbins)==1
    [counts,bins] = hist(q, linspace(min(lim),max(lim),nbins));
else
    [counts,bins] = hist(q, nbins);
end

if nargin==2 || ((nargin>2) && strcmpi(varargin{1},'South'))
   cbPos = 'SouthOutside';
   cbAdd = -[0 0.09 0 0.03];
   cbLoc = 'left';
   axAdd = [0 -.03 0 .03];
   axLim = [min(lim)-eps, max(lim)+eps, -1-eps, max(1,1.1*max(counts))+eps];
   hbar  = @(x,y) bar(x,y,'hist');
elseif strcmpi(varargin{1},'east')
   cbPos = 'EastOutside';
   cbAdd = [.06 .1 -.03 -.2];
   axAdd = [.09 .1 0    -.2];
   cbLoc = 'left';
   axLim = [-1-eps, max(1,1.1*max(counts))+eps, min(lim)-eps, max(lim)+eps];
   hbar  = @(x,y) barh(x,y,'hist');
elseif strcmpi(varargin{1},'west')
   cbPos = 'WestOutside';
   cbAdd = -[.06 -.1 .03 .2];
   axAdd = -[.12 -.1  0  .2];
   cbLoc = 'right';
   axLim = [-max(1,1.1*max(counts))-eps, 1+eps, min(lim)-eps, max(lim)+eps];
   hbar  = @(x,y) barh(x,-y,'hist');
else
   error('Color bar location not supported');
end 

ha = gca;
caxis(lim)
hcb = colorbar(cbPos);
apos = get(ha,'Position');
pos=get(hcb,'Position');
set(hcb,'Position',max(pos + cbAdd,.01), 'YAxisLocation', cbLoc);
   
set(ha,'Position',apos);
hax = axes('Position',max(pos + axAdd,.01));
hbar(bins,counts)
axis tight off, axis(axLim);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','none','EdgeColor',[.2 .2 .2]);
set(gcf, 'CurrentAxes', ha)
end