function varargout = plotBlockAndNeighbors(G, CG, block, varargin)
%Plot a coarse block and its neighbors to current axes (reversed Z axis).
%
% Different colors and levels of transparency are used to distinguish the
% blocks.  The block itself is plotted in blue color and the neighbors in
% cyan, magenta, yellow, white, and green (with a cyclic repeat).  Faults
% are plotted using red patches.
%
% SYNOPSIS:
%       plotBlockAndNeighbors(G, CG, block)
%       plotBlockAndNeighbors(G, CG, block, 'pn1', 'pv1', ...)
%   h = plotBlockAndNeighbors(...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   CG      - Coarse grid data structure.
%
%   block   - The coarse block to be plotted.
%
%   'pn'/pv - List of property names/property values.  OPTIONAL.
%             This list will be passed directly on to function PATCH
%             meaning all properties supported by PATCH are valid.
%
%             As a special case, function 'plotBlockAndNeighbors' supports
%             an option 'PlotFaults', whose argument is either TRUE or
%             FALSE, which indicates whether or not the faults should be
%             added to the graphical output.
%
% RETURNS:
%   h  - Handle to resulting patch object.  The patch object is added
%        directly to the current AXES object (GCA).
%        OPTIONAL.  Only returned if specifically requested.
%
% NOTES:
%   Function 'plotBlockAndNeighbors' is implemented in terms of the
%   high-level functions 'plotGrid' and 'plotFaces', which again use the
%   low-level function PATCH.  If a separate axes is needed for the
%   graphical output, callers should employ function newplot prior to
%   calling 'plotBlockAndNeighbors'. The function relies on a specific set
%   of values for the properties 'FaceColor' and 'FaceAlpha'.
%
% SEE ALSO:
%   plotFaces, plotGrid, patch, newplot.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


opt = struct('PlotFaults', true);
opt = merge_options(opt, varargin{:});

ix = find(strcmpi(varargin(1 : 2 : end), 'plotfaults'));
if ~isempty(ix),
   varargin = varargin([1 : 2*(ix-1), 2*ix + 1 : end]);
end

% Find partitioning
p = CG.partition;

% Set correct alpha values for block, fault, neighbors and faults in
% neighbors
if opt.PlotFaults,
   FaceAlpha = [0.5 0.5 0.05 0.2];
else
   FaceAlpha = [1.0 1.0 0.1 1.0];
end

% Plot the block
h = plotGrid(G, find(p==block),...
             'FaceColor','b','FaceAlpha',FaceAlpha(1), varargin{:});

% Find all faults and plot them
if opt.PlotFaults && isfield(G.faces, 'tag'),
   pB = [0; p];
   pN = pB(G.faces.neighbors+1);
   f  = find(any(pN==block,2));
   plotFaces(G, f(G.faces.tag(f)>0),...
                 'FaceColor','r','FaceAlpha', FaceAlpha(2), varargin{:});
end

% Find neighbors
nb = reshape(CG.faces.neighbors(any(CG.faces.neighbors==block,2),:),[],1);
nb(nb==block) =[];
nb = unique(nb(nb>0));

% Plot the neighbors
col = ['g','c','m','y','w'];
for j=1:numel(nb)

   % plot neighbouring block #j
   plotGrid(G,find(p==nb(j)), 'FaceColor',col(mod(j,numel(col))+1),...
            'FaceAlpha', FaceAlpha(3), varargin{:});

   % find all faults and plot them
   if opt.PlotFaults && isfield(G.faces, 'tag'),
      pN = pB(G.faces.neighbors+1);
      f  = find(any(pN==nb(j),2));
      plotFaces(G, f(G.faces.tag(f)>0), ...
               'FaceColor','r','FaceAlpha',FaceAlpha(4),varargin{:});
   end
end
axis tight off; view(3);

if nargout > 0, varargout{1} = h; end
