function varargout = plotBlockAndNeighbors(CG, block, varargin)
%Plot a coarse block and its neighbours to current axes (reversed Z axis).
%
% Different colours and levels of transparency are used to distinguish the
% blocks.  The block itself is plotted in blue colour and the neighbours
% using colours from a brightened COLORCUBE colour map.  Faults are plotted
% using gray patches (RGB = REPMAT(0.7, [1, 3])) with red edge colours.
%
% SYNOPSIS:
%       plotBlockAndNeighbors(CG, block)
%       plotBlockAndNeighbors(CG, block, 'pn1', 'pv1', ...)
%   h = plotBlockAndNeighbors(...)
%
% PARAMETERS:
%   CG      - Coarse grid data structure.
%
%   block   - The coarse block to be plotted.
%
% KEYWORD ARGUMENTS:
%
%   PlotFaults - Two-element `logical` vector, the entries of which specify
%                whether or not fault faces should be added to the
%                graphical output of the 'block' and its neighbours,
%                respectively.
%
%                DEFAULT: `PlotFaults = TRUE([1,2])` (attach fault faces
%                to both the 'block' and all of its neighbours).
%
%   Alpha      - `(2 + max(find(PlotFaults)))`-element numeric vector,
%                values in [0,1], specifying scalar transparency
%                (`AlphaData`) values for the `block`, its neighbours, and
%                the fault faces of the 'block' and its neighbours,
%                respectively.
%
%                DEFAULT: `Alpha = ONES([1,4])` (no transparency in any of
%                the final objects--all objects drawn opaquely).
%
%   'Any'   -    Additional keyword arguments will be passed directly on to
%                function `patch` meaning all properties supported by
%                `patch` are valid.
%
%
% RETURNS:
%   h - Handle to resulting patch objects.  The patch objects are added
%       directly to the current `axes` object (`gca`).
%       OPTIONAL.  Only returned if specifically requested.
%
% EXAMPLE:
%   % Plot a block and its neighbours from a coarse partitioning of the
%   % "model 3" synthetic geometry
%   require coarsegrid  %  Make "coarse block" concept meaningful
%
%   % Generate geometry
%   G = processGRDECL(makeModel3([100, 60, 15]));
%
%   % Partition geometry
%   p = partitionUI(G, [ 5, 5, 3 ]);
%   p = compressPartition(processPartition(G, p));
%
%   % Generate coarse grid
%   CG = generateCoarseGrid(G, p);
%
%   % Plot selected block (# 37) and its neighbours
%   plotBlockAndNeighbors(CG, 37, 'Alpha', repmat(0.75, [1, 4]))
%   view(-145, 26)
%
% NOTES:
%   Function `plotBlockAndNeighbors` is implemented in terms of plotting
%   function `plotFaces` which in turn uses the built-in function `patch`.
%   If a separate axes is needed for the graphical output, callers should
%   employ function `newplot` prior to calling `plotBlockAndNeighbors`.  This
%   function relies on a specific set of values for the properties
%   `FaceColor` and `FaceAlpha`.
%
% SEE ALSO:
%   `plotFaces`, `patch`, `newplot`

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

opt = struct('PlotFaults', true([1, 2]), ...
             'Alpha'     , ones([1, 4]));
[opt, varargin] = merge_options(opt, varargin{:});

check_input(opt);

% Initialize parameters
p = CG.partition;
G = CG.parent;
h = [];

if any(opt.PlotFaults) && isfield(G.faces, 'tag'),
   findFF = findFaultFaces(G, p);
end

% Plot the block
plot_fault = @(faces, alpha) ...
    plotFaces(G, faces, 'FaceColor', [ 0.7, 0.7, 0.7 ], ...
              'EdgeColor', 'r', 'FaceAlpha', alpha, varargin{:});

plot_block = @(faces, colour, alpha) ...
    plotFaces(G, faces, 'FaceColor', colour, ...
              'FaceAlpha', alpha, varargin{:});

bdryFaces = @(blk) boundaryFaces(G, p == blk);

blockFaces = false([G.faces.num, 1]);
blockFaces(bdryFaces(block)) = true;

% Find all faults and plot them
if opt.PlotFaults(1) && isfield(G.faces, 'tag'),
   faultfaces = findFF(block);

   hf = plot_fault(faultfaces, opt.Alpha(3));
   h  = [hf ; h];

   blockFaces(faultfaces) = false;
end

hf = plot_block(blockFaces, 'b', opt.Alpha(1));
h  = [hf ; h];

%-----------------------------------------------------------------------

% Plot the neighbors
cN = CG.faces.neighbors;
nb = reshape(cN(any(cN == block, 2), :), [], 1);  clear cN
nb(nb == block) = [];
nb = unique(nb(nb > 0));

m = max(8, numel(nb));
col = (2*colorcube(m) + 1) ./ 3;
for j = 1 : numel(nb),
   blockFaces = false([G.faces.num, 1]);
   blockFaces(bdryFaces(nb(j))) = true;

   % Find all faults and plot them
   if opt.PlotFaults(2) && isfield(G.faces, 'tag'),
      faultfaces = findFF(nb(j));

      hf = plot_fault(faultfaces, opt.Alpha(4));
      h  = [hf ; h];                                            %#ok<AGROW>

      blockFaces(faultfaces) = false;
   end

   hf = plot_block(blockFaces, col(j,:), opt.Alpha(2));
   h  = [hf ; h];                                               %#ok<AGROW>
end

axis tight off; view(3);

if nargout > 0, varargout{1} = h; end

%--------------------------------------------------------------------------

function f = findFaultFaces(G, p)
p = [ 0 ; reshape(p, [], 1) ];
N = p(G.faces.neighbors + 1);

f0 = @(b) find(any(N == b, 2));

f1 = @(i) i(G.faces.tag(i) > 0);

f  = @(b) f1(f0(b));

%--------------------------------------------------------------------------

function check_input(opt)
p = find(opt.PlotFaults, 1, 'last');
if isempty(p),
   % Caller doesn't want faults.  Specify position zero to avoid ASSERT
   % failing with a diagnostic that's hard to comprehend:
   %
   %    Error using >=
   %    Matrix dimensions must agree.
   %
   % That results from what would effectively be
   %
   %    NUMEL(opt.Alpha) >= 2 + []
   %
   % otherwise.
   p = 0;
end

assert (numel(opt.Alpha) >= 2 + p, ...
       ['Option ''Alpha'' must specify enough AlphaData values ', ...
        'to cover the requested plotting modes.']);
