function G = emptyGrid(varargin)
%Make empty grid.
%
% SYNOPSIS:
%   G = emptyGrid()
%   G = emptyGrid(numcells)
%   G = emptyGrid(numcells, nodes)
%
% PARAMETERS:
%   numcells - Number of cells.  Used to allocate cells.facePos.  Treated
%              as zero if unspecified.
%
%   nodes    - Node coordinates.  N-by-d array.  Treated as ZEROS([0,3]) if
%              unspecified.
%
% RETURNS:
%   G - Grid structure as described by grid_structure.

   numcells = 0;             if nargin > 0, numcells = varargin{1}; end
   nodes    = zeros(0, 3);   if nargin > 1, nodes    = varargin{2}; end


   %% Empty shell:
   cells = struct('num',       double(numcells), ...
                  'facePos',   ones(numcells+1, 1), ...
                  'faces',     zeros(0, 2));
   faces = struct('num',       double(0), ...
                  'nodePos',   1, ...
                  'nodes',     zeros(0, 1), ...
                  'neighbors', zeros(0, 2));
   nodes = struct('num',       double(size(nodes, 1)),...
                  'coords',    nodes);

   G     = struct('cells',     cells, ...
                  'faces',     faces, ...
                  'nodes',     nodes, ...
                  'type',      {mfilename} );
