function G = emptyGrid(varargin)
%Make empty grid.
%
% SYNOPSIS:
%   G = emptyGrid()
%   G = emptyGrid(numcells)
%   G = emptyGrid(numcells, nodes)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.

   numcells = 0;             if nargin > 0, numcells = varargin{1}; end
   nodes    = zeros(0, 3);   if nargin > 1, nodes    = varargin{2}; end


   %% Empty shell:
   cells = struct('num',       double(numcells), ...
                  'facePos',   ones(numcells+1, 1, 'int32'), ...
                  'faces',     zeros(0, 2, 'int32'));
   faces = struct('num',       double(0), ...
                  'nodePos',   int32(1), ...
                  'nodes',     zeros(0, 1, 'int32'), ...
                  'neighbors', zeros(0, 2, 'int32'));
   nodes = struct('num',       double(size(nodes, 1)),...
                  'coords',    nodes);

   G     = struct('cells',     cells, ...
                  'faces',     faces, ...
                  'nodes',     nodes, ...
                  'type',      {mfilename} );
