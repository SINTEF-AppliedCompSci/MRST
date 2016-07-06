function writeMRSTGrid(G, fn)
%Serialise an MRST Grid Structure to a file.
%
% SYNOPSIS:
%   writeMRSTGrid(G, fn)
%
% PARAMETERS:
%   G  - Grid structure.
%
%   fn - File name.  This will be passed to function FOPEN using mode 'wt'.
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   writeDeck, fopen.

% $Date: 2012-10-10 20:45:49 +0200 (Wed, 10 Oct 2012) $
% $Revision: 10060 $

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

   d = fileparts(fn);
   if ~isdir(d),
      [success, msg, id] = mkdir(d);
      if ~success,
         error(id, 'Failed to create directory ''%s'': %s', d, msg);
      end

      clear success msg id
   end
   clear d

   [fid, msg] = fopen(fn, 'wt');

   if fid < 0, error('Failed to open file ''%s'': %s', fn, msg); end

   if isfield(G, 'nodes');
      nn  = G.nodes.num;
      nfn = size(G.faces.nodes, 1);
   else
      nn  = 1;
      nfn = 1;
   end
   [nc, nf] = deal(G.cells.num, G.faces.num);
   d        = size(G.cells.centroids, 2);
   nhf      = size(G.cells.faces, 1);

   rptfmt = @(s, n) [repmat([' ', s], [1, n]), '\n'];

   fprintf(fid, rptfmt('%d', 6)        , d, nc, nf, nn, nfn, nhf);
   fprintf(fid, '%d\n'                 , size(G.cells.faces, 2) > 1);
   fprintf(fid, '%d\n'                 , isfield(G.cells, 'indexMap'));

   cartDims = zeros([1, G.griddim]);
   if isfield(G, 'cartDims'),
      cartDims = G.cartDims;
   end
   fprintf(fid, rptfmt('%d', G.griddim), cartDims);

   if isfield(G, 'nodes'),
      % Nodes
      fprintf(fid, rptfmt('%.18e', d), G.nodes.coords .');
      fprintf(fid, '%d\n', G.faces.nodePos - 1);
      fprintf(fid, '%d\n', G.faces.nodes - 1);
   else
      fprintf(fid, rptfmt('%.18e', d), zeros([1, d]));
      fprintf(fid, '%d\n', [1; repmat(2, [G.faces.num - 1, 1])] - 1);
      fprintf(fid, '%d\n', 0);
   end

   % Faces
   fprintf(fid, '%d %d\n', G.faces.neighbors .' - 1);

   assert (all(isfield(G.faces, {'areas', 'centroids', 'normals'})), ...
           'Grid must include geometry.');

   fprintf(fid, '%.18e\n'         , G.faces.areas       );
   fprintf(fid, rptfmt('%.18e', d), G.faces.centroids .');
   fprintf(fid, rptfmt('%.18e', d), G.faces.normals   .');

   % Cells
   fprintf(fid, '%d\n', G.cells.facePos - 1);
   ncol = 1 + (size(G.cells.faces, 2) > 1);
   fprintf(fid, rptfmt('%d', ncol), G.cells.faces(:, 1:ncol).' - 1);
   if isfield(G.cells, 'indexMap'),
      fprintf(fid, '%d\n', G.cells.indexMap - 1);
   end

   assert (all(isfield(G.cells, {'volumes', 'centroids'})), ...
           'Grid must include geometry.');

   fprintf(fid, '%.18e\n'         , G.cells.volumes     );
   fprintf(fid, rptfmt('%.18e', d), G.cells.centroids .');

   fclose(fid);
end
