function p = callMetisMatrix(A, n, varargin)
%Partition connectivity graph whilst accounting for connection strengths
%
% SYNOPSIS:
%   p = callMetisMatrix(A, n)
%
% PARAMETERS:
%   A - Symmetric connectivity graph represented as a sparse matrix.  The
%       magnitude (absolute value) of the off-diagonal non-zero entries is
%       used to measure the strengths of the individual connections.
%
%   n - Number of blocks into which the vertices of the connectivity graph
%       'A' should be partitioned.  Must be an integer value exceeding one.
%
%   x - Any additional arguments will be interpreted as options, and will
%       require metis5 to be installed. The function will then call
%       gpmetis.
%
% RETURNS:
%   p - A SIZE(A,1)-by-1 partition vector with semantics similar to those
%       implied by other partitioning functions such as 'partitionUI'.
%
% NOTE:
%   This function is a simple wrapper around the 'kmetis' or 'gpmetis'
%   utilities from the METIS package
%   (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview). 
%   The 'kmetis'/ 'gpmetis' utility is invoked through 'system' on a
%   temporary file whose name is constructed by the 'tempname' function.
%
%   If the METIS function is not avilable in the search path used by the
%   'system' function, then function 'callMetisMatrix' will fail and an
%   appropriate diagnostic message will be issued to the command window.
%   To inform the system about where to find METIS, we use a global
%   variable 'METISPATH' that can be set in your 'startup_user' function.
%   If you, for instance, use Linux and have installed METIS in
%   /usr/local/bin, you add the following line to your 'startup_user' file:
%      global METISPATH; METISPATH = fullfile('/usr','local','bin');
%
% EXAMPLE:
%   % Partition a 10-by-10-by-3 Cartesian grid into 10 blocks according to
%   % the transmissibility field derived from isotropic permeabilities.
%   %
%   G = computeGeometry(cartGrid([10, 10, 3]));
%   K = convertFrom(10 * logNormLayers(G.cartDims), milli*darcy);
%   rock = struct('perm', K);
%
%   T = computeTrans(G, rock);
%   t = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ T, [G.faces.num, 1]);
%
%   i = ~ any(G.faces.neighbors == 0, 2);
%   N = double(G.faces.neighbors(i, :));
%   A = sparse(N, fliplr(N), [ t(i), t(i) ]);
%
%   p = callMetisMatrix(A, 10);
%
% SEE ALSO:
%   `callMetis`, `partitionUI`, `processPartition`, `system`, `tempname`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   global METISPATH

   if n <= 1
      error(msgid('NumBlocks:Invalid'), ...
            'Number of blocks (''n'') must be greater than one.');
   end

   name = tempname();
   [fp, msg] = fopen(name,'w');
   if fp < 0
      error(msgid('File:Open'), ...
            'Failed to open intermediate file ''%s'': %s.', name, msg);
   end

   cleanup = onCleanup(@() fclose(fp));

   vertexWeights = false;
   if nargin < 3
       opts = '';
       binname = 'kmetis';
   else
       if isnumeric(varargin{1})
           weights = varargin{1};
           varargin = varargin(2:end);
           vertexWeights = true;
       end
       opts = sprintf('%s ', varargin{:});
       binname = 'gpmetis';
   end

   if ispc()
       binname = [binname, '.exe'];
   end

   if ~isempty(METISPATH) && isdir(METISPATH) && ...
         exist(fullfile(METISPATH, binname), 'file') == 2

       binname = fullfile(METISPATH, binname);
   end

   [st, attr, id] = fileattrib(binname);
   if ~ st
       error(id, 'METIS binary ''%s'': %s', binname, attr);
   end
   if ~ attr.UserExecute
       error('METIS binary ''%s'' is not executable', binname);
   end
   vertnum = size(A, 1);

   [i, j, v] = find(A);
   t = [i, j, v];                        clear i j v
   t = sortrows(t(t(:,1) ~= t(:,2), :));   % Exclude diagonal
   p = cumsum([ 1 ; accumarray(t(:,1), 1, [vertnum, 1]) ]);

   edgenum = size(t, 1) / 2;  % Exclude diagonal

   modifier = norm(t(:,end), inf);
   t(:,end) = ceil(10e3 * abs(t(:,end)) / modifier);
   t = t(:, 2:end) .';

   fprintf(fp, '%d %d %d1\n', vertnum, edgenum, vertexWeights);
   for i = 1 : vertnum
      if vertexWeights
          fprintf(fp, ' %d ', weights(i));
      end
      fprintf(fp, ' %d %d', t(:, p(i) : p(i + 1) - 1));
      fprintf(fp, '\n');
   end

   command = sprintf('"%s" %s "%s" %d', binname, opts, name, n);
   result  = sprintf('%s.part.%d'  , name, n);

   [stat, output] = system(command);

   if stat ~= 0
      delete(name);  delete(result);
      error('Partition:Failure', ...
           ['Failed to run METIS Graph Partitioner.  Output from ', ...
            'METIS is\n<%s>\nExit Code: %d (%s)'], output, stat, ...
            strtrim(mrstTranslateExitCode(stat)));
   end

   delete(name);

   p = load(result) + 1;
   delete(result);
end
