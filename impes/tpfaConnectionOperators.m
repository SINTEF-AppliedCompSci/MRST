function OP = tpfaConnectionOperators(G, W, np)
%Compute static (i.e., geometry/topology dependent) connection operators
%
% SYNOPSIS:
%   OP = tpfaConnectionOperators(G, wells, np)
%
% PARAMETERS:
%   G  - Grid structure.
%
%   wells -
%        Well data structure as defined by function 'addWell'.  Pass an
%        empty array (e.g., wells=[]) if the current simulation model does
%        not include wells.
%
%   np - Number of fluid phases (and components) involved in the
%        subsequent simulation run.  Typically two or three.
%
% RETURNS:
%   OP - Static connection operators for use in TPFA-based IMPES
%        discretisation of the black-oil system.  This is (currently) a
%        simple data structure, but its fields are intentionally not
%        documented as they may change without prior notice.
%
%        Callers should treat this an (opaque) implementational detail of
%        the 'impesTPFA' solver and related functions.
%
% BUGS:
%   The semantics of return value 'OP' are *strongly* tied to the internal
%   structure of the 'impesTPFA' solver.  There is little use in having an
%   'OP' structure in absence of a TPFA discretisation.
%
% SEE ALSO:
%   impesTPFA, addWell.

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


   N      = double(G.faces.neighbors);
   nf     = diff(G.cells.facePos);
   cellno = rldecode(1 : G.cells.num, nf, 2) .';
   t      = G.faces.neighbors(G.cells.faces(:,1), 1) == cellno;
   sgn    = 2*t - 1;
   i      = ~any(N == 0, 2);
   f2hf   = accumarray([double(G.cells.faces(:,1)), double(2 - t)], ...
                       (1 : numel(cellno)) .', [G.faces.num, 2]);

   Nix    = @(f, c, col) sub2ind(size(N), double(f), col(f, c));
   ocolix = @(f, c     ) 1 + double(N(f,1) == c); % Other column index
   scolix = @(f, c     ) 1 + double(N(f,1) ~= c); % Same column index

   connno = double(G.cells.faces(:,1));           % Connection numbers
   active = i(G.cells.faces(:,1));
   other  = N(Nix(G.cells.faces(:,1), cellno, ocolix));

   gpot = gravity_potential(G, cellno, reshape(gravity, [], 1));

   [I, J] = blockDiagIndex(rldecode(np, numel(nf)), nf);

   if ~isempty(W),
      nperf   = cellfun('prodofsize', { W.cells });
      wcellno = vertcat(W.cells);

      nwc     = numel(wcellno);
      nw      = numel(nperf);
      wother  = rldecode(G.cells.num + (1 : nw), nperf, 2) .';

      wconnno = G.faces.num + (1 : nwc).';
      wsgn    = repmat(-1, [nwc, 1]);

      wgpot   = zeros([nwc, 1]);

      Iw = bsxfun(@plus, (1 : np) .', reshape((wcellno - 1) .* np, 1, []));
      Jw = rldecode(J(end) + (1 : nw), np * nperf, 2);

      assert (numel(Iw) == numel(Jw),                            ...
             ['Counting error whilst defining multi-sys index ', ...
              'expressions for wells.']);

      I  = [I; reshape(Iw, [], 1)];
      J  = [J; reshape(Jw, [], 1)];
   else
      [wcellno, wother, wconnno, wsgn, wgpot] = deal([]);
   end

   OP = struct('active' , active , ...  % Does half-face support flow?
               'cellno' , cellno , ...  % Standard RLDECODE incantation
               'connno' , connno , ...  % Connection # (= cells.faces(:,1))
               'other'  , other  , ...  % Outside cell (wrt. cellno)
               'sgn'    , sgn    , ...  % Outward flux sign (wrt. cellno)
        ...
               'wcellno', wcellno, ...  % VERTCAT(W.cells)
               'wconnno', wconnno, ...  % Perfs numbered from faces.num+1
               'wother' , wother , ...  % Outside DOF (== well number)
               'wsgn'   , wsgn   , ...  % Outward flux sign (wrt. wcellno)
        ...
               'I', I   , 'J', J , ...  % (I,J) for RHS in solve_multiple
        ...
               'gpot'   , gpot   , ...  % (x_{ij} - x_i) * gravity (per hf)
               'wgpot'  , wgpot  , ...  % (x_{wi} - x_i) * gravity (== 0)
        ...
               'hf'     , @(f, c) f2hf(Nix(f, c, scolix)));  % (f,c) -> hf
end

%--------------------------------------------------------------------------

function gpot = gravity_potential(G, cellNo, g)
   g = g(1 : G.griddim);

   if norm(g) > 0,

      gpot = G.faces.centroids(G.cells.faces(:,1), :) - ...
             G.cells.centroids(cellNo            , :);

      gpot = gpot * g;

   else

      gpot = sparse(size(G.cells.faces, 1), 1);

   end
end
