function [fmob, fz, dfmob, Af] = ...
      tpfaUpwindStateVars(G, press, z, rho, cmob, dcmob, bc, wells, OP, Ac)
%Compute upwind mobilities and component volumes
%
% SYNOPSIS:
%   [fmob, fz, dfmob] = ...
%     tpfaUpwindStateVars(G, press, z, rho, cmob, dcmob, bc, wells, OP)
%
%   [fmob, fz, dfmob, Af] = ...
%     tpfaUpwindStateVars(G, press, z, rho, cmob, dcmob, bc, wells, OP, Ac)
%
% PARAMETERS:
%   G  - Grid structure.
%
%   press -
%        Cell pressure.  One (positive) scalar value for each cell in the
%        grid 'G'.
%
%   z  - Component volume densities at surface conditions.  One
%        non-negative value for each fluid component in each cell in the
%        grid 'G'.
%
%   rho -
%        Phase densities.  One positive scalar value for each fluid phase
%        (at reservoir conditions) for each cell in the grid 'G'.
%
%   cmob -
%        Phase mobilities.  One non-negative scalar value for fluid phase
%        (at reservoir conditions) for each cell in the grid 'G'.
%
%   dcmob -
%        Phase mobility derivatives (wrt. saturation).  One scalar value
%        for each fluid phase (at reservoir conditions) differentiated with
%        respect to each phase saturation (at reservoir conditions) for
%        each cell in the grid 'G'.
%
%   bc - Boundary condition structure as defined by function 'addBC'.  If
%        ISEMPTY(bc), the upwind calculation will assume no-flow boundary
%        conditions on all boundary faces.
%
%   wells -
%        Well data structure as defined by function 'addWell'.  Specify an
%        empty array (e.g., wells=[]) if the simulation model does not
%        include well perforations.
%
%   OP - TPFA connection operators as defined by function
%        'tpfaConnectionOperators'.  This operator must be consistent with
%        the values of 'bc' and 'wells'.
%
%   Ac - Fluid composition matrix per cell.  Transforms saturations at
%        reservoir conditions to surface volumes per pore volume at
%        surface conditions.  Mathematically, Ac=Rc*inv(Bc).  OPTIONAL.
%        Only needed if computing upwind composition matrices per
%        interface.
%
% RETURNS:
%   fmob -
%        Upwind phase mobilities.  One non-negative scalar value for each
%        fluid phase (at reservoir conditions) for each connection (i.e.,
%        grid face and, possibly, well perforation) defined in the
%        simulation model specified by (G,bc,wells).
%
%   fz - Upwind component volume densities.  One non-negative scalar value
%        for fluid component (at surface conditions) for each connection
%        defined in the simulation model specified by (G,bc,wells).
%
%   dfmob -
%        Upwind phase mobility derivatives.  One scalar value for each
%        fluid phase (at reservoir conditions) differentiated with respect
%        to each phase saturation (at reservoir conditions) for each
%        connection defined in the simulation model specified by the
%        (G,bc,wells) tuple.
%
%   Af - Upwind fluid composition matrix per interface.  Mathematically,
%        Af=Rf*inv(Bf).  OPTIONAL.  Only computed if specifically requested
%        *and* if the caller supplies input values 'Ac'.  If requested and
%        no input values are present, an error condition occurs.
%
% NOTE:
%   The upwind values are simply taken from the appropriate cell values by
%   considering differences in hydraulic head for each fluid phase (at
%   reservoir conditions).

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


   error(nargchk(9, 10, nargin, 'struct'));

   N = tpfaExtendedConnections(G, bc, wells);

   % Pre-allocate result arrays.
   nconn = size(N, 1);
   np    = size(z, 2);

   fmob  = zeros([nconn, np * 1 ]);
   fz    = zeros([nconn, np * 1 ]);
   dfmob = zeros([nconn, np * np]);

   % Determine upstream directions based on hydraulic head differences.
   %
   [col, i, Neigh] = upwind_directions(N, G.faces.neighbors, ...
                                       press, rho, bc, OP);

   if nargin == 10,
   % Build phase upwind approximation to component-phase density matrix A
   % for each connection: column i*np+k represents the component densities
   % in phase k of cell i.  For each face ij and each phase k, use the
   % column in the upwind cell (i or j) as upwind approximation on the face
   % ij.
   cix       = zeros(nconn, np);
   ix        = sub2ind(size(Neigh), repmat((1:size(Neigh,1))', [1,np]), col);
   cix(i,:)  = Neigh(ix);
   cix(~i,:) = repmat(sum(N(~i, :), 2), [1, np]);
   Af        = expandBlockDiagMatrix(Ac, np, cix);
   elseif nargout == 4,
      % Caller requested 'Af', but failed to supply 'Ac'.  This is an
      % error.
      error('Cannot compute ''Af'' in absence of ''Ac''.');
   end

   % Compute upstream mobilities and phase volumes.  Note that mobility
   % derivatives are only needed for time-step estimation purposes and need
   % therefore only be computed on internal faces.
   %
   rix = (1 : size(Neigh, 1)) .';
   cix = @(c) Neigh(sub2ind(size(Neigh), rix, c));

   for p = 1 : np,
      c = cix(col(:, p));

      fmob(i, p) = cmob(c, p);
      fz  (i, p) = z   (c, p);

      dfmob(i, p : np : end) = dcmob(c, p : np : end);
   end

   c = sum(N(~i, :), 2);
   fmob(~i,:)   = cmob(c, :);
   fz  (~i,:)   = z   (c, :);
end

%--------------------------------------------------------------------------

function [ix, i, N] = upwind_directions(N, neigh, press, rho, bc, OP)
   % N == G.faces.neighbors, with cells repeated for boundary faces...
   phi = repmat(press, [1, size(rho, 2)]);

   cellno = [ OP.cellno(OP.active) ; OP.wcellno ];
   connno = [ OP.connno(OP.active) ; OP.wconnno ];
   gpot   = [ OP.gpot(  OP.active) ; OP.wgpot   ];
   sgn    = [ OP.sgn(   OP.active) ; OP.wsgn    ];

   if ~isempty(bc),
      f = double(bc.face);
      c = double(sum(neigh(f,:), 2));
      i = OP.hf(f, c);

      cellno = [OP.cellno(i) ; cellno];
      connno = [OP.connno(i) ; connno];
      gpot   = [OP.gpot(  i) ; gpot  ];
      sgn    = [OP.sgn(   i) ; sgn   ];

      clear f c i
   end

   i  = all(N ~= 0, 2);
   j  = cumsum(double(i));

   nc = OP.cellno(end);

   pick   = i(connno);
   connno = j(connno(pick));
   cellno = cellno  (pick) ;
   gpot   = gpot    (pick) ;

   G  = sparse(connno, cellno, sgn .* gpot, max(connno), nc);

   dphi       = zeros([size(N,1), size(phi,2)]);
   dphi(i, :) = phi(N(i,1), :) - phi(N(i,2), :) + G*rho(1 : nc, :);

   ix = 1 + (dphi < 0);
   ix = ix(i, :);
   N  = N (i, :);
end

%--------------------------------------------------------------------------

function B = expandBlockDiagMatrix(A, m, f)
%  A - block diagonal matrix
%  m - block size
%  f - upwind cell(i.e. block number) for each new block in B
%
%  Expand A to B by choosing columns from A: block i in B consists of first
%  column from block f(i,1) in A, second column from block f(i,2) in A and
%  so on.

   % compute column numbers to extract from A (in sequence).
   F = m*(f-1) + repmat(1:m,[size(f, 1), 1]);
   F = reshape(F', [], 1);

   % Extract blocks by single indexing, A(k);
   M     = repmat(m, [size(A, 2)/m, 1]);
   [i,j] = blockDiagIndex(M, M);
   %i     = int64(i);
   %j     = int64(j);
   k     = i+(j-1)*size(A,1);

   % Position of each column in A(k)
   pos   = cumsum([1; rldecode(M,m)]);

   % Extract correct block for each entry in f
   ix    = mcolon(pos(F), pos(F+1)-1)';

   % Build matrix
   M     = repmat(m, [size(f,1),1]);
   [I,J] = blockDiagIndex(M, M);
   B     = sparse(I, J, A(k(ix)));
end
