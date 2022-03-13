function [Bv, Phi] = basisMatrixMixed(G, CG, CS)
%Form mixed versions of the multiscale basis function matrices.
%
% SYNOPSIS:
%   [Bv, Phi] = basisMatrixMixed(G, CG, CS)
%
% PARAMETERS:
%   G, CG - Grid and coarse grid data structures, respectively.
%   CS    - Linear system data structure on coarse grid as defined by
%           function 'generateCoarseSystem'.
%
% RETURNS:
%   Bv  - Flux-basis function matrix, B*\Psi.
%   Phi - Pressure-basis function matrix, \Phi.
%
% SEE ALSO:
%   `basisMatrixHybrid`, `evalBasisFunc`, `solveIncompFlowMS`.

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


   % 1) Enumerate the (active) face unknowns according to the
   %    'CS.activeFaces' field.
   %
   actF        = CS.activeFaces;
   renum       = zeros([CG.faces.num, 1]);
   renum(actF) = 1 : numel(actF);

   % 2) Apply mixed basis function expansion to each (active) flux basis
   %    function.  This process forms the SPARSE input triplet [i,j,v] from
   %    which the hybrid flux basis matrix Bv is finally formed.
   %
   %    Each element of the cell array CS.basis(actF) is a 'V' tuple as
   %    defined by function 'evalBasisFunc'.  We need to renumber the
   %    coarse face flux unknowns, so pass the face enumeration array along
   %    with the dynamic 'V' tuples.  Finally, as the flow direction across
   %    face 'f' is from block CG.faces.neighbors(f,1) to block
   %    CG.faces.neighbors(f,2), the expansion process needs access to the
   %    block topology in order to determine the direction of the flow.
   %
   [i,j,v] = cellfun(@(x) expand_mixed(x, renum, CG.faces.neighbors), ...
                     CS.basis(actF), 'UniformOutput', false);
   Bv = sparse(vertcat(i{:}), vertcat(j{:}), vertcat(v{:}), ...
               size(G.cells.faces,1), numel(actF));

   % 3) Apply mixed basis function expansion to each (active) pressure
   %    basis function.  This process forms the SPARSE input triplet
   %    [i,j,v] from which the hybrid pressure basis matrix Phi is finally
   %    formed.
   %
   %    Each element of the cell array CS.basisP(actF) is a 'P' tuple as
   %    defined by function 'evalBasisFunc'.
   %
   [i,j,v] = cellfun(@(x) expand_mixed(x, renum, CG.faces.neighbors), ...
                     CS.basisP(actF), 'UniformOutput', false);
   Phi = sparse(vertcat(i{:}), vertcat(j{:}), vertcat(v{:}), ...
                G.cells.num, numel(actF));
end

%--------------------------------------------------------------------------

function [i, j, v] = expand_mixed(c, fno, n)
   % Expand a basis function tuple (either flux or pressure) 'c' into a
   % SPARSE matrix input triplet [i,j,v].  The values of a mixed basis
   % function are entered directly into the corresponding column (active
   % coarse face number) of the matrix.  We recall that the input tuple,
   % 'c', is a cell array of the form
   %
   %    1  2  3  4  5  6  7
   %   {i, s, f, b, n, m, o}
   %
   % with symbols as defined in 'evalBasisFunc'.  Items 5-7 are ignored in
   % this part of the implementation, while item 1 directly provides the
   % desired output vector 'i' (no further processing required).
   %
   % Our task, then, is to convert item 2 into the output vector 'v', and
   % to build output vector 'j' from scratch.  The 'j' vector is
   % constructed by simply repeating the active coarse face index a number
   % of times equal to the number of values in 'i' (or 'v').  The output
   % 'v' needs to take into account the flow direction induced by the
   % coarse grid topology.  Specifically, the basis functions are generated
   % on the assumption that the flow across coarse face 'f' is from coarse
   % block n(f,1) to coarse block n(f,2).

   % 1) Extract 'i' values directly from input tuple.
   i = c{1};

   % 2) Repeat active coarse face number (fno(c{3})) a number of times
   %    equal to the number of entries in 'i'.
   j = rldecode(fno(c{3}), numel(i));

   % 3) Reverse flow direction if this block (c{4}(1)) is the *inflow*
   %    block for the basis function across coarse face c{3}.
   s = 2*double(n(c{3},1) == c{4}(1)) - 1;
   v = s .* c{2};
end
