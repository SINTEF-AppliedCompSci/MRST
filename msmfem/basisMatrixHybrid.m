function [Bv, Phi] = basisMatrixHybrid(G, CG, CS)
%Form hybrid versions of the multiscale basis function matrices.
%
% SYNOPSIS:
%   [Bv, Phi] = basisMatrixHybrid(G, CG, CS)
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
%   `basisMatrixMixed`, `evalBasisFunc`, `solveIncompFlowMS`.

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


   % 1) Enumerate the (active) cell face unknowns according to the
   %    'CS.activeCellFaces' field.  As the basis functions are stored in
   %    the CS.activeFaces elements of the various 'CS.basis*' fields, we
   %    also create an alias for these coarse faces.
   %
   actB        = CS.activeCellFaces;
   actF        = CS.activeFaces;
   renum       = zeros([size(CG.cells.faces,1), 1]);
   renum(actB) = 1 : numel(actB);

   % 2) Map (face,block) pairs to coarse half-faces, 'j' (assuming j > 0).
   %    Standard algorithm.
   %
   blkno   = rldecode(1 : CG.cells.num, diff(CG.cells.facePos), 2) .';
   colno   = @(f,b) 1 + (CG.faces.neighbors(f,1) ~= b);
   chf_map = accumarray([CG.cells.faces(:,1)      ,          ...
                         colno(CG.cells.faces(:,1), blkno)], ...
                        renum, [CG.faces.num, 2]);
   f2hf_c  = @(f,c) chf_map(sub2ind([CG.faces.num, 2], f, colno(f, c)));

   % 3) Apply hybrid basis function splitting to each (active) flux basis
   %    function.  This process forms the SPARSE input triplet [i,j,v] from
   %    which the hybrid flux basis matrix Bv is finally formed.
   %
   %    Each element of the cell array CS.basis(actF) is a 'V' tuple as
   %    defined by function 'evalBasisFunc'.  We also need to renumber the
   %    block-to-face connections, so pass the static 'f2hf_c' matrix along
   %    with the dynamic 'V' tuples.
   %
   [i, j, v] = cellfun(@(x) split_for_hybrid(x, f2hf_c), ...
                       CS.basis(actF), 'UniformOutput', false);
   Bv = sparse(vertcat(i{:}), vertcat(j{:}), vertcat(v{:}), ...
               size(G.cells.faces,1), numel(actB));

   % 4) Apply hybrid basis function splitting to each (active) pressure
   %    basis function.  This process forms the SPARSE input triplet
   %    [i,j,v] from which the hybrid pressure basis matrix Phi is finally
   %    formed.
   %
   %    Each element of the cell array CS.basisP(actF) is a 'P' tuple
   %    as defined by function 'evalBasisFunc'.
   %
   [i, j, v] = cellfun(@(x) split_for_hybrid(x, f2hf_c), ...
                       CS.basisP(actF), 'UniformOutput', false);
   Phi = sparse(vertcat(i{:}), vertcat(j{:}), vertcat(v{:}), ...
                G.cells.num, numel(actB));
end

%--------------------------------------------------------------------------

function [i, j, v] = split_for_hybrid(c, hfno)
   % Split a basis function tuple (either flux or pressure) 'c' into a
   % SPARSE matrix input triplet [i,j,v].  A hybrid basis function is split
   % into one component for each of (at most) two coarse half-faces.  We
   % recall that the input tuple, 'c', is a cell array of the form
   %
   %    1  2  3  4  5  6  7
   %   {i, s, f, b, n, m, o}
   %
   % with symbols as defined in 'evalBasisFunc'.  Items 6 and 7 are ignored
   % in this part of the implementation, while item 1 directly provides the
   % desired output vector 'i' (no further processing required).
   %
   % Our task, then, is to convert item 2 into the output vector 'v', and
   % to build output vector 'j' from scratch.  The output 'v' is the
   % simplest to define.  We simply negate the v-values (item 2)
   % corresponding to the second block (b(2)) if applicable.  If
   % numel(b)==1, then v==c{2} is returned unchanged.  The 'j' vector
   % contains at most two distinct values, each corresponding to a specific
   % coarse half-face identified by a pair (b,f) of coarse blocks and face.

   % 1) Extract 'i' values directly from input tuple.
   i = c{1};

   % 2) Build 'j' values from scratch by
   %    i   - Defining pairs of coarse blocks and face.
   b = reshape(double(c{4}), [], 1);
   f = repmat (double(c{3}), [numel(b), 1]);

   %    ii  - Extracting coarse half-faces identified by these pairs.
   j = reshape(hfno(f, b), [], 1);  assert (all(j > 0));

   %    iii - Expanding these half-face indices a number of times equal to
   %          the number of entities (item 5) in each of the coarse blocks
   %          identified by b.
   j = rldecode(j, reshape(c{5}, [], 1));

   % 3) Negate b(2) values if numel(b) > 1.  Else somewhat expensive NOP.
   s = [1; -1];
   v = rldecode(s(1 : numel(b)), reshape(c{5}, [], 1)) .* c{2};
end
