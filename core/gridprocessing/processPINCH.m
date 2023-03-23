function nnc = processPINCH(grdecl, G)
%Establish vertical non-neighbouring across pinched-out layers
%
% SYNOPSIS:
%   nnc = processPINCH(grdecl, G)
%
% DESCRIPTION:
%   This function establishes non-neighbouring connections across
%   pinched-out layers that may nevertheless have a non-zero thickness.  In
%   other words, this function creates connections that do not otherwise
%   correspond to geometric interfaces.
%
%   This function is used within function `processGRDECL` to implement the
%   processing of keyword `PINCH` in an `ECLIPSE` input deck.  Note that we
%   currently do not support the complete feature set that may be input
%   through `PINCH`.  Specifically, we only support the `TOPBOT`
%   transmissibility option for item 4 of the keyword.  Any other setting
%   will be reset to `TOPBOT` and a warning will be issued.
%
% PARAMETERS:
%   grdecl - Raw pillar grid structure, as defined (e.g.,) by function
%            `readGRDECL`, with fields `COORDS`, `ZCORN` and, possibly,
%            `ACTNUM`.
%
%   G      - Grid structure--typically created by function `processGRDECL`.
%
% RETURNS:
%   nnc - Structure describing explicit, additional non-neighbouring
%         connections.  Contains the following fields.
%
%           cells - Active cells connected across an NNC.  An m-by-2 array
%           of active cell numbers in the format of `G.faces.neighbors`.
%           Each row of `nnc.cells` represents a single non-neighbouring
%           connection.
%
%           faces - Partially redundant connection information.  An m-by-2
%           array of interface numbers that represent those interfaces that
%           would otherwise be geometrically connected if the NNC were in
%           fact a geometric connection (single interface).
%
%         If the pillar grid contains transmissibility multipliers in the
%         field `MULTZ`, the `nnc` structure will furthermore contain a
%         field `mult` of multipliers defined by item `5` of the `PINCH`
%         keyword.  The `mult` field is an m-by-1 array of non-negative
%         scalars, the i'th of which is the multiplier of the i'th explicit
%         non-neighbouring connection (i.e., row `i` of `nnc.cells`).
%
% NOTE:
%   At present, we only support generating NNCs across explicitly
%   deactivated layers--i.e., layers/cells for which `ACTNUM==0`.  Therefore,
%   function `processPINCH` will fail unless the raw pillar grid structure
%   contains an `ACTNUM` field.
%
% SEE ALSO:
%   `processGRDECL`

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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


   assert (all(isfield(grdecl, {'ACTNUM', 'PINCH'})), ...
          ['Corner-point specification must contain ''ACTNUM'' and ', ...
           '''PINCH'' in order to use function ''%s''.'], mfilename);

   assert ((numel(grdecl.cartDims) == 3) && ...
           (numel(G.cartDims)      == 3) && ...
           all(G.cartDims(:) == grdecl.cartDims(:)), ...
           'Grid dimensions must match dimensions of CP specification.');

   assert ((size(G.cells.faces, 2)  >= 2) && ...
           (max(G.cells.faces(:,2)) >= 6), ...
          ['Grid must include one-sided tags to identify ', ...
           'cardinal directions.']);

   skip = grdecl.ACTNUM == 0;

   if ~ any(skip),
      % No (explicit) candidates for PINCH processing.
      nnc = [];
   else
      % Pinch connects across inactive layers assuming thickness is less
      % than PINCH(1)

      % Switch to column-based indexing to ease subsequent processing
      act = zeros([prod(G.cartDims), 1]);
      act(G.cells.indexMap) = 1 : G.cells.num;

      skip = column_shape(skip, G.cartDims);
      act  = column_shape(act,  G.cartDims);

      % CUMSUM(skip) captures the number of intervening deactivated layers.
      % Consequently, DIFF(x) > 0 identifies those connections that could
      % potentially cross pinched layers.
      x = cumsum(double(skip));
      t = [rldecode(1 : size(x,2), size(x,1), 2) .', x(:), act(:)];
      t = t(act ~= 0, :);
      p = cumsum([1 ; accumarray(t(:,1), 1)]);

      ilo = mcolon(p(1 : end - 1)    , p(2 : end) - 2);
      ihi = mcolon(p(1 : end - 1) + 1, p(2 : end) - 1);

      i   = t(ihi, 2) > t(ilo, 2);
      top = t(ilo(i), 3);  % Top cell
      bot = t(ihi(i), 3);  % Bottom cell

      nnc = identify_conn(G, top, bot);
   end

   if ~ isempty(nnc),
      thickness = calc_thickness(grdecl.ZCORN, G.cartDims, ...
                                 double(G.cells.indexMap), nnc.cells);

      i = ~ (thickness > grdecl.PINCH{1});

      if any(i),
         nnc = structfun(@(a) a(i,:), nnc, 'UniformOutput', false);
      else
         nnc = [];
      end
   end

   if ~ strcmp(grdecl.PINCH{4}, 'TOPBOT'),
      warning(msgid('Trans:Unsupported'), ...
             ['Only ''TOPBOT'' transmissibility option ', ...
              'supported in MRST (''%s'' ignored).'], grdecl.PINCH{4});

      grdecl.PINCH{4} = 'TOPBOT';
   end

   if ~isempty(nnc) && isfield(grdecl, 'MULTZ') && ...
         strcmp(grdecl.PINCH{4}, 'TOPBOT'),  % MULTZ is ignored for 'ALL'.

      nnc.mult = define_multiplier(grdecl, nnc, double(G.cells.indexMap));

      assert (numel(nnc.mult) == size(nnc.cells, 1), ...
             ['Internal error in ''%s''. Number of multipliers must ', ...
              'match number of explicit pinch NNCs.'], mfilename);
   end
end

%--------------------------------------------------------------------------

function q = column_shape(q, sz)
   q = permute(reshape(q, sz), [numel(sz), 1 : (numel(sz) - 1)]);
   q = reshape(q, sz(end), []);
end

%--------------------------------------------------------------------------

function nnc = identify_conn(G, top, bot)
   ctop   = false([G.cells.num, 1]);  ctop(top) = true;
   cbot   = false([G.cells.num, 1]);  cbot(bot) = true;
   rtop   = zeros([G.cells.num, 1]);  rtop(top) = 1 : numel(top);
   rbot   = zeros([G.cells.num, 1]);  rbot(bot) = 1 : numel(bot);
   cellno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';

   % Upper face is bottom face of upper cell
   i   = ctop(cellno) & (G.cells.faces(:,2) == 6);
   fup = zeros([numel(top), 1]);
   fup(rtop(cellno(i))) = G.cells.faces(i, 1);

   % Lower face is upper face of lower cell
   i   = cbot(cellno) & (G.cells.faces(:,2) == 5);
   flo = zeros([numel(bot), 1]);
   flo(rbot(cellno(i))) = G.cells.faces(i, 1);

   if any(fup == 0) || any(flo == 0),
      error('Each cell must have exactly one TOP and one BOTTOM contact.');
   end

   % Only consider NNCs that aren't already established during regular
   % processing.
   include = fup ~= flo;

   nnc = struct('faces', [fup(include), flo(include)], ...
                'cells', [top(include), bot(include)]);
end

%--------------------------------------------------------------------------

function thick = calc_thickness(zcorn, cdims, imap, cells)
   [t{1:3}] = ind2sub(cdims, imap(cells(:,1)));
   [b{1:3}] = ind2sub(cdims, imap(cells(:,2)));

   assert (all(t{1} == b{1}) && all(t{2} == b{2}), ...
           'Unexpected column mismatch between top and bottom cells');

   [i, j] = ndgrid(1 : -1 : 0);
   repeat = [ 1, numel(i) ];

   gapno  = repmat((1 : size(cells, 1)) .', repeat);

   kt = repmat(2*t{3} - 0, repeat);
   kb = repmat(2*b{3} - 1, repeat);
   i  = bsxfun(@minus, 2 * t{1}, reshape(i, 1, []));
   j  = bsxfun(@minus, 2 * t{2}, reshape(j, 1, []));

   zmin = zcorn(sub2ind(2 * cdims, i(:), j(:), kt(:)));
   zmax = zcorn(sub2ind(2 * cdims, i(:), j(:), kb(:)));

   % Pinched-out layer thickness in corner-point geometry is defined as the
   % maximum difference in corner-point depths between the active cells on
   % either side of the pinched-out layer (or layers), for each of the
   % (four) pillars bounding the column.
   %
   thick = accumarray(gapno(:), zmax - zmin, [], @max);
end

%--------------------------------------------------------------------------

function mult = define_multiplier(grdecl, nnc, imap)
   % Two cases:
   % 1) ('TOP'): Use MZ in TOP (active) cell
   % 2) ('ALL'): Use minimum MZ in TOP:(BOTTOM-1)

   if strcmp(grdecl.PINCH{5}, 'TOP'),
      % Case 1)
      % Trivial.  Simply extract the MZ value corresponding to the global
      % index of compressed cell 'nnc.cells(:,1)'.

      mult = reshape(grdecl.MULTZ(imap(nnc.cells(:,1))), [], 1);
   else
      % Case 2)
      % Identify global Cartesian indices of TOP:(BOTTOM - 1) and compute
      % the required minimum value.

      assert (strcmp(grdecl.PINCH{5}, 'ALL'), ...
              'MULTZ handling option in PINCH must be ''TOP'' or ''ALL''.');

      [t{1:3}] = ind2sub(double(grdecl.cartDims), imap(nnc.cells(:,1)));
      [b{1:3}] = ind2sub(double(grdecl.cartDims), imap(nnc.cells(:,2)));

      assert (all(b{3} > t{3} + 1), ...
             ['Internal error in ''%s'': Pinch NNCs appear to be ', ...
              'neighbouring connections?'], mfilename);

      n  = b{3} - t{3};  % TOP:(BOTTOM-1) inclusive.
      ij = rldecode([t{1}, t{2}], n);
      k  = mcolon(t{3}, b{3} - 1) .';
      ix = sub2ind(double(grdecl.cartDims), ij(:,1), ij(:,2), k);

      % Extract minimum MZ value in column.
      mult = accumarray(rldecode(1 : size(nnc.cells, 1), n, 2) .', ...
                        grdecl.MULTZ(ix), [size(nnc.cells, 1), 1], @min);
   end
end
