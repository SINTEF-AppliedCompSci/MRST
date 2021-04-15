function grdecl_new = cutGrdecl(grdecl, ind, varargin)
%Extract logically Cartesian subset of a corner-point description
%
% SYNOPSIS:
%   grdecl = cutGrdecl(grdecl, ind)
%   grdecl = cutGrdecl(grdecl, ind, 'pn1', pv1, ...)
%
% PARAMETERS:
%   grdecl  - Corner-point description as defined by, e.g., function
%             'readGRDECL'.
%
%   ind     - Upper and lower bounds, represented as a 3-by-2 integer array
%             with the first column being lower bounds and the second
%             column being upper bounds, on the logically Cartesian subset.
%             Interpreted in tensor product reference of uncompressed cells.
%
% OPTIONAL ARGUMENTS:
%   'lefthanded_numbering' - Whether or not the corner-point description
%                            assumes a left-handed coordinate system. 
%                            Default value: lefthanded_numbering = false
%                            (assume right-handed coordinate system).
%
% RETURNS:
%   grdecl - Corner-point description corresponding to the uncompressed
%            subset of cell indices represented by 'ind'.
%
% SEE ALSO:
%   `readGRDECL`, `simpleGrdecl`.

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

   % All corner-point descriptions *must* have these fields.
   required = { 'cartDims', 'COORD', 'ZCORN' };

   assert (isstruct(grdecl) && all(isfield(grdecl, required)), ...
          ['Input parameter ''grdecl'' must be a valid corner-point ', ...
           'description as defined, e.g., by function ''readGRDECL''.']);

   opt = struct('lefthanded_numbering', false);
   opt = merge_options(opt, varargin{:});

   cellfelt = {'PERMX', 'PERMY', 'PERMZ', 'PORO', 'ACTNUM', 'SATNUM'};
   cellfelt = reshape(cellfelt(isfield(grdecl, cellfelt)), 1, []);

   nfld       = numel(required) + numel(cellfelt);
   grdecl_new = cell2struct(repmat({ [] }, [1, nfld]), ...
                            [ required, cellfelt ], 2);

   grdecl_new.cartDims = reshape(diff(ind, [], 2) + 1, 1, []);

   extract = subset(grdecl.cartDims, ind);
   for f = cellfelt,
      grdecl_new.(f{1}) = extract(grdecl.(f{1}));
   end

   xyz     = reshape(grdecl.COORD, [6, grdecl.cartDims(1:2) + 1]);
   xyz_new = xyz(         :               , ...
                 ind(1,1) : (ind(1,2) + 1), ...
                 ind(2,1) : (ind(2,2) + 1));

   if opt.lefthanded_numbering,
      xyz_new([2, 5], :, :) = - xyz_new([2, 5], :, :);
   end

   grdecl_new.COORD = reshape(xyz_new, [], 1);

   extract = subset(2 * grdecl.cartDims, ...
                    bsxfun(@plus, 2 * (ind - 1), [1, 2]));

   grdecl_new.ZCORN = extract(grdecl.ZCORN);
end

%--------------------------------------------------------------------------

function e = subset(cartDims, ind)
   assert (numel(cartDims) == size(ind, 1), ...
          ['Index subset must have same number of dimensions ', ...
           '(rows) as there are elements in Cartesian size array.']);

   ix = arrayfun(@colon, ind(:,1), ind(:,2), 'UniformOutput', false);

   [ijk{1:numel(ix)}] = ndgrid(ix{:});
   ijk = cellfun(@(v) reshape(v, [], 1), ijk, 'UniformOutput', false);

   i = sub2ind(reshape(cartDims, 1, []), ijk{:});

   e = @(a) reshape(a(i), [], 1);
end
