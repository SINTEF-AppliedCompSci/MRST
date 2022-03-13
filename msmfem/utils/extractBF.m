function bf = extractBF(basis, sz, cg, varargin)
%Form matrix of resevoir basis function values (mixed/hybrid)
%
% SYNOPSIS:
%   bf = extractBF(basis, m, CG)
%   bf = extractBF(basis, m, CG, 'pn1', 'pv1', ...)
%
% PARAMETERS:
%   basis   - Cell array of packed basis function values as defined by
%             function 'evalBasisFunc'.  May be either of the flux or
%             pressure basis functions (fields 'CS.basis' or 'CS.basisP').
%             The cell array is assumed to be restricted to only those
%             basis functions which will be entered into the resulting
%             matrix.
%
%   m       - Number of rows in the resulting (sparse) basis function value
%             matrix.  This number differs if we're constructing flux
%             values (m==number of fine-scale half-contacts) or pressure
%             values (m==number of fine-scale cells).
%
%   CG      - Coarse grid data structure as defined by function
%             'generateCoarseGrid'.  Only used when forming the 'mixed'
%             edition of the basis function value matrix.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Type -- Which kind of basis function value matrix to
%                         form.  String.  Default value: Type = 'hybrid'.
%                         The supported values are {'mixed', 'hybrid'}.
%
% RETURNS:
%   bf - Matrix, size m-by-NUMEL(basis), of respective basis function
%        values.  The values of the j'th basis function (i.e., basis{j}) is
%        stored in bf(:,j).
%
% NOTE:
%   This function only expands the packed storage format represented by the
%   cell array 'basis'.  This impacts the return values of 'extractBF' when
%   used to derived flux basis function values from the 'CS.basis' field.
%   As function 'evalBasisFunc' returns flux basis function values in the
%   form 'B*v', these are the values which are entered into the 'bf'
%   matrix.  It is the caller's responsibility to compute the required
%   matrix product, e.g., S.BI * bf, in order to derive the correct flux
%   basis function values.
%
% SEE ALSO:
%   `generateCoarseSystem`, `evalBasisFunc`, `generateCoarseGrid`.

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


   opt = struct('Type', 'hybrid');
   opt = merge_options(opt, varargin{:});

   if strcmpi(opt.Type, 'hybrid'),
      bf_type = @hybrid;
   else
      assert (strcmpi(opt.Type, 'mixed'));
      bf_type = @(c) mixed(c, cg.faces.neighbors);
   end

   [i, v] = cellfun(bf_type, basis, 'UniformOutput', false);

   n = cumsum([0; reshape(cellfun(@numel, i), [], 1)]);
   j = zeros([n(end), 1]);
   j(n(1:end-1) + 1) = 1;

   bf = sparse(vertcat(i{:}), cumsum(j), vertcat(v{:}), sz, numel(i));
end


function [i, v] = hybrid(c)
   i = c{1};

   s = [1; -1];
   v = rldecode(s(1 : numel(c{4})), reshape(c{5}, [], 1)) .* c{2};
end


function [i, v] = mixed(c, n)
   i = c{1};

   s = 2*double(n(c{3},1) == c{4}(1)) - 1;
   v = s .* c{2};
end
