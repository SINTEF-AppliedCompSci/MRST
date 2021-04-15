function varargout = gridLogicalIndices(G, varargin)
% Given grid G and optional subset of cells, find logical indices.
%
% SYNOPSIS:
%   ijk = gridLogicalIndices(G)
%   ijk = gridLogicalIndices(G, c)
%   [i,j,k] = gridLogicalIndices(G);
%
% PARAMETERS:
%   G    - Grid structure
%
% OPTIONAL PARAMETERS:
%   c    - cells for which logical indices are to be computed. Defaults to
%          all cells `1:G.cells.num`;
%
% RETURNS:
%  
%   ijk   - If one output paramter is requested: Cell array of size
%           `numel(G.cartDims)` where `ijk{1}` contains logical indices
%           corresponding to `G.cartDims(1)` and so on.
%
%  
%   i,j,k - If multiple outputs are requested: Corresponding to logical
%           indices for different `cartDims`.
%
% EXAMPLE:
%   G = cartGrid([2, 3]);
%   ijk = gridLogicalIndices(G, c);

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


    assert (isfield(G, 'cartDims'), ...
           ['Logical indices requires that the grid ', ...
            'has a valid ''cartDims'' field.']);

    mrstNargInCheck(1, 2, nargin);

    if nargin == 2 && isnumeric(varargin{1}),
        gcell = reshape(G.cells.indexMap(varargin{1}), [], 1);
    else
        gcell = reshape(G.cells.indexMap, [], 1);
    end

    [ijk{1:numel(G.cartDims)}] = ind2sub(G.cartDims, double(gcell));
    if nargout < 2
        varargout{1} = ijk;
    else
        % User requested multiple outputs, give flat arrays
        assert(nargout <= G.griddim, ...
               'Max number of output arguments exceeds grid dimension!');

        nout = nargout;
        [varargout{1:nout}] = ijk{1:nout};
    end
end
