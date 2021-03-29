function rock = makeShaleRock(G, perm, poro, varargin)
%Create rock structure from given permeabilty and porosity values
%
% SYNOPSIS:
%   rock = makeRock(G, perm, poro)
%   rock = makeRock(G, perm, poro, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid for which to construct a rock property structure.
%
%   perm    - Permeability field.  Supported input is:
%               - A scalar.  Interpreted as uniform, scalar (i.e.,
%                 homogenous and isotropic) permeability field repeated for
%                 all active cells in the grid G.
%
%               - A row vector of 1/2/3 columns in two space dimensions or
%                 1/3/6 in columns three space dimensions.  The row vector
%                 will be repetead for each active cell in the grid G and
%                 therefore interpreted as a uniform (i.e., a homogeneous),
%                 possibly anisotropic, permeability field.
%
%               - A matrix with column count as above, but with
%                 `G.cells.num` rows in total.  This input will be treated
%                 as per-cell values, resulting in heterogeneous
%                 permeability.
%
%   poro    - Porosity field.  Can be either a single, scalar value or a
%             column vector with one entry per cell.  Non-positive values
%             will result in a warning.
%
% OPTIONAL PARAMETERS:
%   'ntg' - Net-to-gross factor.  Either a single scalar value that is
%           repeated for all active cells, or a column vector with one
%           entry per cell.  `ntg` acts as a multiplicative factor on
%           porosity when calculating pore volumes. Typically in the range
%           [0 .. 1].    
%
% RETURNS:
%   rock - Valid rock with properties for each active cell in the grid.
%
% EXAMPLE:
%   G = computeGeometry(cartGrid([10, 20, 30], [1, 1, 1]))
%
%   r1 = makeRock(G, 100*milli*darcy, 0.3)
%   r2 = makeRock(G, [ 100, 100, 10 ]*milli*darcy, 0.3)
%
%   K   = logNormLayers(G.cartDims, repmat([400, 0.1, 20], [1, 2]));
%   phi = gaussianField(G.cartDims, [0.2, 0.4]);
%   ntg = rand([G.cells.num, 1]);
%   r3  = makeRock(G, convertFrom(K(:), milli*darcy), phi(:), 'ntg', ntg)
%
%   plotCellData(G, poreVolume(G, r3)), view(3), axis tight
%
% SEE ALSO:
%   `computeTrans`, `poreVolume`, `permTensor`.

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('ntg', []);
    opt = merge_options(opt, varargin{:});
    
    nc = G.cells.num;
    
    perm = expandToCell(perm, nc);
    nt = size(perm, 2);
    
    if G.griddim == 1
        assert (nt == 1, ...
               ['Permeability must have 1 column for ', ...
                'scalar values in 1D.\n', ...
                'You supplied %d components.'], nt);
    elseif G.griddim == 2
        assert (sum(nt == [1, 2, 3]) == 1, ...
               ['Permeability must have 1/2/3 columns for ', ...
                'scalar/diagonal/full tensor respectively in 2D.\n', ...
                'You supplied %d components.'], nt);
    else
        assert (sum(nt == [1, 3, 6]) == 1, ...
               ['Permeability must have 1/3/6 columns for ', ...
                'scalar/diagonal/full tensor respectively in 3D.\n', ...
                'You supplied %d components.'], nt);
    end
    
    poro = expandToCell(poro, nc);
    assert(size(poro, 2) == 1, 'Porosity must be single column');
    lowporo = poro <= 0;
    if any(lowporo)
        warning(['Zero or negative porosity found in cells: ', ...
                num2str(find(lowporo)')]);
    end
    rock = struct('perm', perm, 'poro', poro);
    
    if ~isempty(opt.ntg)
        rock.ntg = expandToCell(opt.ntg, nc);
        assert(size(rock.ntg, 2) == 1, ...
               'Net-to-gross must be single column');
    end

    % To indentify the matrix cell with value 1
    rock.isMatrix = expandToCell(1, nc);


end

%--------------------------------------------------------------------------

function vals = expandToCell(vals, nc)
    if size(vals, 1) == 1
        vals = repmat(vals, nc, 1);
    end
    assert(size(vals, 1) == nc, ...
        ['Supplied values for rock must be either one value per ', ...
        'active cell, or one value for all cells.']);
end
