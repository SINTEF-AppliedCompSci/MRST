function rock = makeRock(G, perm, poro, varargin)
%Construct a rock struct for given permeabilty and porosity
%
% SYNOPSIS:
%   rock = makeRock(1*milli*darcy, 0.5);
%   rock = makeRock(K, poro, 'ntg', nettogross);
%
% PARAMETERS:
%
%   G    - Valid grid structure. This is the grid the user wants to make a
%          rock struct for.
%
%   perm - Permeability field. This can either be:
%           - A scalar. Interpreted as uniform, scalar permeability.
%           - A single row of with (1)/2/3 in 2D columns and (1)/3/6 in 3D.
%           Interpreted as uniform anisotropic permeability.
%           - A matrix with column count as above, but with G.cells.num
%           rows in total. This will be interpreted as cell wise variables,
%           resulting in heterogeneous permeability.
%
%   poro - Porosity. Can be either a single, scalar value or a column
%          vector with one entry per cell. Non-positive values will result
%          in a warning.
%
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   ntg  - Net-to-gross factor. Either a single scalar value, or a column
%          vector with one entry per cell. NTG acts as a multiplicative
%          modifier for the porosity when pore volume is calculated.
%
% RETURNS:
%   rock - Valid rock with properties for each active cell in the grid.
%
% SEE ALSO:
%   computeTrans, poreVolume, permTensor

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

    opt = struct('ntg', []);
    opt = merge_options(opt, varargin{:});
    
    nc = G.cells.num;
    
    perm = expandToCell(perm, nc);
    nt = size(perm, 2);
    
    if G.griddim == 2
        assert(nt == 1 || nt == 2 || nt == 3,...
            ['Permeability must have 1/2/3 columns for scalar/diagonal/full', ...
            ' tensor respectively in 2D. You supplied', num2str(nt)]);
    else
        assert(nt == 1 || nt == 3 || nt == 6,...
            ['Permeability must have 1/3/6 columns for scalar/diagonal/full', ...
            ' tensor respectively in 3D. You supplied', num2str(nt)]);
    end
    
    poro = expandToCell(poro, nc);
    assert(size(poro, 2) == 1, 'Porosity must be single column');
    lowporo = poro <= 0;
    if any(lowporo)
        warning(['Zero or negative porosity found in cells: ', num2str(find(lowporo)')]);
    end
    rock = struct('perm', perm, 'poro', poro);
    
    if ~isempty(opt.ntg)
        rock.ntg = expandToCell(opt.ntg, nc);
        assert(size(rock.ntg, 2) == 1, 'Net-to-gross must be single column');
    end
end

function vals = expandToCell(vals, nc)
    if size(vals, 1) == 1
        vals = repmat(vals, nc, 1);
    end
    assert(size(vals, 1) == nc, ...
        ['Supplied values for rock must be either one value per ', ...
        'active cell, or one value for all cells.']);
end