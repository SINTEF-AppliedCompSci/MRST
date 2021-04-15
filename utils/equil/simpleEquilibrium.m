function s = simpleEquilibrium(G, contacts, vector)
%Routine for creating simple initial equilibriums
%
% SYNOPSIS:
%   s = simpleEquilibrium(G, contacts)
%   s = simpleEquilibrium(G, contacts, vector)
%
% DESCRIPTION:
%   This routine generates cell-wise values that are distributed according
%   to a list of different contacts. The most common usage of this routine
%   is to create initial saturations where the problem is at hydrostatic
%   equilibrium, i.e. the initial phases are distributed according to
%   density.
%
%   As the underlying algorithm makes certain simplifications without
%   taking fluid physics into account, it should be used with care and only
%   for problems with relatively structured grids and no capillary
%   pressure.
%
% REQUIRED PARAMETERS:
%   G        - The grid structure for which the properties are to be
%            calculated. Note that we require the grid to have geometry
%            information computed by a call to computeGeometry before usage.
%
%   contacts - A list of up to N contacts. Each contact represents the
%            interface between two properties and N contacts must be given
%            for N + 1 properties. Each contact is given as the distance
%            along the vector (third argument) where the contact occurs. If
%            no third argument is given, the contacts will be assumed to be
%            along the z-coordinate, which is also the default gravity
%            direction in MRST.
%
%   vector   - Optional argument. Two modes are supported:
%            If given as a vector, this vector will be the direction that
%            the contacts are parametrized along. For instance, the default
%            [0, 0, 1] interprets the contacts along the z-direction of the
%            grid. If [1, 0, 1] is given, contacts are interpreted along
%            the 45 degree angle between the x and z axes.
% 
%            If a single value is given, it is interpreted as a index into
%            the cell centroids, i.e. vector = 2 will use the second column
%            of the cell coordinates.
%
%
% RETURNS:
%   s       - A G.cell.num x (N + 1) matrix of properties, distributed
%           according to the contacts. It will not necessarily be
%           guaranteed to be the true equilibrium of the grid is not
%           Cartesian or the fluids physics include capillary forces.
%
% EXAMPLE:
%   G = cartGrid([10, 1, 10], [1, 1, 1]);
%   G = computeGeometry(G);
% 
%   % First contact at .37, second at .8 for a total of three phases present
%   contacts = [.37, .8];
% 
%   % Compute equilibrium saturations
%   s = simpleEquilibrium(G, contacts);
%
% SEE ALSO:
%   `initEclipseState`, `test_simpleEquilibrium`

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

    if nargin < 3
        vector = [0, 0, 1];
    end
    if numel(vector) > 1
        vector = vector./norm(vector);
    else
        assert(mod(vector, 1) == 0, ...
            'Single input for third argument must be a integer!')
    end
    
    assert(all(diff(contacts) >= 0), ...
        'Contacts must be monotonically increasing!')
    assert(isfield(G.cells, 'centroids'), ...
        ['Geometry information in grid required. Consider calling', ...
        ' ''computeGeometry'' beforehand.']);
    nPh = numel(contacts) + 1;
    s = zeros(G.cells.num, nPh);
    
    cellNo = rldecode(1 : G.cells.num, ...
                    diff(G.cells.facePos), 2) .';
    fz = proj(G.faces.centroids(G.cells.faces(:, 1), :), vector);
    cz = proj(G.cells.centroids, vector);
    top = accumarray(cellNo, fz, [], @max);
    bottom = accumarray(cellNo, fz, [], @min);
    
    width = top - bottom;
    bad = width == 0;
    
    contacts = [-inf, contacts, inf];
    for i = 2:numel(contacts)
        c_prev = contacts(i-1);
        c_next = contacts(i);
        
        distTop = min(c_next, top);
        distBottom = max(bottom, c_prev);
        
        s(:, i-1) = max(distTop - distBottom, 0)./width;
        if any(bad)
            % Cells with zero thickness
            s(bad, i-1) = cz(bad) <= c_next & cz(bad) > c_prev;
        end
    end
    % Reverse to be consistent with MRST's convention of heavier fluids
    % coming first.
    s = s(:, end:-1:1);
end

function x = proj(x, vec)
    if numel(vec) == 1
        x = x(:, vec);
    else
        x = sum(x.*repmat(vec, size(x, 1), 1), 2);
    end
end