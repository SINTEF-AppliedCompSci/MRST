function u = formatField(u, dim, type)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    switch type
      case {'displacement', 'u'}
        u = reshape(u, dim, [])';
      case {'stress'}        
        uu = reshape(u, dim*dim, [])';
        nc = size(uu, 1);
        switch dim
          case 2
            vdim = 3; % voigt dimension
            u = NaN(nc, vdim);
            u(:, 1) = uu(:, 1);
            u(:, 2) = uu(:, 4);
            u(:, 3) = 0.5*(uu(:, 2) + uu(:, 3));
          case 3
            vdim = 6; % voigt dimension
            u = NaN(nc, vdim);
            u(:, 1) = uu(:, 1);
            u(:, 2) = uu(:, 5);
            u(:, 3) = uu(:, 9);
            u(:, 4) = 0.5*(uu(:, 6) + uu(:, 8));
            u(:, 5) = 0.5*(uu(:, 3) + uu(:, 7));
            u(:, 6) = 0.5*(uu(:, 2) + uu(:, 4));
          otherwise
            error('dimension not valid (accepted : 2D or 3D)');
        end
      otherwise
        error('type not recognized');
    end
end
