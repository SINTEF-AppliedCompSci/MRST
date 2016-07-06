function [flag_v, flag_g] = hybridUpwind(Gi, vT)
% Hybrid upwinding - two-phase only at the moment
    nPh = numel(Gi);
    assert(nPh == 2, ...
        'Hybrid upwinding currently only supported for two-phase problems.');
    
    % Total velocity for viscous flow, density sorting for buoyancy terms
    flag_v = repmat(vT > 0, 1, nPh);
    dp = Gi{2} - Gi{1};
    d_diff = dp <= 0;
    flag_g = [d_diff, ~d_diff];
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
