function [flag_v, flag_g] = getSaturationUpwind(upwtype, state, G, vT, T, K, upstr)
% Get upwind flags for viscous and gravity parts of flux function.
    switch lower(upwtype)
        case 'potential'
            % Standard PPU upwind
            flag = multiphaseUpwindIndices(G, vT, T, K, upstr);
            [flag_v, flag_g] = deal(flag);
        case {'hybrid', 'hybrid_combined'}
            % Hybrid upwinding
            [flag_v, flag_g] = hybridUpwind(G, vT, T, K, upstr);
        case 'static'
            % Use whatever flags are kept in state, likely from the
            % pressure solver. Debug only.
            [flag_v, flag_g] = deal(state.upstreamFlag);
        otherwise
            error(['Unknown upwind type: ''', model.upwindType, '''']);
    end
end

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
