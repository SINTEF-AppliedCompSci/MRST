function states = convertReservoirFluxesToSurface(model, states)
%Compute surface fluxes from reservoir fluxes
%
% SYNOPSIS:
%   states = convertReservoirFluxesToSurface(model, states)
%
% DESCRIPTION:
%   This function, given states with .bfactors and .flux, will compute the
%   surface/standard condition fluxes and place them under the field
%   surfaceFlux. To ensure bfactors and fluxes are added to states during a
%   simulation, enable the flag "model.extraStateOutput" before simulating.
%
% REQUIRED PARAMETERS:
%   model   - ReservoirModel subclass used to produce the states.
%
%   states  - States with valid fields .flux and .bfactors.
%
%
% RETURNS:
%   states  - States with additional field surfaceFlux.
%

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
    if numel(states) > 1
        assert(iscell(states));
    else
        states = {states};
    end
    
    if ~isfield(states{1}, 'bfactor')
        error(['Missing bfactors in state - were the states produced from',...
               ' a model with the ''''extraStateOutput'''' flag enabled']);
    end
    
    upstr = model.operators.faceUpstr;
    intx  = model.operators.internalConn;
    
    for i = 1:numel(states)
        state = states{i};
        
        [m, n] = size(state.flux);
        state.surfaceFlux = zeros(m, n);
        for j = 1:n
            up = state.upstreamFlag(:, j);
            b = state.bfactor(:, j);
            v = state.flux(intx, j);
            
            state.surfaceFlux(intx, j) = upstr(up, b).*v;
        end
        states{i} = state;
    end
end
