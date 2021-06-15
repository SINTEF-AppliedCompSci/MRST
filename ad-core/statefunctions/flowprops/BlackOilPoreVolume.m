classdef BlackOilPoreVolume < PoreVolume
    % Effective pore-volume after pressure-multiplier
    properties
    end
    
    methods
        function gp = BlackOilPoreVolume(model, varargin)
            gp@PoreVolume(model, varargin{:});
            gp = gp.dependsOn({'pressure'}, 'state');
            assert(isfield(model.fluid, 'pvMultR'), 'pvMultR missing from fluid.');
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for rock-compressibility
            pv = evaluateOnDomain@PoreVolume(prop, model, state);
            p = model.getProp(state, 'pressure');
            pvMult = prop.evaluateFluid(model, 'pvMultR', p);
            pv = pv.*pvMult;
        end
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
