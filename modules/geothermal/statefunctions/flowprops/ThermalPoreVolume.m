classdef ThermalPoreVolume < StateFunction
%State function for pressure- and temperature dependent pore volume
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = ThermalPoreVolume(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isfield(model.fluid, 'pvMultR')
                gp = gp.dependsOn('Temperature');
                gp = gp.dependsOn({'pressure'}, 'state');
            end
            gp.label = '\phi';
        end
        
        %-----------------------------------------------------------------%
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for possible
            % rock-compressibility-dilatation as a function of p and T
            f = model.fluid;
            pv = model.operators.pv;
            if isfield(f, 'pvMultR')
                [p, T] = model.getProps(state, 'pressure', 'Temperature');
                pvMult = prop.evaluateFluid(model, 'pvMultR', p, T);
                pv     = pv.*pvMult;
            end
        end
    end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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