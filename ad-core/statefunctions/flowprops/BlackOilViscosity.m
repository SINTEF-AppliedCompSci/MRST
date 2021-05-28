classdef BlackOilViscosity < Viscosity
    % Black-oil style viscosity functions that account for rs and Rv
    properties
        useSaturatedFlag = true;
        disgas = false;
        vapoil = false;
    end
    
    methods
        function mu = BlackOilViscosity(model, varargin)
            mu@Viscosity(model, varargin{:});
            assert(isa(model, 'ThreePhaseBlackOilModel'), ...
                ['Model must be derived from the black-oil model. ', ...
                'Did you want the regular Viscosity class instead?'])
            mu.disgas = model.disgas;
            mu.vapoil = model.vapoil;
            if mu.disgas
                mu = mu.dependsOn('rs', 'state');
            end
            if mu.vapoil
                mu = mu.dependsOn('rv', 'state');
            end
            if mu.useSaturatedFlag
                mu = mu.dependsOn('s', 'state');
            end
        end

        function mu = evaluatePhaseViscosity(prop, model, state, name, p)
            if prop.disgas && strcmp(name, 'O')
                % Oileic phase with dissolved gas component
                rs = model.getProp(state, 'rs');
                if prop.useSaturatedFlag
                    sG = model.getProp(state, 'sg');
                    flag = sG > 0;
                else
                    flag = false(numelValue(rs), 1);
                end
                mu = prop.evaluateFluid(model, 'muO', p, rs, flag);
            elseif prop.vapoil && strcmp(name, 'G')
                % Gaseous phase with vaporized oil component
                rv = model.getProp(state, 'rv');
                if prop.useSaturatedFlag
                    sO = model.getProp(state, 'so');
                    flag = sO > 0;
                else
                    flag = false(numelValue(rv), 1);
                end
                mu = prop.evaluateFluid(model, 'muG', p, rv, flag);
            else
                % Can use base class directly
                mu = evaluatePhaseViscosity@Viscosity(prop, model, state, name, p);
            end
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
