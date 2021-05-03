classdef GasComponent < ImmiscibleComponent
    % A black-oil component description of the gas type, which modifies the
    % immiscible implementation to account for parts of the gas component
    % dissolving into the oileic phase (if disgas is present) and that the
    % gas density can be modified by vaporized oil (if vapoil is present)

    properties
        disgas
        vapoil
    end
    
    methods
        function c = GasComponent(name, gasIndex, disgas, vapoil)
            c@ImmiscibleComponent(name, gasIndex);
            c.disgas = disgas;
            c.vapoil = vapoil;
            c = c.functionDependsOn('getComponentDensity', 'ShrinkageFactors', 'PVTPropertyFunctions');
            if disgas
                c = c.functionDependsOn('getComponentDensity', 'rs', 'state');
            end
        end

        function c = getPhaseComposition(component, model, state, varargin)
            c = getPhaseComposition@ImmiscibleComponent(component, model, state, varargin{:});
        end

        function c = getComponentDensity(component, model, state, varargin)
            if component.disgas || component.vapoil % Check for black-oil behavior
                phasenames = model.getPhaseNames();
                nph = numel(phasenames);
                c = cell(nph, 1);
                pvt = model.PVTPropertyFunctions;
                gix = (phasenames == 'G');
                rho = pvt.get(model, state, 'Density', true);
                b = pvt.get(model, state, 'ShrinkageFactors', true);
                rhoS = pvt.get(model, state, 'SurfaceDensity', true);
                rhoGS = rhoS{gix};
                if component.vapoil % Component density is not phase density
                    bG = b{gix};
                    c{gix} = rhoGS.*bG;
                else
                    c{gix} = rho{gix};
                end
                if component.disgas % There is mass of gas in oileic phase
                    oix = (phasenames == 'O');
                    bO = b{oix};
                    rs = model.getProp(state, 'rs');
                    c{oix} = rs.*rhoGS.*bO;
                end
            else
                c = getComponentDensity@ImmiscibleComponent(component, model, state, varargin{:});
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
