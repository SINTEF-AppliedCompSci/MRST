classdef OilComponent < ImmiscibleComponent
    % A black-oil component description of the oil type, which modifies the
    % immiscible implementation to account for parts of the oil component
    % vaporizing into the gaseous phase (if vapoil is present) and that the
    % oil density can be modified by dissolved gas (if disgas is present)
    properties
        disgas
        vapoil
    end
    
    methods
        function c = OilComponent(name, gasIndex, disgas, vapoil)
            c@ImmiscibleComponent(name, gasIndex);
            c.disgas = disgas;
            c.vapoil = vapoil;
            c = c.functionDependsOn('getComponentDensity', ...
                                    {'ShrinkageFactors', 'SurfaceDensity'},...
                                    'PVTPropertyFunctions');
            if vapoil
                c = c.functionDependsOn('getComponentDensity', 'rv', 'state');
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
                rhoS = pvt.get(model, state, 'SurfaceDensity', true);
                rho = pvt.get(model, state, 'Density', true);
                b = pvt.get(model, state, 'ShrinkageFactors', true);
                oix = (phasenames == 'O');
                rhoOS = rhoS{oix};
                if component.disgas % Component density is not phase density
                    bO = b{oix};
                    c{oix} = rhoOS.*bO;
                else
                    c{oix} = rho{oix};
                end
                if component.vapoil % There is mass of oil in gaseous phase
                    gix = (phasenames == 'G');
                    bG = b{gix};
                    rv = model.getProp(state, 'rv');
                    c{gix} = rv.*rhoOS.*bG;
                end
            else
                c = getComponentDensity@ImmiscibleComponent(component, model, state, varargin);
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
