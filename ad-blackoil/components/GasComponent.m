classdef GasComponent < ImmiscibleComponent
    properties
        disgas
        vapoil
    end
    
    methods
        function c = GasComponent(name, gasIndex, disgas, vapoil)
            c@ImmiscibleComponent(name, gasIndex);
            c.disgas = disgas;
            c.vapoil = vapoil;
            c = c.dependsOn('ShrinkageFactors', 'PVTPropertyFunctions');
            if disgas
                c = c.dependsOn('rs', 'state');
            end
        end
        
        function c = getComponentDensity(component, model, state)
            c = getComponentDensity@ImmiscibleComponent(component, model, state);
            if component.disgas || component.vapoil
                phasenames = model.getPhaseNames();
                gix = phasenames == 'G';
                oix = phasenames == 'O';
                b = model.getProps(state, 'ShrinkageFactors');
                rhoS = model.getSurfaceDensities();
                if size(rhoS, 1) > 1
                    reg = model.FlowPropertyFunctions.ShrinkageFactors.regions;
                else
                    reg = 1;
                end
                rhoGS = rhoS(reg, gix);
                if component.vapoil
                    bG = b{gix};
                    c{gix} = rhoGS.*bG;
                end
                if component.disgas
                    bO = b{oix};
                    rs = model.getProp(state, 'rs');
                    c{oix} = rs.*rhoGS.*bO;
                end
            end
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
