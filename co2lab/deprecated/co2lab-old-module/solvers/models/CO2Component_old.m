classdef CO2Component_old < GenericComponent

    methods
        function c = CO2Component_old(name, disgas, gasIndex)
            c@GenericComponent(name, gasIndex);
            c = c.functionDependsOn('getComponentDensity', ...
                                    {'ShrinkageFactors', 'SurfaceDensity'}, ...
                                    'PVTPropertyFunctions');
            if disgas
                c = c.functionDependsOn('getComponentDensity', 'rs', 'state');
            end
        end
        
        function c = getPhaseComposition(component, model, state, varargin)
            % @@ This is modeled after the GasComponent in the black oil module.  But should 
            %    it not rather call the GenericComponent version, to account for presence in 
            %    multiple phases?
            c = getPhaseComposition@ImmiscibleComponent(component, model, state, varargin{:});
        end
        
        function c = getComponentDensity(component, model, state, varargin)
            if model.disgas
                % establish (empty) cell array with one entry per phase
                phasenames = model.getPhaseNames();
                nph = numel(phasenames);
                c = cell(nph, 1);
                
                gix = (phasenames == 'G');
                pvt = model.PVTPropertyFunctions;
                
                rho = pvt.get(model, state, 'Density', true);
                b = pvt.get(model, state, 'ShrinkageFactors', true);
                rhoS = pvt.get(model, state, 'SurfaceDensity', true);
                rhoGS = rhoS{gix};
                
                % density of CO2 in CO2 phase
                c{gix} = rho{gix};
                
                % density of CO2 in brine phase
                if model.disgas
                    wix = (phasenames == 'W');
                    bW = b{wix};
                    rs = model.getProp(state, 'rs');
                    c{wix} = rs .* rhoGS .* bW;
                end
            else
                % density of CO2 in CO2 phase (which is the only phase with CO2)
                c = getComponentDensity@ImmiscibleComponent(component, model, state, varargin{:});
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
