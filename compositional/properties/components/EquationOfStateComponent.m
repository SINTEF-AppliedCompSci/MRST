classdef EquationOfStateComponent < GenericComponent
    properties
        componentIndex % Global component numbering
        surfacePhaseMassFractions % Mass fraction for each phase
        surfacePhaseDensityPure % Density
        pressure % pressure
        T % Temp
    end
    
    methods
        function c = EquationOfStateComponent(name, p, T, cindex, surfaceMassFractions, density, mw)
            c@GenericComponent(name);
            c.componentIndex = cindex;
            c.pressure = p;
            c.T = T;
            c.surfacePhaseMassFractions = surfaceMassFractions;
            c.surfacePhaseDensityPure = density;
            c.molarMass = mw;
            c = c.functionDependsOn('getPhaseComposition', {'ComponentPhaseMassFractions'}, 'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentDensity', {'Density','ComponentPhaseMassFractions'}, 'PVTPropertyFunctions');
        end
        
        function c = getComponentDensity(component, model, state, extra)
            c = component.getPhaseComposition(model, state);
            if nargin < 4
                rho = model.getProps(state, 'Density');
            else
                rho = extra.rho;
            end
            for ph = 1:numel(c)
                if ~isempty(c{ph})
                    c{ph} = rho{ph}.*c{ph};
                end
            end
        end
        
        function c = getPhaseComposition(component, model, state, varargin)
            massFractions = model.getProps(state, 'ComponentPhaseMassFractions');
            ix = component.componentIndex;
            nph = size(massFractions, 2);
            c = cell(1, nph);
            for ph = 1:nph
                mf = massFractions{ix, ph};
                if ~isempty(mf)
                    c{ph} = mf;
                end
            end
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            for i = 1:nph
                c{i} = component.surfacePhaseMassFractions(i);
            end
        end

        function c = getPurePhaseDensitySurface(component, model, state, pressure, temperature)
            % Surface density, for a pure component
            rho = component.surfacePhaseDensityPure;
            % Ideal gas scaling
            scale = (pressure./component.pressure).*(component.T./temperature);
            rho = bsxfun(@times, rho, scale);
            c = arrayfun(@(x) x, rho, 'UniformOutput', false);
        end

        function c = getPhaseComponentFractionInjection(component, model, state, force)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            comp_i = vertcat(force.components);
            comp_i = model.EOSModel.getMassFraction(comp_i);
            index = component.componentIndex;
            ci = comp_i(:, index);
            if any(ci ~= 0)
                c{model.getLiquidIndex()} = ci;
                c{model.getVaporIndex()} = ci;
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
