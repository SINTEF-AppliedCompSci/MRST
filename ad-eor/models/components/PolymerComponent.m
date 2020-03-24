classdef PolymerComponent < GenericComponent
    properties

    end
    
    methods
        function c = PolymerComponent()
            c@GenericComponent('polymer');
        end
        
        
        function c = getComponentDensity(component, model, state, varargin)
            cp = model.getProp(state, 'polymer');
            rho = model.getProps(state, 'Density');
            nph = numel(rho);
            c = cell(1, nph);
            c{1} = cp.*rho{1};
        end
        
        function c = getComponentMass(component, model, state, varargin)
            f = model.fluid;
            cp = model.getProp(state, 'polymer');
            rho = model.getProps(state, 'Density');
            pv = model.getProp(state, 'PoreVolume');
            
            nph = numel(rho);
            c = cell(1, nph);
            
            rhoW = rho{1};
            sw = model.getProp(state, 'sW');
            % In mobile water 
            acc = (1-f.dps).*sw.*cp.*rhoW;
            % Adsorbed part
            poro = model.rock.poro;
            ads = model.getProp(state, 'PolymerAdsorption');
            adsorbed = f.rhoWS.*f.rhoR.*((1-poro)./poro).*ads;
            
            c{1} = pv.*(adsorbed + acc);
        end
        
        function cmob = getComponentMobility(component, model, state, varargin)
        % We use a Todd-Longstaff model. It implies that the mobility of the
        % polymer is a non-linear function of the polymer concentration.
            mass = component.getComponentDensity(model, state, varargin{:});
            [mob, rho, c] = model.getProps(state, 'Mobility', 'Density', 'polymer');
            wIx = 1;
            mobW = mob{wIx};
            rhoW = rho{wIx};
            fluid  = model.fluid;
            mixpar = fluid.mixPar;
            cpbar  = c/fluid.cpmax;
            a  = fluid.muWMult(fluid.cpmax).^(1-mixpar);
            mobP = c.*rhoW.*mobW./(a+(1-a)*cpbar);
            
            nphase = numel(mass);
            cmob = cell(1, nphase);
            cmob{wIx} = mobP;
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            % Polymer does not enter into any phase stream
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
        end

        function c = getPhaseComponentFractionInjection(component, model, state, force)
            c = cell(model.getNumberOfPhases(), 1);
            if isfield(force, 'compi')
                comp_i = vertcat(force.compi);
            else
                comp_i = vertcat(force.sat);
            end
            wIx = 1;
            cp = vertcat(force.cp);
            ci = comp_i(:, wIx).*cp;
            if any(ci ~= 0)
                c{wIx} = ci;
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
