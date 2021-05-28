classdef SurfactantComponent < ConcentrationComponent
    properties

    end

    methods
        function c = SurfactantComponent()
            % this surfactant model is transport in water phase
            wIx = 1;
            c@ConcentrationComponent('surfactant', wIx);
            
            c = c.functionDependsOn('getComponentDensity', 'surfactant',        'state');
            c = c.functionDependsOn('getComponentDensity', 'ShrinkageFactors',  'PVTPropertyFunctions');

            c = c.functionDependsOn('getComponentMass', {'s', 'surfactant'},    'state');
            c = c.functionDependsOn('getComponentMass', 'ShrinkageFactors',     'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMass', 'PoreVolume',           'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMass', 'SurfactantAdsorption', 'FlowPropertyFunctions');
            
            c = c.functionDependsOn('getComponentMobility', 'surfactant',       'state');
            c = c.functionDependsOn('getComponentMobility', 'ShrinkageFactors', 'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMobility', 'Mobility',         'FlowPropertyFunctions');
        end

        function c = getComponentDensity(component, model, state, varargin)
            cs = model.getProp(state, 'surfactant');
            b = model.getProps(state, 'ShrinkageFactors');
            nph = numel(b);
            c = cell(1, nph);
            % water phase index
            wIx = 1;
            c{wIx} = cs .* b{wIx};
        end


        function c = getComponentMass(component, model, state, varargin)
             f = model.fluid;
             cs = model.getProp(state, 'surfactant');
             pv = model.getProp(state, 'PoreVolume');
             b = model.getProps(state, 'ShrinkageFactors');

             nph = model.getNumberOfPhases;
             c = cell(1, nph);

             wIx = 1;
             bW = b{wIx};
             sw = model.getProp(state, 'sW');
             % In mobile water
             acc = sw.*cs.*bW;
             % Adsorbed part
             poro = model.rock.poro;
             ads = model.getProp(state, 'SurfactantAdsorption');
             adsorbed = f.rhoRSft .* ((1-poro)./poro) .* ads;

             c{wIx} = pv.*(adsorbed + acc);
        end

        function cmob = getComponentMobility(component, model, state, varargin)
             % The mobility of surfactant is a nonlinear function of
             % surfactant concentration, which affects both the water
             % viscosity and the relative permeabilities (by changing the
             % interfacial tension between oil and water).
             [mob, b, c] = model.getProps(state, 'Mobility', 'ShrinkageFactors', 'surfactant');
             wIx = 1;
             mobW = mob{wIx};
             bW = b{wIx};
             mobS = c.*bW.*mobW;

             nphase = model.getNumberOfPhases;
             cmob = cell(1, nphase);
             cmob{wIx} = mobS;
        end

        function c = getInjectionMassFraction(component, model, force)
            c = vertcat(force.cs)./model.fluid.rhoWS;
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
