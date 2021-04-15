classdef PolymerComponent < ConcentrationComponent
    properties

    end

    methods
        function c = PolymerComponent()
            % this polymer model is transport in the water phase
            wIx = 1;
            c@ConcentrationComponent('polymer', wIx);
            
            c = c.functionDependsOn('getComponentDensity', 'polymer',             'state');
            c = c.functionDependsOn('getComponentDensity', 'ShrinkageFactors',    'PVTPropertyFunctions');

            c = c.functionDependsOn('getComponentMass', {'s', 'polymer'},         'state');
            c = c.functionDependsOn('getComponentMass', 'ShrinkageFactors',       'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMass', 'PoreVolume',             'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMass', 'PolymerAdsorption',      'FlowPropertyFunctions');
            
            c = c.functionDependsOn('getComponentMobility', 'polymer',            'state');
            c = c.functionDependsOn('getComponentMobility', 'ShrinkageFactors',   'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMobility', 'Mobility',           'FlowPropertyFunctions');
            c = c.functionDependsOn('getComponentMobility', 'PolymerEffViscMult', 'PVTPropertyFunctions');
            c = c.functionDependsOn('getComponentMobility', 'PolymerViscMult',    'PVTPropertyFunctions');
        end


        function c = getComponentDensity(component, model, state, varargin)
            cp = model.getProp(state, 'polymer');
            b = model.getProps(state, 'ShrinkageFactors');
            nph = numel(b);
            c = cell(1, nph);
            % water phase index
            c{1} = cp .* b{1};
        end

        function c = getComponentMass(component, model, state, varargin)
             f = model.fluid;
             cp = model.getProp(state, 'polymer');
             pv = model.getProp(state, 'PoreVolume');
             b = model.getProps(state, 'ShrinkageFactors');

             nph = model.getNumberOfPhases;
             c = cell(1, nph);

             wIx = 1;
             bW = b{wIx};
             sw = model.getProp(state, 'sW');
             % In mobile water
             acc = (1-f.dps).*sw.*cp.*bW;
             % Adsorbed part
             poro = model.rock.poro;
             ads = model.getProp(state, 'PolymerAdsorption');
             adsorbed = f.rhoR .* ((1-poro)./poro) .* ads;

             c{wIx} = pv.*(adsorbed + acc);
        end

        function cmob = getComponentMobility(component, model, state, varargin)
             % We use a Todd-Longstaff model. It implies that the mobility of the
             % polymer is a non-linear function of the polymer concentration.
             [mob, b, c, effviscmult, pviscmult] = ...
                 model.getProps(state, 'Mobility', 'ShrinkageFactors', ...
                     'polymer', 'PolymerEffViscMult', 'PolymerViscMult');
             wIx = 1;
             mobW = mob{wIx};
             bW = b{wIx};
             mobP = c.*bW.*mobW .* effviscmult ./ pviscmult;

             nphase = model.getNumberOfPhases;
             cmob = cell(1, nphase);
             cmob{wIx} = mobP;
         end

        function c = getInjectionMassFraction(component, model, force)
            c = vertcat(force.cp)./model.fluid.rhoWS;
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
