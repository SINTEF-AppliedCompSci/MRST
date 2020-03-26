classdef SurfactantComponent < GenericComponent
    properties

    end

    methods
        function c = SurfactantComponent()
            c@GenericComponent('surfactant');
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
             adsorbed = f.rhoR .* ((1-poro)./poro) .* ads;

             c{wIx} = pv.*(adsorbed + acc);
        end

        function cmob = getComponentMobility(component, model, state, varargin)
             % We use a Todd-Longstaff model. It implies that the mobility of the
             % polymer is a non-linear function of the polymer concentration.
             [mob, b, c] = model.getProps(state, 'Mobility', 'ShrinkageFactors', 'surfactant');
             wIx = 1;
             mobW = mob{wIx};
             bW = b{wIx};
             mobS = c.*bW.*mobW;

             nphase = model.getNumberOfPhases;
             cmob = cell(1, nphase);
             cmob{wIx} = mobS;
        end


        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            % Polymer enters into the water stream
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
            c{1} = 1;
        end

        % FIXME: the tricky part here is that it is not composition, it is
        % a concentration, we need to check how this value is used
        function c = getPhaseComponentFractionInjection(component, model, state, force)
            c = cell(model.getNumberOfPhases(), 1);
            if isfield(force, 'compi')
                comp_i = vertcat(force.compi);
            else
                comp_i = vertcat(force.sat);
            end
            wIx = 1;
            cs = vertcat(force.cs);
            ci = comp_i(:, wIx) .* cs./model.fluid.rhoWS;
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
