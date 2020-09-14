classdef GenericOverallCompositionModel < OverallCompositionCompositionalModel & GenericReservoirModel
    properties
        
    end
    
    methods
        function model = GenericOverallCompositionModel(varargin)
            model = model@OverallCompositionCompositionalModel(varargin{:});
            model.OutputStateFunctions = {'ComponentTotalMass', 'Density'};
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@ReservoirModel(model, state0, state, dt, drivingForces, varargin{:});
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            % Discretize
            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Assemble equations and add in sources
            [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ComponentPhaseMassFractions');
            comps = cellfun(@(x, y) {x, y}, X(:, 1+model.water), X(:, 2+model.water), 'UniformOutput', false);
            
            
            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                             pressures, sat, mob, rho, ...
                                                             {}, comps, ...
                                                             drivingForces);
            
            % Add sources
            eqs = model.insertSources(eqs, src);
            % Assemble equations
            for i = 1:numel(eqs)
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            % Finally assemble
            eqs = [eqs, weqs];
            names = [names, wnames];
            types = [types, wtypes];
        end
        
        function names = getComponentNames(model)
            names = cellfun(@(x) x.name, model.Components, 'UniformOutput', false);
        end

        function [state, report] = updateState(model, state, problem, dx, forces)
            if isfield(state, 'FractionalDerivatives')
                state = rmfield(state, 'FractionalDerivatives');
            end
            [state, report] = updateState@OverallCompositionCompositionalModel(model, state, problem, dx, forces);
            if ~isempty(model.FacilityModel)
                state = model.FacilityModel.applyWellLimits(state);
            end
        end
        
        function model = validateModel(model, varargin)
            % Validate model.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.validateModel`
            if isempty(model.FacilityModel) || ~isa(model.FacilityModel, 'GenericFacilityModel')
                model.FacilityModel = GenericFacilityModel(model);
            end
            if isempty(model.Components)
                names_hc = model.EOSModel.CompositionalMixture.names;
                if model.water
                    names = [names_hc, 'water']; % Put water last
                else
                    names = names_hc;
                end
                nc = numel(names);
                model.Components = cell(1, nc);
                p = model.FacilityModel.pressure;
                T = model.FacilityModel.T;
                for ci = 1:nc
                    name = names{ci};
                    switch name
                        case 'water'
                            c = ImmiscibleComponent('water', 1);
                        otherwise
                            c = getEOSComponent(model, p, T, name, ci);
                    end
                    model.Components{ci} = c;
                end
            end
            model = validateModel@OverallCompositionCompositionalModel(model, varargin{:});
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@OverallCompositionCompositionalModel(model, varargin{:});
            model.FlowDiscretization.GravityPotentialDifference.saturationWeighting = true;
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if model.outputFluxes
                f = model.getProp(state, 'PhaseFlux');
                nph = numel(f);
                state.flux = zeros(model.G.faces.num, nph);
                state.flux(model.operators.internalConn, :) = [f{:}];
            end
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            % Get primary variables from state, before a possible
            % initialization as AD.'
            [p, sW, z] = model.getProps(state, 'pressure', 'sW', 'z');
            z_tol = model.EOSModel.minimumComposition;
            z = ensureMinimumFraction(z, z_tol);
            z = expandMatrixToCell(z);
            cnames = model.EOSModel.getComponentNames();
            names = [{'pressure'}, cnames(2:end)];
            vars = [p, z(2:end)];
            if model.water
                names = [names, {'sW'}];
                vars = [vars, {sW}];
            end
            origin = cell(1, numel(names));
            [origin{:}] = deal(class(model));
            if ~isempty(model.FacilityModel)
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
            end
        end

        function state = initStateAD(model, state, vars, names, origin)
            isP = strcmp(names, 'pressure');
            isAD = any(cellfun(@(x) isa(x, 'ADI'), vars));
            state = model.setProp(state, 'pressure', vars{isP});
            removed = isP;
            
            cnames = model.EOSModel.getComponentNames();
            ncomp = numel(cnames);
            z = cell(1, ncomp);
            z_end = 1;
            for i = 1:ncomp
                name = cnames{i};
                sub = strcmp(names, name);
                if any(sub)
                    z{i} = vars{sub};
                    z_end = z_end - z{i};
                    removed(sub) = true;
                else
                    fill = i;
                end
            end
            z{fill} = z_end;
            state = model.setProp(state, 'components', z);
            if isAD
                [state.x, state.y, state.L, state.FractionalDerivatives] = ...
                    model.EOSModel.getPhaseFractionAsADI(state, state.pressure, state.T, state.components);
            end
            if ~isempty(model.FacilityModel)
                % Select facility model variables and pass them off to attached
                % class.
                fm = class(model.FacilityModel);
                isF = strcmp(origin, fm);
                state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                removed = removed | isF;
            end
            if model.water
                isWater = strcmp(names, 'sW');
                sW = vars{isWater};
                removed(isWater) = true;
                offset = 1;
            else
                offset = 0;
            end
            % Set up state with remaining variables
            state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
            
            % Now that props have been set up, we can compute the
            % saturations from the mole fractions.
            if isAD
                % We must get the version with derivatives
                Z = model.getProps(state, 'PhaseCompressibilityFactors');
                Z_L = Z{offset+1};
                Z_V = Z{offset+2};
            else
                % Already stored in state - no derivatives needed
                Z_L = state.Z_L;
                Z_V = state.Z_V;
            end
            
            L = state.L;
            volL = L.*Z_L;
            volV = (1-L).*Z_V;
            volT = volL + volV;
            sL = volL./volT;
            sV = volV./volT;
            
            if model.water
                [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
                void = 1 - sW;
                sL = sL.*void;
                sV = sV.*void;
                [sL, sV] = model.setMinimumTwoPhaseSaturations(state, sW, sL, sV, pureLiquid, pureVapor, twoPhase);

                s = {sW, sL, sV};
            else
                s = {sL, sV};
            end
            state = model.setProp(state, 's', s);
        end

        function forces = validateDrivingForces(model, forces, varargin)
            forces = validateDrivingForces@OverallCompositionCompositionalModel(model, forces, varargin{:});
            forces = validateCompositionalForces(model, forces, varargin{:});
        end
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
