classdef GenericBlackOilModel < ThreePhaseBlackOilModel & GenericReservoirModel
    properties
    end

    methods
        function model = GenericBlackOilModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});
            model.OutputStateFunctions = {'ComponentTotalMass'};
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            % Get primary variables from state, before a possible
            % initialization as AD.
            phases = model.getPhaseNames();
            nph = numel(phases);
            if model.oil
                ix = find(phases == 'O');
            else
                ix = nph;
            end
            phases(ix) = [];

            snames = arrayfun(@(x) ['s', x], phases, 'UniformOutput', false);
            s = cell(1, nph-1);
            [p, s{:}, rs, rv] = model.getProps(state, ...
                'pressure', snames{:}, 'rs', 'rv');

            if model.disgas || model.vapoil
                % In this case, gas saturation is replaced with rs/rv in cells
                % where free gas is not present
                assert(model.oil, 'Cannot have disgas/vapoil without oil phase.');
                assert(model.gas, 'Cannot have disgas/vapoil without gas phase.');
                % X is either Rs, Rv or Sg, depending on each cell's saturation status
                if model.water
                    sW = s{phases == 'W'};
                else
                    sW = 0;
                end
                isG = phases == 'G';
                sG = s{isG};
                st  = model.getCellStatusVO(state,  1-sW-sG, sW, sG);
                x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
                s{isG} = x;
                snames{isG} = 'x';
            end
            vars = [p, s];
            names = ['pressure', snames];
            origin = cell(1, numel(vars));
            [origin{:}] = deal(class(model));

            if not(isempty(model.FacilityModel))
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
            end
        end

        function state = initStateAD(model, state, vars, names, origin)
            removed = false(size(vars));
            if model.disgas || model.vapoil
                % Black-oil specific variable switching
                if model.water
                    isw = strcmpi(names, 'sw');
                    sW = vars{isw};
                    removed = removed | isw;
                else
                    sW = 0;
                end

                isx = strcmpi(names, 'x');
                x = vars{isx};
                sG = model.getProps(state, 'sg');
                st  = model.getCellStatusVO(state, 1-sW-sG, sW, sG);
                sG = st{2}.*(1-sW) + st{3}.*x;
                sO = st{1}.*(1-sW) + ~st{1}.*(1 - sW - sG);
                if model.water
                    sat = {sW, sO, sG};
                else
                    sat = {sO, sG};
                end
                removed(isx) = true;
            else
                % Without variable switching
                phases = model.getPhaseNames();
                nph = numel(phases);
                sat = cell(1, nph);
                fill = ones(model.G.cells.num, 1);
                removed_sat = false(1, nph);
                for i = 1:numel(phases)
                    sub = strcmpi(names, ['s', phases(i)]);
                    if any(sub)
                        fill = fill - vars{sub};
                        removed = removed | sub;
                        removed_sat(i) = true;
                        sat{i} = vars{sub};
                    end
                end
                if any(~removed_sat)
                    sat{~removed_sat} = fill;
                end
            end
            state = model.setProp(state, 's', sat);

            if not(isempty(model.FacilityModel))
                % Select facility model variables and pass them off to attached
                % class.
                fm = class(model.FacilityModel);
                isF = strcmp(origin, fm);
                state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                removed = removed | isF;
            end

            % Set up state with remaining variables
            state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
            % Account for dissolution changing variables
            if model.disgas
                rsSat = model.getProp(state, 'RsMax');
                rs = ~st{1}.*rsSat + st{1}.*x;
                % rs = rs.*(value(sO) > 0);
                state = model.setProp(state, 'rs', rs);
            end

            if model.vapoil
                rvSat = model.getProp(state, 'RvMax');
                rv = ~st{2}.*rvSat + st{2}.*x;
                % rv = rv.*(value(sG) > 0);
                state = model.setProp(state, 'rv', rv);
                % No rv, no so -> zero on diagonal in matrix
                bad_oil = value(sO) == 0 & value(rv) == 0;
                if any(bad_oil)
                    sO(bad_oil) = 1 - sW(bad_oil) - value(sG(bad_oil));
                    state = model.setProp(state, 'sO', sO);
                end
            end
        end
        
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Treat source or bc terms
            if ~isempty(drivingForces.bc) || ~isempty(drivingForces.src)
                [pressures, sat, mob, rho, rs, rv] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'Rs', 'Rv');
                dissolved = model.getDissolutionMatrix(rs, rv);
                eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                                 pressures, sat, mob, rho, ...
                                                                 dissolved, {}, ...
                                                                 drivingForces);
            end

            % Add aquifer contributions if any.
            if ~isempty(model.AquiferModel)
                eqs = addAquifersContribution(model.AquiferModel, eqs, names, state, dt);
            end

            % Add sources
            eqs = model.insertSources(eqs, src);
            % Assemble equations
            for i = 1:numel(eqs)
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            eqs = [eqs, weqs];
            names = [names, wnames];
            types = [types, wtypes];
        end

        function names = getComponentNames(model)
            names = cellfun(@(x) x.name, model.Components, 'UniformOutput', false);
        end

        function [state, report] = updateState(model, state, problem, dx, forces)
            [state, report] = updateState@ThreePhaseBlackOilModel(model, state, problem, dx, forces);
            if ~isempty(model.FacilityModel)
                state = model.FacilityModel.applyWellLimits(state);
            end
        end

        function model = validateModel(model, varargin)
            % Validate model.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.validateModel`
            assert(isa(model.FacilityModel, 'GenericFacilityModel'), ...
                'Generic model can only be used with GenericFacilityModel.')
            if isempty(model.Components)
                nph = model.getNumberOfPhases();
                model.Components = cell(1, nph);
                names = model.getPhaseNames();
                disgas = model.disgas;
                vapoil = model.vapoil;
                for ph = 1:nph
                    switch names(ph)
                        case 'W'
                            c = ImmiscibleComponent('water', ph);
                        case 'O'
                            if disgas || vapoil
                                c = OilComponent('oil', ph, disgas, vapoil);
                            else
                                c = ImmiscibleComponent('oil', ph);
                            end
                        case 'G'
                            if disgas || vapoil
                                c = GasComponent('gas', ph, disgas, vapoil);
                            else
                                c = ImmiscibleComponent('gas', ph);
                            end
                        otherwise
                            error('Unknown phase');
                    end
                    model.Components{ph} = c;
                end
            end
            model = validateModel@ThreePhaseBlackOilModel(model, varargin{:});
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if model.outputFluxes
                state_flow = model.FlowDiscretization.buildFlowState(model, state, state0, dt);
                f = model.getProp(state_flow, 'PhaseFlux');
                nph = numel(f);
                state.flux = zeros(model.G.faces.num, nph);
                state.flux(model.operators.internalConn, :) = [f{:}];

                if ~isempty(drivingForces.bc)
                    [p, s, mob, rho, b] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ShrinkageFactors');
                    sat = num2cell(s, [1, size(s, 1)]);               
                    [~, ~, ~, fRes] = getBoundaryConditionFluxesAD(model, p, sat, mob, rho, b, drivingForces.bc);
                    idx = model.getActivePhases();
                    fWOG = cell(3, 1);
                    fWOG(idx) = fRes;

                    state = model.storeBoundaryFluxes(state, fWOG{1}, fWOG{2}, fWOG{3}, drivingForces);
                end
            end
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
