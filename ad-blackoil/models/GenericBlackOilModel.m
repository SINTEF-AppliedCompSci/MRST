classdef GenericBlackOilModel < ThreePhaseBlackOilModel & GenericReservoirModel
    properties

    end

    methods
        function model = GenericBlackOilModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});
            model.OutputStateFunctions = {'ComponentTotalMass', 'Density'};
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = getEquations@PhysicalModel(model, state0, state, dt, drivingForces, varargin{:});
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
                    [p, s, mob, r, b] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ShrinkageFactors');
                    sat = expandMatrixToCell(s);
                    rho = expandMatrixToCell(r);

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
