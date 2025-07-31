classdef BiochemistryModel < GenericOverallCompositionModel
    % Bio-chemical model for compositional mixture with Hydrogen (H2)
    %
    % SYNOPSIS:
    %   model = BiochemistryModel(G, rock, fluid)
    %   model = BiochemistryModel(G, rock, fluid, compFluid)
    %   model = BiochemistryModel(..., 'pn1', vn1, ...)
    %
    % DESCRIPTION:
    %   This model forms the basis for simulation of bio-chemical systems within
    %   compositional models. The model couples compositional model to bio-chemical
    %   reactions model. Currently supports microbial growth and decay of Monod type.
    %
    % REQUIRED PARAMETERS:
    %   G         - Simulation grid
    %   rock      - Rock structure
    %   fluid     - Fluid model
    %
    % OPTIONAL PARAMETERS:
    %   compFluid - Compositional fluid mixture (defaults to H2O mixture if empty)
    %
    % SEE ALSO:
    %   `ReservoirModel`, `ThreePhaseCompositionalModel`

    properties
        % Boolean indicating if biochemical effects are considered
        bacteriamodel = true;

        % Formulation type for bacterial transport
        bacterialFormulation = 'bacterialmodel';

        % Compositional fluid mixture (default: Methanogenesis with H2-CO2-H2O-CH4)
        compFluid

        % Physical quantities and bounds
        Y_H2 = 3.90875e11;      % [1/mole(H2)]
        gammak = [];             % Stoichiometric coefficients
        mol_diff = [];           % Molecular diffusion coefficients
        alphaH2  = 3.6e-7;       % H2 Monod coefficient
        alphaCO2 = 1.98e-6;      % CO2 Monod coefficient
        Psigrowthmax = 1.338e-4; % Maximum growth rate [1/s]
        b_bact = 2.35148e-6;     % Decay rate [1/s]
        Db = 1e-8*meter/second;  % Bacterial diffusion coefficient
        nbactMax = 1e9;          % Maximum bacterial concentration [1/m^3]

        % Model configuration flags
        bDiffusionEffect = false;    % Enable bacterial diffusion
        moleculardiffusion = false;  % Enable molecular diffusion

        % Metabolic reaction type
        metabolicReaction = 'MethanogenicArchae';
    end

    methods
        %-----------------------------------------------------------------%
        function model = BiochemistryModel(G, rock, fluid, compFluid, includeWater, backend, varargin)
            % Class constructor

            % Initialize parent model
            model = model@GenericOverallCompositionModel(G, rock, fluid, compFluid, ...
                'water', includeWater, 'AutoDiffBackend', backend);
            model = merge_options(model, varargin{:});

            % Set up operators
            model = model.setupOperators();

            % Phase validation
            model.gas = true;
            if ~includeWater
                assert(model.oil, 'Model requires a liquid phase');
            end

            % Set up molecular diffusion if enabled
            if model.moleculardiffusion
                namecp = compFluid.names();
                indices = struct('H2', find(strcmp(namecp, 'H2')), ...
                    'C1', find(strcmp(namecp, 'C1')), ...
                    'CO2', find(strcmp(namecp, 'CO2')), ...
                    'H2O', find(strcmp(namecp, 'H2O')), ...
                    'N2', find(strcmp(namecp, 'N2')), ...
                    'C2', find(strcmp(namecp, 'C2')), ...
                    'C3', find(strcmp(namecp, 'C3')), ...
                    'NC4', find(strcmp(namecp, 'NC4')));

                % Diffusion coefficients
                coeffs = struct('H2',  [4.5e-9, 6.1e-5], ...
                    'C1',  [2.6e-9, 1.6e-5], ...
                    'H2O', [2.3e-9, 1.5e-5], ...
                    'CO2', [1.9e-9, 1.4e-5], ...
                    'N2',  [2.1e-9, 1.8e-5], ...
                    'C2',  [3.2e-9, 2.5e-5], ...
                    'C3',  [2.8e-9, 2.2e-5], ...
                    'NC4', [2.4e-9, 1.9e-5]);

                % Assign diffusion coefficients
                fields = fieldnames(indices);
                model.mol_diff = zeros(numel(fields), 2);
                for i = 1:numel(fields)
                    comp = fields{i};
                    if isfield(coeffs, comp) && ~isempty(indices.(comp))
                        model.mol_diff(indices.(comp),:) = coeffs.(comp);
                    end
                end
            end

            % Set up compositional fluid and EOS model
            if isempty(compFluid)
                if strcmp(model.metabolicReaction, 'MethanogenicArchae')
                    compFluid = TableCompositionalMixture(...
                        {'Hydrogen', 'Water', 'Nitrogen', 'CarbonDioxide', 'Methane'}, ...
                        {'H2', 'H2O', 'N2', 'CO2', 'C1'});
                    model.gammak = [-4.0, 2.0, 0.0, -1.0, 1.0];
                    eos = SoreideWhitsonEquationOfStateModel([], compFluid, 'pr');
                    model.EOSModel = eos;
                else
                    warning('Only MethanogenicArchae is currently implemented');
                end
            else
                ncomp = compFluid.getNumberOfComponents();
                model.gammak = zeros(1, ncomp);
                if strcmp(model.metabolicReaction, 'MethanogenicArchae')
                    namecp = compFluid.names();
                    indH2  = find(strcmp(namecp, 'H2'));
                    indH2O = find(strcmp(namecp, 'H2O'));
                    indCO2 = find(strcmp(namecp, 'CO2'));
                    indC1  = find(strcmp(namecp, 'C1'));
                    model.gammak(indH2)  = -4.0;
                    model.gammak(indH2O) = 2.0;
                    model.gammak(indCO2) = -1.0;
                    model.gammak(indC1)  = 1.0;
                end
            end

            model.compFluid = compFluid;
            eos = SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
            model.EOSModel = eos;

            % Validate bacterial formulation
            assert(any(strcmpi(model.bacterialFormulation, {'bacterialmodel'})), ...
                'BiochemistryModel currently supports only one micro-organism');

            % Set output state functions
            model.OutputStateFunctions = {'ComponentTotalMass', 'Density'};

            % Set up flow discretization if bacteria model is active
            if model.bacteriamodel
                model.FlowDiscretization = BiochemicalFlowDiscretization(model);
            end
        end

        %-----------------------------------------------------------------%
        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@GenericOverallCompositionModel(model);
        end

        %-----------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
            % Set up operators, accounting for possible dynamic properties

            if nargin < 3, rock = model.rock; end
            if nargin < 2, G = model.G; end

            % Handle dynamic properties
            drock = rock;
            if model.dynamicFlowTrans()
                nbact0 = 0;
                drock.perm = rock.perm(1*barsa(), nbact0);
            end

            if model.dynamicFlowPv()
                if ~model.dynamicFlowTrans()
                    nbact0 = 0;
                end
                drock.poro = rock.poro(1*barsa(), nbact0);
            end

            % Set up operators using parent class
            model = setupOperators@ReservoirModel(model, G, drock, varargin{:});
            model.rock = rock;
        end

        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            % Validate model before simulation

            if model.bacteriamodel
                if isempty(model.FacilityModel) || ...
                        ~isa(model.FacilityModel, 'BiochemistryGenericFacilityModel')
                    model.FacilityModel = BiochemistryGenericFacilityModel(model);
                end
            else
                if isempty(model.FacilityModel) || ...
                        ~isa(model.FacilityModel, 'GenericFacilityModel')
                    model.FacilityModel = GenericFacilityModel(model);
                end
            end
            model = validateModel@GenericOverallCompositionModel(model, varargin{:});
        end

        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)
            % Set up state function groupings for biochemical simulation

            model = setupStateFunctionGroupings@GenericOverallCompositionModel(model, varargin{:});

            % Get property containers
            fluxprops = model.FlowDiscretization;
            pvtprops = model.PVTPropertyFunctions;
            flowprops = model.FlowPropertyFunctions;

            % Add biochemical properties if enabled
            if model.bacteriamodel
                flowprops = flowprops.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
                flowprops = flowprops.setStateFunction('PsiDecayRate', DecayBactRateSRC(model));
                flowprops = flowprops.setStateFunction('BactConvRate', BactConvertionRate(model));
            end

            % Set up pore volume calculation
            pvt = pvtprops.getRegionPVT(model);
            if isfield(model.fluid, 'pvMultR')
                pv = DynamicFlowPoreVolume(model, pvt);
            else
                pv = PoreVolume(model, pvt);
            end
            pvtprops = pvtprops.setStateFunction('PoreVolume', pv);

            % Update model properties
            model.PVTPropertyFunctions = pvtprops;
            model.FlowPropertyFunctions = flowprops;
            model.FlowDiscretization = fluxprops;
        end

        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            % Validate state before simulation

            state = validateState@ThreePhaseCompositionalModel(model, state);

            if model.bacteriamodel && ~isfield(state, 'nbact')
                % Initialize bacterial concentration if not present
                nbact0 = 1e6; % Initial concentration [1/m^3]
                state.nbact = repmat(nbact0, model.G.cells.num, 1);
            end
        end

        %-----------------------------------------------------------------%
        function [vars, names, origin] = getPrimaryVariables(model, state)
            % Get primary variables from state

            [p, z] = model.getProps(state, 'pressure', 'z');
            z_tol = model.EOSModel.minimumComposition;
            z = ensureMinimumFraction(z, z_tol);
            z = expandMatrixToCell(z);

            % Get component names
            cnames = model.EOSModel.getComponentNames();

            % Handle non-EoS phases
            extra = model.getNonEoSPhaseNames();
            ne = numel(extra);
            enames = cell(1, ne);
            evars = cell(1, ne);
            for i = 1:ne
                sn = ['s', extra(i)];
                enames{i} = sn;
                evars{i} = model.getProp(state, sn);
            end

            if model.bacteriamodel
                % Include bacterial concentration if enabled
                nbact = model.getProps(state, 'bacteriamodel');
                names = [{'pressure'}, cnames(2:end), {'nbact'}, enames];
                vars = [p, z(2:end), nbact, evars];
            else
                names = [{'pressure'}, cnames(2:end), enames];
                vars = [p, z(2:end), evars];
            end

            origin = repmat({class(model)}, 1, numel(names));

            % Add facility model variables if present
            if ~isempty(model.FacilityModel)
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
            end
        end

        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            % Get model equations including biochemical terms

            % Get component conservation equations
            [eqs, flux, names, types] = ...
                model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);

            % Add component sources
            src = model.FacilityModel.getComponentSources(state);

            % Get phase properties
            [pressures, sat, mob, rho, X] = model.getProps(state, ...
                'PhasePressures', 's', 'Mobility', 'Density', 'ComponentPhaseMassFractions');

            % Prepare component fractions
            comps = cellfun(@(x, y) {x, y}, ...
                X(:, model.getLiquidIndex), ...
                X(:, model.getVaporIndex), ...
                'UniformOutput', false);

            % Add boundary conditions and sources
            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                pressures, sat, mob, rho, {}, comps, drivingForces);

            % Insert sources
            eqs = model.insertSources(eqs, src);

            % Add bacterial conversion terms if enabled
            if model.bacteriamodel
                cnames = model.EOSModel.getComponentNames();
                ncomp = numel(cnames);
                src_rate = model.FacilityModel.getProps(state, 'BactConvRate');
                for i = 1:ncomp
                    if ~isempty(src_rate{i})
                        eqs{i} = eqs{i} - src_rate{i};
                    end
                end
            end

            % Accumulate and divide fluxes
            for i = 1:numel(eqs)
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end

            % Add bacterial conservation equation if enabled
            if model.bacteriamodel
                [beqs, bflux, bnames, btypes] = ...
                    model.FlowDiscretization.bacteriaConservationEquation(model, state, state0, dt);

                % Add bacterial sources
                fd = model.FlowDiscretization;
                src_growthdecay = model.FacilityModel.getBacteriaSources(fd, state, state0, dt);
                beqs{1} = beqs{1} - src_growthdecay;

                % Handle bacterial diffusion
                if any(model.bDiffusionEffect > 0)
                    beqs{1} = model.operators.AccDiv(beqs{1}, bflux{1});
                end
            else
                [beqs, bnames, btypes] = deal([]);
            end

            % Get well equations
            [weqs, wnames, wtypes, state] = ...
                model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);

            % Combine all equations
            eqs   = [eqs, beqs, weqs];
            names = [names, bnames, wnames];
            types = [types, btypes, wtypes];
        end

        %-----------------------------------------------------------------%
        function forces = validateDrivingForces(model, forces, varargin)
            forces = validateDrivingForces@GenericOverallCompositionModel(model, forces, varargin{:});
        end

        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)
            % Initialize AD state with biochemical variables

            if model.bacteriamodel
                % Handle pressure and bacterial concentration
                isP = strcmp(names, 'pressure');
                isB = strcmp(names, 'nbact');
                isAD = any(cellfun(@(x) isa(x, 'ADI'), vars));

                state = model.setProp(state, 'pressure', vars{isP});
                state = model.setProp(state, 'nbact', vars{isB});
                removed = isP | isB;

                % Handle components
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

                % Compute phase fractions if AD variables
                if isAD
                    [state.x, state.y, state.L, state.FractionalDerivatives] = ...
                        model.EOSModel.getPhaseFractionAsADI(state, state.pressure, state.T, state.components);
                end

                % Handle facility model variables
                if ~isempty(model.FacilityModel)
                    fm = class(model.FacilityModel);
                    isF = strcmp(origin, fm);
                    state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                    removed = removed | isF;
                end

                % Handle saturations
                nph = model.getNumberOfPhases();
                phnames = model.getPhaseNames();
                s = cell(1, nph);
                extra = model.getNonEoSPhaseNames();
                ne = numel(extra);
                void = 1;

                for i = 1:ne
                    sn = ['s', extra(i)];
                    isVar = strcmp(names, sn);
                    si = vars{isVar};
                    removed(isVar) = true;
                    void = void - si;
                    s{phnames == extra(i)} = si;
                end

                % Initialize remaining variables
                state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));

                % Compute saturations from mole fractions
                li = model.getLiquidIndex();
                vi = model.getVaporIndex();

                if isAD
                    Z = model.getProps(state, 'PhaseCompressibilityFactors');
                    Z_L = Z{li};
                    Z_V = Z{vi};
                else
                    Z_L = state.Z_L;
                    Z_V = state.Z_V;
                end

                L = state.L;
                propmodel = model.EOSModel.PropertyModel;

                if isempty(propmodel.volumeShift)
                    volL = L.*Z_L;
                    volV = (1-L).*Z_V;
                else
                    volL = L./propmodel.computeMolarDensity(model.EOSModel, state.pressure, state.x, Z_L, state.T, true);
                    volV = (1-L)./propmodel.computeMolarDensity(model.EOSModel, state.pressure, state.y, Z_V, state.T, false);
                end

                volT = volL + volV;
                sL = volL./volT;
                sV = volV./volT;

                [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
                sL = sL.*void;
                sV = sV.*void;

                [s{li}, s{vi}] = model.setMinimumTwoPhaseSaturations(state, 1 - void, sL, sV, pureLiquid, pureVapor, twoPhase);
                state = model.setProp(state, 's', s);
            else
                state = initStateAD@GenericOverallCompositionModel(model, state, vars, names, origin);
            end
        end

        %-----------------------------------------------------------------%
        function [v_eqs, tolerances, names] = getConvergenceValues(model, problem, varargin)
            % Get values for convergence check with adjusted tolerances

            [v_eqs, tolerances, names] = ...
                getConvergenceValues@ReservoirModel(model, problem, varargin{:});

            if model.bacteriamodel
                % Adjust tolerance for bacterial equation
                bacteriaIndex = find(strcmp(names, 'bacteria (cell)'));
                tolerances(bacteriaIndex) = 5.0e-2;

                % Scale equations
                scale = model.getEquationScaling(problem.equations, ...
                    problem.equationNames, problem.state, problem.dt);
                ix = ~cellfun(@isempty, scale);
                v_eqs(ix) = cellfun(@(scale, x) norm(scale.*value(x), inf), ...
                    scale(ix), problem.equations(ix));

                % Optionally adjust tolerance based on iteration
                iter = problem.iterationNo;
                maxIter = model.EOSNonLinearSolver.LinearSolver.maxIterations;
                if (v_eqs(bacteriaIndex) > tolerances(bacteriaIndex) && (iter+5>maxIter))
                    % tolerances(bacteriaIndex) = 0.75*tolerances(bacteriaIndex).*iter;
                end
            end
        end

        %-----------------------------------------------------------------%
        function scale = getEquationScaling(model, eqs, names, state0, dt)
            % Get scaling for residual equations

            scale = cell(1, numel(eqs));
            cnames = model.getComponentNames();

            % Get component and bacterial masses
            if model.bacteriamodel
                [cmass, chemistry] = model.getProps(state0, ...
                    'ComponentTotalMass', 'BacterialMass');
                cmass = value(cmass);
                chemistry = value(chemistry);
            else
                cmass = model.getProps(state0, 'ComponentTotalMass');
                cmass = value(cmass);
            end

            if ~iscell(cmass), cmass = {cmass}; end

            % Compute total mass for scaling
            ncomp = model.getNumberOfComponents();
            mass = 0;
            for i = 1:ncomp
                mass = mass + cmass{i};
            end

            % Set component scaling
            scaleMass = dt./mass;
            for n = cnames
                ix = strcmpi(n{1}, names);
                if ~any(ix), continue; end
                scale{ix} = scaleMass;
            end

            % Set bacterial scaling if enabled
            if model.bacteriamodel
                ix = strcmpi(names, 'bacteria');
                if any(ix)
                    scaleChemistry = dt./max(chemistry, dt);
                    %scaleChemistry = filloutliers(scaleChemistry, "nearest", "mean");
                    scale{ix} = scaleChemistry;
                end
            end
        end

        %-----------------------------------------------------------------%
        function scaling = getScalingFactorsCPR(model, problem, names, solver)
            % Get scaling factors for CPR preconditioner
            scaling = model.getEquationScaling(problem.equations, ...
                problem.equationNames, problem.state, problem.dt);
        end

        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get field name for variables including biochemical terms

            switch lower(name)
                case {'nbact', 'bacteriamodel'}
                    index = ':';
                    fn = 'nbact';
                otherwise
                    [fn, index] = getVariableField@OverallCompositionCompositionalModel(model, name, varargin{:});
            end
        end

        %-----------------------------------------------------------------%
        function names = getComponentNames(model)
            % Get names of fluid components
            names = getComponentNames@GenericOverallCompositionModel(model);
        end

        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            % Update state after convergence
            [state, report] = updateAfterConvergence@GenericOverallCompositionModel(model, state0, state, dt, drivingForces);
        end

        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dz, drivingForces)
            % Update state with bounds checking

            [state, report] = updateState@GenericOverallCompositionModel(model, state, problem, dz, drivingForces);

            if model.bacteriamodel
                % Apply bounds to bacterial concentration and saturations
                state = model.capProperty(state, 'nbact', 1.e-3, 120);
                state = model.capProperty(state, 's', 1.0e-8, 1);
                state.components = ensureMinimumFraction(state.components, model.EOSModel.minimumComposition);
            end
        end

        %-----------------------------------------------------------------%
        function isDynamic = dynamicFlowTrans(model)
            % Check if fluid flow transmissibility is dynamic
            isDynamic = isa(model.rock.perm, 'function_handle');
        end

        %-----------------------------------------------------------------%
        function isDynamic = dynamicFlowPv(model)
            % Check if pore volume is dynamic
            isDynamic = isa(model.rock.poro, 'function_handle');
        end

        %-----------------------------------------------------------------%
        function state = computeBactPopulation(model, state)
            % Compute bacterial population using quadratic growth equation

            % Extract properties
            s = model.getProps(state, 's');
            x = model.getProps(state, 'x');
            nbact = model.getProps(state, 'nbact');

            % Identify component indices
            namecp = model.EOSModel.getComponentNames();
            idx_H2 = find(strcmp(namecp, 'H2'));
            idx_CO2 = find(strcmp(namecp, 'CO2'));

            % Extract values
            L_ix = model.getLiquidIndex();
            xH2 = x(:, idx_H2);
            xCO2 = x(:, idx_CO2);
            sL = s(:, L_ix);

            % Model parameters
            aH2 = model.alphaH2;
            aCO2 = model.alphaCO2;
            PsigrowthMax = model.Psigrowthmax;
            nbMax = model.nbactMax;
            bbact = model.b_bact;
            dt = 288; % Time step (5 seconds)

            % Compute coefficients for quadratic equation
            A = PsigrowthMax .* (xH2 ./ (aH2 + xH2)) .* (xCO2 ./ (aCO2 + xCO2));
            B = bbact ./ nbMax;

            % Solve quadratic equation: nbact_new - (1 + dt*A)*nbact + dt*B*nbact^2 = 0
            a_quad = dt * B;
            b_quad = -(1 + dt * A);
            c_quad = nbact;

            % Ensure non-negative discriminant
            discriminant = max(b_quad.^2 - 4 * a_quad .* c_quad, 0);

            % Compute new bacterial concentration
            nbact_new = (-b_quad + sqrt(discriminant)) ./ (2 * a_quad);
            nbact_new = max(nbact_new, 0);

            % Update state with bounds
            state = model.setProp(state, 'nbact', nbact_new);
            state = model.capProperty(state, 'nbact', 0.8, 1.0e12);
        end
    end
end

% Helper function for saturation capping
function state = capSaturation(model, state, name, minvalue, maxvalue)
% Cap saturation values between specified bounds

v = model.getProp(state, name);

if iscell(v)
    % Handle cell arrays
    for i = 1:numel(v)
        value = v{i};
        value = max(minvalue, value);
        if nargin > 4
            value = min(value, maxvalue);
        end
        v{i} = value;
    end
else
    % Handle numeric arrays
    v = max(minvalue, v);
    if nargin > 4
        v = min(v, maxvalue);
    end
end

state = model.setProp(state, name, v);
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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