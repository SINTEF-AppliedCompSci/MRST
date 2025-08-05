classdef BiochemistryModel <  GenericOverallCompositionModel 
%Bio-chemical model for compositional mixture with Hydrogen (H2)
%
% SYNOPSIS:
%   model = BioCompositionalModel(G, rock, fluid)
%   model = BioCompositionalModel(G, rock, fluid, compFluid)
%   model = BioCompositionalModel(..., 'pn1', vn1, ...)
%
% DESCRIPTION:mrs
%   This model forms the basis for simulation of bio-chemical systems within compositional 
%   models. The model couple compositional modelto  bio-chemical reactions model. We use only 
%   microbial growth and decay of MOnod type.
%
% REQUIRED PARAMETERS:
%   G         - Simulation goprid.
%
%   rock      - Valid rock used for the model.
%
%   fluid     - Fluid model used for the model.
%
%   compFluid - Compositional fluid mixture. Optional, defaults to H2O
%               mixture if left empty.
%
% OPTIONAL PARAMETERS:
%   'property' - Set property to the specified value.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   `ReservoirModel`, `ThreePhaseCompotitionalModel`
%

    properties

        % Boolean indicating if we are considering bio-chemical
        % effects (for debugging). Default = true.
        bactrial = true; 
        %  Default is 'nbact'.
        bacterialFormulation = 'bacterialmodel';

        % Compositional fluid mixutre. Default is a `Methanogenesis`
        % with H2-CO2-H2O-CH4
        compFluid 
        % Physical quantities and bounds
        Y_H2 = 3.90875e11; % in 1/mole(H2)
        gammak   = [];
        mol_diff = [];
        alphaH2  = 3.6e-7;
        alphaCO2 = 1.98e-6;

        Psigrowthmax = 1.338e-4; % 1/s        
        b_bact       = 2.35148e-6; % 1/s
        
        Db = 10^(-8)*meter/second
        bDiffusionEffect = false;
        moleculardiffusion = false;
        
        nbactMax = 1e9; % 1/m^3
        
        bacteriamodel = true;
        metabolicReaction = 'MethanogenicArchae';

    end
    
    methods
        %-----------------------------------------------------------------%
        function model = BiochemistryModel(G, rock, fluid, compFluid, includeWater, backend, varargin)
        % Class constructor. Required arguments are G, rock and fluid.
            model = model@GenericOverallCompositionModel(G, rock, fluid, compFluid, 'water', includeWater, 'AutoDiffBackend', backend);
            model = merge_options(model, varargin{:});
            % set up operators later
            model = model.setupOperators();
            % Check if we have a water phase
            model.gas = true;
            if ~includeWater
                assert(model.oil, ...
                'we need a liquid phase');
            end

            namecp = compFluid.names();
            if model.moleculardiffusion
                indices = struct('H2', find(strcmp(namecp, 'H2')), ... 
                    'C1', find(strcmp(namecp, 'C1')), ... 
                    'CO2', find(strcmp(namecp, 'CO2')), ...
                    'H2O', find(strcmp(namecp, 'H2O')), ...
                    'N2', find(strcmp(namecp, 'N2')), ...
                    'C2', find(strcmp(namecp, 'C2')), ...
                    'C3', find(strcmp(namecp, 'C3')), ...
                    'NC4', find(strcmp(namecp, 'NC4')));
                
                fields = fieldnames(indices);
                nfields=numel(fields);
                model.mol_diff=zeros(nfields,2);

                coeffs = struct(...
                    'H2',  [4.5e-9, 6.1e-5], ...
                    'C1', [2.6e-9, 1.6e-5], ...
                    'H2O', [2.3e-9, 1.5e-5], ...
                    'CO2', [1.9e-9, 1.4e-5], ... 
                    'N2',  [2.1e-9, 1.8e-5], ... 
                    'C2',  [3.2e-9, 2.5e-5], ... 
                    'C3',  [2.8e-9, 2.2e-5], ... 
                    'NC4', [2.4e-9, 1.9e-5]);
                
                for i = 1:nfields
                    comp = fields{i};
                    indComp = indices.(comp);
                    
                    if ~isempty(indComp) && isfield(coeffs, comp)
                      model.mol_diff(indices.(comp),:)= coeffs.(comp);
                    end
                 end
            end

            % Set compositinal fluid
            if isempty(compFluid)
            
                % Default is Methanogenesis
                %compFluid = TableCompositionalMixture({'Hydrogen', 'Water','Nitrogen', 'CarbonDioxide', 'Methane'}, ...
                %   {'H2', 'Water', 'N2', 'CO2', 'CH4'});
                if strcmp(model.metabolicReaction,'MethanogenicArchae')
                    compFluid = TableCompositionalMixture({'Hydrogen', 'Water','Nitrogen', 'CarbonDioxide', 'Methane'}, ...
                    {'H2', 'H2O', 'N2', 'CO2', 'C1'});
                    model.gammak = [-4.0, 2.0, 0.0, -1.0, 1.0];
                    eos =SoreideWhitsonEquationOfStateModel([], compFluid, 'pr');
                    model.EOSModel = eos;
               else
                  warning('MethanogenicArchae is the default now; other reaction not yet implemented');
                end
            else
                ncomp=compFluid.getNumberOfComponents();
                model.gammak=zeros(1,ncomp);
                if strcmp(model.metabolicReaction,'MethanogenicArchae')
                    indH2=find(strcmp(namecp,'H2'));
                    indH2O= find(strcmp(namecp,'H2O'));
                    indCO2=find(strcmp(namecp,'CO2'));
                    indC1= find(strcmp(namecp,'C1'));
                    model.gammak(indH2)  = -4.0;
                    model.gammak(indH2O) = 2.0;
                    model.gammak(indCO2) = -1.0;
                    model.gammak(indC1)  = 1.0;
                end                 
            end
            model.compFluid = compFluid;
                                   
            model.compFluid = compFluid;                
            eos = SoreideWhitsonEquationOfStateModel([], compFluid, 'sw');
            model.EOSModel = eos;
            % Set nonlinear tolerance
%             model.nonlinearTolerance = 1e-3;
            % Check that we have a valid thermal formulation
            assert(any(strcmpi(model.bacterialFormulation, {'bacterialmodel'})), ...
                'BioChemsitryModel supports currently only one micro-organism');
            % Set output state functions
            model.OutputStateFunctions = {'ComponentTotalMass', 'Density'};
            if model.bacteriamodel
                model.FlowDiscretization = BiochemicalFlowDiscretization(model);
            end
        end
        
        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@GenericOverallCompositionModel(model);
        end
        
        function model = setupOperators(model, G, rock, varargin)
            % Set up operators, potentially accounting for dynamic
            % transmissibilites

            % Set rock and grid from model if not provided
            if nargin < 3, rock = model.rock; end
            if nargin < 2, G = model.G;       end

            drock = rock;
            if model.dynamicFlowTrans()
                % Assign dummy transmissibilities to appease
                % model.setupOperators
                drock = rock;
                nbact0 = 0;
                drock.perm = rock.perm(1*barsa(),nbact0);
            end

            if model.dynamicFlowPv()
                % Assign dummy transmissibilities to appease
                % model.setupOperators
                if ~model.dynamicFlowTrans()                
                    drock = rock;
                    nbact0 = 0;
                end
                drock.poro = rock.poro(1*barsa(),nbact0);
            end
            % Let reservoir model set up operators
            model = setupOperators@ReservoirModel(model, G, drock, varargin{:});
            model.rock = rock;
        end

        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
            % Validate model to see if it is ready for simulation
            % Check that we have a facility model
            if model.bacteriamodel
                if isempty(model.FacilityModel) ...
                        || ~isa(model.FacilityModel, 'BiochemistryGenericFacilityModel')
                    model.FacilityModel = BiochemistryGenericFacilityModel(model);
                end
            else
                if isempty(model.FacilityModel) ...
                        || ~isa(model.FacilityModel, 'GenericFacilityModel')
                    model.FacilityModel = GenericFacilityModel(model);
                end
            end
            model = validateModel@GenericOverallCompositionModel(model, varargin{:});
        end


        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)
            % Set up state function groupings for geothermal simulation


            %model.PVTPropertyFunctions = []; % make sure this ir reset
            model = setupStateFunctionGroupings@GenericOverallCompositionModel(model, varargin{:});
            fluxprops = model.FlowDiscretization;
            pvtprops  = model.PVTPropertyFunctions;
            flowprops = model.FlowPropertyFunctions;
            if model.bacteriamodel
                flowprops = flowprops.setStateFunction('PsiGrowthRate', GrowthBactRateSRC(model));
                flowprops = flowprops.setStateFunction('PsiDecayRate', DecayBactRateSRC(model));
                flowprops = flowprops.setStateFunction('BactConvRate', BactConvertionRate(model));
            end

            pvt = pvtprops.getRegionPVT(model);
            if isfield(model.fluid, 'pvMultR')
                % Check for multiplier
                pv = DynamicFlowPoreVolume(model, pvt);
            else
                pv = PoreVolume(model, pvt);
            end

            pvtprops = pvtprops.setStateFunction('PoreVolume', pv);
            model.PVTPropertyFunctions  = pvtprops;
            model.FlowPropertyFunctions = flowprops;
            model.FlowDiscretization    = fluxprops;

        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            % Validate state and check if it is ready for simulation

            % Let parent model do it's thing
            state = validateState@ThreePhaseCompositionalModel(model, state);
            if model.bacteriamodel
                if ~isfield(state, 'nbact')
                    % Set temperature if it is not given
                    nbact0 = 10^6;
                    state.nbact = repmat(nbact0, model.G.cells.num, 1);
                end
            end
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            % Get primary variables from state, before a possible
            % initialization as AD.
            [p, z] = model.getProps(state, 'pressure', 'z');
            z_tol = model.EOSModel.minimumComposition;
            z = ensureMinimumFraction(z, z_tol);
            z = expandMatrixToCell(z);
            cnames = model.EOSModel.getComponentNames();
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
                nbact = model.getProps(state, 'bacteriamodel'); 
                names = [{'pressure'}, cnames(2:end), {'nbact'}, enames];
                vars = [p, z(2:end), nbact, evars];
                origin = cell(1, numel(names));
                [origin{:}] = [deal(class(model))];
            else
                names = [{'pressure'}, cnames(2:end), enames];
                vars = [p, z(2:end), evars];
                origin = cell(1, numel(names));
                [origin{:}] = deal(class(model));
            end

            if ~isempty(model.FacilityModel)
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
            end
        end

        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            % Discretize
            % state = capSaturation(model,state, 's', 1.0e-8, 1-1.0e-8);
            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Assemble equations and add in sources
            [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ComponentPhaseMassFractions');
            comps = cellfun(@(x, y) {x, y}, X(:, model.getLiquidIndex), X(:, model.getVaporIndex), 'UniformOutput', false);


            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                pressures, sat, mob, rho, ...
                {}, comps, ...
                drivingForces);

            % Add sources
            eqs = model.insertSources(eqs, src);
            % Assemble equations

            if model.bacteriamodel
                cnames = model.EOSModel.getComponentNames();
                ncomp = numel(cnames);
                src_rate = model.FacilityModel.getProps(state, 'BactConvRate');
                for i = 1:ncomp
                    if ~isempty(src_rate{i})
                        eqs{i} = eqs{i} -src_rate{i};
                    end
                end
            end

            for i = 1:numel(eqs)
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            if model.bacteriamodel
                [beqs, bflux, bnames, btypes] = model.FlowDiscretization.bacteriaConservationEquation(model, state, state0, dt);
                fd = model.FlowDiscretization;
                src_growthdecay = model.FacilityModel.getBacteriaSources(fd, state, state0, dt);
                beqs{1} = beqs{1} - src_growthdecay;
                %  treat diffusion separately
                if any(model.bDiffusionEffect > 0)
                    beqs{1} = model.operators.AccDiv(beqs{1}, bflux{1});
                end
            else
                [beqs, bnames, btypes] = deal([]);
            end
            % Concatenate
            eqs   = [eqs, beqs];
            names = [names, bnames];
            types = [types, btypes];


                
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            % Concatenate
            eqs   = [eqs  , weqs  ];
            names = [names, wnames];
            types = [types, wtypes];

        end
       
        function forces = validateDrivingForces(model, forces, varargin)
            forces = validateDrivingForces@GenericOverallCompositionModel(model, forces, varargin{:});
        end 

        function state = initStateAD(model, state, vars, names, origin)        
            if model.bacteriamodel

                isP = strcmp(names, 'pressure');
                isB = strcmp(names, 'nbact');
                isAD = any(cellfun(@(x) isa(x, 'ADI'), vars));
                state = model.setProp(state, 'pressure', vars{isP});
                state = model.setProp(state, 'nbact', vars{isB});

                removed = isP | isB;

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
                li = model.getLiquidIndex();
                vi = model.getVaporIndex();
                % Set up state with remaining variables
                state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));

                % Now that props have been set up, we can compute the
                % saturations from the mole fractions.
                if isAD
                    % We must get the version with derivatives
                    Z = model.getProps(state, 'PhaseCompressibilityFactors');
                    Z_L = Z{li};
                    Z_V = Z{vi};
                else
                    % Already stored in state - no derivatives needed
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
            % Get values for convergence check
            [v_eqs, tolerances, names] = getConvergenceValues@ReservoirModel(model, problem, varargin{:});
            bacteriaIndex = find(strcmp(names, 'bacteria (cell)'));
            tolerances(bacteriaIndex) = 5.0e-2;
            if model.bacteriamodel
                scale = model.getEquationScaling(problem.equations, problem.equationNames, problem.state, problem.dt);
                ix    = ~cellfun(@isempty, scale);
                v_eqs(ix) = cellfun(@(scale, x) norm(scale.*value(x), inf), scale(ix), problem.equations(ix));
                % Reduce Tolerance
                iter = problem.iterationNo;
                maxIter = model.EOSNonLinearSolver.LinearSolver.maxIterations;
                if (v_eqs(bacteriaIndex) > tolerances(bacteriaIndex) && (iter+5>maxIter))
                     % tolerances(bacteriaIndex) = 0.75*tolerances(bacteriaIndex).*iter;
                end
            end
        end
        
        function scale = getEquationScaling(model, eqs, names, state0, dt)
            % Get scaling for the residual equations to determine convergence

            scale = cell(1, numel(eqs));
            cnames = model.getComponentNames();

            if model.bacteriamodel
                [cmass, chemistry] = model.getProps(state0, 'ComponentTotalMass', ...
                    'BacterialMass');
                cmass = value(cmass);
                chemistry = value(chemistry);
            else
                cmass= model.getProps(state0, 'ComponentTotalMass');
                cmass = value(cmass);
            end

            if ~iscell(cmass), cmass = {cmass}; end
            ncomp = model.getNumberOfComponents();
            mass = 0;
            for i = 1:ncomp
                mass = mass + cmass{i};
            end

            scaleMass = dt./mass;
            for n = cnames
                ix = strcmpi(n{1}, names);
                if ~any(ix), continue; end
                scale{ix} = scaleMass;
            end
            if model.bacteriamodel
                ix = strcmpi(names, 'bacteria');
                if any(ix)
                    scaleChemistry = dt./max(chemistry, dt);
                    scaleChemistry = filloutliers(scaleChemistry, "nearest","mean");
                    scale{ix} = scaleChemistry;
                end
            end

        end

        function scaling = getScalingFactorsCPR(model, problem, names, solver) %#ok

            scaling = model.getEquationScaling(problem.equations, problem.equationNames, problem.state, problem.dt);

        end
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
                case {'nbact', 'bacteriamodel'} %Bacteria model
                    index = ':';
                    fn = 'nbact';
                    %========================================
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@OverallCompositionCompositionalModel(model, name, varargin{:});
            end
        end

        function names = getComponentNames(model)
            % Get names of the fluid components
            names  = getComponentNames@GenericOverallCompositionModel(model);

        end


        function  [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@GenericOverallCompositionModel(model, state0, state, dt, drivingForces);
        end

        function [state, report] = updateState(model, state, problem, dz, drivingForces)
            [state, report] = updateState@GenericOverallCompositionModel(model, state, problem, dz, drivingForces);
            if model.bacteriamodel
                state = model.capProperty(state, 'nbact', 1.e-3, 120);
                                   
                state = model.capProperty(state, 's', 1.0e-8, 1);
                state.components = ensureMinimumFraction(state.components, model.EOSModel.minimumComposition);
            end
        end


        function isDynamic = dynamicFlowTrans(model)
            % Get boolean indicating if the fluid flow transmissibility is
            % dynamically calculated
            isDynamic = isa(model.rock.perm, 'function_handle');

        end

        function isDynamic = dynamicFlowPv(model)
            % Get boolean indicating if the fluid flow porevolume is
            % dynamically calculated

            isDynamic = isa(model.rock.poro, 'function_handle');

        end

        function state = computeBactPopulation(model, state)
            % COMPUTEBACTPOPULATION Computes bacterial population using a quadratic equation.
            % This function estimates the bacterial concentration by solving the microbial
            % growth equation as a quadratic polynomial for a single time step.
            %
            % INPUTS:
            %   model - The simulation model containing bacterial properties.
            %   state - The current simulation state.
            %
            % OUTPUT:
            %   state - Updated simulation state with computed bacterial population.

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

            % Compute coefficients for the quadratic equation
            A = PsigrowthMax .* (xH2 ./ (aH2 + xH2)) .* (xCO2 ./ (aCO2 + xCO2));
            B = bbact ./ nbMax;

            % Quadratic equation: nbact_new - (1 + dt * A) * nbact + dt * B * nbact^2 = 0
            a_quad = dt * B;
            b_quad = -(1 + dt * A);
            c_quad = nbact;

            % Solve using quadratic formula
            discriminant = b_quad.^2 - 4 * a_quad .* c_quad;

            % Ensure discriminant is non-negative
            discriminant = max(discriminant, 0);

            nbact_new = (-b_quad + sqrt(discriminant)) ./ (2 * a_quad);

            % Ensure non-negative values
            nbact_new = max(nbact_new, 0);

            % Update state with computed bacterial population
            state = model.setProp(state, 'nbact', nbact_new);
            state = model.capProperty(state, 'nbact', 08, 1.0e12);

        end

    end
end
function state = capSaturation(model, state, name, minvalue, maxvalue)
% Assume 'state' is a cell array and 'model.getProp' and 'model.setProp' work with cell content
v = model.getProp(state, name); % Get the property, assuming 'v' is a cell array

if iscell(v)
    % Apply operations to each cell
    for i = 1:numel(v)
        % Extract the value from the cell
        value = v{i};

        % Perform the operations
        value = max(minvalue, value);
        if nargin > 4
            value = min(value, maxvalue);
        end

        % Update the cell with the processed value
        v{i} = value;
    end
else
    v = max(minvalue, v);
    if nargin > 4
        v = min(v, maxvalue);
    end
end
    % Update the state property
    state = model.setProp(state, name, v);
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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