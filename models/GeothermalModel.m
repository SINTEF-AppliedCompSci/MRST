classdef GeothermalModel < ReservoirModel & GenericReservoirModel
%Geothermal model for H20 (liquid and vapor) and optionally NaCl
%
% SYNOPSIS:
%   model = GeothermalModel(G, rock, fluid)
%   model = GeothermalModel(G, rock, fluid, compFluid)
%   model = GeothermalModel(..., 'pn1', vn1, ...)
%
% DESCRIPTION:
%   This model forms the basis for simulation of geothermal systems with
%   H2O and (optionally) NaCl. The model also has support for phase changes
%   where water can exist in liquid and vapor form. This is currently under
%   development.
%
% REQUIRED PARAMETERS:
%   G         - Simulation grid.
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

        % Boolean indicating if we are considering thermal
        % effects (for debugging). Default = true.
        thermal = true; 
        % Choice of thermal primary variable. Possible values are 'temperature'
        % and 'enthaply. Default is 'temperaure'. Note: Two-phase liquid-vapor
        % flow requires enthalpy formulation.
        thermalFormulation = 'temperature';
        % Compositional fluid mixutre. Default is a `CompositionalBrineFluid`
        % with H2O as the only component
        compFluid 
        % Physical quantities and bounds
        geothermalGradient = 30*Kelvin/(kilo*meter);
        minimumTemperature = -inf;
        maximumTemperature =  inf;
        % Update limits
        dTMaxRel = 0.2;
        dTMaxAbs = Inf;
        dhMaxRel = 0.2;
        dhMaxAbs = Inf;
        dxMaxAbs = 0.1;

    end
    
    methods
        %-----------------------------------------------------------------%
        function model = GeothermalModel(G, rock, fluid, varargin)
        % Class constructor. Required arguments are G, rock and fluid.
        
            % Check if compositional fluid is provided
            compFluid = [];
            if nargin > 3 && isa(varargin{1}, 'CompositionalMixture')
                compFluid = varargin{1};
                varargin  = varargin(2:end);
            end
            % Call parent constructor, but whithout rock because we will
            % set up operators later
            model = model@ReservoirModel(G, [], fluid);
            model.rock = rock;
            model = model.setupOperators();
            % Check if we have a gas phase
            model.water = true;
            if isfield(model.fluid, 'rhoG')
                model.gas = true;
            end
            % Set compositinal fluid
            if isempty(compFluid)
                % Default is a one-component system with H20 only
                compFluid = CompositionalBrineFluid({'H2O'}           , ...
                                                    18.015281*gram/mol, ...
                                                    0                 );
            end
            model.compFluid = compFluid;
            % Set nonlinear tolerance
            model.nonlinearTolerance = 1e-3;
            % Merge in additional optjions
            model = merge_options(model, varargin{:});
            % Check that we have a valid thermal formulation
            assert(any(strcmpi(model.thermalFormulation, {'temperature', 'enthalpy'})), ...
                'GeothermalModel currently only supports temperature formulation');
            % Set output state functions
            model.OutputStateFunctions = {'ComponentTotalMass', 'TotalThermalEnergy'};
            
        end
        
        %-----------------------------------------------------------------%
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
                drock.perm = rock.perm(1*barsa, 273.15 + 20*Kelvin);
            end
            % Let reservoir model set up operators
            model = setupOperators@ReservoirModel(model, G, drock);
            model.rock = rock;
            
            pv  = model.operators.pv;
            vol = model.G.cells.volumes;
            % Compute heat transmissibilities if we are not using dynamic
            if ~model.dynamicHeatTransRock()
                % Compute static heat transmissibility
                lambdaR = model.rock.lambdaR.*(vol - pv)./vol;
                r       = struct('perm', lambdaR);
                Thr     = getFaceTransmissibility(model.G, r);
                % Assign to model.operators
                model.operators.Thr     = Thr(model.operators.internalConn);
                model.operators.Thr_all = Thr;
            end
            
            if ~model.dynamicHeatTransFluid()
                % Compute static heat transmissibility
                lambdaF = repmat(model.fluid.lambdaF, model.G.cells.num, 1).*pv./vol;
                r       = struct('perm', lambdaF);
                Thf     = getFaceTransmissibility(model.G, r);
                % Assign to model.operators
                model.operators.Thf     = Thf(model.operators.internalConn);
                model.operators.Thf_all = Thf;
            end
            
            model.operators.vol = model.G.cells.volumes;
            
            % Compute hash
            model.operators.hashRock = obj2hash(rock);
            
        end

        %-----------------------------------------------------------------%
        function model = validateModel(model, varargin)
        % Validate model to see if it is ready for simulation
        
            % Check that we have a facility model
            if isempty(model.FacilityModel) ...
                    || ~isa(model.FacilityModel, 'GeothermalGenericFacilityModel')
                model.FacilityModel = GeothermalGenericFacilityModel(model);
            end
            % Set up components
            if isempty(model.Components)
                names = model.getComponentNames();
                nc = numel(names);
                for i = 1:nc
                    name        = names{i};
                    molarMass   = model.compFluid.molarMass(i);
                    diffusivity = model.compFluid.molecularDiffusivity(i);
                    c = BrineComponent(name, molarMass, diffusivity, i);
                    model.Components{i} = c;
                end
            end
            % Call parent model validation
            model = validateModel@ReservoirModel(model, varargin{:});
            
        end
        
        %-----------------------------------------------------------------%
        function ctrl = validateDrivingForces(model, ctrl, index)
        % Validate dirving forces before simulation
            
            ctrl = validateDrivingForces@ReservoirModel(model, ctrl, index);
            if isfield(ctrl, 'bc') && ~isempty(ctrl.bc)
                % Check that bcs have thermal fields
                assert( ...
                    (isfield(ctrl.bc, 'T') ...
                        && numel(ctrl.bc.T) == numel(ctrl.bc.face)) || ...
                    (isfield(ctrl.bc, 'Hflux') ...
                        && numel(ctrl.bc.Hflux) == numel(ctrl.bc.face)),  ...
                    'bc must have a given temperature (T) or heat flux (Hflux)');
            end
            if isfield(ctrl, 'W') && ~isempty(ctrl.W)
                % Check that wells have a prescribed temperature
                assert(all(arrayfun(@(w) isfield(w, 'T'), ctrl.W))  , ...
                       'All wells must have a temperature field (T)');
            end
            
        end
        
        % --------------------------------------------------------------------%
        function forces = getValidDrivingForces(model)
        % Get valid forces. This class adds support for wells, bc and
        % src. Wells can also be operated in groups
         
            forces = getValidDrivingForces@ReservoirModel(model);
            % Support for group control
            forces.groups = [];
            
        end
        
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)
        % Set up state function groupings for geothermal simulation
            
            pvt = model.PVTPropertyFunctions;
            fp  = model.FlowPropertyFunctions;
            fd  = model.FlowDiscretization;
            
            % Set up state functions for parent
            model = setupStateFunctionGroupings@ReservoirModel(model, varargin{:});
            
            % PVT properties
            if isempty(pvt)
                pvt = model.PVTPropertyFunctions;
                % Density and mass
                rho   = ThermalDensity(model);
                b     = DensityDerivedShrinkageFactors(model);
                rhoR  = RockDensity(model);
                massR = RockMass(model);
                pvt   = pvt.setStateFunction('Density'         , rho  );
                pvt   = pvt.setStateFunction('ShrinkageFactors', b    );
                pvt   = pvt.setStateFunction('RockDensity'     , rhoR );
                pvt   = pvt.setStateFunction('RockMass'        , massR); 
                % Viscosity
                mu  = ThermalViscosity(model);
                pvt = pvt.setStateFunction('Viscosity', mu);
                % Enthalpy
                h   = PhaseEnthalpy(model);
                pvt = pvt.setStateFunction('PhaseEnthalpy', h);
                % Components
                massf = ComponentPhaseMassFractionsBrine(model);
                molf  = ComponentPhaseMoleFractionsBrine(model);
                pvt   = pvt.setStateFunction('ComponentPhaseMassFractions', massf);
                pvt   = pvt.setStateFunction('ComponentPhaseMoleFractions', molf );
                % Replace
                model.PVTPropertyFunctions = pvt;
            end
            % Flow properties
            if isempty(fp)
                fp = model.FlowPropertyFunctions;
                % Internal energy
                u  = PhaseInternalEnergy(model);
                uR = RockInternalEnergy(model);
                fp = fp.setStateFunction('PhaseInternalEnergy', u);
                fp = fp.setStateFunction('RockInternalEnergy' , uR);
                % Thermal energy
                energyPh = PhaseThermalEnergy(model);
                energyR  = RockThermalEnergy(model);
                energy   = TotalThermalEnergy(model);
                fp = fp.setStateFunction('PhaseThermalEnergy', energyPh);
                fp = fp.setStateFunction('RockThermalEnergy' , energyR );
                fp = fp.setStateFunction('TotalThermalEnergy', energy  );
                % Replace
                model.FlowPropertyFunctions = fp;
            end
            % Flux discretization
            if isempty(fd)
                fd = model.FlowDiscretization;
                if isempty(fd) || ~isa(fd, 'GeothermalFlowDiscretization')
                    warning('Assuming default flux discretization');
                    fd = GeothermalFlowDiscretization(model);
                    model.FlowDiscretization = fd;
                end
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
        % Validate state and check if it is ready for simulation
    
            % Let parent model do it's thing
            state = validateState@ReservoirModel(model, state);
            % Set component overall mass fractions if they are not given
            ncomp = model.getNumberOfComponents();
            if isfield(state, 'components')
                assert(all(size(state.components) == [model.G.cells.num, ncomp]));
            else
                state.components = zeros(model.G.cells.num, ncomp);
                state.components(:,1) = 1;
            end
            if ~isfield(state, 'T')
                % Set temperature if it is not given
                T0 = (273.15 + 20)*Kelvin;
                state.T = repmat(T0, model.G.cells.num, 1);
            end
            if ~model.thermal
                return
            end
            if isfield(state, 'enthalphy')
                warning(['Enthalpy given in initial state. I''ll assume ', ...
                         'this is consistent with the EOS and carry on'  ]);
                return
            end
            if strcmpi(model.thermalFormulation, 'enthalpy')
                % Set enthalpy from temperature
                [p, T] = model.getProps(state, 'pressure'   , ...
                                               'temperature');
                if model.getNumberOfPhases() == 2
                    h = model.fluid.h(p, T);
                else
                    rho = model.fluid.rhoW(p, T);
                    u   = model.fluid.uW(p, T);
                    h = u + p./rho;
                end
                state.enthalpy = h;
            end
            % Flash computation to ensure that we are in thermodynamic
            % equilibrium
            state = computeFlashGeothermal(model, state);

        end
        
        %-----------------------------------------------------------------%
        function [vararg, control] = getDrivingForces(model, control)
        % Get driving forces in expanded format.
    
            [vararg, control] = getDrivingForces@ReservoirModel(model, control);
            ix = find(strcmpi(vararg, 'bc'));
            if ~isempty(ix) && ~isempty(vararg{ix+1})
                bc = vararg{ix+1};
                if ~isfield(bc, 'components')
                    assert(numel(model.compFluid.names) == 1, ...
                        ['Model has more than one component - please ', ...
                         'provide component field to bc (bc.cmomponents)']);
                    bc.components = ones(numel(bc.face),1);
                    vararg{ix+1}  = bc;
                end
            end
            ix = find(strcmpi(vararg, 'W'));
            if ~isempty(ix) && ~isempty(vararg{ix+1})
                W = vararg{ix+1};
                if isempty(W(1).components)
                    assert(numel(model.compFluid.names) == 1, ...
                        ['Model has more than one component - ...    ', ...
                         'please provide component field to well (W.components)']);
                    [W.components] = deal(1);
                    vararg{ix+1} = W;
                end
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [vars, names, origin] = getPrimaryVariables(model, state)
        % Get primary variables from state
   
            [p, x] = model.getProps(state, 'pressure', 'components');
            x = ensureMinimumFraction(x, 1e-8);
            x = expandMatrixToCell(x);
            cnames = model.getComponentNames();
            vars   = [p, x(:, 2:end)];
            names  = [{'pressure'}, cnames(2:end)];
            if model.thermal
                thName = model.thermalFormulation;
                thVar  = model.getProps(state, thName);
                vars   = [vars , thVar ];
                names  = [names, thName];
            end
            origin = repmat({class(model)}, 1, numel(cnames) + model.thermal);
            if ~isempty(model.FacilityModel)
                [v, n, o] = model.FacilityModel.getPrimaryVariables(state);
                vars = [vars, v];
                names = [names, n];
                origin = [origin, o];
            end
            
        end
        
        %-----------------------------------------------------------------%
        function state = initStateAD(model, state, vars, names, origin)
        % Initialize AD state from double state
    
            isP = strcmp(names, 'pressure');
            state = model.setProp(state, 'pressure', vars{isP});
            removed = isP;
            
            cnames = model.getComponentNames();
            ncomp = numel(cnames);
            x = cell(1, ncomp);
            x_end = ones(model.G.cells.num, 1);
            for i = 1:ncomp
                name = cnames{i};
                sub = strcmp(names, name);
                if any(sub)
                    x{i} = vars{sub};
                    x_end = x_end - x{i};
                    removed(sub) = true;
                else
                    fill = i;
                end
            end
            x{fill} = x_end;
            state = model.setProp(state, 'components', x);
            if ~isempty(model.FacilityModel)
                % Select facility model variables and pass them off to attached
                % class.
                fm = class(model.FacilityModel);
                isF = strcmp(origin, fm);
                state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
                removed = removed | isF;
            end
            % Set up state with remaining variables
            state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
            state = computeFlashGeothermal(model, state);

        end

        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
        % Get equations from AD states

            [eqs, flux, names, types] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Add sources
            eqs = model.insertSources(eqs, src);
            % Assemble equations
            for i = 1:numel(eqs)
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
            end
            if model.thermal
                [eeqs, eflux, enames, etypes] = model.FlowDiscretization.energyConservationEquation(model, state, state0, dt);
                % Add sources
                esrc = model.FacilityModel.getEnergySources(state);
                eeqs = model.insertSources(eeqs, esrc);
                % Assemble equations
                eeqs{1} = model.operators.AccDiv(eeqs{1}, eflux{1});
            else
                [eeqs, enames, etypes] = deal([]);
            end
            % Get facility equations
            [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
            % Concatenate
            eqs   = [eqs  , eeqs  , weqs  ];
            names = [names, enames, wnames];
            types = [types, etypes, wtypes];
            % Add in boundary conditions
            eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, drivingForces);
            
        end
        
        %-----------------------------------------------------------------%
        function [v_eqs, tolerances, names] = getConvergenceValues(model, problem, varargin)
        % Get values for convergence check
        
            [v_eqs, tolerances, names] = getConvergenceValues@ReservoirModel(model, problem, varargin{:});
            scale = model.getEquationScaling(problem.equations, problem.equationNames, problem.state, problem.dt);
            ix    = ~cellfun(@isempty, scale);
            v_eqs(ix) = cellfun(@(scale, x) norm(scale.*value(x), inf), scale(ix), problem.equations(ix));
            
        end

        %-----------------------------------------------------------------%
        function scale = getEquationScaling(model, eqs, names, state0, dt)
        % Get scaling for the residual equations to determine convergence
            
            scale = cell(1, numel(eqs));
            cnames = model.getComponentNames();
            [cmass, energy] = model.getProps(state0, 'ComponentTotalMass', ...
                                                     'TotalThermalEnergy');
            cmass = value(cmass); energy = value(energy);
            if ~iscell(cmass), cmass = {cmass}; end
            ncomp = model.getNumberOfComponents();
            mass = 0;
            for i = 1:ncomp
                mass = mass + cmass{i};
            end
            scaleMass = dt./mass;
            for n = cnames
               ix = strcmpi(n{1}, names);
               if ~any(ix)
                   continue
               end
               scale{ix} = scaleMass;
            end
            ix = strcmpi(names, 'energy');
            if any(ix)
                scaleEnergy = dt./energy;
                scale{ix} = scaleEnergy;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function massFraction = getMassFraction(model, moleFraction)
        % Convert molar fraction to mass fraction
        
            if iscell(moleFraction)
                ncomp = numel(moleFraction);
                mass = cell(1, ncomp);
                totMass = 0;
                for i = 1:ncomp
                    mi = moleFraction{i};
                    if ~isempty(mi)
                        mass{i} = model.Components{i}.molarMass.*moleFraction{i};
                        totMass = totMass + mass{i};
                    end
                end
                massFraction = cell(size(mass));
                for i = 1:ncomp
                    if ~isempty(mass{i})
                        massFraction{i} = mass{i}./totMass;
                    end
                end
            else
                molarMass = reshape(cellfun(@(c) c.molarMass, model.Components), 1, []);
                mass = bsxfun(@times, moleFraction, molarMass);
                massFraction = bsxfun(@rdivide, mass, sum(mass, 2));
            end
            
        end
        
        %-----------------------------------------------------------------%
        function moleFraction = getMoleFraction(model, massFraction)
        % Convert molar fraction to mass fraction
        
            if iscell(massFraction)
                assert(~isa(massFraction{1}, 'ADI'));
                massFraction = horzcat(massFraction{:});
            end
            model = model.validateModel();
            mass = massFraction;
            molarMass = reshape(cellfun(@(c) c.molarMass, model.Components), 1, []);
            mole = bsxfun(@rdivide, mass, molarMass);
            moleFraction = bsxfun(@rdivide, mole, sum(mole, 2));
            
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
        % Map known variable by name to field and column index in `state`

            % Temprature and enthalpy are both read from state even though
            % only one of them is a primary variable
            switch lower(name)
                case {'x', 'molefractions', 'components'}
                    % Liquid phase mole fraction
                    fn = 'components';
                    index = ':';
                case {'enthalpy', 'h'}
                    % Enthalpy
                    fn    = 'enthalpy';
                    index = ':';
                case {'temperature', 't'}
                    % Temperature
                    fn    = 'T';
                    index = ':';
                otherwise
                    cnames = model.getComponentNames();
                    sub = strcmpi(cnames, name);
                    if any(sub)
                        fn    = 'components';
                        index = find(sub);
                    else
                        % This will throw an error for us
                        [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
                    end
            end
            
        end
        
        %-----------------------------------------------------------------%
        function names = getComponentNames(model)
        % Get names of the fluid components
        
            names  = getComponentNames@ReservoirModel(model);
            names  = horzcat(names,  model.compFluid.names);
            
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dz, drivingForces)
        % Update the state based on increments of the primary values
            
            state0 = state;
            vars0 = problem.primaryVariables;
            vars = vars0;
            switch model.thermalFormulation
                case 'temperature'
                    % Temperature
                    state = model.updateStateFromIncrement(state, dz, problem, 'temperature', model.dTMaxRel, model.dTMaxAbs);
                    state = model.capProperty(state, 'temperature', model.minimumTemperature, model.maximumTemperature);
                    [vars, removed] = model.stripVars(vars, 'temperature');
                case 'enthalpy'
                    state = model.updateStateFromIncrement(state, dz, problem, 'enthalpy', model.dhMaxRel, model.dhMaxAbs);
                    [vars, removed] = model.stripVars(vars, 'enthalpy');
            end
            % Components
            cnames = model.getComponentNames();
            ncomp = numel(cnames);
            ok = false(ncomp, 1);
            x = state.components;
            rm = 0;
            for i = 1:ncomp
                name = lower(cnames{i});
                cix = strcmpi(vars0, name);
                if any(cix)
                    x0 = x(:, i);
                    dx = dz{cix};
                    if isfinite(model.dxMaxAbs)
                        dx = sign(dx).*min(abs(dx), model.dxMaxAbs);
                    end
                    x(:, i) = min(max(x0 + dx, 0), 1);
                    ok(i) = true;
                    [vars, ix] = model.stripVars(vars, {name});
                    removed(~removed) = removed(~removed) | ix;
                    rm = rm - (x(:, i) - x0);
                end
            end
            if any(ok)
                % We had components as active variables somehow
                assert(nnz(~ok) == 1)
                x(:, ~ok) = min(max(x(:, ~ok) + rm, 0), 1);
                x = bsxfun(@rdivide, x, sum(x, 2));
                state.components = x;
            else
                state.dx = zeros(1, ncomp);
            end

            % Parent class handles almost everything for us
            problem.primaryVariables = vars;
            dz(removed) = [];
            [state, report] = updateState@ReservoirModel(model, state, problem, dz, drivingForces);
            if ~isempty(model.FacilityModel)
                state = model.FacilityModel.applyWellLimits(state);
            end
            if problem.iterationNo == 1
                state.switched = false(model.G.cells.num, 1);
                state.switchCount = zeros(model.G.cells.num, 1);
            end
            state.components = ensureMinimumFraction(state.components, 1e-8);
            state = computeFlashGeothermal(model, state, state0);
            
        end
        
        %-----------------------------------------------------------------%
        function [eqs, state, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, forces)
        % Add sources from BCs and source terms
            
            [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures'             , ...
                                                                  's'                          , ...
                                                                  'Mobility'                   , ...
                                                                  'Density'                    , ...
                                                                  'ComponentPhaseMassFractions');
            comps = cellfun(@(x) {x}, X, 'UniformOutput', false);
            if ~iscell(sat)
                sat = expandMatrixToCell(sat);
            end
            if model.dynamicFlowTrans
                prop = model.FlowDiscretization.getStateFunction('Transmissibility');
                T    = prop.evaluateOnDomain(model, state, true);
                model.operators.T_all = value(T);
            end
            
            if ~isempty(forces.bc)
                forces.bc = getBCProperties(forces.bc, model, state);
            end
            [eqs, state, src] = addBoundaryConditionsAndSources@ReservoirModel(model, eqs, names, types, state, ...
                                                                      pressures, sat, mob, rho, {}, comps, forces);
            if ~model.thermal
                return
            end
            eix = strcmpi(names, 'energy');
            if ~isempty(src.bc.sourceCells)
                src = getHeatFluxBoundary(model, src, forces);
                q = src.bc.mapping*src.bc.heatFlux;
                eqs{eix}(src.bc.sourceCells) = eqs{eix}(src.bc.sourceCells) - q;
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
            
            if isempty(force)
                return
            end
            cnames = model.getComponentNames();
            sub = strcmpi(cnames, cname);
            if any(sub)
                cells = src.sourceCells;
                x_bc = model.getProp(force, 'components');
                mf_bc = model.getMassFraction(x_bc);    
                massFractions = {mf_bc};
                qC = zeros(size(cells));
                nph = model.getNumberOfPhases;
                for ph = 1:nph
                    q_ph = src.phaseMass{ph};
                    inj = q_ph > 0;
                    qC = qC + ~inj.*component{1}(cells).*q_ph ...
                            +  inj.*massFractions{1}(:, sub).*q_ph;
                end
                if ~isempty(src.mapping)
                    qC = src.mapping*qC;
                end
                eq(cells) = eq(cells) - qC;
                src.components{end+1} = qC;
            else
                [eq, src] = addComponentContributions@ReservoirModel(model, cname, eq, component, src, force);
            end
            
        end
        
        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        % Final update to the state after convergence has been achieved

            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if model.extraStateOutput
                rho = model.getProps(state, 'Density');
                state.rho = horzcat(rho{:});
            end
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
                [heatFluxAdv, heatFluxCond] = model.getProps(state, 'AdvectiveHeatFlux', 'ConductiveHeatFlux');
                state.heatFluxAdv  = zeros(model.G.faces.num, nph);
                state.heatFluxCond = zeros(model.G.faces.num, 1);
                state.heatFluxAdv(model.operators.internalConn,:) = horzcat(heatFluxAdv{:});
                state.heatFluxCond(model.operators.internalConn)  = heatFluxCond;
            end
            
        end
        
        % ----------------------------------------------------------------%
        function [model, state] = updateForChangedControls(model, state, forces)
        % Update model and state when controls/drivingForces has changed
         
            state0 = state;
            [model, state] = updateForChangedControls@ReservoirModel(model, state, forces);
            if isfield(state0, 'wellSol') && ~isempty(state0.wellSol)
                inactive = ~vertcat(state.wellSol.status);
                if any(inactive) % this safeguard test required by octave
                   [state.wellSol(inactive).T] = deal(state0.wellSol(inactive).T);
                end
            end
            
        end
        
        % ----------------------------------------------------------------%
        function scaling = getScalingFactorsCPR(model, problem, names, solver) %#ok
            scaling = model.getEquationScaling(problem.equations, problem.equationNames, problem.state, problem.dt);
        end
        
        %-----------------------------------------------------------------%
        function isDynamic = dynamicFlowTrans(model)
        % Get boolean indicating if the fluid flow transmissibility is
        % dynamically calculated
        
           isDynamic = isa(model.rock.perm, 'function_handle');
           
        end
        
        %-----------------------------------------------------------------%
        function isDynamic = dynamicHeatTransFluid(model)
        % Get boolean indicating if the fluid heat transmissibility is
        % dynamically calculated
        
           isDynamic = isa(model.fluid.lambdaF, 'function_handle');
           
        end
        
        %-----------------------------------------------------------------%
        function isDynamic = dynamicHeatTransRock(model)
        % Get boolean indicating if the rock heat transmissibility is
        % dynamically calculated
        
           isDynamic = isa(model.rock.lambdaR, 'function_handle');
           
        end
        
    end
    
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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