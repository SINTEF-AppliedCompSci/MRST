classdef ThreePhaseCompositionalModel < ReservoirModel
    % Base class for compositional models
    %
    % SYNOPSIS:
    %   model = ThreePhaseCompositionalModel(G, rock, fluid, compFluid)
    %
    % DESCRIPTION:
    %   This is the base class for several compositional models in MRST. It
    %   contains common functionality and is not intended for direct use.
    %
    % PARAMETERS:
    %   G         - Grid structure
    %   rock      - Rock structure for the reservoir
    %   fluid     - The flow fluid, containing relative permeabilities,
    %               surface densities and flow properties for the
    %               aqueous/water phase (if present)
    %   compFluid - CompositionalMixture instance describing the species
    %               present.
    %
    % RETURNS:
    %  model - Initialized class instance
    %
    % SEE ALSO:
    %   `NaturalVariablesCompositionalModel`, `OverallCompositionCompositionalModel`

    properties
        EOSModel % EquationOfState model used for PVT and phase behavior
        EOSNonLinearSolver % NonLinearSolver used to solve EOS problems that appear
        dzMaxAbs % Maximum allowable change in any mole fraction
        incTolPressure % Relative increment tolerance for pressure
        incTolComposition % Increment tolerance for composition
        useIncTolComposition % If true, use increment tolerance for composition. Otherwise, use mass-balance.
        fugacityTolerance % Tolerance for fugacity equality (in units 1/barsa)
        liquidPhase = 'O';
        vaporPhase = 'G';
    end
    
    methods
        function model = ThreePhaseCompositionalModel(G, rock, fluid, compFluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            
            if isa(compFluid, 'CompositionalMixture')
                model.EOSModel = EquationOfStateModel(G, compFluid);
            elseif isa(compFluid, 'EquationOfStateModel')
                model.EOSModel = compFluid;
            end
            
            model.EOSModel.verbose = false;
            model.EOSNonLinearSolver = getDefaultFlashNonLinearSolver();
            
            model.nonlinearTolerance = 1e-3;
            model.incTolPressure = 1e-3;
            model.useIncTolComposition = false;
            model.incTolComposition = 1e-3;
            
            model.fugacityTolerance = 1e-3;
            
            model.water = true;
            model.oil = true;
            model.gas = true;
            
            model.dzMaxAbs = 0.1;
            model.dsMaxAbs = 0.1;
            
            model.minimumPressure = 0;
            model.dpMaxRel = 0.25;
            model = merge_options_relaxed(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, varargin)
            error('Superclass not intended for direct use!');
        end
        
        function names = getComponentNames(model)
            names = getComponentNames@ReservoirModel(model);
            names = horzcat(names, model.EOSModel.getComponentNames());
        end
        
        function isEoS = getEoSComponentMask(model)
            if isempty(model.Components)
                % Old, non-generic version
                nc = model.EOSModel.getNumberOfComponents();
                isEoS = true(1, nc + model.water);
                if model.water
                    isEoS(1) = false;
                end
            else
                isEoS = cellfun(@(x) isa(x, 'EquationOfStateComponent'), model.Components);
            end
        end
        
        function ix = getLiquidIndex(model)
            ph = model.getPhaseNames();
            ix = find(ph == model.liquidPhase);
        end

        function ix = getVaporIndex(model)
            ph = model.getPhaseNames();
            ix = find(ph == model.vaporPhase);
        end
        
        function ix = getEoSPhaseIndices(model)
            ix = [model.getLiquidIndex, model.getVaporIndex];
        end
        
        function [ph, varargout] = getNonEoSPhaseNames(model)
            varargout = cell(1, nargout-1);
            [ph, varargout{1:end}] = model.getPhaseNames();
            act = ph ~= model.liquidPhase & ph ~= model.vaporPhase;
            ph = ph(act);
            for i = 1:nargout-1
                varargout{i} = varargout{i}(act);
            end
        end

        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
                case {'z', 'components'}
                    % Overall mole fraction
                    fn = 'components';
                    index = ':';
                case {'x', 'liquidmf', 'liquidmolefractions'}
                    % Liquid phase mole fraction
                    fn = 'x';
                    index = ':';
                case {'y', 'vapormf', 'vapormolefractions'}
                    % Vapor phase mole fraction
                    fn = 'y';
                    index = ':';
                case {'l', 'liquid'}
                    % Liquid mole fraction
                    fn = 'L';
                    index = 1;
                case {'z_l'}
                    % Liquid compressibility
                    fn = 'Z_L';
                    index = 1;
                case {'z_v'}
                    % Vapor compressibility
                    fn = 'Z_V';
                    index = 1;
                case {'k', 'equilibriumconstants'}
                    fn = 'K';
                    index = ':';
                case 'sl'
                    [fn, index] = getVariableField@ReservoirModel(model, ['s', model.liquidPhase], varargin{:});
                case 'sv'
                    [fn, index] = getVariableField@ReservoirModel(model, ['s', model.vaporPhase], varargin{:});
                otherwise
                    names = model.EOSModel.getComponentNames();
                    sub = strcmpi(names, name);
                    if any(sub)
                        fn = 'components';
                        index = find(sub);
                    else
                        % This will throw an error for us
                        [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
                    end
            end
        end

        function state = validateState(model, state)
            state = validateState@ReservoirModel(model, state);
            ncell = model.G.cells.num;
            ncomp = model.EOSModel.getNumberOfComponents();
            model.checkProperty(state, 'Components', [ncell, ncomp], [1, 2]);
            assert(all(max(model.getProp(state, 'Components')) <= 1), ...
                'Initial mole fractions are larger than unity.')
            T = model.getProp(state, 'Temperature');
            if numel(T) == 1
                % Expand single temperature to all grid cells
                fn = model.getVariableField('Temperature');
                state = rmfield(state, fn);
                state = model.setProp(state, 'Temperature', repmat(T, ncell, 1));
            end
            model.checkProperty(state, 'Temperature', [ncell, 1], [1, 2]);
            if ~isfield(state, 'x') || ~isfield(state, 'K')
                state.components = ensureMinimumFraction(state.components, model.EOSModel.minimumComposition);
                state = model.computeFlash(state, inf);
            end
            if isfield(state, 'wellSol') && ~isfield(state.wellSol, 'components')
                for i = 1:numel(state.wellSol)
                    state.wellSol(i).components = [];
                end
            end
        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
            if ~isfield(state, 'dpRel')
                % Set dpRel if it is not present (typically at the very
                % first iteration)
                if nonlinsolver.minIterations > 0
                    state.dpRel = inf(model.G.cells.num, 1);
                else
                    state.dpRel = zeros(model.G.cells.num, 1);
                end
            end
            [state, report] = stepFunction@ReservoirModel(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p0 = state.pressure;
            [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
            range = max(p0) - min(p0);
            tol = model.incTolPressure;
            if isinf(tol)
                tol = 1e-3;
            end
            range = max(range, mean(p0)*tol);
            state.dpRel = (state.pressure - p0)./range;
            state.dpAbs = state.pressure - p0;
        end
        
        function state = computeFlash(model, state, dt, iteration)
            % Flash a state with the model's EOS.
            if nargin < 4
                iteration = 1;
            end
            if nargin < 3
                dt = inf;
            end
            state0 = state;
            if iteration == 1
                state.eos.iterations = 0;
                state.eos.cellcount = 0;
                state = model.EOSModel.validateState(state);
            end
            [state, report] = model.EOSNonLinearSolver.solveTimestep(state, dt, model.EOSModel);

            if ~isempty(report)
                state.eos.iterations = state.eos.iterations + report.StepReports{1}.Iterations;

                if ~report.StepReports{1}.Converged
                    state = model.EOSModel.updateAfterConvergence(state0, state, dt, struct());
                    disp   ('********************************');
                    disp   ('*    Flash did not converge    *');
                    disp   ('********************************');
                    fprintf('* Final residuals after %d iterations:\n', report.StepReports{1}.Iterations);
                    for i = 1:numel(report.StepReports{1}.NonlinearReport{end}.Residuals)
                        fprintf('* %s: %1.4g \n', model.EOSModel.CompositionalMixture.names{i},....
                                report.StepReports{1}.NonlinearReport{end}.Residuals(i))
                    end
                end
            end
            L = state.L;
            Z_L = state.Z_L;
            Z_V = state.Z_V;
            sL = L.*Z_L./(L.*Z_L + (1-L).*Z_V);
            void = 1;
            nph = model.getNumberOfPhases();
            if nph > 2
                ph = getNonEoSPhaseNames(model);
                for i = 1:numel(ph)
                    void = void - model.getProp(state, ['s', ph(i)]);
                end
            end
            l = model.liquidPhase;
            v = model.vaporPhase;
            state = model.setProp(state, ['s', l], void.*sL);
            state = model.setProp(state, ['s', v], void.*(1-sL));

            assert(all(all(state.s >= 0)), 'Negative saturations after flash.');
        end
    
        
        function [v_eqs, tolerances, names] = getConvergenceValues(model, problem, varargin)
            [v_wells, tol_wells, names_wells, is_well] = ...
                model.FacilityModel.getFacilityConvergenceValues(problem);

            v_eqs = norm(problem, inf);
            % Check components
            [v_comp, tol_comp, names_comp, is_comp] = model.getComponentConvergenceValues(problem);
            % Check fugacity
            [v_f, tol_f, names_f, is_f] = model.getFugacityConvergenceValues(problem);
            % Remaining values use nonlinear tolerances
            rest = ~(is_well | is_f | is_comp);
            v_rest = v_eqs(rest);
            tol_rest = repmat(model.nonlinearTolerance, size(v_rest));
            names_rest = problem.equationNames(rest);
            % Define tolerances, values and names
            tolerances = [tol_comp, tol_rest, tol_f, tol_wells];
            v_eqs = [v_comp, v_rest, v_f, v_wells];
            names = horzcat(names_comp, names_rest, names_f, names_wells);
            % Pressure tolerances
            if any(strcmpi(problem.primaryVariables, 'pressure'))
                % Check if dpRel is present (has been set in stepFunction,
                % but we may have come here from somewhere else)
                if isfield(problem.state, 'dpRel')
                    dp = norm(problem.state.dpRel, inf);
                else
                    dp = inf;
                end
                if isfinite(model.incTolPressure)
                    isp = strcmpi(names, 'pressure');
                    if any(isp)
                        v_eqs(isp) = dp;
                        tolerances(isp) = model.incTolPressure;
                        names{isp} = 'dPressure';
                    else
                        v_eqs = [dp, v_eqs];
                        tolerances = [model.incTolPressure, tolerances];
                        names = ['deltaP', names];
                    end
                elseif isa(model, 'PressureNaturalVariablesModel')
                    pRes = norm(problem.b, inf);
                    if isempty(pRes)
                        pRes = inf;
                    end
                    v_eqs = [pRes, v_eqs];
                    tolerances = [model.nonlinearTolerance, tolerances];
                    names = ['Pressure', names];
                end
            end
        end
        
        function state = storeDensities(model, state, rhoW, rhoO, rhoG)
            % Densities
            isActive = model.getActivePhases();

            state.rho = zeros(model.G.cells.num, sum(isActive));
            rho = {value(rhoW), value(rhoO), value(rhoG)};
            state = model.setPhaseData(state, rho, 'rho');
        end
        
        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ReservoirModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            names = model.getComponentNames();
            ncomp = numel(names);
            
            
            z_well = model.getProp(well.W, 'components');
            mf_well = model.EOSModel.getMassFraction(z_well);
            l = model.getLiquidIndex();
            v = model.getVaporIndex();
            cqLs = qMass{l};
            cqVs = qMass{v};
            ncell = numelValue(cqLs);

            N = numel(compSrc);
            compSrc = [compSrc, cell(1, ncomp)];
            
            wellSol.components = zeros(ncell, ncomp);
            for cNo = 1:ncomp
                Z_well = mf_well(cNo);
                X_res = packed.components{cNo}{1};
                Y_res = packed.components{cNo}{2};

                injO = cqLs > 0;
                injG = cqVs > 0;
                % Account for both phases.
                q_i = (cqLs.*injO + cqVs.*injG).*Z_well ...
                       + ~injO.*X_res.*cqLs + ~injG.*Y_res.*cqVs;
                compSrc{N+cNo} = q_i;
                wellSol.components(:, cNo) = value(q_i);
            end
        end

        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
            if isempty(force)
                return
            end
            cnames = model.getComponentNames();
            sub = strcmpi(cnames, cname);
            if any(sub)
                cells = src.sourceCells;
                if isfield(force, 'componentMass')
                    qC = force.componentMass(:, sub);
                else
                    if isfield(force, 'xM')
                        assert(isfield(force, 'yM'))
                        massFractions = {force.xM, force.yM};
                    elseif isfield(force, 'x')
                        assert(isfield(force, 'y'))
                        [x_bc, y_bc] = model.getProps(force, 'x', 'y');
                        massFractions = {model.EOSModel.getMassFraction(x_bc), ...
                                         model.EOSModel.getMassFraction(y_bc)};
                    else
                        z_bc = model.getProp(force, 'components');
                        mf_bc = model.EOSModel.getMassFraction(z_bc);
                        
                        massFractions = {mf_bc, mf_bc};
                    end
                    qC = zeros(size(cells));
                    for ph = 1:2
                        if ph == 1
                            ix = model.getLiquidIndex();
                        else
                            ix = model.getVaporIndex();
                        end
                        q_ph = src.phaseMass{ix};
                        inj = q_ph > 0;

                        qC = qC + ~inj.*component{ph}(cells).*q_ph ...
                                +  inj.*massFractions{ph}(:, sub).*q_ph;
                    end
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
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if isfield(state, 'wellSol') && isfield(state.wellSol, 'components')
                names = model.getComponentNames();
                for i = 1:numel(names)
                    names{i} = names{i}(isstrprop(names{i}, 'alphanum'));
                end
                isGeneric = isa(model, 'GenericReservoirModel');
                ncomp = numel(names);
                for i = 1:numel(state.wellSol)
                    for j = 1:ncomp
                        if state.wellSol(i).status
                            if isGeneric
                                state.wellSol(i).(names{j}) = sum(state.wellSol(i).ComponentTotalFlux(:, j));
                            else
                                state.wellSol(i).(names{j}) = sum(state.wellSol(i).components(:, j));
                            end
                        end
                    end
                end
            end
            tmpflds = {'dz', 'eos', 'dpRel', 'dpAbs', 'switched', 'switchCount'};
            for fld = tmpflds
                f = fld{1};
                if isfield(state, f)
                    state = rmfield(state, f);
                end
            end
        end
        
        function state = setFlag(model, varargin)
            state = model.EOSModel.setFlag(varargin{:});
        end
        
        function [isLiquid, isVapor, is2ph] = getFlag(model, state)
            [isLiquid, isVapor, is2ph] = model.EOSModel.getFlag(state);
        end
        
        function is2ph = getTwoPhaseFlag(model, state)
            is2ph = model.EOSModel.getTwoPhaseFlag(state);
        end
        
        function m = PropertyModel(model)
            m = model.EOSModel.PropertyModel;
        end
        
        function scaling = getScalingFactorsCPR(model, problem, names, solver) %#ok
            % Get scaling factors for CPR reduction in `CPRSolverAD`
            %
            % PARAMETERS:
            %   model   - Class instance
            %   problem - `LinearizedProblemAD` which is intended for CPR
            %             preconditioning.
            %   names   - The names of the equations for which the factors are
            %             to be obtained.
            %   solver  - The `LinearSolverAD` class requesting the scaling
            %             factors.
            %
            % RETURNS:
            %   scaling - Cell array with either a scalar scaling factor for
            %             each equation, or a vectogetComponentScalingr of equal length to that 
            %             equation.
            %
            % SEE ALSO
            %   `CPRSolverAD`
            scaling = getScalingFactorsCPR@ReservoirModel(model, problem, names, solver);
        end

        function [sL, sV] = setMinimumTwoPhaseSaturations(model, state, s_rem, sL, sV, pureLiquid, pureVapor, twoPhase)
            % Set a minumum phase saturation value for the EOS-governed
            % phases. This may be required for numerical stability, as the
            % component-in-phase conservation equations degenerate when the
            % saturations are zero. The minimum value is taken from the
            % configured value i nthe EOS.
            stol = model.EOSModel.minimumSaturation;
            if any(value(s_rem) > 0)
                if nargin > 7
                    mustHaveVapor = pureVapor | twoPhase;
                    mustHaveLiquid = pureLiquid | twoPhase;
                else
                    mustHaveVapor = pureVapor;
                    mustHaveLiquid = pureLiquid;
                end
                sT = sum(state.s, 2);
                if any(pureVapor)
                    sV(pureVapor) = sT(pureVapor) - s_rem(pureVapor);
                end
                
                if any(mustHaveVapor)
                    if isa(sV, 'ADI')
                        sV.val(mustHaveVapor) = max(sV.val(mustHaveVapor), stol);
                    else
                        sV(mustHaveVapor) = max(sV(mustHaveVapor), stol);
                    end
                end
                
                if any(pureLiquid)
                    sL(pureLiquid) = sT(pureLiquid) - s_rem(pureLiquid);
                end
                
                if any(mustHaveLiquid)
                    if isa(sL, 'ADI')
                        sL.val(mustHaveLiquid) = max(sL.val(mustHaveLiquid), stol);
                    else
                        sL(mustHaveLiquid) = max(sL(mustHaveLiquid), stol);
                    end
                end
            end
        end

        function model = validateModel(model, varargin)
            model = validateModel@ReservoirModel(model, varargin{:});
            % Use matching AD backends for EOS and for flow model
            model.EOSModel.AutoDiffBackend = model.AutoDiffBackend;
            isGeneric = isa(model, 'GenericReservoirModel');
            if isGeneric
                phases = model.getPhaseNames;
                assert(any(phases == model.liquidPhase));
                assert(any(phases == model.vaporPhase));
            else
                assert(strcmp(model.liquidPhase, 'O'))
                assert(strcmp(model.vaporPhase, 'G'))
                assert(model.gas, 'Gaseous phase must be present for compositional');
                assert(model.oil, 'Oileic phase must be present for compositional');
            end
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@ReservoirModel(model, varargin{:});
            % Compositional specializations
            pvtprops = model.PVTPropertyFunctions;
            pvtprops = pvtprops.setStateFunction('ShrinkageFactors', DensityDerivedShrinkageFactors(model));
            pvtprops = pvtprops.setStateFunction('Density', CompositionalDensity(model));
            
            pvtprops = pvtprops.setStateFunction('ComponentPhaseMassFractions', ComponentPhaseMassFractionsLV(model));
            pvtprops = pvtprops.setStateFunction('ComponentPhaseMoleFractions', ComponentPhaseMoleFractionsLV(model));

            pvtprops = pvtprops.setStateFunction('PhaseMixingCoefficients', PhaseMixingCoefficientsLV(model));
            pvtprops = pvtprops.setStateFunction('Fugacity', FugacityLV(model));
            pvtprops = pvtprops.setStateFunction('PhaseCompressibilityFactors', PhaseCompressibilityFactorsLV(model));
            pvtprops = pvtprops.setStateFunction('Viscosity', CompositionalViscosityLV(model));

            model.PVTPropertyFunctions = pvtprops;
            
            fp = model.FlowPropertyFunctions;
            cmass = fp.getStateFunction('ComponentTotalMass');
            if isempty(cmass.getMinimumDerivatives())
                ncomp = model.getNumberOfComponents();
                mv = model.EOSModel.minimumComposition;
                md = repmat(mv/100, 1, ncomp);
                md(1) = mv/barsa; % Pressure alignment - assumed 1.
                cmass = cmass.setMinimumDerivatives(md);
                fp = fp.setStateFunction('ComponentTotalMass', cmass);
            end
            model.FlowPropertyFunctions = fp;
        end
    end
    
    methods (Access=protected)
       function [v_comp, tol_comp, names_comp, isComponent] = getComponentConvergenceValues(model, problem)
            % Check convergence criterion of components
            isComponent = false(size(problem.equations));
            cnames = model.EOSModel.getComponentNames();
            names = problem.equationNames;
            
            % Get state, timestep length
            state = problem.state;
            dt = problem.dt;
            rho = value(model.getProp(state, 'Density'));
            s = value(model.getProp(state, 's'));
            pv = value(model.getProp(state, 'PoreVolume'));
            mass = bsxfun(@times, pv, rho.*s);
            ncomp = numel(cnames);
            for i = 1:ncomp
                f = strcmpi(names, cnames{i});
                if any(f)
                    isComponent(f) = true;
                end
            end
            [shortnames, snames] = model.getNonEoSPhaseNames();
            if model.useIncTolComposition
                if problem.iterationNo == 1
                    v_comp = inf(1, ncomp);
                else
                    v_comp = problem.state.dz;
                end
                tol_comp = model.incTolComposition;
                names_comp = model.getComponentNames();
                names_comp = names_comp(~ismember(names_comp, snames));
                names_comp = cellfun(@(x) ['d', x], names_comp, 'UniformOutput', false);
            else
                phix = model.getEoSPhaseIndices();
                tol_comp = model.nonlinearTolerance;
                massT = sum(mass(:, phix), 2);
                scale = dt./massT;
                if numel(snames)
                    maxEOSSat = sum(s(:, phix), 2);
                    scale(maxEOSSat < 1e-4) = 0;
                end
                v_comp = cellfun(@(x) norm(scale.*value(x), inf), problem.equations(isComponent));
                names_comp = names(isComponent);
            end
            tol_comp = repmat(tol_comp, size(v_comp));
            phnames = model.getPhaseNames();
            for i = 1:numel(shortnames)
                sn = snames{i};
                isSimplePhase = strcmpi(names, sn);
                if any(isSimplePhase)
                    rhoi = rho(:, phnames == shortnames(i));
                    scale_s = dt./(pv.*rhoi);                
                    w_simple = value(problem.equations{isSimplePhase}).*scale_s;
                    v_comp = [v_comp, norm(w_simple, inf)]; %#ok
                    tol_comp = [tol_comp, model.toleranceCNV]; %#ok
                    isComponent(isSimplePhase) = true;
                    names_comp = [names_comp, sn]; %#ok
                end
            end
       end
       
       function [v_f, tol_f, names_f, isFugacity] = getFugacityConvergenceValues(model, problem)
            % Check fugacity constraints if present
            isFugacity = strcmpi(problem.types, 'fugacity');
            phix = model.getEoSPhaseIndices();
            state = problem.state;
            s = value(model.getProp(state, 's'));
            % Scale by amount of EoS phases (equations are meaningless if
            % none are present)
            scale = sum(s(state.flag == 0, phix), 2);

            v_f = cellfun(@(x) norm(value(x).*scale, inf), problem.equations(isFugacity));
            tol_f = repmat(model.fugacityTolerance, size(v_f));
            names_f = problem.equationNames(isFugacity);
       end
       
       function dz = computeChange(model, dz, s_hc)
           dz = bsxfun(@times, abs(dz), s_hc);
           tol = model.toleranceCNV/10;
           dz(s_hc < tol, :) = 0;
           dz = max(dz, [], 1);
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
