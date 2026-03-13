classdef ECPAEquationOfStateModel < PhysicalModel
    % Equation of state model. Implements electrolyte cubic plus association
    % equation of state with successive substitution and newton solvers.
    properties
        ECPACompositionalMixture % ECPACompositionalMixture
        method = 'ssi' % Type of method
        maxSSI = inf;
        PropertyModel % Model to be used for property evaluations
        selectGibbsMinimum = true; % Use minimum Gibbs energy to select Z
        minimumComposition = 1e-8; % Minimum composition value (for numerical stability)
        minimumSaturation  = 1e-8; % Minimum total EOS phase saturation
        equilibriumConstantFunctions = {};
        extraOutput = 0;
    end
    
    properties (Access = private)
        eosType
    end
    
    methods
        function model = ECPAEquationOfStateModel(G, fluid, eosname)
            if nargin < 3
                eosname = 'eCPA';
            end
            model = model@PhysicalModel(G);
            assert(isa(fluid, 'ECPACompositionalMixture'));
            model.ECPACompositionalMixture = fluid;
            model.nonlinearTolerance = 1e-4;
            model.PropertyModel = ECPACompositionalPropertyModel(fluid);
            model = model.setType(eosname);
        end
        
        function n = getNumberOfComponents(model)
            n = model.ECPACompositionalMixture.getNumberOfComponents();
        end
        
        function n = getNumberOfIons(model)
            n = model.ECPACompositionalMixture.getNumberOfIons();
        end
        
        function n = getNumberOfMolecules(model)
            n = model.ECPACompositionalMixture.getNumberOfMolecules();
        end
        
        function n = getComponentNames(model)
            n = model.ECPACompositionalMixture.names;
        end
        
        function v = solveCubicEOS(model, cubic)
            % Use vectorized cubic solver
            v = cubicPositive(cubic(:,2:end));
            % Avoid picking any complex roots. At least one root will be
            % real. We replace these invalid values with NaN so that when
            % min or max is used, they will ignore any previously complex
            % values.
            mtol = 0;
            v(abs(imag(v)) > mtol) = nan;
            v(v <= 0) = nan;
            v = real(v);
        end
        
        function model = setType(model, arg)
            if ischar(arg)
                switch(lower(arg))
                    case {'cpa', 'ecpa','e-cpa','electrolyte-cubic-plus-association'}
                        %electrolyte-cubic-plus-association
                        model.eosType = 6;
                    otherwise
                        error('Invalid string ''%s''.\n Valid choices are:\n PR: Peng-Robinson\n', arg);
                end
            end
        end
        
        function t = shortname(model)
            switch model.eosType
                case 1
                    t = 'pr';
                case 2
                    t = 'srk';
                case 3
                    t = 'zj';
                case 4
                    t = 'rk';
                case 5
                    t = 'prcorr';
                case 6
                    t = 'ecpa';
                otherwise
                    t = model.eosType;
            end
        end
        
        function [Z, vm, XA, XC] = computeCompressibilityZ(model, p, T, xy, A, B, Ai, Si, Tr, isLiquid)
            if iscell(xy)
                xy = cellfun(@value, xy, 'UniformOutput', false);
                xy = [xy{:}];
            end
            if iscell(Ai)
                Ai = cellfun(@value, Ai, 'UniformOutput', false);
                Ai = [Ai{:}];
            end
            if iscell(Si)
                Si = cellfun(@value, Si, 'UniformOutput', false);
                Si = [Si{:}];
            end
            if iscell(Tr)
                Tr = cellfun(@value, Tr, 'UniformOutput', false);
                Tr = [Tr{:}];
            end
            p = value(p); T = value(T); A = value(A); B = value(B);
            
            R = 8.314462618;
            Bi = model.ECPACompositionalMixture.b;
            epsilonAB = model.ECPACompositionalMixture.epsilonAB;
            betaAB = model.ECPACompositionalMixture.betaAB;
            alpha0 = model.ECPACompositionalMixture.alpha0;
            rBorn = model.ECPACompositionalMixture.rborn;
            rBorn = rBorn(~isnan(rBorn));
            mu0 = model.ECPACompositionalMixture.mu0;
            mu0 = mu0(~isnan(mu0));
            d = model.ECPACompositionalMixture.d;
            d = d(~isnan(d));
            hyd = model.ECPACompositionalMixture.hyd;
            charge = model.ECPACompositionalMixture.charge;
            nmole = model.getNumberOfMolecules();
            
            if isLiquid
                vm = B + 5*1e-6;
            else
                vm = B + 10000*1e-6;
            end
            nc = numel(T);
            active = true(nc, 1);
            xy_loc = xy(active, :);
            T_loc = T(active);
            B_loc = B(active);
            A_loc = A(active);
            vm_loc = vm(active);
            Tr_loc = Tr(active, :);
            p_loc = p(active);
            Ai_loc = Ai(active, :);
            Si_loc = Si(active, :);
            XA = ones(nc, nmole);
            XC = ones(nc, 1);
            XA_loc = XA;
            XC_loc = XC;
            for i = 1:10000
                if ~isempty(mu0)
                    [XA_loc, XC_loc, gxy] = model.pressureA(T_loc, xy_loc(:,1:nmole), B_loc, epsilonAB, betaAB, Bi, vm_loc, Tr_loc);
                    assoic = 2 .* xy_loc(:,1) .* (2 - XA_loc(:,1) - XC_loc(:,1)) ...
                        + sum(xy_loc(:,2:nmole) .* (1 - XA_loc(:,2:nmole)), 2);
                    assoic = 0.5 .* R .* T_loc .* gxy .* assoic;
                    if ~isempty(d)
                        IFP = model.infiniteFrequencyPermittivity(xy_loc, alpha0, vm_loc);
                        g = 3.22 - 0.94 .* Tr_loc(:,1);
                        theta = model.pressureAion(xy_loc, hyd);
                        SP = model.staticPermittivity(T_loc, xy_loc, theta, g, mu0, vm_loc, IFP);
                        Kappa = model.debyeLength(T_loc, xy_loc, charge, vm_loc, SP);
                        Kai = model.kaiFunction(d, Kappa, false);
                        Sigama = model.sigamaFunction(d, Kappa, false);
                        DIFPV = model.derivativeIFPV(IFP);
                        DSPV = model.derivativeSPV(DIFPV, SP, IFP);
                        PDH = model.pressureDH(T_loc, xy_loc(:,nmole+1:end), charge(nmole+1:end), Kai, Sigama, Kappa, SP, DSPV);
                        PB = model.pressureB(xy_loc(:,nmole+1:end), charge(nmole+1:end), SP, DSPV, vm_loc, rBorn);
                    else
                        PDH = 0; PB = 0;
                    end
                else
                    assoic = 0; PDH = 0; PB = 0;
                end
                cubic = [(p_loc - PDH - PB), (assoic - R .* T_loc), ...
                    - ((p_loc - PDH - PB) .* B_loc .^ 2 + R .* T_loc .* B_loc - A_loc),...
                    - A_loc .* B_loc - B_loc .^ 2 .* assoic];
                cubic = bsxfun(@rdivide, cubic, cubic(:,1));
                if model.selectGibbsMinimum
                    v0 = model.solveCubicEOS(cubic);
                    bad = bsxfun(@lt, v0, B_loc);
                    v0(bad) = nan;
                    
                    candidates = isfinite(v0);
                    numRoots = sum(candidates, 2);
                    single = numRoots == 1;
                    multiple = ~single;
                    Vm = max(v0, [], 2);
                    if any(multiple)
                        V_max = max(v0(multiple, :), [], 2);
                        V_min = min(v0(multiple, :), [], 2);
                        xi = xy_loc(multiple, :);
                        [~, phi] = model.computeFugacity(T_loc(multiple), p_loc(multiple), xi, A_loc(multiple), B_loc(multiple), Ai_loc(multiple, :), V_max, Si_loc(multiple, :), XA_loc(multiple, :), XC_loc(multiple), Tr_loc(multiple, :));
                        g_max = sum(phi .* xi, 2);
                        [~, phi] = model.computeFugacity(T_loc(multiple), p_loc(multiple), xi, A_loc(multiple), B_loc(multiple), Ai_loc(multiple, :), V_min, Si_loc(multiple, :), XA_loc(multiple, :), XC_loc(multiple), Tr_loc(multiple, :));
                        g_min = sum(phi .* xi, 2);
                        vi = V_max;
                        smallest = g_min < g_max;
                        vi(smallest) = V_min(smallest);
                        Vm(multiple) = vi;
                    end
                else
                    if isLiquid
                        Vm = computeLiquidZ(model, cubic, B_loc);
                    else
                        Vm = computeVaporZ(model, cubic, B_loc);
                    end
                end
                R_norm = sum(abs(vm_loc-Vm)./vm_loc,2);
                done = R_norm < 1e-6;
                keep = ~done;
                replace = active;
                replace(active) = done;
                XC(replace) = XC_loc(done);
                XA(replace, :) = XA_loc(done, :);
                vm(replace) = vm_loc(done);
                active(active) = ~done;
                
                if all(done)
                    break
                end
                xy_loc = xy_loc(keep, :);
                T_loc = T_loc(keep);
                B_loc = B_loc(keep);
                A_loc = A_loc(keep);
                vm_loc = Vm(keep);
                Tr_loc = Tr_loc(keep, :);
                p_loc = p_loc(keep);
                Ai_loc = Ai_loc(keep, :);
                Si_loc = Si_loc(keep, :);
            end
            if ~all(done)
                warning('Volume solver did not converge');
            end
            Z = p.*vm./(R.*T);
        end
        
        function v = computeLiquidZ(model, cubic, B)
            % Pick smallest v for liquid phase (least energy)
            v = model.solveCubicEOS(cubic);
            bad = bsxfun(@lt, v, B);
            v(bad) = nan;
            v = min(v, [], 2);
            model.checkv(v);
        end
        
        function v = computeVaporZ(model, cubic, B)
            % Pick largest v for vapor phase(most energy)
            v = model.solveCubicEOS(cubic);
            bad = bsxfun(@lt, v, B);
            v(bad) = nan;
            v = max(v, [], 2);
            model.checkv(v);
        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, iteration, varargin) %#ok
            % Compute a single step of the solution process for a given
            % equation of state (thermodynamic flash calculation);
            timer = tic();
            T = state.T;
            P = state.pressure;
            nc = numel(P);
            assert(all(P >= 0))
            
            z = state.components;
            ncomp = model.getNumberOfComponents();
            nmole = model.getNumberOfMolecules();
            
            % Basic assertions
            assert(all(sum(z, 2) > 0.999), ...
                'Molar fractions must sum up to unity.')
            assert(iteration > 0);
            [t_flash, t_stability] = deal(0);
            % If we get here, we are using a proper flash
            if iteration == 1
                % Book-keeping
                state.eos.itCount = zeros(nc, 1);
                state.eos.converged = false(nc, 1);
            end
            hasKvalues = ~isempty(model.equilibriumConstantFunctions);
            if hasKvalues
                state.K = model.evaluateEquilibriumConstants(P, T, z);
            end
            L0 = state.L;
            K0 = state.K;
            flash_method = lower(model.method);
            % Only apply calculations for cells that have not converged yet
            if iteration == 1
                x0 = state.x;
                y0 = state.y;
                % Cells which are two-phase can be flashed directly. For
                % single-phase cells, we perform a stability test in order
                % to check if the single-phase state is still stable.
                twoPhase = model.getTwoPhaseFlag(state);
                % NaN means that L was not initialized. We then start with
                % a phase stability test to get good estimates for the
                % K-values.
                twoPhase(isnan(L0)) = false;
                initSingle = ~twoPhase;
                stable = initSingle;
                % Call stability test.
                ts = tic();
                [stable(initSingle), x0(initSingle, :), y0(initSingle, :), K0(initSingle, :)] = ...
                    model.performPhaseStabilityTest(state.pressure(initSingle), state.T(initSingle),...
                    state.components(initSingle, :), state.K(initSingle, :));
                t_stability = toc(ts);
                [A_L, A_V, B_L, B_V, Si_L, Si_V, Ai, Tr] = model.getMixtureFugacityCoefficients(T, x0, y0);
                % Solve EOS for each phase
                Z0_L = model.computeCompressibilityZ(P, T, x0, A_L, B_L, Ai, Si_L, Tr, true);
                Z0_V = model.computeCompressibilityZ(P, T, y0, A_V, B_V, Ai, Si_V, Tr, false);
                L0(~stable) = model.solveRachfordRice(L0(~stable), K0(~stable, :), z(~stable, :));
                L0(stable) = model.singlePhaseLabel(P(stable), T(stable), z(stable, :));
                active = ~stable;
                % Flag stable cells as converged
                state.eos.converged(stable) = true;
            else
                Z0_L = state.Z_L;
                Z0_V = state.Z_V;
                x0 = state.x;
                y0 = state.y;
                active = ~state.eos.converged;
            end
            if hasKvalues
                % All cells have been solved already - we are done.
                active = active | false;
                values = zeros(1, ncomp);
                resConv = true(1, ncomp);
            elseif any(active)
                state.eos.itCount(active) = state.eos.itCount(active) + 1;
                K_init = K0(active, :);
                K = K_init;
                L = L0(active);
                z = z(active, :);
                P = P(active);
                T = T(active);
                switch flash_method
                    case 'newton'
                        % Newton-based solver for equilibrium
                        updatefn = @model.newtonCompositionUpdate;
                    case 'ssi'
                        % Successive substitution solver for equilibrium
                        if iteration > model.maxSSI
                            % Remaining cells are above SSI threshold,
                            % switch to Newton
                            updatefn = @model.newtonCompositionUpdate;
                            flash_method = 'newton';
                        else
                            % Just do SSI as asked
                            updatefn = @model.substitutionCompositionUpdate;
                        end
                    otherwise
                        error('Unknown flash method %s. Try ''newton'' or ''ssi''.', model.method);
                end
                timer_flash = tic();
                [x, y, K, Z_L, Z_V, L, equilvals] = updatefn(P, T, z, K, L);
                t_flash = toc(timer_flash);
                % Wei
                if iteration == 1
                    L = model.solveRachfordRice(L, K, z);                    
                end
                singlePhase = L == 0 | L == 1;
                % Single phase cells are converged
                equilvals(singlePhase, :) = 0;
                values = max(equilvals(:,1:nmole), [], 1);
                conv = max(equilvals(:,1:nmole), [], 2) <= model.nonlinearTolerance;
                conv = conv & iteration > nonlinsolve.minIterations;
                resConv = values <= model.nonlinearTolerance & iteration > nonlinsolve.minIterations;
                % Insert back the local values into global arrays
                state.eos.converged(active) = conv;
                % Insert updated values in active cells
                % Just single-phase state
                K(singlePhase, :) = K_init(singlePhase, :);
                x(singlePhase, :) = z(singlePhase, :);
                y(singlePhase, :) = y(singlePhase, :);
                Z_L = value(Z_L);
                Z_V = value(Z_V);
                liquid = L == 1; Z_V(liquid) = Z_L(liquid);
                vapor = L == 0; Z_L(vapor) = Z_V(vapor);
                
                L0(active) = L;
                Z0_L(active) = Z_L;
                Z0_V(active) = Z_V;
                K0(active, :) = K;
                x0(active, :) = x;
                y0(active, :) = y;
            else
                values = zeros(1, ncomp);
                resConv = true(1, ncomp);
            end
            state.K = K0;
            state.L = L0;
            state.x = x0;
            state.y = y0;
            state.Z_L = Z0_L;
            state.Z_V = Z0_V;
            
            failure = false;
            failureMsg = '';
            values_converged = values <= model.nonlinearTolerance;
            if model.verbose
                printConvergenceReport(model.ECPACompositionalMixture.names, values, values_converged, iteration);
            end
            report = model.makeStepReport(...
                'Failure',      failure, ...
                'FailureMsg',   failureMsg, ...
                'Converged',    all(values_converged), ...
                'ResidualsConverged', resConv, ...
                'Residuals',    values);
            report.ActiveCells = sum(active);
            if model.extraOutput
                report.ActiveFlag = active;
            end
            report.Method = flash_method;
            report.TotalTime = toc(timer);
            report.FlashTime = t_flash;
            report.StabilityTime = t_stability;
        end
                
        function L = singlePhaseLabel(eos, p, T, z)
            % Li's method for phase labeling
            Vc = eos.ECPACompositionalMixture.Vcrit;
            Tc = eos.ECPACompositionalMixture.Tcrit;
            n_m = sum(~isnan(Tc(1,:)));
            Vz = bsxfun(@times, Vc(:,1:n_m), z(:,1:n_m));
            
            Tc_est = sum(bsxfun(@times, Vz, Tc(:,1:n_m)), 2)./sum(Vz, 2);
            L = double(T < Tc_est);
        end
        
        function [x, y, K, Z_L, Z_V, L, values] = substitutionCompositionUpdate(model, P, T, z, K, L)
            % Determine overall liquid fraction
            L = model.solveRachfordRice(L, K, z);
            % Compute liquid component fraction
            [x, sx] = model.computeLiquid(L, K, z);
            % Vapor component fraction
            y = model.computeVapor(L, K, z);
            [Z_L, Z_V, f_L, f_V] = model.getCompressibilityAndFugacity(P, T, x, y, z, [], []);
            % Compute fugacity ratios
            f_r = bsxfun(@times, sx, f_L./f_V);
            % Update equilibrium constant estimates based on fugacity ratio
            mc = model.minimumComposition;
            % Ignore tiny compositions for the purpose of updating K-values
            % and estimating convergence.
            f_r(z <= mc) = 1;
            values = abs(f_r - 1);
            K = K.*abs(f_r);
            K(~isfinite(K)) = 1e-5;
            
            nmole = model.getNumberOfMolecules();
            K(:,nmole+1:end) = min(max(K(:,nmole+1:end), 1e-5), 1);
        end
        
        function [stable, x, y, K] = performPhaseStabilityTest(model, P, T, z, K)
            if nargin < 5
                K = [];
            end
            if isempty(z)
                stable = [];
                [x, y] = deal(zeros(0, size(z, 2)));
            else
                if isempty(model.equilibriumConstantFunctions)
                    % Use fugacity-based stability test
                    z = ensureMinimumFraction(z, model.minimumComposition);
                    [stable, x, y, K] = eCPAphaseStabilityTest(model, z, P, T, K);
                else
                    % We have K-values, use that instead
                    K = model.evaluateEquilibriumConstants(P, T, z);
                    hasStabFn = isprop(model.PropertyModel, 'checkStabilityFunction') && ...
                        ~isempty(model.PropertyModel.checkStabilityFunction);
                    if hasStabFn
                        [stable, L] = model.PropertyModel.checkStabilityFunction(P, T, expandMatrixToCell(z));
                    else
                        L = model.solveRachfordRice([], K, z);
                        L_tol = 1e-10;
                        stable = abs(L - 1) <= L_tol | L <= L_tol;
                    end
                    x = model.computeLiquid(L, K, z);
                    y = model.computeVapor(L, K, z);
                end
            end
        end
                
        function K = evaluateEquilibriumConstants(model, P, T, z)
            assert(~isempty(model.equilibriumConstantFunctions), ...
                'Equilibrium constants not defined. This function is only valid when using K-value formulation.');
            isCellInput = iscell(z);
            if ~isCellInput
                z = expandMatrixToCell(z);
            end
            K = cellfun(@(fn) fn(P, T, z), model.equilibriumConstantFunctions, 'UniformOutput', false);
            
            mv = 1e8;
            for i = 1:numel(K)
                K{i}(K{i} > mv) = mv;
                K{i}(K{i} < 1/mv) = 1/mv;
            end
            if ~isCellInput
                K = value(K);
            end
        end
        
        function [x, y, K, Z_L, Z_V, L, vals] = newtonCompositionUpdate(model, P, T, z, K, L)
            [x, y, K, Z_L, Z_V, L, vals] = eCPAnewtonFugacityEquilibrium(model, P, T, z, K, L);
        end
        
        function state = validateState(model, state)
            n_cell = size(state.pressure, 1);
            if ~isfield(state, 'L')
                % Add empty L to indicate that we do not know. The flash
                % routine will then perform the stability test for all
                % cells.
                state.L = nan(n_cell, 1);
            end
            
            if ~isfield(state, 'K')
                if isfield(state, 'x') && isfield(state, 'y')
                    K = state.y./state.x;
                else
                    K = eCPAestimateEquilibriumWilson(model, state.pressure, state.T);
                end
                state.K = K;
            end
            
            if ~isfield(state, 'Z_V')
                state.Z_V = ones(n_cell, 1);
            end
            if ~isfield(state, 'Z_L')
                state.Z_L = ones(n_cell, 1);
            end
            
            if ~isfield(state, 'x')
                state.x = state.components;
            end
            if ~isfield(state, 'y')
                state.y = state.components;
            end
        end
        
        function [Pr, Tr] = getReducedPT(model, P, T, useCell)
            if nargin < 4
                useCell = true;
            end
            Tc = model.ECPACompositionalMixture.Tcrit;
            Pc = model.ECPACompositionalMixture.Pcrit;
            if useCell
                n = model.CompositionalMixture.getNumberOfComponents();
                Tr = cell(1, n);
                Pr = cell(1, n);
                
                Ti = 1./Tc;
                Pi = 1./Pc;
                for i = 1:n
                    Tr{i} = T.*Ti(i);
                    Pr{i} = P.*Pi(i);
                end
            else
                Tr = bsxfun(@rdivide, T, Tc);
                Pr = bsxfun(@rdivide, P, Pc);
            end
        end
                
        function [A_L, A_V, B_L, B_V, Si_L, Si_V, Ai, Tr] = getMixtureFugacityCoefficients(model, T, x, y)
            % Calculate intermediate values for fugacity computation           
            [A_ij, Ai, Tr] = model.getMixingParameters(T, iscell(x));
            [Si_L, A_L, B_L] = model.getPhaseMixCoefficients(x, A_ij);
            [Si_V, A_V, B_V] = model.getPhaseMixCoefficients(y, A_ij);
        end
        
        function [Z_L, Z_V, f_L, f_V] = getCompressibilityAndFugacity(model, P, T, x, y, z, Z_L, Z_V, varargin)
            [A_L, A_V, B_L, B_V, Si_L, Si_V, Ai, Tr] = model.getMixtureFugacityCoefficients(T, x, y);
            R = 8.314462618;
            Bi = model.ECPACompositionalMixture.b;
            epsilonAB = model.ECPACompositionalMixture.epsilonAB;
            betaAB = model.ECPACompositionalMixture.betaAB;
            P_L = P; P_V = P;
            if isempty(Z_L)
                [~, v_L, XA_L, XC_L] = model.computeCompressibilityZ(P_L, T, x, A_L, B_L, Ai, Si_L, Tr, true);
            else
                v_L = Z_L .* R .* value(T) ./ value(P_L);
                [XA_L, XC_L] = model.pressureA(value(T), value(x), value(B_L), epsilonAB, betaAB, Bi, v_L, value(Tr));          
            end
            if isempty(Z_V)
                [~, v_V, XA_V ,XC_V] = model.computeCompressibilityZ(P_V, T, y, A_V, B_V, Ai, Si_V, Tr, false);
            else
                v_V = Z_V .* R .* value(T) ./ value(P_V);
                [XA_V, XC_V] = model.pressureA(value(T), value(y), value(B_V), epsilonAB, betaAB, Bi, v_V, value(Tr));  
            end
            [v_L, XA_L, XC_L] = model.setvDerivatives(T, P_L, x, A_L, B_L, v_L, XA_L, XC_L, Tr, varargin{:});
            [v_V, XA_V, XC_V] = model.setvDerivatives(T, P_V, y, A_V, B_V, v_V, XA_V, XC_V, Tr, varargin{:});
            Z_V = P_V .* v_V ./ (R .* T);
            Z_L = P_L .* v_L ./ (R .* T);
            f_L = model.computeFugacity(T, P_L, x, A_L, B_L, Ai, v_L, Si_L, XA_L, XC_L, Tr);
            f_V = model.computeFugacity(T, P_V, y, A_V, B_V, Ai, v_V, Si_V, XA_V, XC_V, Tr);
        end
        
        function [A_ij, Ai, Tr] = getMixingParameters(model, T, useCell)
            if nargin < 3
                useCell = true;
            end
            
            nc = numel(value(T));
            ncomp = model.getNumberOfComponents();
            nmole = model.getNumberOfMolecules();
            Tc = model.ECPACompositionalMixture.Tcrit;
            a0 = model.ECPACompositionalMixture.a0;
            c1 = model.ECPACompositionalMixture.c1;
            
            bic = model.ECPACompositionalMixture.getBinaryInteraction();
            if useCell
                Tr = cell(1, ncomp);
                Ai = cell(1, ncomp);
                Ti = 1./Tc;
                for i = 1:ncomp
                    if i <= nmole
                        Tr{i} = T.*Ti(i);
                        Ai{i} = a0(i).*(1+ c1(i).* (1- Tr{i}.^0.5)).^2;
                    else
                        Tr{i} = NaN(nc, 1);
                        Ai{i} = zeros(nc, 1);
                    end
                end
                % aij
                A_ij = cell(ncomp, ncomp);
                for i = 1:ncomp
                    for j = i:ncomp
                        A_ij{i, j} = (Ai{i}.*Ai{j}).^0.5.*(1 - bic(i, j));
                        A_ij{j, i} = A_ij{i, j};
                    end
                end
            else
                Tr = bsxfun(@rdivide, T, Tc);
                Ai = bsxfun(@times,a0, (1+bsxfun(@times, c1, 1- Tr.^0.5 )).^2);
                Ai(isnan(Ai)) = 0;
                A_ij = cell(ncomp, 1);
                for j = 1:ncomp
                    A_ij{j} = bsxfun(@times, bsxfun(@times, Ai, Ai(:, j)).^(1/2), 1 - bic(j, :));
                end
            end
        end
        
        function [Si, A, B] = getPhaseMixCoefficients(model, x, A_ij)
            [A, B] = deal(0);
            Bi = model.ECPACompositionalMixture.b;
            if iscell(x)
                ncomp = numel(x);
                Si = cell(1, ncomp);
                [Si{:}] = deal(0);
                % A = sum_ij x_i A_ij x_j
                % S_i = sum_j A_ij x_j
                % Note: A_ij is symmetric and we exploit that when we
                % evaluate the sums
                xjAij = cell(ncomp, ncomp);
                for i = 1:ncomp
                    B = B + x{i}.*Bi(i);
                    for j = 1:ncomp
                        result = A_ij{i, j}.*x{j};
                        xjAij{i, j} = result;
                        Si{i} = Si{i} + result;
                    end
                end
                for i = 1:ncomp
                    for j = i:ncomp
                        tmp = x{j}.*xjAij{j, i};
                        if i == j
                            A = A + tmp;
                        else
                            A = A + 2*tmp;
                        end
                    end
                end
            else
                ncomp = size(x, 2);
                Si = zeros(size(x));
                B = sum(x.*Bi, 2);
                A = 0;
                for i = 1:ncomp
                    A = A + sum(bsxfun(@times, bsxfun(@times, A_ij{i}, x(:, i)), x), 2);
                    Si = Si + bsxfun(@times, A_ij{i}, x(:, i));
                end
            end
        end
              
        function [f, LOGphi] = computeFugacity(model, T, P, xy, A, B, Ai, vm, Si, XAx, XCx, Tr)
            Bi = model.ECPACompositionalMixture.b;
            alpha0 = model.ECPACompositionalMixture.alpha0;
            rBorn = model.ECPACompositionalMixture.rborn;
            rBorn = rBorn(~isnan(rBorn));
            mu0 = model.ECPACompositionalMixture.mu0;
            mu0 = mu0(~isnan(mu0));
            d = model.ECPACompositionalMixture.d;
            d = d(~isnan(d));
            charge = model.ECPACompositionalMixture.charge;
            stoichiometric = flip(abs(charge));
            t_sto = sum(stoichiometric);
            hyd = model.ECPACompositionalMixture.hyd;
            ncomp = numel(Bi);
            nion = numel(d);
            nmole = ncomp - nion;
            nc = numel(T);
            
            FuSRK = model.srkFu(T, P, A, B, vm, Ai, Bi, Si);
            if ~isempty(mu0)
                FuAsso = model.assocWaterFu(xy, B, vm, Bi, XAx, XCx);
                if ~isempty(d)
                    IFP = model.infiniteFrequencyPermittivity(xy, alpha0, vm);
                    if iscell(xy)
                        g = 3.22 - 0.94 .* Tr{1};
                    else
                        g = 3.22 - 0.94 .* Tr(:,1);
                    end
                    theta = model.pressureAion(xy, hyd);
                    SP = model.staticPermittivity(T, xy, theta, g, mu0, vm, IFP);
                    Kappa = model.debyeLength(T, xy, charge, vm, SP);
                    Sigama = model.sigamaFunction(d, Kappa, iscell(xy));
                    DIFPn = model.derivativeIFPn(IFP, vm, alpha0, iscell(xy));
                    DSPns = model.derivativeSPns(DIFPn, SP, IFP, T, vm, g, mu0, theta, iscell(xy));
                    MuDHs = model.chemicalPotentialDHs(T, xy, charge, vm, Sigama, Kappa, SP, DSPns);
                    MuBs = model.chemicalPotentialBs(T, xy, charge, SP, DSPns, rBorn);
                    Kai = model.kaiFunction(d, Kappa, iscell(xy));
                    DSPni = model.derivativeSPni(DIFPn, SP, IFP);
                    MuDHi = model.chemicalPotentialDHi(T, xy, charge, vm, Sigama, Kappa, SP, DSPni, Kai);
                    MuBi = model.chemicalPotentialBi(T, xy, charge, SP, DSPni, rBorn);
                else
                    if iscell(xy)
                        [MuDHs, MuBs] = deal(cell(1, nmole));
                        for i = 1 : nmole
                            MuDHs{i} = zeros(nc,1);
                            MuBs{i} = zeros(nc,1);
                        end
                        [MuDHi, MuBi] = deal(cell(1, nion));
                        for i = 1 : nion
                            MuDHi{i} = zeros(nc,1);
                            MuBi{i} = zeros(nc,1);
                        end
                    else
                        MuDHs(1:nmole) = 0; MuDHi(1:nion) = 0; MuBs(1:nmole) = 0; MuBi(1:nion) = 0;
                    end
                end
            else
                if iscell(xy)
                    FuAsso = cell(1, ncomp);
                    for i = 1 : ncomp
                        FuAsso{i} = zeros(nc,1);
                    end
                    [MuDHs, MuBs] = deal(cell(1, nmole));
                    for i = 1 : nmole
                        MuDHs{i} = zeros(nc,1);
                        MuBs{i} = zeros(nc,1);
                    end
                    [MuDHi, MuBi] = deal(cell(1, nion));
                    for i = 1 : nion
                        MuDHi{i} = zeros(nc,1);
                        MuBi{i} = zeros(nc,1);
                    end
                else
                    FuAsso(1:ncomp) = 0; MuDHs(1:nmole) = 0; MuDHi(1:nion) = 0; MuBs(1:nmole) = 0; MuBi(1:nion) = 0;
                end
            end
            if iscell(xy)
                [f, phi, LOGphi] = deal(cell(1, ncomp));
                for i = 1 : ncomp
                    if i<=nmole
                        phi{i} = exp(FuSRK{i} + FuAsso{i} + MuDHs{i} + MuBs{i});
                    else
                        phi{i} = exp(FuSRK{i} + FuAsso{i} + MuDHi{i-nmole} + MuBi{i-nmole}./5);
                    end
                end
                if nion > 0
                    phi_ion = (phi{nmole+1}.^stoichiometric(1).*phi{nmole+2}.^stoichiometric(2)).^(1/t_sto);
                    [phi{nmole+1}, phi{nmole+2}] = deal(phi_ion);
                end
                for i = 1:ncomp
                    f{i} = phi{i} .* P .* xy{i};
                    LOGphi{i} = log(phi{i});
                end
            else
                phi = exp(FuSRK + FuAsso + [MuDHs, MuDHi] + [MuBs, MuBi]);
                if nion>0
                    phi_ion = (phi(:,nmole+1).^stoichiometric(1).*phi(:,nmole+2).^stoichiometric(2)).^(1/t_sto);
                    [phi(:,nmole+1), phi(:,nmole+2)] = deal(phi_ion);
                end
                f = phi .* bsxfun(@times, P, xy);
                LOGphi = log(phi);
            end
        end
        
        function [eqs, f_L, f_V, Z_L, Z_V] = equationsEquilibrium(model, P, T, x, y, z, L, Z_L, Z_V)
            useFugacity = isempty(model.equilibriumConstantFunctions);
            ncomp = numel(x);
            if useFugacity
                [Z_L, Z_V, f_L, f_V] = model.getCompressibilityAndFugacity(P, T, x, y, z, Z_L, Z_V);
            else
                [Z_L, Z_V] = model.getCompressibilityAndFugacity(P, T, x, y, z, Z_L, Z_V);
                f_L = cell(1, ncomp);
                f_V = cell(1, ncomp);
                for i = 1:ncomp
                    K = model.equilibriumConstantFunctions{i}(P, T, z);
                    f_L{i} = K.*x{i};
                    f_V{i} = y{i};
                end
            end
            eqs = cell(1, 2*ncomp + 1);
            isLiq = value(L) == 1;
            isVap = value(L) == 0;
            isPure = isLiq | isVap;
            eqs{end} = zeros(numelValue(L), 1);
            nmole = model.getNumberOfMolecules();
            for i = 1:ncomp
                eqs{i} = z{i} - L.*x{i} - (1-L).*y{i};
                if i <= nmole
                    eqs{i+ncomp} = (f_V{i} - f_L{i});
                else
                    scale = value(f_V{i})./value(f_L{i});
                    eqs{i+ncomp} = (f_V{i}- f_L{i}.*scale );
                end
                eqs{end} = eqs{end} - ~isLiq.*y{i} + ~isVap.*x{i};
                if any(isPure)
                    % Treat pure phase behavior
                    eqs{i}(isPure) = z{i}(isPure) - x{i}(isPure);
                    eqs{i+ncomp}(isPure) = z{i}(isPure) - y{i}(isPure);
                    if any(isLiq)
                        eqs{end}(isLiq) = L(isLiq) - 1;
                    end
                    if any(isVap)
                        eqs{end}(isVap) = L(isVap);
                    end
                end
            end
        end
        
        function [XAxy, XCxy, gxy] = pressureA(model, T, xy, B, epsilonAB, betaAB, Bi, vm, Tr)
            nmole = model.getNumberOfMolecules();
            nc = numel(T);
            
            rouxy = bsxfun(@rdivide, 1, vm);
            etaxy = 0.25 * rouxy .* B;
            gxy = bsxfun(@rdivide, 1, 1 - bsxfun(@times, 1.9, etaxy));
            dexy(:,1) = bsxfun(@times, gxy .* (exp(bsxfun(@rdivide, epsilonAB(1), T)) - 1), Bi(1) * betaAB(1));
            Xcxy = (-1+(1+8.*rouxy.*xy(:,1).*dexy(:,1)).^0.5)./(4*rouxy.*xy(:,1).*dexy(:,1));
            
            bxy = 1;
            XAxy = ones(nc, nmole);
            ap = model.ECPACompositionalMixture.getAssociationParameter();
            for i = 2 : nmole
                dexy(:,i) = dexy(:,1) .* ap(1,i);
                XAxy(:,i) = 1./(1+2.*rouxy.*xy(:,1).*Xcxy.*dexy(:,i));
                bxy = bxy+rouxy.*xy(:,i).*XAxy(:,i).*dexy(:,i);
            end
            XCxy=(-bxy+(bxy.^2+8*bxy.*rouxy.*xy(:,1).*dexy(:,1)).^0.5)./(4*bxy.*rouxy.*xy(:,1).*dexy(:,1));
            XAxy(:,1) = 1./(1+2*rouxy.*xy(:,1).*XCxy.*dexy(:,1));
        end
        
        function  IFP = infiniteFrequencyPermittivity(model, xy, alpha0, v_m)
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            if iscell(xy)
                ncomp = model.getNumberOfComponents();
                a0 = 0;
                for i = 1:ncomp
                    a0 = a0 + xy{i} .* alpha0(i);
                end
                MP = NA .* a0 ./ (3 .* VP .* v_m);
                IFP = (2 .* MP + 1) ./ (1 - MP);
            else
                MP = NA .* sum(bsxfun(@times, xy, alpha0), 2) ./ (3 .* VP .* v_m);
                IFP = (2 .* MP + 1) ./ (1 - MP);
            end
        end
        
        function theta = pressureAion(model, xy, hyd)
            if iscell(xy)
                a0 = 0;
                for i = 1:numel(hyd)
                    a0 = a0 + xy{i} .* hyd(i);
                end
                b = 1 - xy{1} + a0;
                theta = (-b + (b.^2 + 4*xy{1}).^0.5) ./ (2 * xy{1});
            else
                b = 1 - xy(:,1) + sum(bsxfun(@times, xy, hyd), 2);
                theta = (-b + sqrt(b.^2 + 4 .* xy(:,1))) ./ (2 .* xy(:,1));
            end
        end
        
        function  SP = staticPermittivity(model, T, xy, theta, g, mu0, v_m, IFP)
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            kB = 1.380649e-23;
            if iscell(xy)
                VDM = IFP + NA .* ((IFP + 2) ./ 3).^2 .* ...
                    xy{1}.*theta .* g .* mu0.^2 ./ (VP .* kB .* T .* v_m);
                SP=(VDM + (VDM.^2 + 8 .* IFP.^2).^0.5) ./ 4;
            else
                VDM = IFP + NA .* ((IFP + 2) ./ 3).^2 .* ...
                    sum(bsxfun(@times, bsxfun(@times, xy(:,1), theta), bsxfun(@times, g, mu0.^2)) ,2) ...
                    ./ (VP .* kB .* T .* v_m);
                SP=(VDM + (VDM.^2 + 8 .* IFP.^2).^0.5) ./ 4;
            end
        end
        
        function  Kappa = debyeLength(model, T, xy, charge, v_m, SP)
            EC = 1.602176634e-19;
            NA = 6.02214076e23;
            kB = 1.380649e-23;
            VP = 8.8541878128e-12;
            if iscell(xy)
                ncomp = model.getNumberOfComponents();
                nion = model.getNumberOfIons();
                nmole = ncomp - nion;
                a0 = 0;
                for i = nmole+1:ncomp
                    a0 = a0 + xy{i} .* charge(i).^2;
                end
                Kappa = (EC.^2 .* NA .* a0 ./ (kB .* T .* VP .* SP .* v_m)).^0.5;
            else
                Kappa = (EC.^2 .* NA .* sum(bsxfun(@times, xy, charge.^2), 2) ...
                    ./ (kB .* T .* VP .* SP .* v_m)).^0.5;
            end
        end
        
        function  Kai = kaiFunction(model, d, Kappa, iscell)
            if iscell
                ncomp = numel(d);
                Kai = cell(1, ncomp);
                for i = 1 : ncomp
                    Kai{i} = (log(1 + Kappa .* d(i)) - Kappa .* d(i) ...
                        + 0.5*(Kappa .* d(i)).^2) ./ (d(i).^3);
                end
            else
                Kai = (log(1 + bsxfun(@times, Kappa, d)) - bsxfun(@times, Kappa, d) ...
                    + 0.5 .* bsxfun(@times, Kappa, d) .^ 2) ./ (d .^ 3);
            end
        end
        
        function  Ksi = ksiFunction(model, d, Kappa, iscell)
            if iscell
                ncomp = numel(d);
                Ksi = cell(1, ncomp);
                for i = 1 : ncomp
                    Ksi{i} =(3*Kappa.^2+2*Kappa.^3 .*d) ./(1+Kappa .*d).^2;
                end
            else
                Ksi = bsxfun(@rdivide, 3*Kappa.^2+bsxfun(@times, 2*Kappa.^3, d), (1+bsxfun(@times, Kappa, d)).^2);
            end
        end
        
        function  Sigama = sigamaFunction(model, d, Kappa, iscell)
            if iscell
                ncomp = numel(d);
                Sigama = cell(1, ncomp);
                for i = 1 : ncomp
                    Sigama{i} =  Kappa.^2 ./ (1 + Kappa .* d(i));
                end
            else
                Sigama = bsxfun(@rdivide, Kappa .^2, 1 + bsxfun(@times, Kappa, d));
            end
        end
        
        function  DIFPV = derivativeIFPV(model, IFP)
            DIFPV = bsxfun(@rdivide, -(IFP + 2) .* (IFP - 1), 3);
        end
        
        function  DSPV = derivativeSPV(model, DIFPV, SP, IFP)
            DSPV = ((IFP + 2) .* SP) .^ 2 ./ (2 .* SP .^ 2 + IFP .^ 2) ...
                .* ((4 .* IFP + (2 + 4 .* SP - IFP) .* SP) ./ ((IFP + 2) .^ 3 .* SP) .* DIFPV ...
                - (2 .* SP + IFP) .* (SP - IFP) ./ (SP .* (IFP + 2) .^ 2));
        end
        
        function  PDH = pressureDH(model, T, xy, charge, Kai, Sigama, Kappa, SP, DSPV)
            kB = 1.380649e-23;
            PDH = kB .* T .* sum(bsxfun(@times, xy, charge .^ 2) .* Kai, 2) ...
                ./ (4 .* pi .* sum(bsxfun(@times, xy, charge .^ 2), 2)) ...
                - kB .* T .* sum(bsxfun(@times, xy, charge .^ 2) .* Sigama, 2)...
                ./ (8 .* pi .* sum(bsxfun(@times, xy, charge .^ 2), 2)) ...
                .* (Kappa + Kappa ./ SP .* DSPV);
        end
        
        function  PB = pressureB(model, xy, charge, SP, DSPV, v_m, rBorn)
            EC = 1.602176634e-19;
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            PB = bsxfun(@rdivide, NA .* EC .^ 2, (8 .* pi .* VP .* SP .^ 2)) ...
                .* sum(bsxfun(@rdivide, bsxfun(@times, xy, charge .^ 2), rBorn), 2) .* DSPV ./ v_m;
        end
        
        function Fu = srkFu(model, T, P, A, B, vm, Ai, Bi, Si)
            R = 8.314462618;
            if iscell(Ai)
                ncomp = model.getNumberOfComponents();
                Fu = deal(cell(1, ncomp));
                f1 = 1./(vm-B) - A./(R.*T.*B.*(vm+B));
                f2 = - log(P .* (vm - B) ./ (R .* T));
                f3 = log((vm+B)./vm);
                for i = 1:ncomp
                    Fu{i} = f2 + Bi(i).*f1 - A./(R.*T.*B).*(2.*Si{i}./A - Bi(i)./B).*f3;
                end
            else
                Fu = bsxfun(@times, (bsxfun(@rdivide, 1, (vm - B)) ...
                    - A ./ (B .* R .* T .* (vm + B))), Bi) ...
                    - log(P .* (vm - B) ./ (R .* T)) ...
                    - A./(B.*R.*T).*(2.*Si./A - Bi./B) .*log((vm + B) ./ vm);
            end
        end
        
        function Fu = assocWaterFu(model, xy, B, vm, Bi, XAx, XCx)
            ncomp = model.getNumberOfComponents();
            nion = model.getNumberOfIons();
            nmole = ncomp - nion;
            if iscell(xy)
                XAx = expandMatrixToCell(XAx);
                assox = deal(cell(1, numel(xy)));
                assox{1} = 2 .* log(XAx{1} .* XCx);
                Fu = deal(cell(1, ncomp));
                f1 = 1.9./(4*vm-1.9*B);
                f2 = 2 .* xy{1} .* (2 - XAx{1} - XCx);
                for i = 2:nmole
                    assox{i} = log(XAx{i});
                    f2 = f2 + xy{i} .* (1 - XAx{i});
                end
                for i = 1:ncomp
                    if i > nmole
                        assox{i} = 0;
                    end
                    Fu{i} = assox{i}-0.5.*f1.*Bi(i).*f2;
                end
            else
                assox = zeros(size(value(xy)));
                assox(:,1) = 2 .* log(XAx(:,1) .* XCx(:,1));
                assox(:,2:nmole) = log(XAx(:,2:nmole));
                asso = 2 .* xy(:,1) .* (2 - XAx(:,1) - XCx(:,1))...
                    + sum(xy(:,2:nmole) .* (1 - XAx(:,2:nmole)), 2);
                Fu = assox - 0.5 .* 1.9 .* bsxfun(@rdivide, bsxfun(@times, Bi, asso), 4 .* vm - 1.9 .* B);
            end
        end
        
        function  DIFPn = derivativeIFPn(model, IFP, v_m, alpha0, iscell)
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            if iscell
                ncomp = numel(alpha0);
                DIFPn = cell(1, ncomp);
                for i = 1 : ncomp
                    DIFPn{i} = ((IFP + 2)./3).^2.*NA./(VP.*v_m) .* alpha0(i);
                end
            else
                DIFPn = bsxfun(@times, ((IFP + 2) ./ 3) .^ 2 .* NA ./ (VP .* v_m), alpha0);
            end
        end
        
        function  DSPns = derivativeSPns(model, DIFPn, SP, IFP, T, v_m, g, mu0, theta, iscell)
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            kB = 1.380649e-23;
            ncomp = model.getNumberOfComponents();
            nion = model.getNumberOfIons();
            nmole = ncomp - nion;
            if iscell
                DSPns = cell(1, nmole);
                for i = 1 : nmole
                    DSPns{i} = ((IFP + 2) .* SP).^2 ./ (2*SP.^2 + IFP.^2) .* ...
                        ((4*IFP + (2 + 4*SP - IFP) .* SP) ./ ((IFP + 2).^3 .* SP) ...
                        .* DIFPn{i} + NA .* (g .* mu0 .^ 2 .* theta) ...
                        ./ (9 .* VP .* kB .* T .* v_m));
                end
            else
                DSPns = bsxfun(@times, ((IFP + 2) .* SP) .^ 2 ./ (2 .* SP .^ 2 + IFP .^ 2), ...
                    (bsxfun(@times, (4 .* IFP + (2 + 4 .* SP - IFP) .* SP) ...
                    ./ ((IFP + 2) .^ 3 .* SP), DIFPn(:,1:nmole)) + NA .* (g .* mu0 .^ 2 .* theta) ...
                    ./ (9 .* VP .* kB .* T .* v_m) ));
            end
        end
        
        function  MuDHs = chemicalPotentialDHs(model, T, xy, charge, v_m, Sigama, Kappa, SP, DSPns)
            kB = 1.380649e-23;R = 8.314462618;
            ncomp = model.getNumberOfComponents();
            nion = model.getNumberOfIons();
            nmole = ncomp - nion;
            if iscell(xy)
                MuDHs = cell(1, nmole);
                d1 = 0; d2 = 0;
                for i = nmole+1 : ncomp
                    d1 = d1 + xy{i} .* charge(i).^2;
                    d2 = d2 + xy{i} .* charge(i).^2 .* Sigama{i-nmole};
                end
                for i = 1 : nmole
                    MuDHs{i} = kB .* T .* v_m ./ (8 .* pi .*  d1) ...
                        .* d2 .* Kappa ./ SP .* DSPns{i} ./ (R .* T);
                end
            else
                MuDHs = bsxfun(@times, kB .* T .* v_m ./ (8 .* pi .* sum(bsxfun(@times, xy(:,1+nmole:end), charge(:,1+nmole:end) .^ 2), 2)) ...
                    .* sum(bsxfun(@times, xy(:,1+nmole:end) .* Sigama, charge(:,1+nmole:end) .^ 2), 2) .* Kappa ./ SP, DSPns)./ (R .* T);
            end
        end
        
        function  MuBs = chemicalPotentialBs(model, T, xy, charge, SP, DSPns, rBorn)
            EC = 1.602176634e-19;R = 8.314462618;
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            ncomp = model.getNumberOfComponents();
            nion = model.getNumberOfIons();
            nmole = ncomp - nion;
            if iscell(xy)
                MuBs = cell(1, nmole);
                d1 = 0;
                for i = nmole+1 : ncomp
                    d1 = d1 + xy{i} .* charge(i).^2 ./ rBorn(i-nmole);
                end
                for i = 1 : nmole
                    MuBs{i} = - NA*EC^2 ./ (8*pi*VP .* SP.^2) .* d1 .* DSPns{i} ./ (R*T);
                end
            else
                MuBs = -bsxfun(@times, NA .* EC .^ 2 ./ (8 .* pi .* VP .* SP .^ 2) ...
                    .* sum(bsxfun(@times, xy(:,1+nmole:end), charge(:,1+nmole:end) .^2 ./ rBorn), 2), DSPns)./ (R .* T);
            end
        end
        
        function  DSPni = derivativeSPni(model, DIFPn, SP, IFP)
            ncomp = model.getNumberOfComponents();
            nion = model.getNumberOfIons();
            nmole = ncomp - nion;
            if iscell(DIFPn)
                DSPni = cell(1, nion);
                for i = 1 : nion
                    DSPni{i} = ((IFP + 2) .* SP).^2 ./ (2*SP.^2 + IFP.^2) .* ...
                        (4*IFP + (2 + 4*SP - IFP) .* SP) ./ ((IFP + 2).^3 .* SP) ...
                        .* DIFPn{i+nmole};
                end
            else
                DSPni = bsxfun(@times, ((IFP + 2) .* SP) .^ 2 ./ (2 .* SP .^ 2 + IFP .^ 2), ...
                    bsxfun(@times, (4 .* IFP + (2 + 4 .* SP - IFP) .* SP) ./ ((IFP + 2) .^ 3 .* SP), DIFPn(:,nmole+1:end)));
            end
        end
        
        function  MuDHi = chemicalPotentialDHi(model, T, xy, charge, v_m, Sigama, Kappa, SP, DSPni, Kai)
            kB = 1.380649e-23;R = 8.314462618;
            ncomp = model.getNumberOfComponents();
            nion = model.getNumberOfIons();
            nmole = ncomp - nion;
            if iscell(xy)
                MuDHi = cell(1, nion);
                d1 = 0; d2 = 0;d3=0;
                for i = nmole+1 : ncomp
                    d1 = d1 + xy{i} .* charge(i).^2 .* Kai{i-nmole};
                    d2 = d2 + xy{i} .* charge(i).^2;
                    d1 = d1 + xy{i} .* charge(i).^2 .* Sigama{i-nmole};
                end
                for i = 1 : nion
                    MuDHi{i} = (kB .* T .* v_m .* d1 ./ (4*pi*d2.^2) .* charge(i+nmole)^2 ...
                        - kB .* T .* v_m .* charge(i+nmole).^2 ./ (4*pi .* d2) ...
                        .* (Kai{i} + 0.5*Kappa .* d3 ./ d2) ...
                        +  kB*T .* v_m ./ (8*pi .* d2) .* d3 .* Kappa ./ SP .* DSPni{i})./(R*T);
                end
            else
                MuDHi = kB*T .* v_m .* sum(bsxfun(@times, xy(:,nmole+1:end), charge(nmole+1:end) .^ 2) .* Kai, 2) ...
                    ./ bsxfun(@rdivide,(4*pi* (sum(bsxfun(@times, xy(:,nmole+1:end), charge(nmole+1:end).^2), 2)).^2), charge(nmole+1:end).^2) ...
                    - bsxfun(@times, kB .* T .* v_m, charge(nmole+1:end) .^ 2) ./ (4 .* pi .* sum(bsxfun(@times, xy(:,nmole+1:end), charge(nmole+1:end) .^ 2), 2)) ...
                    .* (Kai + 0.5 .* Kappa .* sum(bsxfun(@times, xy(:,nmole+1:end), charge(nmole+1:end) .^ 2) .* Sigama, 2) ./ sum(bsxfun(@times, xy(:,nmole+1:end), charge(nmole+1:end) .^ 2), 2)) ...
                    + bsxfun(@times, kB .* T .* v_m ./ (8 .* pi .* sum(bsxfun(@times, xy(:,nmole+1:end), charge(nmole+1:end) .^ 2), 2)) ...
                    .* sum(bsxfun(@times, xy(:,nmole+1:end), charge(nmole+1:end) .^ 2) .* Sigama, 2) .* Kappa ./ SP, DSPni);
                MuDHi = MuDHi ./ (R .* T);
            end
        end
        
        function  MuBi = chemicalPotentialBi(model, T, xy, charge, SP, DSPni, rBorn)
            EC = 1.602176634e-19;
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            R = 8.314462618;
            ncomp = model.getNumberOfComponents();
            nion = model.getNumberOfIons();
            nmole = ncomp - nion;
            if iscell(xy)
                MuBi = cell(1, nion);
                d1 = 0;
                for i = nmole+1 : ncomp
                    d1 = d1 + xy{i} .* charge(i).^2 ./ rBorn(i-nmole);
                end
                for i = 1 : nion
                    MuBi{i} = (NA*EC^2.*charge(i+nmole).^ 2 ./ (8*pi*VP*rBorn(i)) .* (1./SP - 1) ...
                        - NA*EC^2 ./ (8*pi*VP.*SP.^2) .* d1 .* DSPni{i}) ./ (R .* T);
                end
            else
                MuBi = bsxfun(@times, NA .* EC .^ 2 .* charge(nmole+1:end) .^ 2 ./ (8 .* pi .* VP .* rBorn), (1 ./ SP - 1)) ...
                    - bsxfun(@times, NA .* EC .^2 ./ (8 .* pi .* VP .* SP .^ 2) .* sum(bsxfun(@times, xy(:,nmole+1:end), charge(nmole+1:end) .^ 2 ./ rBorn), 2), DSPni);
                MuBi = MuBi ./ (R .* T);
            end
        end
        
        function [v, XA, XC]= setvDerivatives(model, T, P, x, A, B, v, XA, XC, Tr, cellJacMap)
            S = getSampleAD(A, B);
            if isnumeric(S)
                % This function does nothing if A, B are both doubles.
                return;
            end
            
            ncomp = model.getNumberOfComponents();
            nmole = model.getNumberOfMolecules();
            XA = expandMatrixToCell(XA);
            
            backend = model.AutoDiffBackend;
            isDiagonal = isa(backend, 'DiagonalAutoDiffBackend');
            if isnumeric(v)
                v = backend.convertToAD(v, S);
                XC = backend.convertToAD(XC, S);
            end
            for k = 1:nmole
                XA{k} = backend.convertToAD(XA{k}, S);
            end
            if nargin < 11
                if isDiagonal
                    cellJacMap = cell(v.offsets(end)-1, 1);
                else
                    cellJacMap = cell(numel(v.jac), 1);
                end
            end
            
            am = value(A);
            bm = value(B);
            vm = value(v);
            xm = value(x);
            if iscell(xm)
                xm = [xm{:}];
            end
            T = value(T);
            Tr = value(Tr);
            XAm = value(XA);
            XCm = value(XC);
            
            alpha0 = model.ECPACompositionalMixture.alpha0;
            rBorn = model.ECPACompositionalMixture.rborn;
            rBorn = rBorn(~isnan(rBorn));
            mu0 = model.ECPACompositionalMixture.mu0;
            mu0 = mu0(~isnan(mu0));
            d = model.ECPACompositionalMixture.d;
            d = d(~isnan(d));
            charge = model.ECPACompositionalMixture.charge;
            hyd = model.ECPACompositionalMixture.hyd;
            b = model.ECPACompositionalMixture.b;
            epsilonAB = model.ECPACompositionalMixture.epsilonAB;
            betaAB = model.ECPACompositionalMixture.betaAB;

            
            IFP = model.infiniteFrequencyPermittivity(xm, alpha0, vm);
            g = 3.22 - 0.94 .* Tr(:,1);
            theta = model.pressureAion(x, hyd); % ADI
            thetam = value(theta);
            SP = model.staticPermittivity(T, xm, thetam, g, mu0, vm, IFP);
            Kappa = model.debyeLength(T, xm, charge, vm, SP);
            Sigama = model.sigamaFunction(d, Kappa, false);
            Ksi = model.ksiFunction(d, Kappa, false);
            Kai = model.kaiFunction(d, Kappa, false);
            rou = 1./vm;
            eta = 0.25 *rou .* bm;
            ga = 1./(1 - 1.9*eta);
            
            de(:,1) = ga .* (exp(epsilonAB(1)./ T) - 1).*b(1) * betaAB(1);
            ap = model.ECPACompositionalMixture.getAssociationParameter();
            for i = 2 : nmole
                de(:,i) = de(:,1) .* ap(1,i);
            end
            
            dx = cell(1,ncomp);
            dXA = cell(1, nmole);
            if isDiagonal
                offset = v.offsets;
                rowMajor = backend.rowMajor;
                if rowMajor
                    T = T';
                    xm = xm';
                    ga = ga';
                    de = de';
                    vm = vm';
                    am = am';
                    bm = bm';
                    XAm = XAm';
                    XCm = XCm';
                    IFP = IFP';
                    g = g';
                    SP = SP';
                    thetam = thetam';
                    charge = charge';
                    Kappa = Kappa';
                    Sigama = Sigama';
                    Ksi = Ksi';
                    Kai = Kai';
                    rBorn = rBorn';
                    Tr = Tr';
                end
                for i = 1:numel(v.jac)
                    act = offset(i):offset(i+1)-1;
                    map = cellJacMap(act);
                    map = map(~cellfun(@isempty, map));
                    if isempty(map)
                        if v.jac{i}.dim(1) ~= numel(vm) && numel(v.jac{i}.subset) ~= numel(vm)
                            continue
                        end
                        % All cells
                        dA = model.getDiagonalJac(A, i);
                        dB = model.getDiagonalJac(B, i);
                        dtheta = model.getDiagonalJac(theta, i);
                        for j = 1 : ncomp
                            dx{j} = model.getDiagonalJac(x{j}, i);
                        end
                        [row, column] = size(dA);
                        if isa(P, 'ADI')
                            dP = P.jac{i}.diagonal;
                        else
                            dP = zeros(row, column);
                        end
                        
                        [psrk, s2] = model.setSRKDerivatives(T, vm, am, bm, dA, dB);
                        if ~isempty(d)
                            [pdh, d40, d23, d26, d27, d3, d8] = model.setDHDerivatives(T, xm, ...
                                vm, SP, IFP, dx, alpha0,thetam,dtheta,...
                                g,mu0,charge,Sigama,Kappa,Kai,Ksi);
                            [pb, b9] = model.setBornDerivatives(xm, vm, ...
                                SP, charge, rBorn, dx, d23, d26, d27, d3, d8);
                        else
                            if rowMajor
                                d40 = zeros(1,column);
                                b9 = zeros(1,column);
                                pdh = zeros(row,column);
                                pb = zeros(row,column);
                            else
                                d40 = zeros(row,1);
                                b9 = zeros(row,1);
                                pdh = zeros(row,column);
                                pb = zeros(row,column);
                            end
                        end
                        
                        if ~isempty(mu0)
                            [pa, pa1, a11, a1, a2] = model.setAssoc1Derivatives(T, xm, ...
                                vm, XAm, XCm, ga, de, dx, dB);
                            
                            if isempty(dP) && ~isempty(psrk)
                                dP = zeros(row, column);
                            elseif isempty(psrk)
                                dP = zeros(row, column);
                            end
                            
                            dv_loc = (dP-(psrk+pa+pdh+pb)) ./ (s2+a11+d40+b9);
                            if rowMajor
                                active = true(column, 1);
                                activeRow = active';
                                T_loc = T(activeRow);
                                Tr_loc = Tr(:,activeRow);
                                xm_loc = xm(:,activeRow);
                                a1_loc = a1(activeRow);
                            else
                                active = true(row, 1);
                                T_loc = T(active);
                                Tr_loc = Tr(active, :);
                                xm_loc = xm(active, :);
                                a1_loc = a1(active);
                            end
                            B_loc = B(active);
                            x_loc = cellfun(@(x) x(active), x, 'UniformOutput', false);
                            v_loc = v(active);
                            XA_loc = cellfun(@(x) x(active), XA, 'UniformOutput', false);
                            XC_loc = XC(active);
                            for itea = 1:10000
                                % t_outer = tic;

                                if v_loc.jac{i}.isZero
                                    v_loc.jac{i} = v_loc.jac{i}.expandZero();
                                end
                                v_loc.jac{i}.diagonal = dv_loc;
                                if all(v_loc.jac{i}.subset == 0)
                                    v_loc.jac{i}.subset = [];
                                end
                                [XA_loc, XC_loc] = setXDerivatives(model, T_loc, x_loc, B_loc, epsilonAB, betaAB, b, v_loc, Tr_loc, XA_loc, XC_loc, i);
                                for j = 1 : nmole
                                    dXA{j} = XA_loc{j}.jac{i}.diagonal;
                                end
                                dXC = XC_loc.jac{i}.diagonal;
                                pa3 = model.setAssoc2Derivatives(T_loc, ...
                                    xm_loc, a1_loc, dXA, dXC);
                                
                                dv_next = (dP-(psrk+pa1+pa3+pdh+pb)) ./ (s2+a2+d40+b9);
                                if rowMajor
                                    nz = all(dv_next, 2);
                                    R_norm = max(abs((dv_next(nz, :) - dv_loc(nz, :)) ./ dv_loc(nz, :)), [],1)';
                                else
                                    nz = all(dv_next, 1);
                                    R_norm = max(abs((dv_next(:,nz) - dv_loc(:,nz)) ./ dv_loc(:,nz)), [],2);
                                end
                                converged = R_norm < 1e-6;
                                done = converged;
                                if isempty(R_norm)
                                    done = active;
                                end
                                keep = ~done;
                                replace = active;
                                replace(active) = done;
                                v(replace) = v_loc(done);
                                XC(replace) = XC_loc(done);
                                for j = 1 : nmole
                                    XA{j}(replace) = XA_loc{j}(done);
                                end
                                active(active) = ~done;
                                if all(done)
                                    break
                                end
                                if rowMajor
                                    keepRow = keep';
                                    psrk = psrk(:, keepRow);
                                    pa1 = pa1(:, keepRow);
                                    pdh = pdh(:, keepRow);
                                    pb = pb(:, keepRow);
                                    dP = dP(:, keepRow);
                                    s2 = s2(keepRow);
                                    a2 = a2(keepRow);
                                    d40 = d40(keepRow);
                                    b9 = b9(keepRow);
                                    T_loc = T_loc(keepRow);
                                    Tr_loc = Tr_loc(:, keepRow);
                                    dv_loc = dv_next(:, keepRow);
                                    xm_loc = xm_loc(:, keepRow);
                                    a1_loc = a1_loc(keepRow);
                                else
                                    psrk = psrk(keep, :);
                                    pa1 = pa1(keep, :);
                                    pdh = pdh(keep, :);
                                    pb = pb(keep, :);
                                    dP = dP(keep, :);
                                    s2 = s2(keep);
                                    a2 = a2(keep);
                                    d40 = d40(keep);
                                    b9 = b9(keep);
                                    T_loc = T_loc(keep);
                                    Tr_loc = Tr_loc(keep, :);
                                    dv_loc = dv_next(keep, :);
                                    xm_loc = xm_loc(keep, :);
                                    a1_loc = a1_loc(keep);
                                end
                                v_loc = v_loc(keep);
                                B_loc = B_loc(keep);
                                x_loc = cellfun(@(x) x(keep, :), x_loc , 'UniformOutput', false);
                                XA_loc = cellfun(@(x) x(keep, :), XA_loc , 'UniformOutput', false);
                                XC_loc = XC_loc(keep);

                                % time(itea) = toc(t_outer);
                            end
                            if ~all(done)
                                warning('Volume derivative solver did not converge');
                            end
                        else
                            dv = (dP-psrk)./s2;
                            if any(any(dv~=0))
                                if v.jac{i}.isZero
                                    v.jac{i} = v.jac{i}.expandZero();
                                end
                                v.jac{i}.diagonal = dv;
                                if all(v.jac{i}.subset == 0)
                                    v.jac{i}.subset = [];
                                end
                            end
                        end
                    else
                        % Subset of cells
                        map = map{1};
                        vmi = vm(map);
                        bmi = bm(map);
                        Ti = T(map);
                        ami = am(map);
                        SPi = SP(map);
                        IFPi = IFP(map);
                        thetami = thetam(map);
                        gi = g(map);
                        Kappai = Kappa(map);
                        XCmi = XCm(map);
                        gai = ga(map);
                        if rowMajor
                            xmi = xm(:, map);
                            Sigamai = Sigama(:, map);
                            Kaii = Kai(:, map);
                            Ksii = Ksi(:, map);
                            XAmi = XAm(:, map);
                            dei = de(:, map);
                            Tri = Tr(:, map);
                        else
                            xmi = xm(map, :);
                            Sigamai = Sigama(map, :);
                            Kaii = Kai(map, :);
                            Ksii = Ksi(map, :);
                            XAmi = XAm(map, :);
                            dei = de(map, :);
                            Tri = Tr(map, :);
                        end
                        Bi = B(map);
                        Ai = A(map);
                        thetai = theta(map);
                        vi = v(map);
                        XCi = XC(map);
                        XAi = cellfun(@(x) x(map), XA , 'UniformOutput', false);
                        xi = cellfun(@(x) x(map), x , 'UniformOutput', false);
                        Pi = P(map);
                        
                        dA = Ai.jac{i}.getDiagonal;
                        dB = Bi.jac{i}.getDiagonal;
                        
                        [row, column] = size(dA);
                        dtheta = thetai.jac{i}.getDiagonal;
                        for j = 1 : ncomp
                            dx{j} = xi{j}.jac{i}.getDiagonal;
                        end
                        if isa(Pi, 'ADI')
                            dP = Pi.jac{i}.getDiagonal;
                        else
                            dP = zeros(row, column);
                        end
                        
                        [psrk, s2] = model.setSRKDerivatives(Ti, vmi, ami, bmi, dA, dB);
                        if ~isempty(d)
                            [pdh, d40, d23, d26, d27, d3, d8] = model.setDHDerivatives(Ti, xmi, ...
                                vmi, SPi, IFPi, dx, alpha0,thetami,dtheta,...
                                gi,mu0,charge,Sigamai,Kappai,Kaii,Ksii);
                            [pb, b9] = model.setBornDerivatives(xmi, vmi, ...
                                SPi, charge, rBorn, dx, d23, d26, d27, d3, d8);
                        else
                            if rowMajor
                                d40 = zeros(1,column);
                                b9 = zeros(1,column);
                                pdh = zeros(row,column);
                                pb = zeros(row,column);
                            else
                                d40 = zeros(row,1);
                                b9 = zeros(row,1);
                                pdh = zeros(row,column);
                                pb = zeros(row,column);
                            end
                        end
                        
                        if ~isempty(mu0)
                            [pa, pa1, a11, a1, a2] = model.setAssoc1Derivatives(Ti, xmi, ...
                                vmi, XAmi, XCmi, gai, dei, dx, dB);
                            
                            if isempty(dP) && ~isempty(psrk)
                                dP = zeros(row, column);
                            end
                            dv_loc = (dP-(psrk+pa+pdh+pb)) ./ (s2+a11+d40+b9);
                            if rowMajor
                                active = true(column, 1);
                                activeRow = active';
                                T_loc = Ti(activeRow);
                                Tr_loc = Tri(:,activeRow);
                                xm_loc = xmi(:,activeRow);
                                a1_loc = a1(activeRow);
                            else
                                active = true(row, 1);
                                T_loc = Ti(active);
                                Tr_loc = Tri(active, :);
                                xm_loc = xmi(active, :);
                                a1_loc = a1(active);
                            end
                            B_loc = Bi(active);
                            x_loc = cellfun(@(x) x(active), xi, 'UniformOutput', false);
                            v_loc = vi(active);
                            XA_loc = cellfun(@(x) x(active), XAi, 'UniformOutput', false);
                            XC_loc = XCi(active);
                            for itea = 1:10000
                                if v_loc.jac{i}.isZero
                                    v_loc.jac{i} = v_loc.jac{i}.expandZero();
                                end
                                v_loc.jac{i}.diagonal = dv_loc;
                                if all(v_loc.jac{i}.subset == 0)
                                    v_loc.jac{i}.subset = [];
                                end
                                [XA_loc, XC_loc] = setXDerivatives(model, T_loc, x_loc, B_loc, epsilonAB, betaAB, b, v_loc, Tr_loc, XA_loc, XC_loc, i);
                                for j = 1 : nmole
                                    dXA{j} = XA_loc{j}.jac{i}.diagonal;
                                end
                                dXC = XC_loc.jac{i}.diagonal;
                                pa3 = model.setAssoc2Derivatives(T_loc, ...
                                    xm_loc, a1_loc, dXA, dXC);
                                
                                dv_next = (dP-(psrk+pa1+pa3+pdh+pb)) ./ (s2+a2+d40+b9);
                                if rowMajor
                                    nz = all(dv_next, 2);
                                    R_norm = max(abs((dv_next(nz, :) - dv_loc(nz, :)) ./ dv_loc(nz, :)), [],1)';
                                else
                                    nz = all(dv_next, 1);
                                    R_norm = max(abs((dv_next(:,nz) - dv_loc(:,nz)) ./ dv_loc(:,nz)), [],2);
                                end
                                converged = R_norm < 1e-6;
                                done = converged;
                                if isempty(R_norm)
                                    done = active;
                                end
                                keep = ~done;
                                replace = active;
                                replace(active) = done;
                                vi(replace) = v_loc(done);
                                XCi(replace) = XC_loc(done);
                                for j = 1 : nmole
                                    XAi{j}(replace) = XA_loc{j}(done);
                                end
                                active(active) = ~done;
                                if all(done)
                                    v(map) = vi;
                                    XC(map) = XCi;
                                    for j = 1 : nmole
                                        XA{j}(map) = XAi{j};
                                    end
                                    break
                                end
                                if rowMajor
                                    keepRow = keep';
                                    psrk = psrk(:, keepRow);
                                    pa1 = pa1(:, keepRow);
                                    pdh = pdh(:, keepRow);
                                    pb = pb(:, keepRow);
                                    dP = dP(:, keepRow);
                                    s2 = s2(keepRow);
                                    a2 = a2(keepRow);
                                    d40 = d40(keepRow);
                                    b9 = b9(keepRow);
                                    T_loc = T_loc(keepRow);
                                    Tr_loc = Tr_loc(:, keepRow);
                                    dv_loc = dv_next(:, keepRow);
                                    xm_loc = xm_loc(:, keepRow);
                                    a1_loc = a1_loc(keepRow);
                                else
                                    psrk = psrk(keep, :);
                                    pa1 = pa1(keep, :);
                                    pdh = pdh(keep, :);
                                    pb = pb(keep, :);
                                    dP = dP(keep, :);
                                    s2 = s2(keep);
                                    a2 = a2(keep);
                                    d40 = d40(keep);
                                    b9 = b9(keep);
                                    T_loc = T_loc(keep);
                                    Tr_loc = Tr_loc(keep, :);
                                    dv_loc = dv_next(keep, :);
                                    xm_loc = xm_loc(keep, :);
                                    a1_loc = a1_loc(keep);
                                end
                                v_loc = v_loc(keep);
                                B_loc = B_loc(keep);
                                x_loc = cellfun(@(x) x(keep, :), x_loc , 'UniformOutput', false);
                                XA_loc = cellfun(@(x) x(keep, :), XA_loc , 'UniformOutput', false);
                                XC_loc = XC_loc(keep);
                            end
                            if ~all(done)
                                warning('Volume derivative solver did not converge');
                            end
                        else
                            dv = (dP-psrk)./s2;
                            if any(any(dv~=0))
                                if v.jac{i}.isZero
                                    v.jac{i} = v.jac{i}.expandZero();
                                end
                                v.jac{i} = v.jac{i}.setDiagonal(map, dv);
                                if ~isempty(v.jac{i}.subset)
                                    v.jac{i}.subset(map) = (1:numel(map))';
                                end
                            end
                        end                    
                    end
                end
            else
                for i = 1:numel(v.jac)
                    if isempty(cellJacMap{i})
                        [n, m] = size(v.jac{i});
                        if n == 0 || m ~= numel(vm)
                            % There are either no derivatives present or the
                            % derivatives are not of the right dimension
                            continue
                        end
                        dA = model.getJac(A, i);
                        dB = model.getJac(B, i);
                        dP = model.getJac(P, i);
                        dtheta = model.getJac(theta, i);
                        for j = 1 : ncomp
                            dx{j} = model.getJac(x{j}, i);
                        end
                        
                        [psrk, s2] = model.setSRKDerivatives(T, vm, am, bm, dA, dB);
                        if ~isempty(d)
                            [pdh, d40, d23, d26, d27, d3, d8] = model.setDHDerivatives(T, xm, ...
                                vm, SP, IFP, dx, alpha0,thetam,dtheta,...
                                g, mu0, charge, Sigama, Kappa, Kai, Ksi);
                            [pb, b9] = model.setBornDerivatives(xm, vm, ...
                                SP, charge, rBorn, dx, d23, d26, d27, d3, d8);                         
                        else
                            d40 = zeros(m,1); 
                            pdh = zeros(m,1);
                            b9 = zeros(m,1); 
                            pb = zeros(m,1);
                        end
                        
                        if ~isempty(mu0)
                            [pa, pa1, a11, a1, a2] = model.setAssoc1Derivatives(T, xm, ...
                                vm, XAm, XCm, ga, de, dx, dB);
                            
                            dv_loc = (dP-(psrk+pa+pdh+pb)) ./ (s2+a11+d40+b9);
                            active = true(m, 1);
                            B_loc = model.activeADI(B, i, active);
                            T_loc = model.activeADI(T, i, active);
                            v_loc = model.activeADI(v, i, active);
                            XC_loc = model.activeADI(XC, i, active);
                            Tr_loc = model.activeADI(Tr, i, active);
                            xm_loc = model.activeADI(xm, i, active);
                            a1_loc = model.activeADI(a1, i, active);
                            XA_loc =cell(1,nmole);x_loc =cell(1,ncomp);
                            for j = 1 : ncomp
                                if j<=nmole
                                    XA_loc{j} = model.activeADI(XA{j}, i, active);
                                end
                                x_loc{j} = model.activeADI(x{j}, i, active);
                            end
                            for itea = 1:10000
                                n = numel(dv_loc);
                                v_loc.jac = {sparse((1:n)', (1:n)', dv_loc, n, n)};
                                [XA_loc, XC_loc] = setXDerivatives(model, T_loc, x_loc, B_loc, epsilonAB, betaAB, b, v_loc, Tr_loc, XA_loc, XC_loc, 1);
                                for j = 1 : nmole
                                    dXA{j} = model.getJac(XA_loc{j},1);
                                end
                                dXC = model.getJac(XC_loc,1);
                                pa3 = model.setAssoc2Derivatives(T_loc, ...
                                    xm_loc, a1_loc, dXA, dXC);
                                
                                dv_next = (dP-(psrk+pa1+pa3+pdh+pb))./ (s2+a2+d40+b9);
                                R_norm = abs((dv_next - dv_loc) ./ dv_loc);
                                R_zero = dv_next - dv_loc;
                                converged = R_norm < 1e-6;
                                trivial = R_zero == 0;
                                done = trivial | converged;
                                keep = ~done;
                                replace = active;
                                replace(active) = done;
                                v = model.replaceADI(v, v_loc, i, replace, done);
                                XC = model.replaceADI(XC, XC_loc, i, replace, done);
                                for j = 1 : nmole
                                    XA{j} = model.replaceADI(XA{j},XA_loc{j},i, replace, done);
                                end
                                active(active) = ~done;
                                if all(done)
                                    break
                                end
                                psrk = psrk(keep);
                                pa1 = pa1(keep);
                                pdh = pdh(keep);
                                pb = pb(keep);
                                if numel(dP) ~= 1
                                    dP = dP(keep);
                                end
                                s2 = s2(keep);
                                a2 = a2(keep);
                                d40 = d40(keep);
                                b9 = b9(keep);
                                a1_loc = a1_loc(keep);
                                T_loc = T_loc(keep);
                                Tr_loc = Tr_loc(keep, :);
                                dv_loc = dv_next(keep);
                                xm_loc = xm_loc(keep,:);
                                v_loc = model.keepADI(v_loc, keep);
                                B_loc = model.keepADI(B_loc, keep);
                                XC_loc = model.keepADI(XC_loc, keep);
                                for j = 1 : ncomp
                                    if j<=nmole
                                        XA_loc{j} = model.keepADI(XA_loc{j}, keep);
                                    end
                                    x_loc{j} = model.keepADI(x_loc{j}, keep);
                                end
                            end
                            if ~all(done)
                                warning('Volume derivative solver did not converge');
                            end
                        else
                            dv = (dP-psrk) ./ s2;
                            v.jac{i} = sparse((1:n)', (1:m)', dv, n, m);
                        end
                    else
                        map = cellJacMap{i};
                        vmi = vm(map);
                        bmi = bm(map);
                        Ti = T(map);
                        ami = am(map);
                        SPi = SP(map);
                        IFPi = IFP(map);
                        thetami = thetam(map);
                        xmi = xm(map, :);
                        gi = g(map);
                        Sigamai = Sigama(map, :);
                        Kappai = Kappa(map);
                        Kaii = Kai(map, :);
                        Ksii = Ksi(map, :);
                        XAmi = XAm(map, :);
                        XCmi = XCm(map);
                        gai = ga(map);
                        dei = de(map, :);
                        Bi = B(map);
                        Ai = A(map);
                        Pi = P(map);
                        thetai = theta(map);
                        vi = v(map);
                        XCi = XC(map);
                        XAi = cellfun(@(x) x(map), XA , 'UniformOutput', false);
                        xi = cellfun(@(x) x(map), x , 'UniformOutput', false);
                        Tri = Tr(map, :);
                        [n, m] = size(v.jac{i});
                        if n == 0 || m ~= numel(vmi)
                            % There are either no derivatives present or the
                            % derivatives are not of the right dimension
                            continue
                        end
                        dA = model.getJac(Ai, i);
                        dB = model.getJac(Bi, i);
                        dP = model.getJac(Pi, i);
                        dtheta = model.getJac(thetai, i);
                        for j = 1 : ncomp
                            dx{j} = model.getJac(xi{j}, i);
                        end
                        
                        [psrk, s2] = model.setSRKDerivatives(Ti, vmi, ami, bmi, dA, dB);
                        if ~isempty(d)
                            [pdh, d40, d23, d26, d27, d3, d8] = model.setDHDerivatives(Ti, xmi, ...
                                vmi, SPi, IFPi, dx, alpha0,thetami,dtheta,...
                                gi,mu0,charge,Sigamai,Kappai,Kaii,Ksii);
                            [pb, b9] = model.setBornDerivatives(xmi, vmi, ...
                                SPi, charge, rBorn, dx, d23, d26, d27, d3, d8);
                        else
                            d40 = zeros(m,1);
                            pdh = zeros(m,1);
                            b9 = zeros(m,1);
                            pb = zeros(m,1);
                        end
                        
                        if ~isempty(mu0)
                            [pa, pa1, a11, a1, a2] = model.setAssoc1Derivatives(Ti, xmi, ...
                                vmi, XAmi, XCmi, gai, dei, dx, dB);
                                                        
                            dv_loc = (dP-(psrk+pa+pdh+pb)) ./ (s2+a11+d40+b9);
                            active = true(m, 1);
                            B_loc = model.activeADI(Bi, i, active);
                            T_loc = model.activeADI(Ti, i, active);
                            v_loc = model.activeADI(vi, i, active);
                            XC_loc = model.activeADI(XCi, i, active);
                            Tr_loc = model.activeADI(Tri, i, active);
                            xm_loc = model.activeADI(xmi, i, active);
                            a1_loc = model.activeADI(a1, i, active);
                            XA_loc =cell(1,nmole);x_loc =cell(1,ncomp);
                            for j = 1 : ncomp
                                if j<=nmole
                                    XA_loc{j} = model.activeADI(XAi{j}, i, active);
                                end
                                x_loc{j} = model.activeADI(xi{j}, i, active);
                            end
                            for itea = 1:10000
                                n = numel(dv_loc);
                                v_loc.jac = {sparse((1:n)', (1:n)', dv_loc, n, n)};
                                [XA_loc, XC_loc] = setXDerivatives(model, T_loc, x_loc, B_loc, epsilonAB, betaAB, b, v_loc, Tr_loc, XA_loc, XC_loc, 1);
                                for j = 1 : nmole
                                    dXA{j} = model.getJac(XA_loc{j},1);
                                end
                                dXC = model.getJac(XC_loc,1);
                                pa3 = model.setAssoc2Derivatives(T_loc, ...
                                    xm_loc, a1_loc, dXA, dXC);
                                
                                dv_next = (dP-(psrk+pa1+pa3+pdh+pb))./ (s2+a2+d40+b9);
                                R_norm = abs((dv_next - dv_loc) ./ dv_loc);
                                R_zero = dv_next - dv_loc;
                                converged = R_norm < 1e-6;
                                trivial = R_zero == 0;
                                done = trivial | converged;
                                keep = ~done;
                                replace = active;
                                replace(active) = done;
                                v(map) = model.replaceADI(v(map), v_loc, i, replace, done);
                                XC(map) = model.replaceADI(XC(map), XC_loc, i, replace, done);
                                for j = 1 : nmole
                                    XA{j}(map) = model.replaceADI(XA{j}(map),XA_loc{j},i, replace, done);
                                end
                                active(active) = ~done;
                                if all(done)
                                    break
                                end
                                psrk = psrk(keep);
                                pa1 = pa1(keep);
                                pdh = pdh(keep);
                                pb = pb(keep);
                                if numel(dP) ~= 1
                                    dP = dP(keep);
                                end
                                s2 = s2(keep);
                                a2 = a2(keep);
                                d40 = d40(keep);
                                b9 = b9(keep);
                                a1_loc = a1_loc(keep);
                                T_loc = T_loc(keep);
                                Tr_loc = Tr_loc(keep, :);
                                dv_loc = dv_next(keep);
                                xm_loc = xm_loc(keep,:);
                                v_loc = model.keepADI(v_loc, keep);
                                B_loc = model.keepADI(B_loc, keep);
                                XC_loc = model.keepADI(XC_loc, keep);
                                for j = 1 : ncomp
                                    if j<=nmole
                                        XA_loc{j} = model.keepADI(XA_loc{j}, keep);
                                    end
                                    x_loc{j} = model.keepADI(x_loc{j}, keep);
                                end
                            end
                            if ~all(done)
                                warning('Volume derivative solver did not converge');
                            end
                        else
                            dv = (dP-psrk) ./ s2;
                            v.jac{i} = sparse((1:n)', (1:m)', dv, n, m);
                        end
                    end
                end
            end
        end
        
        function [XAxy, XCxy] = setXDerivatives(model, T, xy, B, epsilonAB, betaAB, Bi, vm, Tr, XA, XC, i)
            nmole = model.getNumberOfMolecules();

            isDiagonal = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
            if isDiagonal
                rowMajor = model.AutoDiffBackend.rowMajor;
            else
                rowMajor = false;
            end
            if rowMajor
                T = T';
            end
            
            rouxy = 1./vm;
            etaxy = 0.25 * rouxy .* B;
            gxy = 1 ./(1 - 1.9*etaxy);
            dexy = cell(1, nmole);
            dexy{1} = gxy .* (exp(epsilonAB(1)./T) - 1).* Bi(1) .* betaAB(1);
            Xcxy = (-1+(1+8.*rouxy.*xy{1}.*dexy{1}).^0.5)./(4*rouxy.*xy{1}.*dexy{1});
            ap = model.ECPACompositionalMixture.getAssociationParameter();
            bxy = 1;
            XAxy = cell(1, nmole);
            for j = 2 : nmole
                dexy{j} = dexy{1} .* ap(1,j);
                XAxy{j} = 1./(1+2.*rouxy.*xy{1}.*Xcxy.*dexy{j});
                bxy = bxy+rouxy.*xy{j}.*XAxy{j}.*dexy{j};
            end
            XCxy=(-bxy+(bxy.^2+8*bxy.*rouxy.*xy{1}.*dexy{1}).^0.5)./(4*bxy.*rouxy.*xy{1}.*dexy{1});
            XAxy{1} = 1./(1+2*rouxy.*xy{1}.*XCxy.*dexy{1});
        end
        
        function [psrk, s2] = setSRKDerivatives(model, T, vm, am, bm, dA, dB)
            R = 8.314462618;
            s0 = vm.*(vm+bm);
            s1 = R.*T./(vm-bm).^2;
            s2 = am.*(2*vm+bm)./s0.^2 - s1;
            s3 = am.*vm./s0.^2 + s1;
            s4 = -1./ s0;
            psrk = s3 .* dB + s4 .* dA;
        end
        
        function [pdh, d40, d23, d26, d27, d3, d8] = setDHDerivatives(model, ...
                T, xm, vm, SP, IFP, dx, alpha0, thetam, dtheta, g, mu0, ...
                charge, Sigama, Kappa, Kai, Ksi)
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            kB = 1.380649e-23;
            ncomp = model.getNumberOfComponents();
            nmole = model.getNumberOfMolecules();
            isDiagonal = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
            if isDiagonal
                rowMajor = model.AutoDiffBackend.rowMajor;
            else
                rowMajor = false;
            end
            
            d0 = 2*SP.^2 + IFP.^2;
            d1 = 4*IFP + SP.*(4*SP - IFP + 2);
            d2 = 4*SP.^2.*IFP - SP.*IFP.^2 - 2*SP - 4*IFP;
            d3 = -(d0 + d2) ./ (3*vm.*d0);
            d4 = 0;
            for j = 1 : ncomp
                d4 = d4 + dx{j}.* alpha0(j);
            end
            d12 = 0;d14 = 0;d28 = 0;
            if rowMajor
                d5 = (dx{1}.*thetam + xm(1,:).*dtheta).*g.*mu0.^2;
                d9 = sum(xm.*charge.^2, 1);
                d10 = sum(xm(nmole+1:end,:).*charge(nmole+1:end).^2.*Sigama, 1);
                for j = nmole+1 : ncomp
                    d12 = d12 + dx{j} .* (charge(j).^2.* Kai(j-nmole,:));
                    d28 = d28 + dx{j} .* (charge(j).^2.* Sigama(j-nmole,:));
                    d14 = d14 + dx{j} .* charge(j).^2;
                end
                d17 = sum(xm(nmole+1:end,:).*charge(nmole+1:end).^2.*Kai, 1);
                d30 = sum(xm(nmole+1:end,:).*charge(nmole+1:end).^2.*Ksi, 1);
            else
                d5 = (dx{1}.*thetam + xm(:,1).*dtheta).*g.*mu0.^2;
                d9 = sum(xm.*charge.^2, 2);
                d10 = sum(xm(:,nmole+1:end).*charge(nmole+1:end).^2.*Sigama, 2);
                for j = nmole+1 : ncomp
                    d12 = d12 + dx{j}.* charge(j).^2.* Kai(:,j-nmole);
                    d28 = d28 + dx{j}.* charge(j).^2.* Sigama(:,j-nmole);
                    d14 = d14 + dx{j} .* charge(j).^2;
                end
                d17 = sum(xm(:,nmole+1:end).*charge(nmole+1:end).^2.*Kai, 2);
                d30 = sum(xm(:,nmole+1:end).*charge(nmole+1:end).^2.*Ksi, 2);
            end
            d6 = (NA * (IFP + 2).* d1 ./(9*VP.*vm.*d0)) .* d4;
            d7 = (NA * (IFP + 2).^2 .* SP ./(9*VP.*kB.*T.*vm.*d0)) .* d5;
            d8 = d6 + d7;
            
            d11 = -(kB*T.*Kappa.*d10) ./ (8*pi.*vm.*d9);
            d13 = kB*T.*d12 ./ (4*pi.*d9);
            d15 = kB*T.*Kappa.*d10.*d14 ./ (8*pi.*d9.^2);
            d16 = -kB*T.*Kappa.*d10 ./ (8*pi.*d9);
            d18 = -kB*T.*d17.*d14 ./ (4*pi.*d9.^2);
            d19 = d11 + d16.*d3;
            d20 = d13 + d15 + d16.*d8 + d18;
            d21 = -(IFP+2).*(IFP-1)./(3.*vm);
            d22 = ((IFP+2)./3).^2 .*NA./(VP.*vm) .* d4;
            d23 = 2 - d2./d0;
            d24 = -8*SP.*IFP.^3+ IFP.^4 +2*IFP.^2 -2*SP.^2.*IFP.^2 -4*SP.^2-16*SP.*IFP;
            d25 = -8*SP.^4 +4*SP.^3.*IFP+ 8*SP.^2 + 4*SP.^2.*IFP.^2- 4*SP.*IFP -4*IFP.^2;
            d26 = (SP.*d24.*d3+d25.*d21) ./ d0.^2;
            d27 = (SP.*d24.*d8+d25.*d22) ./ d0.^2;
            d29 = -kB*T.*Kappa.*d23.*d28 ./ (24*pi.*d9);
            d31 = -0.5*Kappa .* (1./vm+d3);
            d32 = 0.5*Kappa .* (d14./d9 - d8);
            d33 = -kB*T.*d30.*d23.*d31 ./ (24*pi.*d9);
            d34 = -kB*T.*d30.*d23.*d32 ./ (24*pi.*d9);
            d35 = -kB*T.*Kappa.*d10.*d26 ./ (24*pi.*d9);
            d36 = -kB*T.*Kappa.*d10.*d27 ./ (24*pi.*d9);
            d37 = kB*T.*Kappa.*d10.*d23.*d14 ./ (24*pi.*d9.^2);
            d38 = d33 + d35;
            d39 = d29 + d34 + d36 + d37;
            d40 = d19 + d38;
            pdh = d20 + d39;
        end
        
        function [pb, b9] = setBornDerivatives(model,xm,vm,SP,charge,rBorn,dx,d23,d26,d27,d3,d8)
            NA = 6.02214076e23;
            VP = 8.8541878128e-12;
            EC = 1.602176634e-19;
            ncomp = model.getNumberOfComponents();
            nmole = model.getNumberOfMolecules();
            isDiagonal = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
            if isDiagonal
                rowMajor = model.AutoDiffBackend.rowMajor;
            else
                rowMajor = false;
            end
            
            b0 = d23 - 3;
            if rowMajor
                b1 = sum(xm(nmole+1:end,:).*charge(nmole+1:end).^2./rBorn, 1);
            else
                b1 = sum(xm(:,nmole+1:end).*charge(nmole+1:end).^2./rBorn, 2);
            end
            b2 = -NA*EC.^2.*b1.*b0./(24*pi.*VP.*SP.*vm.^2);
            b3 = 0;
            for j = nmole+1 : ncomp
                b3 = b3 + dx{j} .* (charge(j).^2./ rBorn(j-nmole));
            end
            b4 = NA*EC.^2.*b0.*b3 ./ (24*pi.*VP.*SP.*vm);
            b5 = d26 - b0.*d3;
            b6 = d27 - b0.*d8;
            b7 = NA*EC.^2.*b1.*b5 ./ (24*pi.*VP.*SP.*vm);
            b8 = NA*EC.^2.*b1.*b6 ./ (24*pi.*VP.*SP.*vm);
            b9 = b2+b7;
            pb = b4+b8;
        end
        
        function [pa, pa1, a11, a1, a2] = setAssoc1Derivatives(model,T,xm,vm,XAm,XCm,ga,de,dx,dB)
            R = 8.314462618;
            nmole = model.getNumberOfMolecules();
            
            isDiagonal = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
            if isDiagonal
                rowMajor = model.AutoDiffBackend.rowMajor;
            else
                rowMajor = false;
            end
            if rowMajor
                a0 = 2 * xm(1,:) .* (2 - XAm(1,:) - XCm(1,:)) ...
                    + sum(xm(2:nmole,:) .* (1 - XAm(2:nmole,:)), 1);
                a4 = 2 * dx{1} .* (2 - XAm(1,:) - XCm(1,:));
                for j = 2 : nmole
                    a4 = a4 + dx{j} .* (1-XAm(j,:));
                end
                a6 = 2*(xm(1,:).*XAm(1,:).^2.*(1./XAm(1,:)-1) ...
                    +xm(1,:).*XCm(1,:).^2.*(1./XCm(1,:)-1) );
                for j = 2 : nmole
                    a6 = a6 + xm(j,:).*XAm(j,:).^2.*(1./XAm(j,:)-1);
                end
                dA0 = dx{1} .* (2 * XCm(1,:).*de(1,:));
                dC0 = dx{1} .* (2 * XAm(1,:).*de(1,:));
                for j = 2 : nmole
                    dC0 = dC0 + dx{j} .* XAm(j,:).*de(j,:);
                end
                a9 = 2*(xm(1,:).*XAm(1,:).^2.*dA0 ...
                    +xm(1,:).*XCm(1,:).^2.* dC0 );
                for j = 2 : nmole
                    a9 = a9 + xm(j,:).*XAm(j,:).^2 .* (2*dx{1}.*XCm(1,:).*de(j,:));
                end
            else
                a0 = 2 * xm(:,1) .* (2 - XAm(:,1) - XCm(:,1)) ...
                    + sum(xm(:,2:nmole) .* (1 - XAm(:,2:nmole)), 2);
                a4 =  2 * dx{1} .*  (2 - XAm(:,1) - XCm(:,1));
                for j = 2 : nmole
                    a4 = a4 + dx{j} .*(1 - XAm(:,j));
                end
                a6 = 2*(xm(:,1).*XAm(:,1).^2.*(1./XAm(:,1)-1) ...
                    +xm(:,1).*XCm(:,1).^2.*(1./XCm(:,1)-1) );
                for j = 2 : nmole
                    a6 = a6 + xm(:,j).*XAm(:,j).^2.*(1./XAm(:,j)-1);
                end
                dA0 = dx{1}.* (2 * XCm(:,1).*de(:,1));
                dC0 = dx{1}.* (2 * XAm(:,1).*de(:,1));
                for j = 2 : nmole
                    dC0 = dC0 + dx{j}.* (XAm(:,j).*de(:,j));
                end
                a9 = 2*(xm(:,1).*XAm(:,1).^2.*dA0 ...
                    +xm(:,1).*XCm(:,1).^2.* dC0 );
                for j = 2 : nmole
                    a9 = a9 + xm(:,j).*XAm(:,j).^2 .* dx{1}.* (2 *XCm(:,1).*de(:,j));
                end
            end
            a1 = ga ./ vm;
            a2 = 0.5 * R .* T .* a1.^2 .* a0;
            a3 = - 1.9 / 8 * R .* T .* a1.^2 .* a0;
            a5 = -0.5 * R .* T .* a1 .* a4;
            pa1 = a3.*dB + a5;
            a7 = 0.5 * R .* T .* a1.^2 .* a6;
            a8 = -1.9/8 * R .* T .* a1.^2 .* a6;
            
            a10 = -0.5 * R .* T .* a1./vm .* a9;
            pa2 = a8.*dB + a10;
            a11 = a2 + a7;
            pa = pa1 + pa2;
        end
        
        function pa3 = setAssoc2Derivatives(model,T_loc,xm_loc,a1_loc,dXA,dXC)
            R = 8.314462618;
            nmole = model.getNumberOfMolecules();
            isDiagonal = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
            if isDiagonal
                rowMajor = model.AutoDiffBackend.rowMajor;
            else
                rowMajor = false;
            end
            
            if rowMajor
                pa3 = 2*xm_loc(1,:).*(dXA{1}+dXC);
                for j = 2 : nmole
                    if ~isempty(dXA{j})
                        pa3 = pa3 + xm_loc(j,:).*dXA{j};
                    end
                end
            else
                pa3 = 2*xm_loc(:,1).*(dXA{1}+dXC);
                for j = 2 : nmole
                    if ~isempty(dXA{j})
                        pa3 = pa3 + xm_loc(:,j).*dXA{j};
                    end
                end
            end
            pa3 = 0.5 * R .* T_loc .* a1_loc .* pa3;
        end
                
        function kij = getMoleculeBIP(model, Tr, names)
            ncomp = numel(names);
            kij = cell(ncomp,ncomp);
            for i = 1 : ncomp-1
                switch(lower(names{i}))
                    case {'water', 'h2o'}
                        for j = i + 1 : ncomp
                            switch(lower(names{j}))
                                case {'carbondioxide', 'co2'}
                                    if iscell(Tr)
                                        kij{i,j} = - 0.404265 .* Tr{j} .^ 2 ...
                                            + 1.263 .* Tr{j} - 0.76809;
                                    else
                                        kij{i,j} = - 0.404265 .* Tr(:,j) .^ 2 ...
                                            + 1.263 .* Tr(:,j) - 0.76809;
                                    end
                                    kij{j,i} = kij{i,j};
                                case {'hydrogensulfide', 'h2s'}
                                    if iscell(Tr)
                                        kij{i,j} = - 0.40299 .* Tr{j} .^ 2 ...
                                            + 0.825348 .* Tr{j} - 0.25978;
                                    else
                                        kij{i,j} = - 0.40299 .* Tr(:,j) .^ 2 ...
                                            + 0.825348 .* Tr(:,j) - 0.25978;
                                    end
                                    kij{j,i} = kij{i,j};
                                case {'sulfurdioxide', 'so2'}
                                    if iscell(Tr)
                                        kij{i,j} = -0.1394 .* Tr{j} + 0.0468;
                                    else
                                        kij{i,j} = -0.1394 .* Tr(:,j) + 0.0468;
                                    end
                                    kij{j,i} = kij{i,j};
                                case {'methane', 'ch4'}
                                    if iscell(Tr)
                                        kij{i,j} = - 0.364714 .* Tr{j} .^ 2 ...
                                            + 1.73109 .* Tr{j} - 1.817027;
                                    else
                                        kij{i,j} = - 0.364714 .* Tr(:,j) .^ 2 ...
                                            + 1.73109 .* Tr(:,j) - 1.817027;
                                    end
                                    kij{j,i} = kij{i,j};
                                case {'nitrogen', 'n2'}
                                    if iscell(Tr)
                                        kij{i,j} = 0.0578675 .* Tr{j} .^ 4 ...
                                            - 0.844736 .* Tr{j} .^ 3 ...
                                            + 4.2075455 .* Tr{j} .^2 ...
                                            - 8.30813 .* Tr{j} + 5.275865;
                                    else
                                        kij{i,j} = 0.0578675 .* Tr(:,j) .^ 4 ...
                                            - 0.844736 .* Tr(:,j) .^ 3 ...
                                            + 4.2075455 .* Tr(:,j) .^2 ...
                                            - 8.30813 .* Tr(:,j) + 5.275865;
                                    end
                                    kij{j,i} = kij{i,j};
                                case {'oxygen', 'o2'}
                                    if iscell(Tr)
                                        kij{i,j} = 0.0896695 .* Tr{j} .^ 3 ...
                                            - 0.928507 .* Tr{j} .^ 2 ...
                                            + 3.2950445 .* Tr{j} - 3.471575;
                                    else
                                        kij{i,j} = 0.0896695 .* Tr(:,j) .^ 3 ...
                                            - 0.928507 .* Tr(:,j) .^ 2 ...
                                            + 3.2950445 .* Tr(:,j) - 3.471575;
                                    end
                                    kij{j,i} = kij{i,j};
                                case {'argon', 'ar'}
                                    if iscell(Tr)
                                        kij{i,j} = 0.0977033 .* Tr{j} .^ 3 ...
                                            - 0.96074999 .* Tr{j} .^ 2 ...
                                            + 3.3311513 .* Tr{j} - 3.550028;
                                    else
                                        kij{i,j} = 0.0977033 .* Tr(:,j) .^ 3 ...
                                            - 0.96074999 .* Tr(:,j) .^ 2 ...
                                            + 3.3311513 .* Tr(:,j) - 3.550028;
                                    end
                                    kij{j,i} = kij{i,j};
                            end
                        end
                    case {'carbondioxide', 'co2'}
                        for j = i + 1 : ncomp
                            switch(lower(names{j}))
                                case {'methane', 'ch4'}
                                    kij{i,j} = 0.1;
                                    kij{j,i} = kij{i,j};
                                case {'hydrogensulfide', 'h2s'}
                                    kij{i,j} = 0.1;
                                    kij{j,i} = kij{i,j};
                                case {'nitrogen', 'n2'}
                                    kij{i,j} = -0.0712;
                                    kij{j,i} = kij{i,j};
                                case {'oxygen', 'o2'}
                                    kij{i,j} = 0.077;
                                    kij{j,i} = kij{i,j};
                                case {'argon', 'ar'}
                                    kij{i,j} = 0.089;
                                    kij{j,i} = kij{i,j};
                                case {'sulfurdioxide', 'so2'}
                                    kij{i,j} = 0.029;
                                    kij{j,i} = kij{i,j};
                            end
                        end
                    case {'hydrogensulfide', 'h2s'}
                        for j = i + 1 : ncomp
                            switch(lower(names{j}))
                                case {'sulfurdioxide', 'oxygen', 'argon', 'so2', 'o2', 'ar'}
                                    kij{i,j} = 0;
                                    kij{j,i} = kij{i,j};
                                case {'methane', 'ch4'}
                                    kij{i,j} = 0.0852;
                                    kij{j,i} = kij{i,j};
                                case {'nitrogen', 'n2'}
                                    kij{i,j} = 0.1499;
                                    kij{j,i} = kij{i,j};
                            end
                        end
                    case {'sulfurdioxide', 'so2'}
                        for j = i + 1 : ncomp
                            switch(lower(names{j}))
                                case {'oxygen', 'o2'}
                                    kij{i,j} = 0.099;
                                    kij{j,i} = kij{i,j};
                                case {'methane', 'ch4'}
                                    kij{i,j} = 0.047;
                                    kij{j,i} = kij{i,j};
                                case {'nitrogen', 'n2'}
                                    kij{i,j} = -0.056;
                                    kij{j,i} = kij{i,j};
                                case {'argon', 'ar'}
                                    kij{i,j} = 0.086;
                                    kij{j,i} = kij{i,j};
                            end
                        end
                    case {'methane', 'ch4'}
                        for j = i + 1 : ncomp
                            switch(lower(names{j}))
                                case {'oxygen', 'o2'}
                                    kij{i,j} = 0;
                                    kij{j,i} = kij{i,j};
                                case {'nitrogen', 'n2'}
                                    kij{i,j} = 0.0393;
                                    kij{j,i} = kij{i,j};
                                case {'argon', 'ar'}
                                    kij{i,j} = 0.011;
                                    kij{j,i} = kij{i,j};
                            end
                        end
                    case {'nitrogen', 'n2'}
                        for j = i + 1 : ncomp
                            switch(lower(names{j}))
                                case {'oxygen', 'o2'}
                                    kij{i,j} = -0.028;
                                    kij{j,i} = kij{i,j};
                                case {'argon', 'ar'}
                                    kij{i,j} = -0.005;
                                    kij{j,i} = kij{i,j};
                            end
                        end
                    case {'oxygen', 'o2'}
                        for j = i + 1 : ncomp
                            switch(lower(names{j}))
                                case {'argon', 'ar'}
                                    kij{i,j} = 0.008;
                                    kij{j,i} = kij{i,j};
                            end
                        end
                end
            end
        end
        
        function [Uref, Alp, Ta] = getMoleIonBIP(model, names, nmole)
            ncomp = numel(names);
            ion = {names{nmole + 1 : ncomp}};
            Uref = zeros(1,nmole); Alp = zeros(1,nmole); Ta = zeros(1,nmole);
            for i = 1 : nmole
                switch(lower(names{i}))
                    case {'water', 'h2o'}
                        if strcmp(ion,{'Na+','Cl-'})
                            Uref(i) = -365; Alp(i) = 2125; Ta(i) = 356;
                        elseif strcmp(ion,{'K+','Cl-'})
                            Uref(i) = -239; Alp(i) = 1757; Ta(i) = 350;
                        elseif strcmp(ion,{'Ca2+','Cl-'})
                            Uref(i) = -291; Alp(i) = 872; Ta(i) = 285;
                        elseif strcmp(ion,{'Mg2+','Cl-'})
                            Uref(i) = -323; Alp(i) = 808; Ta(i) = 266;
                        elseif strcmp(ion,{'Na+','SO42-'})
                            Uref(i) = 129; Alp(i) = 1225; Ta(i) = 300;
                        end
                    case {'carbondioxide', 'co2'}
                        if strcmp(ion,{'Na+','Cl-'})
                            Uref(i) = 614; Alp(i) = 2701; Ta(i) = 320;
                        elseif strcmp(ion,{'K+','Cl-'})
                            Uref(i) = 371; Alp(i) = 282; Ta(i) = 181;
                        elseif strcmp(ion,{'Ca2+','Cl-'})
                            Uref(i) = 852; Alp(i) = 630; Ta(i) = 200;
                        elseif strcmp(ion,{'Mg2+','Cl-'})
                            Uref(i) = 680; Alp(i) = 773; Ta(i) = 200;
                        elseif strcmp(ion,{'Na+','SO42-'})
                            Uref(i) = 1500; Alp(i) = 191; Ta(i) = 200;
                        end
                    case {'methane', 'ch4'}
                        if strcmp(ion,{'Na+','Cl-'})
                            Uref(i) = 1140; Alp(i) = 1725; Ta(i) = 315;
                        elseif strcmp(ion,{'K+','Cl-'})
                            Uref(i) = 1033; Alp(i) = 19612; Ta(i) = 359;
                        elseif strcmp(ion,{'Ca2+','Cl-'})
                            Uref(i) = 1103; Alp(i) = 7617; Ta(i) = 300;
                        end
                    case {'hydrogensulfide', 'h2s'}
                        if strcmp(ion,{'Na+','Cl-'})
                            Uref(i) = 314; Alp(i) = 1757; Ta(i) = 322;
                        end
                    case {'nitrogen', 'n2'}
                        if strcmp(ion,{'Na+','Cl-'})
                            Uref(i) = 1086; Alp(i) = 11869; Ta(i) = 333;
                        elseif strcmp(ion,{'K+','Cl-'})
                            Uref(i) = 1100; Alp(i) = 1851; Ta(i) = 317;
                        elseif strcmp(ion,{'Na+','SO42-'})
                            Uref(i) = 2189; Alp(i) = 31814; Ta(i) = 680;
                        end
                    case {'oxygen', 'o2'}
                        if strcmp(ion,{'Na+','Cl-'})
                            Uref(i) = 1331; Alp(i) = 29477; Ta(i) = 300;
                        elseif strcmp(ion,{'K+','Cl-'})
                            Uref(i) = 1039; Alp(i) = 31401; Ta(i) = 305;
                        elseif strcmp(ion,{'Ca2+','Cl-'})
                            Uref(i) = 1323; Alp(i) = 34987; Ta(i) = 280;
                        elseif strcmp(ion,{'Mg2+','Cl-'})
                            Uref(i) = 1261; Alp(i) = 4104; Ta(i) = 314;
                        elseif strcmp(ion,{'Na+','SO42-'})
                            Uref(i) = 2249; Alp(i) = 33668; Ta(i) = 297;
                        end
                    case {'argon', 'ar'}
                        if strcmp(ion,{'Na+','Cl-'})
                            Uref(i) = 1281; Alp(i) = 8302; Ta(i) = 325;
                        elseif strcmp(ion,{'K+','Cl-'})
                            Uref(i) = 892; Alp(i) = 33351; Ta(i) = 242;
                        end
                end
            end
        end
        
        function [x, y, L, derivatives] = getPhaseFractionAsADI(model, state, p, T, z)
            % Compute derivatives for values obtained by solving the
            % equilibrium equations (molar fractions in each phase, liquid
            % mole fraction).
            singlePhase = model.getSinglePhase(state);
            twoPhase = ~singlePhase;
            
            x = z;
            y = z;
            L = state.L;
            derivatives = struct('dL', [], 'dx', [], 'dy', []);
            if ~any(twoPhase)
                return
            end
            sample = getSampleAD(p, T, z{:});
            if isfield(state, 'FractionalDerivatives')
                fd = state.FractionalDerivatives;
                [dLdpTz, dxdpTz, dydpTz] = deal(fd.dL, fd.dx, fd.dy);
            else
                [dLdpTz, dxdpTz, dydpTz] = model.getPhaseFractionDerivativesPTZ(state, twoPhase);
            end
            ncomp = model.getNumberOfComponents();
            
            p2ph = p(twoPhase);
            T2ph = T(twoPhase);
            z2ph = cellfun(@(x) x(twoPhase), z, 'UniformOutput', false);
            
            L = model.AutoDiffBackend.convertToAD(L, sample);
            L(twoPhase) = model.addDerivativePTZ(L(twoPhase), dLdpTz, p2ph, T2ph, z2ph);
            for i = 1:ncomp
                x{i}(twoPhase) = state.x(twoPhase, i);
                y{i}(twoPhase) = state.y(twoPhase, i);
                if ~isa(z{i}, 'ADI')
                    x{i} = model.AutoDiffBackend.convertToAD(x{i}, sample);
                    y{i} = model.AutoDiffBackend.convertToAD(y{i}, sample);
                end
                x{i}(twoPhase) = model.addDerivativePTZ(x{i}(twoPhase), dxdpTz{i}, p2ph, T2ph, z2ph);
                y{i}(twoPhase) = model.addDerivativePTZ(y{i}(twoPhase), dydpTz{i}, p2ph, T2ph, z2ph);
            end
            if nargout > 3
                derivatives.dL = dLdpTz;
                derivatives.dx = dxdpTz;
                derivatives.dy = dydpTz;
            end
        end
        
        function [dLdpTz, dxdpTz, dydpTz] = getPhaseFractionDerivativesPTZ(model, state, cells)
            if nargin < 3
                cells = ':';
            end
            % Locally remove AD-context
            L = value(state.L);
            x = value(state.x);
            y = value(state.y);
            p = value(state.pressure);
            T = value(state.T);
            z = value(state.components);
            % Will get derivatives of these
            L = L(cells);
            x = expandMatrixToCell(x(cells, :));
            y = expandMatrixToCell(y(cells, :));
            % With respect to these
            p = p(cells);
            T = T(cells);
            z = expandMatrixToCell(z(cells, :));
            % And we need these
            Z_L = state.Z_L(cells);
            Z_V = state.Z_V(cells);
            % A few constants
            ncell = numel(p);
            ncomp = model.getNumberOfComponents();
            [xAD, yAD, zAD] = deal(cell(1, ncomp));
            
            fullZ = true;
            if fullZ
                [pAD, TAD, zAD{:}] = model.AutoDiffBackend.initVariablesAD(p, T, z{:});
                nprimary = 2 + ncomp;
            else
                [pAD, TAD, zAD{1:end-1}] = model.AutoDiffBackend.initVariablesAD(p, T, z{1:end-1});
                zAD{end} = 1;
                for i = 1:ncomp-1
                    zAD{end} = zAD{end} - zAD{i};
                end
                nprimary = 1 + ncomp;
                nprimary = nprimary + 1;
            end
            [LAD, xAD{:}, yAD{:}] = model.AutoDiffBackend.initVariablesAD(L, x{:}, y{:});
            
            eqsPrim  = model.equationsEquilibrium(pAD, TAD, x, y, zAD, L, Z_L, Z_V);
            eqsSec = model.equationsEquilibrium(p, T, xAD, yAD, z, LAD, Z_L, Z_V);
            ep = combineEquations(eqsPrim{:});
            es = combineEquations(eqsSec{:});
            
            % Let F ble the equilibrium equations (fugacity balance,
            % vapor/liquid balance etc). We can then write by the chain
            % rule (and using partial derivatives instead of total
            % derivatives),
            % DF / Dp = dF / dp +  dF / ds * ds/dp.
            % where p is the primary variables and s the secondary
            % variables. We then obtain
            % ds / dp = -inv(dF / ds)*(DF / Dp)
            dFdp = ep.jac{1};
            dFds = es.jac{1};
            % We really want: dsdp = dFds\dFdp, but right hand side is
            % large so we just use LU factorization since it is much faster
            [dFds_L, dFds_U] = lu(dFds);
            dsdp = -(dFds_U\(dFds_L\dFdp));
            [I, J, V] = find(dsdp);
            % P, T and each component
            dLdpTz = zeros(ncell, nprimary);
            dxdpTz = cell(ncomp, 1);
            [dxdpTz{:}] = deal(zeros(ncell, nprimary));
            dydpTz = dxdpTz;
            
            I = I-1;
            J = J-1;
            cellNo = mod(I, ncell)+1;
            secondaryNo = floor(I/ncell)+1;
            primaryNo = floor(J/ncell)+1;
            
            isdL = secondaryNo == 1;
            makeSub = @(act) sub2ind([ncell, nprimary], cellNo(act), primaryNo(act));
            dLdpTz(makeSub(isdL)) = V(isdL);
            if ~fullZ
                dLdpTz(:, end) = -sum(dLdpTz(:, 3:end), 2);
            end
            
            % Equivialent version - if there are no zeros!
            % dLdpTz = reshape(V(isdL), ncell, nprimary);
            for i = 1:ncomp
                isdx = secondaryNo == 1 + i;
                dxdpTz{i}(makeSub(isdx)) = V(isdx);
                if ~fullZ
                    dxdpTz{i}(:, end) = -sum(dxdpTz{i}(:, 3:end), 2);
                end
                
                isdy = secondaryNo == ncomp + 1 + i;
                dydpTz{i}(makeSub(isdy)) = V(isdy);
                if ~fullZ
                    dydpTz{i}(:, end) = -sum(dydpTz{i}(:, 3:end), 2);
                end
            end
        end
        
        function L = solveRachfordRice(model, L, K, z)
            L = solveRachfordRiceVLE(L, K, z, 'min_z', model.minimumComposition);
        end
        
        function y = computeVapor(model, L, K, z)
            if isstruct(L)
                y = model.computeVapor(L.L, L.K, L.components);
                return
            end
            if iscell(z)
                y = z;
                sv = 0;
                for i = 1:numel(z)
                    bad = ~isfinite(K{i});
                    y{i} = K{i}.*z{i}./(L + (1-L).*K{i});
                    y{i}(bad) = 0;
                    sv = sv + value(y{i});
                end
                y = cellfun(@(k, zi) k.*zi./(L + (1-L).*k), K, z, 'UniformOutput', false);
                y = cellfun(@(x) x./sv, y, 'UniformOutput', false);
                assert(all(cellfun(@(x) all(isfinite(value(x))), y)));
            else
                y = K.*z./bsxfun(@plus, L, bsxfun(@times, 1-L, K));
                y(~isfinite(K)) = z(~isfinite(K));
                y = bsxfun(@rdivide, y, sum(y, 2));
                assert(all(isfinite(y(:))));
            end
        end
        
        function [x, sv] = computeLiquid(model, L, K, z)
            if isstruct(L)
                x = model.computeLiquid(L.L, L.K, L.components);
                return
            end
            if iscell(z)
                x = z;
                sv = 0;
                for i = 1:numel(z)
                    bad = ~isfinite(K{i});
                    x{i} = z{i}./(L + (1-L).*K{i});
                    x{i}(bad) = 0;
                    sv = sv + value(x{i});
                end
                x = cellfun(@(x) x./sv, x, 'UniformOutput', false);
                assert(all(cellfun(@(x) all(isfinite(value(x))), x)));
            else
                tmp = bsxfun(@times, 1-L, K);
                x = z./bsxfun(@plus, L, tmp);
                x(~isfinite(K)) = 0;
                sv = sum(x, 2);
                x = bsxfun(@rdivide, x, sv);
                assert(all(isfinite(x(:))));
            end
        end
        
        function [sL, sV] = computeSaturations(model, rhoL, rhoV, x, y, L, Z_L, Z_V)
            sL = L.*Z_L./(L.*Z_L + (1-L).*Z_V);
            sV = (1-L).*Z_V./(L.*Z_L + (1-L).*Z_V);
        end
        
        function frac = getMoleFraction(model, massfraction)
            % Convert mass fraction to molar fraction
            mm = model.CompositionalMixture.molarMass;
            if iscell(massfraction)
                ncomp = numel(massfraction);
                moles = cell(1, ncomp);
                totMole = 0;
                for i = 1:numel(massfraction)
                    moles{i} = massfraction{i}./mm(i);
                    totMole = totMole + moles{i};
                end
                frac = cellfun(@(x) x./totMole, moles, 'UniformOutput', false);
            else
                moles = bsxfun(@times, massfraction, 1./mm);
                frac = bsxfun(@rdivide, moles, sum(moles, 2));
            end
        end
        
        function frac = getMassFraction(model, molfraction)
            % Convert molar fraction to mass fraction
            if iscell(molfraction)
                ncomp = numel(molfraction);
                mass = cell(1, ncomp);
                totMass = 0;
                for i = 1:ncomp
                    mi = molfraction{i};
                    if ~isempty(mi)
                        mass{i} = model.ECPACompositionalMixture.molarMass(i).*molfraction{i};
                        totMass = totMass + mass{i};
                    end
                end
                frac = cell(size(mass));
                for i = 1:ncomp
                    if ~isempty(mass{i})
                        frac{i} = mass{i}./totMass;
                    end
                end
            else
                mass = bsxfun(@times, molfraction, model.ECPACompositionalMixture.molarMass);
                frac = bsxfun(@rdivide, mass, sum(mass, 2));
            end
        end
        
        function L = estimateSinglePhaseState(model, p, T, z, L, stable)
            z = z(stable, :);
            p = p(stable);
            T = T(stable);
            
            K = eCPAestimateEquilibriumWilson(model, p, T);
            L0 = repmat(0.5, nnz(stable), 1);
            L(stable) = model.solveRachfordRice(L0, K, z);
            
            L(stable & L > 0.5) = 1;
            L(stable & L <= 0.5) = 0;
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
            changed = false;
            state = model.setFlag(state);
            if ~isfield(state, 'x')
                state.x = model.computeLiquid(state);
                changed = true;
            end
            if ~isfield(state, 'y')
                state.y = model.computeVapor(state);
                changed = true;
            end
            if changed
                singlePhase = state.L == 1 | state.L == 0;
                z_1ph = state.components(singlePhase, :);
                state.x(singlePhase, :) = z_1ph;
                state.y(singlePhase, :) = z_1ph;
            end
        end
        
        function state = setFlag(model, state, pureLiquid, pureVapor)
            if nargin < 4
                L = value(state.L);
                pureVapor = L == 0;
                if nargin < 3
                    pureLiquid = L == 1;
                end
            end
            
            pureVapor(pureLiquid) = false;
            state.flag = 1*pureLiquid + 2*pureVapor;
        end
        
        function isSinglePhase = getSinglePhase(model, state)
            if ~isfield(state, 'flag')
                state = model.setFlag(state);
            end
            isSinglePhase = state.flag ~= 0;
        end
        
        function [isLiquid, isVapor, is2ph] = getFlag(model, state)
            if ~isfield(state, 'flag')
                state = model.setFlag(state);
            end
            flag = state.flag;
            isLiquid = flag == 1;
            isVapor = flag == 2;
            is2ph = ~(isLiquid | isVapor);
        end
        
        function is2ph = getTwoPhaseFlag(model, state)
            if ~isfield(state, 'flag')
                state = model.setFlag(state);
            end
            is2ph = state.flag == 0;
        end
    end
    
    methods (Static, Access=protected)
        function x = activeADI(x, ix, subs)
            if islogical(subs)
                subs = find(subs);
            end
            if isa(x, 'ADI')
                x.val = x.val(subs);
                x.jac = {x.jac{ix}};
                x.jac{1} = x.jac{1}(subs,:);
                x.jac{1} = x.jac{1}(:,subs);
            else
                x = x(subs, :);
            end
        end
        
        function x = keepADI(x, subs)
            if islogical(subs)
                subs = find(subs);
            end
            if isa(x, 'ADI')
                x.val = x.val(subs);
                x.jac{1} = x.jac{1}(subs,:);
                x.jac{1} = x.jac{1}(:,subs);
            else
                x = x(subs);
            end
        end
        
        function x = replaceADI(x, y, ix, replace, subs)
            if islogical(subs)
                subs = find(subs);
            end
            if islogical(replace)
                replace = find(replace);
            end
            if ~isempty(subs)
                x.val(replace) = y.val(subs);
                x.jac{ix}(replace,replace) = y.jac{1}(subs,subs);
            end
        end
        
        function dx = getDiagonalJac(x, ix)
            if isa(x, 'GenericAD')
                dx = x.jac{ix}.diagonal;
            else
                dx = 0;
            end
        end
        
        function dx = getJac(x, ix)
            if isa(x, 'ADI')
                dx = diag(x.jac{ix});
            else
                dx = 0;
            end
        end
        
        
        function dx = getJacSub(x, ix, map)
            if isa(x, 'ADI')
                [ii, jj, vv] = find(x.jac{ix});
                if isempty(vv)
                    dx = 0;
                else
                    keep = ii == map(jj);
                    dx = zeros(size(map));
                    dx(jj(keep)) = vv(keep);
                end
            else
                dx = 0;
            end
        end
        
        function D = makeDiagonal(x, vx, ncell)
            D = sparse(vx, vx, x, ncell, ncell);
        end
        
        function checkv(v)
            % Throw error if v takes on unphysical values: Negative, or
            % non-finite values.
            if any(v < 0)
                error('%d negative Z-factors detected...', sum(v < 0));
            end
            if any(~isfinite(v))
                error('%d non-finite Z-factors detected...', sum(~isfinite(v)));
            end
        end
        
        function result = addDerivativePTZ(result, dvdptz, p, T, z)
            if ~isa(result, 'ADI')
                return
            end
            nJac = numel(result.jac);
            for jacNo = 1:nJac
                result.jac{jacNo} = 0*result.jac{jacNo};
            end
            variables = [{p}, {T}, z];
            for varNo = 1:numel(variables)
                v = variables{varNo};
                if isa(v, 'ADI')
                    v_der = dvdptz(:, varNo);
                    D = [];
                    for jacNo = 1:nJac
                        [inc, D] = diagMult(v_der, v.jac{jacNo}, D);
                        result.jac{jacNo} = result.jac{jacNo} + inc;
                    end
                end
            end
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

