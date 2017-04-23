classdef EquationOfStateModel < PhysicalModel
    % Equation of state model. Implements generalized two-parameter cubic
    % equation of state with Newton and successive substitution solvers, as
    % well as standard functions for computing density and viscosity.
    properties
        fluid
        omegaA
        omegaB
        % Use Newton based solver for flash. If set to false, successive
        % substitution if used instead.
        useNewton
        fastDerivatives
        performFlash
    end
    
    properties (Access = private)
        eosA
        eosB
        eosType
    end
    
    methods
        
        function model = EquationOfStateModel(G, fluid, eosname)
            if nargin < 3
                eosname = 'PR';
            end
            model = model@PhysicalModel(G);
            assert(isa(fluid, 'CompositionalFluid'));
            model.fluid = fluid;
            model.nonlinearTolerance = 1e-4;
            model.useNewton = false;
            model.fastDerivatives = true;
            model.performFlash = true;
            
            model = model.setType(eosname);
        end
        
        function Z = solveCubicEOS(model, A, B)
            % Peng Robinson equation of state in form used by Coats
            [E2, E1, E0] = model.getCubicCoefficients(A, B);
            if 1
                % Use vectorized cubic solver
                Z = cubicPositive(E2, E1, E0);
%                 Z = mrstCubic(1, E2, E1, E0);
            else
                % Fallback to Matlab's root solver
                n = numel(A);
                Z = zeros(n, 3);
                for i = 1:n
                    Z(i, :) = roots([1, E2(i), E1(i), E0(i)]);
                end
            end
            % Avoid picking any complex roots. At least one root will be
            % real. We replace these invalid values with NaN so that when
            % min or max is used, they will ignore any previously complex
            % values.
            mtol = 0;
            Z(abs(imag(Z)) > mtol) = nan;
            Z(Z <= 0) = nan;
            Z = real(Z);
        end
        
        function model = setType(model, arg, c1, c2)
            if ischar(arg)
                switch(lower(arg))
                    case {'pr', 'peng-robinson'}
                        % Peng-Robinson
                        model.eosA = 1 + sqrt(2);
                        model.eosB = 1 - sqrt(2);
                        model.omegaA = 0.4572355;
                        model.omegaB = 0.0779691;
                        model.eosType = 1;
                    case {'srk', 'soave-redlich-kwong'}
                        % Soave-Redlich-Kwong
                        model.eosA = 0;
                        model.eosB = 1;
                        model.omegaA = 0.427480;
                        model.omegaB = 0.086640;
                        model.eosType = 2;
                    case {'zj', 'zudkevitch-joffe'}
                        % Zudkevitch-Joffe
                        model.eosA = 0;
                        model.eosB = 1;
                        model.omegaA = 0.427480;
                        model.omegaB = 0.086640;
                        model.eosType = 3;
                        error('Not implemented yet.')
                    case {'rk', 'redlich-kwong'}
                        % Redlich-Kwong
                        model.eosA = 0;
                        model.eosB = 1;
                        model.omegaA = 0.427480;
                        model.omegaB = 0.086640;
                        model.eosType = 4;
                    otherwise
                    error('Invalid string ''%s''.\n Valid choices are:\n PR: Peng-Robinson\n', arg);
                end
            end
            
            if nargin > 2
                model.eosA = c1;
                if nargin > 3
                    model.eosB = c2;
                end
            end
        end
                
        function [m1, m2] = getEOSCoefficients(model)
            m1 = model.eosA;
            m2 = model.eosB;
        end
        
        function Z = computeLiquidZ(model, A, B)
            % Pick smallest Z factors for liquid phase (least energy)
            Z = model.solveCubicEOS(A, B);
            bad = bsxfun(@lt, Z, B);
            Z(bad) = nan;
            Z = min(Z, [], 2);
            checkZ(Z);
        end
        function Z = computeVaporZ(model, A, B)
            % Pick largest Z factors for vapor phase (most energy)
            Z = model.solveCubicEOS(A, B);
            bad = bsxfun(@lt, Z, B);
            Z(bad) = nan;
            Z = max(Z, [], 2);
            checkZ(Z);
        end
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolve, nonlinsolve, iteration, varargin)%#ok
            % Compute a single step of the solution process for a given
            % equation of state (thermodynamic flash calculation);
            T = state.T;
            P = state.pressure;
            assert(all(P >= 0))
            
            z = state.components;
            ncomp = model.fluid.getNumberOfComponents();
            
            % Basic assertions
            assert(all(sum([z{:}], 2) > 0.999), ...
                                'Molar fractions must sum up to unity.')
            assert(iteration > 0);
            
            if iteration == 1
                % Book-keeping
                state = model.setEquilibriumFields(state);
                state.eos.itCount = zeros(size(state.L));
            end

            L0 = state.L;
            K0 = state.K;
            Z0_V = state.Z_V;
            Z0_L = state.Z_L;
            x0 = state.x;
            y0 = state.y;
            
            % Only apply calculations for cells that have not converged yet
            if iteration == 1
                [stable, x0, y0, L_stable, V_stable] = PhaseStabilityTest(model, state.components, state.pressure, state.T);
                L0(L_stable) = 1;
                L0(V_stable) = 0;
                
                acf = model.fluid.acentricFactors;
                [Si_L, Si_V, A_L, A_V, B_L, B_V] = model.getMixtureFugacityCoefficients(P, T, x0, y0, acf);
                % Solve EOS for each phase
                Z0_L = model.computeLiquidZ(A_L, B_L);
                Z0_V = model.computeVaporZ(A_V, B_V);

                active = ~stable;
            else
                active = ~state.eos.converged;
            end
            
            if any(active)
                state.eos.itCount(active) = state.eos.itCount(active) + 1;
                L = L0(active);
                subset = @(x, active) cellfun(@(y) y(active, :), x, 'UniformOutput', false);

                z = subset(z, active);
                K = subset(K0, active);
                P = P(active);
                T = T(active);

                if model.performFlash
                    if model.useNewton
                        % Newton-based solver for equilibrium
                        [x, y, K, Z_L, Z_V, L, equilvals] = model.newtonCompositionUpdate(P, T, z, K, L);
                    else
                        % Successive substitution solver for equilibrium
                        [x, y, K, Z_L, Z_V, L, equilvals] = model.substitutionCompositionUpdate(P, T, z, K, L);
                    end
                else
                    [x, y, K, Z_L, Z_V, L, equilvals] = model.updateWithoutFlash(P, T, z, K, L);
                end
                values = max(equilvals, [], 1);
                conv = max(equilvals, [], 2) <= model.nonlinearTolerance;
                conv = conv & iteration > nonlinsolve.minIterations;
                resConv = values <= model.nonlinearTolerance & iteration > nonlinsolve.minIterations;

                isConverged = all(conv);
                % Insert back the local values into global arrays
                state.eos.converged(active) = conv;
                % Insert updated values in active cells
                [K_next, x_next, y_next] = deal(cell(1, ncomp));
                for i = 1:ncomp
                    K_next{i} = K0{i};
                    x_next{i} = x0{i};
                    y_next{i} = y0{i};

                    K_next{i}(active) = K{i};
                    x_next{i}(active) = x{i};
                    y_next{i}(active) = y{i};
                    L0(active) = L;

                    Z0_L(active ) = Z_L;
                    Z0_V(active ) = Z_V;

                end
            else
                K_next = K0;
                x_next = x0;
                y_next = y0;
                isConverged = true;
                values = zeros(1, ncomp);
                resConv = true(1, ncomp);
            end
            state.K = K_next;
            state.L = L0;
            state.x = x_next;
            state.y = y_next;
            state.Z_L = Z0_L;
            state.Z_V = Z0_V;

            state.mixing.K = state.K;

            failure = false;
            failureMsg = '';
            
            if model.verbose
                printConvergenceReport(model.fluid.names, values, isConverged, iteration);
            end
            report = model.makeStepReport(...
                            'Failure',      failure, ...
                            'FailureMsg',   failureMsg, ...
                            'Converged',    isConverged, ...
                            'ResidualsConverged', resConv, ...
                            'Residuals',    values);
            report.ActiveCells = sum(active);
        end
        
        function [x, y, K, Z_L, Z_V, L, values] = updateWithoutFlash(model, P, T, z, K, L)
            % Determine overall liquid fraction
            L = model.solveRachfordRice(L, K, z);
            % Compute liquid component fraction
            x = model.computeLiquid(L, K, z);
            % Vapor component fraction
            y = model.computeVapor(L, K, z);
            % Reduced values for pressure etc (divided by critical values)
            % [Pr, Tr] = model.getReducedPT(P, T);
            acf = (toCell(model.fluid.acentricFactors));
            [Si_L, Si_V, A_L, A_V, B_L, B_V] = model.getMixtureFugacityCoefficients(P, T, x, y, acf);
            % Solve EOS for each phase
            Z_L = model.computeLiquidZ(A_L, B_L);
            Z_V = model.computeVaporZ(A_V, B_V);
            values = zeros(model.G.cells.num, numel(x));
        end
        
        function [x, y, K, Z_L, Z_V, L, values] = substitutionCompositionUpdate(model, P, T, z, K, L)
            % Determine overall liquid fraction
            L = model.solveRachfordRice(L, K, z);
            % Compute liquid component fraction
            [x, sx] = model.computeLiquid(L, K, z);
            % Vapor component fraction
            y = model.computeVapor(L, K, z);
            % Reduced values for pressure etc (divided by critical values)
            % [Pr, Tr] = model.getReducedPT(P, T);
            acf = (toCell(model.fluid.acentricFactors));
            ncomp = numel(x);
            
            [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = model.getMixtureFugacityCoefficients(P, T, x, y, acf);
            % Solve EOS for each phase
            Z_L = model.computeLiquidZ(A_L, B_L);
            Z_V = model.computeVaporZ(A_V, B_V);
            
            % Compute fugacities
            f_L = model.computeFugacity(P, x, Z_L, A_L, B_L, Si_L, Bi);
            f_V = model.computeFugacity(P, y, Z_V, A_V, B_V, Si_V, Bi);
            
            % Pure liquid / pure vapor is converged
            ok = L == 1 | L == 0;
            f_r = cell(1, ncomp);
            % Compute fugacity ratios
            for i = 1:ncomp
                f_r{i} = sx.*f_L{i}./f_V{i};
                f_r{i}(ok | z{i} == 0) = 1;
            end
            % Update equilibrium constant estimates based on fugacity ratio
            values = abs([f_r{:}] - 1);
            for i = 1:ncomp
                K{i} = max(K{i}.*abs(f_r{i}), 1e-12);
                K{i}(~isfinite(K{i})) = 1;
            end
        end
        
        function [x, y, K, Z_L, Z_V, L, vals] = newtonCompositionUpdate(model, P, T, z, K, L)
            [x, y, K, Z_L, Z_V, L, vals] = newtonFugacityEquilibrium(model, P, T, z, K, L);
        end

        function state = setEquilibriumFields(model, state)
            n_cell = size(state.pressure, 1);
            ncomp = model.fluid.getNumberOfComponents();
            if ~isfield(state, 'L') 
                % Initial guess of 0.5
                state.L = repmat(0.5, n_cell, 1);
            end
            acf = (toCell(model.fluid.acentricFactors));
            % Estimate equilibrium constants using Wilson equation
            K_est = cell(1, ncomp);
            for i = 1:ncomp
                Tr = state.T./model.fluid.Tcrit(i);
                Pr = state.pressure./model.fluid.Pcrit(i);
                K_est{i} = exp(5.37.*(1 + acf{i}).*(1 - 1./Tr))./Pr;
                K_est{i}(~isfinite(K_est{i})) = 1000;
            end

            if ~isfield(state, 'K')
                state.K = K_est;
            else
                pure = state.L == 1 | state.L == 0;
                for i = 1:numel(state.K)
                    state.K{i}(pure) = K_est{i}(pure);
                end
            end
            
            state.Z_V = zeros(n_cell, 1);
            state.Z_L = zeros(n_cell, 1);
            
            if ~isfield(state, 'x')
                state.x = cellfun(@(x) 0*x, state.K, 'uniformoutput', false);
            end
            
            if ~isfield(state, 'y')
                state.y = cellfun(@(x) 0*x, state.K, 'uniformoutput', false);
            end
            
            state.eos.converged = false(n_cell, 1);
        end
        
        function [Pr, Tr] = getReducedPT(model, P, T)
            n = model.fluid.getNumberOfComponents();
            Tr = cell(1, n);
            Pr = cell(1, n);
            
            for i = 1:n
                Tr{i} = T./model.fluid.Tcrit(i);
                Pr{i} = P./model.fluid.Pcrit(i);
            end
        end
        
        function [A_ij, Bi] = getMixingParameters(model, P, T, acf)
            % Calculate intermediate values for fugacity computation
            ncomp = model.fluid.getNumberOfComponents();
            [Pr, Tr] = model.getReducedPT(P, T);
            if ~iscell(acf)
                acf = toCell(acf);
            end
            [Ai, Bi] = deal(cell(1, ncomp));
            [oA, oB] = deal(cell(1, ncomp));
            
            switch model.eosType
                case 1
                    % PR
                    for i = 1:ncomp
                        oA{i} = model.omegaA.*(1 + (0.37464 + 1.54226.*acf{i} - 0.26992.*acf{i}.^2).*(1-Tr{i}.^(1/2))).^2;
                    end
                case 2
                    % SRK
                    for i = 1:ncomp
                        oA{i} = model.omegaA.*(1 + (0.48 + 1.574.*acf{i} - 0.176.*acf{i}.^2).*(1-Tr{i}.^(1/2))).^2;
                    end
                case 3
                    % ZJ
                    error('Not implemented yet.')
                case 4
                    % RK
                    for i = 1:ncomp
                        oA{i} = model.omegaA.*Tr{i}.^(-1/2);
                    end
                otherwise
                    error('Unknown eos type: %d', model.eosType);
            end
            [oB{:}] = deal(model.omegaB);
            bic = model.fluid.getBinaryInteraction();
            A_ij = cell(ncomp, ncomp);

            for i = 1:ncomp
                Ai{i} = oA{i}.*Pr{i}./Tr{i}.^2;
                Bi{i} = oB{i}.*Pr{i}./Tr{i};
            end
            for i = 1:ncomp
                for j = 1:ncomp
                    A_ij{i, j} = (Ai{i}.*Ai{j}).^(1/2).*(1 - bic(i, j));
                end
            end
        end
        
        function [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = getMixtureFugacityCoefficients(model, P, T, x, y, acf)
            % Calculate intermediate values for fugacity computation
            [A_ij, Bi] = model.getMixingParameters(P, T, acf);
            [Si_L, A_L, B_L] = model.getPhaseMixCoefficients(x, A_ij, Bi);
            [Si_V, A_V, B_V] = model.getPhaseMixCoefficients(y, A_ij, Bi);
        end
        
        function eosdata = getPropertiesFastAD(model, P, T, x, y, Z_L, Z_V)
            % Get packed properties (fugacity, Z-factors) with FastAD type
            % to easily get derivatives
            P = double(P);
            T = double(T);
            x = cellfun(@double, x, 'UniformOutput', false);
            y = cellfun(@double, y, 'UniformOutput', false);
            
            [P, x{:}, y{:}] = initVariablesFastAD(P, x{:}, y{:});
            [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = model.getMixtureFugacityCoefficients(P, T, x, y, model.fluid.acentricFactors);
            
            if nargin < 6
                Z_L = model.computeLiquidZ(double(A_L), double(B_L));
                Z_V = model.computeVaporZ(double(A_V), double(B_V));

                Z_L = FastAD(Z_L, 0*P.jac);
                Z_V = FastAD(Z_V, 0*P.jac);
                
                Z_L = model.setZDerivatives(Z_L, A_L, B_L);
                Z_V = model.setZDerivatives(Z_V, A_V, B_V);
            end
            f_L = model.computeFugacity(P, x, Z_L, A_L, B_L, Si_L, Bi);
            f_V = model.computeFugacity(P, y, Z_V, A_V, B_V, Si_V, Bi);
            
            eosdata = struct();
            eosdata.Z_L = Z_L;
            eosdata.Z_V = Z_V;
            eosdata.f_L = f_L;
            eosdata.f_V = f_V;
        end
        
        function [Z_L, Z_V, f_L, f_V] = getProperties(model, P, T, x, y)
            [Si_L, Si_V, A_L, A_V, B_L, B_V, Bi] = model.getMixtureFugacityCoefficients(P, T, x, y, model.fluid.acentricFactors);

            Z_L = model.computeLiquidZ(double(A_L), double(B_L));
            Z_V = model.computeVaporZ(double(A_V), double(B_V));
            s = getSampleAD(P, T, x{:}, y{:});
            if isa(s, 'ADI')
                Z_L = double2ADI(Z_L, s);
                Z_V = double2ADI(Z_V, s);
            end
                
            Z_L = model.setZDerivatives(Z_L, A_L, B_L);
            Z_V = model.setZDerivatives(Z_V, A_V, B_V);
            
            f_L = model.computeFugacity(P, x, Z_L, A_L, B_L, Si_L, Bi);
            f_V = model.computeFugacity(P, y, Z_V, A_V, B_V, Si_V, Bi);
        end
        
        function [Z_L, Z_V, f_L, f_V, rep] = getEquilibriumProperties(model, P, T, x, y, Z_L, Z_V, packed)
            [s, isAD] = getSampleAD(P, x{:}, y{:});
            wantFugacity = nargout > 2;
            
            useFast = isAD && model.fastDerivatives;
            rep = struct('t_compressibility', 0, 't_mixture', 0, 't_fugacity', 0);
            if useFast
                if isempty(packed)
                    packed = model.getPropertiesFastAD(P, T, x, y);
                end
                [Z_L, Z_V, f_L, f_V] = deal(packed.Z_L, packed.Z_V, packed.f_L, packed.f_V);
            else
                [Z_L, Z_V, f_L, f_V] = model.getProperties(P, T, x, y);
            end
            
            if useFast
                % Fit derivatives into normal ADI structure by way of the
                % chain rule
                x0 = x;
                y0 = y;
                P0 = P;
                ncell = model.G.cells.num;

                ncomp = model.fluid.getNumberOfComponents();
                Z_L0 = double2ADI(double(Z_L), s);
                Z_V0 = double2ADI(double(Z_V), s);
                
                hasx = isa(x0{1}, 'ADI');
                hasy = isa(y0{1}, 'ADI');
                hasP = isa(P0, 'ADI');
                
                if wantFugacity
                    [f_L0, f_V0] = deal(cell(1, ncomp));
                    for i = 1:ncomp
                        f_L0{i} = double2ADI(double(f_L{i}), s);
                        f_V0{i} = double2ADI(double(f_V{i}), s);
                    end
                end
                vx = (1:ncell)';
                DM = @(x) makeDiagonal(x, vx, ncell);
                
                xoffset = 1;
                yoffset = ncomp + 1;
                for jacNo = 1:numel(s.jac)
                    if ~all(size(s.jac{jacNo}) == ncell)
                        continue
                    end
                    % Diagonal derivatives we want
                    [dZ_L, dZ_V] = deal(0);
                    if hasP
                        dP = diag(P0.jac{jacNo});
                        dZ_L = Z_L.jac(:, 1).*dP;
                        dZ_V = Z_V.jac(:, 1).*dP;
                    end
                    extractDiag = @(v) cellfun(@(z) diag(z.jac{jacNo}), v, 'UniformOutput', false);
                    
                    if hasx
                        dx = extractDiag(x0);
                        dX = [dx{:}];
                        xIx = (xoffset+1):(xoffset+ncomp);
                        if isa(Z_L, 'FastAD')
                            dZ_L = dZ_L + sum(Z_L.jac(:, xIx).*dX, 2);
                        end
                    end
                    if hasy
                        dy = extractDiag(y0);
                        dY = [dy{:}];
                        yIx = (yoffset+1):(yoffset+ncomp);
                        if isa(Z_V, 'FastAD')
                            dZ_V = dZ_V + sum(Z_V.jac(:, yIx).*dY, 2);
                        end
                    end
                    
                    if wantFugacity
                        for i = 1:ncomp
                            has_f_L = isa(f_L{i}, 'FastAD');
                            has_f_V = isa(f_V{i}, 'FastAD');
                            pi = 1;
                            df_L = 0;
                            df_V = 0;
                            if hasP
                                if has_f_L
                                    df_L = df_L + f_L{i}.jac(:, pi).*dP;
                                end
                                if has_f_V
                                    df_V = df_V + f_V{i}.jac(:, pi).*dP;
                                end
                            end

                            if hasx
                                if has_f_L
                                    df_L = df_L + sum(f_L{i}.jac(:, xIx).*dX, 2);
                                end
                            end
                            if hasy
                                if has_f_V
                                    df_V = df_V + sum(f_V{i}.jac(:, yIx).*dY, 2);
                                end
                            end

                            f_L0{i}.jac{jacNo} = DM(df_L);
                            f_V0{i}.jac{jacNo} = DM(df_V);
                            % Finally make AD
                        end
                    end
                    Z_L0.jac{jacNo} = DM(dZ_L);
                    Z_V0.jac{jacNo} = DM(dZ_V);
                end
                
                Z_L = Z_L0;
                Z_V = Z_V0;
                if wantFugacity
                    f_L = f_L0;
                    f_V = f_V0;
                end
            end
        end
        
        function [Si, A, B] = getPhaseMixCoefficients(model, x, A_ij, Bi)
            ncomp = numel(x);
            [A, B] = deal(0);
            Si = cell(1, ncomp);
            [Si{:}] = deal(0);
            for i = 1:ncomp
                B = B + x{i}.*Bi{i};
                for j = 1:ncomp
                    A_ijx_j = A_ij{i, j}.*x{j};
                    A = A + x{i}.*A_ijx_j;
                    Si{i} = Si{i} + A_ijx_j;
                end
            end
        end
                
        function [E2, E1, E0] = getCubicCoefficients(model, A, B)
            [m1, m2] = model.getEOSCoefficients();
            
            E0 = -(A.*B + m1.*m2.*B.^2.*(B+1));
            E1  = A - (m1 + m2 - m1.*m2).*B.^2 - (m1 + m2).*B;
            E2 = (m1 + m2 - 1).*B - 1;
        end
        
        function f = computeFugacity(model, p, x, Z, A, B, Si, Bi)
            % Compute fugacity based on EOS coefficients
            [m1, m2] = model.getEOSCoefficients();
            ncomp = model.fluid.getNumberOfComponents();
            f = cell(1, ncomp);
            a1 = -log(Z-B);
            b1 = log((Z + m2.*B)./(Z + m1.*B)).*(A./((m1-m2).*B));
            b2 = (Z-1)./B;
            for i = 1:ncomp
                v = a1 + b1.*(2.*Si{i}./A - Bi{i}./B) + Bi{i}.*b2;
                % v = -log(max(Z - B, 0)) + (A./((m1-m2).*B)).*(2.*Si{i}./A - Bi{i}./B).*log((Z + m2.*B)./(Z + m1.*B)) + Bi{i}./B.*(Z-1);
                f{i} = exp(v).*p.*x{i};
            end
        end
        
        
        function [Z_L, Z_V] = getCompressibility(model, state, P, T, x, y)
            [Z_L, Z_V] = getEquilibriumProperties(model, P, T, x, y, state.Z_L, state.Z_V, state.eos.packed);
        end
        
        function [eqs, f_L, f_V, Z_L, Z_V, rep] = equationsEquilibrium(model, P, T, x, y, z, L, Z_L, Z_V, K, packed)
            t1 = tic();
            if nargin < 10
                packed = [];
            end
            [Z_L, Z_V, f_L, f_V, rep] = model.getEquilibriumProperties(P, T, x, y, Z_L, Z_V, packed);
            
            ncomp = numel(x);
            timer = tic();
            eqs = cell(1, 2*ncomp + 1);
            emptyJac = 0*P + 0*z{1};
            
            isLiq = double(L) == 1;
            isVap = double(L) == 0;
            isPure = isLiq | isVap;
            eqs{end} = emptyJac + (isLiq | isVap);
            for i = 1:ncomp
                eqs{i} = z{i} - L.*x{i} - (1-L).*y{i} + emptyJac;
                if model.performFlash 
                    eqs{i+ncomp} = (f_V{i} - f_L{i}) + emptyJac;
                else
                    eqs{i+ncomp} = (x{i} - K{i}.*y{i}) + emptyJac;
                end
                eqs{end} = eqs{end} - ~isLiq.*y{i} + ~isVap.*x{i} + emptyJac;
                
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
            rep.t_assembly = toc(timer);
            rep.t_remainder = toc(t1) - (rep.t_fugacity + rep.t_mixture + rep.t_compressibility + rep.t_assembly);
        end

        function [x, y, LL, packed] = getPhaseFractionAsADI(model, state, pP, TP, z0)
            % Compute derivatives for values obtained by solving the
            % equilibrium equations (molar fractions in each phase, liquid
            % mole fraction).
            zP = z0;
            L = state.L;
            K = state.K;
            zD = cellfun(@double, zP, 'UniformOutput', false);
            % Compute secondary variables as pure doubles
            xD = state.x;
            yD = state.y;
            ncomp = numel(xD);
            ncell = numel(L);
            Z_L = state.Z_L;
            Z_V = state.Z_V;

            xS = cell(1, ncomp);
            yS = cell(1, ncomp);
            zS = zD;
            
            % Secondary as primary variables to get equilibrium jacobian
            % with respect to them.
            [xS{:}, yS{:}, LS] = initVariablesADI(xD{:}, yD{:}, L);

            xP = xD;
            yP = yD;

            pS = double(pP);
            TS = double(TP);
            LP = L;
            packed = state.eos.packed;
            eqsPrim  = model.equationsEquilibrium(pP, TP, xP, yP, zP, LP, Z_L, Z_V, K, packed);
            eqsSec = model.equationsEquilibrium(pS, TS, xS, yS, zS, LS, Z_L, Z_V, K, packed);
            
            ep = cat(eqsPrim{:});
            es = cat(eqsSec{:});
            
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

            if isa(pP, 'ADI')
                primVar = pP;
            else
                primVar = zP{1};
            end
            % We now have a big lumped matrix of all the derivatives, which
            % we can then extract into the regular ADI format of cell
            % arrays
            [I, J, V] = find(dsdp);
            jSz = cellfun(@(x) size(x, 2), primVar.jac);
            V(~isfinite(V)) = 0;

            jOffset = cumsum([0, jSz]);

            x = cell(1, ncomp);
            y = cell(1, ncomp);
            
            jump = numel(L);
            nj = numel(primVar.jac);
            for i = 1:ncomp
                xjac = cell(1, nj);
                yjac = cell(1, nj);


                startIx = jump*(i-1);
                endIx = jump*i;
                
                startIy = startIx + ncomp*jump;
                endIy = endIx + ncomp*jump;

                for j = 1:nj
                    startJ = jOffset(j);
                    endJ = jOffset(j+1);

                    actx = I > startIx & I <= endIx & J > startJ & J <= endJ;
                    acty = I > startIy & I <= endIy & J > startJ & J <= endJ;
                    
                    xjac{j} = sparse(I(actx) - startIx, J(actx) - startJ, V(actx), endIx-startIx, endJ-startJ);
                    yjac{j} = sparse(I(acty) - startIy, J(acty) - startJ, V(acty), endIy-startIy, endJ-startJ);
                end
                x{i} = ADI(xD{i}, xjac);
                y{i} = ADI(yD{i}, yjac);
            end
            
            Ljac = cell(1, nj);
            startI = 2*ncomp*ncell;
            endI = 2*ncomp*ncell + ncell;
            for i = 1:nj
                startJ = jOffset(i);
                endJ = jOffset(i+1);
                actL = I > startI & J > startJ & J <= endJ;
                Ljac{i} = sparse(I(actL) - startI, J(actL) - startJ, V(actL), endI-startI, endJ-startJ);
            end
            LL = ADI(L, Ljac);
        end
        
        function L = solveRachfordRice(model, L, K, z) %#ok
            % Solve Rachford Rice equations to find liquid and vapor
            % distribution.
            ncomp = numel(K);
            tmp1 = warning('query','MATLAB:nearlySingularMatrix');
            tmp2 = warning('query','MATLAB:singularMatrix');
            warning('off','MATLAB:nearlySingularMatrix')
            warning('off','MATLAB:singularMatrix')

            n_L = numel(L);
            L_final = L;
            active = true(n_L, 1);
            maxit = 25;
            tol = 1e-6;
            tolRes = 1e-12;
            for it = 1:maxit
                % fprintf('Iteration %d: %d active\n', it, nnz(active));
                L0 = L;
                L = initVariablesADI(L);
                
                eq = 0;
                for i = 1:ncomp
                    Ki = K{i}(active);
                    zi = z{i}(active);
                    v = ((Ki - 1).*zi)./(1 + (1-L).*(Ki - 1));
                    present = zi > 0;
                    eq = eq + v.*present;
                end
                J = -eq.jac{1};
                r = eq.val;
                vNorm = abs(r);
                % Converged values do not need to be updated
                convResidual = vNorm < tolRes;
                dL = ~convResidual.*mldivide(J, r);
                dL = sign(dL).*min(abs(dL), 0.2);
                L = double(L) + dL;
                L = max(L, 0);
                L = min(L, 1);
                dLNorm = abs(L - L0);
                
                conv = dLNorm < tol | convResidual;
                if it == maxit
                    disp('reached max iterations in Racheford rice')
                    conv = conv | true;
                end
                update = false(n_L, 1);
                update(active) = conv;
                
      
                L_final(update) = L(conv);
                active(update) = false;
                L = L(~conv);
                
                if all(conv)
                    break
                end
                
            end
            
            L = L_final;

            warning(tmp1.state,'MATLAB:nearlySingularMatrix')
            warning(tmp2.state,'MATLAB:singularMatrix')

        end
        
        function y = computeVapor(model, L, K, z)
            if isstruct(L)
                y = model.computeVapor(L.L, L.K, L.components);
                return
            end
            y = cellfun(@(k, zi) k.*zi./(L + (1-L).*k), K, z, 'UniformOutput', false);
            
            sv = 0;
            for i = 1:numel(y)
                sv = sv + double(y{i});
            end
            y = cellfun(@(x) x./sv, y, 'UniformOutput', false);
            
            assert(all(cellfun(@(x) all(isfinite(double(x))), y)));
        end
        
        function [x, sv] = computeLiquid(model, L, K, z)
            if isstruct(L)
                x = model.computeLiquid(L.L, L.K, L.components);
                return
            end
            x = cellfun(@(k, zi) zi./(L + (1-L).*k), K, z, 'UniformOutput', false);
           
            sv = 0;
            for i = 1:numel(x)
                sv = sv + double(x{i});
            end
            x = cellfun(@(x) x./sv, x, 'UniformOutput', false);
            assert(all(cellfun(@(x) all(isfinite(double(x))), x)));
        end
        
        function rho = computeDensity(model, p, x, Z, T, isLiquid)
            ncomp = numel(model.fluid.names);
            M = 0;
            for i = 1:ncomp
                M = M + x{i}.*model.fluid.molarMass(i);
            end
            R = 8.3144598;
            rho = p.*M./(R.*T.*Z);
        end
        
        function rho = computeMolarDensity(model, p, x, Z, T, isLiquid)
            ncomp = numel(model.fluid.names);
            R = 8.3144598;
            rho = p./(R.*T.*Z);
        end
        function [sL, sV] = computeSaturations(model, rhoL, rhoV, x, y, L, Z_L, Z_V)
            sL = L.*Z_L./(L.*Z_L + (1-L).*Z_V);
            sV = 1 - sL;
        end
        
        function mu = computeViscosity(model, P, rho, T, x, isLiquid)
            % Compute viscosity using the Lohrenz, Bray and Clark
            % correlation for hydrocarbon mixtures (LBC viscosity)
            ncomp = numel(x);
            molfactor = 1/gram;
            
            % We first compute an estimated low pressure (surface)
            % viscosity for the mixture
            MW = (molfactor*model.fluid.molarMass).^(1/2);
            [a, b] = deal(0);
            for i = 1:ncomp
                tr = T./model.fluid.Tcrit(i);
                Tc = model.fluid.Tcrit(i)./Rankine();
                Pc = model.fluid.Pcrit(i)./psia();

                mwi = MW(i);
                e_i = (5.4402*Tc.^(1/6))./(mwi.*Pc.^(2/3).*(centi*poise));
                
                large = double(tr) > 1.5;
                % Different estimates based on how far above we are from
                % the critical temp
                mu_i = (~large.*34e-5.*tr.^(0.94) + large.*17.78e-5.*(4.58*tr - 1.67).^0.625)./e_i;
                a = a + x{i}.*mu_i.*mwi;
                b = b + x{i}.*mwi;
            end
            % Final atmospheric / low pressure viscosity is the ratio of a and b
            mu_atm = a./b;
            % Compute critical properties and coefficient for final
            % expression
            [P_pc, T_pc, Vc, mwc] = model.computePseudoCriticalPhaseProperties(x);            
            e_mix = 5.4402*(T_pc./Rankine()).^(1/6)./((molfactor*mwc).^(1/2).*(P_pc./psia()).^(2/3).*(centi*poise));
            % Reduced density via definition of critical density 
            rhor = Vc.*rho./mwc;
            % Final adjusted viscosity at current conditions
            mu = mu_atm + ((0.1023 + 0.023364.*rhor + 0.058533.*rhor.^2 ...
                - 0.040758.*rhor.^3 + 0.0093324.*rhor.^4).^4 - 1e-4)./e_mix;
        end
        
        function [P_pc, T_pc, Vc, mw] = computePseudoCriticalPhaseProperties(model, x)
            % Molar fraction weighted aggregated properties to get
            % pseudocritical properties of mixture.
            ncomp = numel(x);
            [T_pc, P_pc, Vc, mw] = deal(0);
            for i = 1:ncomp
                T_pc = T_pc + model.fluid.Tcrit(i)*x{i};
                P_pc = P_pc + model.fluid.Pcrit(i)*x{i};
                mw = mw + model.fluid.molarMass(i)*x{i};
                Vc = Vc + model.fluid.Vcrit(i)*x{i};
            end
        end
        
        function frac = getMassFraction(model, molfraction)
            % Convert molar fraction to mass fraction
            ncomp = numel(molfraction);
            mass = cell(1, ncomp);
            totMass = 0;
            for i = 1:numel(molfraction)
                mass{i} = model.fluid.molarMass(i).*molfraction{i};
                totMass = totMass + mass{i};
            end
            frac = cellfun(@(x) x./totMass, mass, 'UniformOutput', false);
        end
        
        function Z = setZDerivatives(model, Z, A, B, cellJacMap)
            % Z comes from the solution of a cubic equation of state, so
            % the derivatives are not automatically computed. By
            % differentiating the cubic EOS manually and solving for dZ/du
            % where u is some primary variable, we can still obtain
            % derivatives without making any assumptions other than the EOS
            % being a cubic polynomial
            if nargin < 5
                if isa(Z, 'ADI')
                cellJacMap = cell(numel(Z.jac), 1);
                else
                    cellJacMap = {};
                end
            end
            [E2, E1, E0] = model.getCubicCoefficients(A, B);
            e2 = double(E2);
            e1 = double(E1);
            z = double(Z);
            if isa(Z, 'ADI')
                for i = 1:numel(Z.jac)
                    if isempty(cellJacMap{i})
                        [n, m] = size(Z.jac{i});
                        if n == 0 || m ~= numel(z);
                            % There are either no derivatives present or the
                            % derivatives are not of the right dimension
                            continue
                        end
                        dE2 = getJac(E2, i);
                        dE1 = getJac(E1, i);
                        dE0 = getJac(E0, i);

                        d = -(dE2.*z.^2 + dE1.*z + dE0)./(3*z.^2 + 2*z.*e2 + e1);
                        Z.jac{i} = sparse((1:n)', (1:m)', d, n, m);
                    else
                        map = cellJacMap{i};
                        zi = z(map);
                        e2i = e2(map);
                        e1i = e1(map);
                        
                        [n, m] = size(Z.jac{i});
                        if n == 0 || m ~= numel(zi);
                            % There are either no derivatives present or the
                            % derivatives are not of the right dimension
                            continue
                        end
                        dE2 = getJacSub(E2, i, map);
                        dE1 = getJacSub(E1, i, map);
                        dE0 = getJacSub(E0, i, map);

                        d = -(dE2.*zi.^2 + dE1.*zi + dE0)./(3*zi.^2 + 2*zi.*e2i + e1i);
                        Z.jac{i} = sparse(map, (1:m)', d, n, m);
                    end
                end
            elseif isa(Z, 'FastAD')
                for i = 1:size(Z.jac, 2)
                    dE2 = getJac(E2, i);
                    dE1 = getJac(E1, i);
                    dE0 = getJac(E0, i);
                    
                    Z.jac(:, i) = -(dE2.*z.^2 + dE1.*z + dE0)./(3*z.^2 + 2*z.*e2 + e1);
                end
            end
        end

        
    function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
        
        pure = state.L == 0 | state.L == 1;
        for i = 1:numel(state.K)
            state.K{i}(pure) = 1;
        end
        state.x = model.computeLiquid(state);
        state.y = model.computeVapor(state);
        
        L = state.L;
        for i = 1:numel(state.K)
            x = state.x{i};
            z = state.components{i};
            y = (z - L.*x)./(1-L);
            y(L == 1) = z(L == 1);

            state.y{i} = y;
            state.K{i} = y./x;
            state.K{i}(x == 0) = 1;
            state.K{i}(~isfinite(state.K{i})) = 1;
        end
        % Set derivatives
        state.eos.packed = model.getPropertiesFastAD(state.pressure, state.T, state.x, state.y);
    end

    end
end

function v = toCell(x)
    v = cell(size(x));
    for i = 1:numel(x)
        v{i} = x(i);
    end
end

function dx = getJac(x, ix)
    if isa(x, 'ADI')
        dx = diag(x.jac{ix});
    elseif isa(x, 'FastAD')
        dx = x.jac(:, ix);
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
%             dx = vv;
        end
    else
        dx = 0;
    end
end

function D = makeDiagonal(x, vx, ncell)
    D = sparse(vx, vx, x, ncell, ncell);
end

function checkZ(Z)
   % Give warning if bad values 
   if any(Z < 0)
       warning([num2str(nnz(Z<0)), ' negative Z-factors detected...']);
   end
   if any(~isfinite(Z))
       warning([num2str(nnz(~isfinite(Z))), ' non-finite Z-factors detected...']);
   end
end

function Z = fallbackCubic(A, E2, E1, E0)
    n = numel(A);
    Z = zeros(n, 3);
    for i = 1:n
        Z(i, :) = roots([1, E2(i), E1(i), E0(i)]);
    end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
