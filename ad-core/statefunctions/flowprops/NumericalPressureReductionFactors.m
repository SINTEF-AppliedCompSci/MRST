classdef NumericalPressureReductionFactors < StateFunction
    properties
        useDiagonalReduction
        includeDerivatives = true;
    end
    methods
        
        function prf = NumericalPressureReductionFactors(model, varargin)
            prf@StateFunction(model, varargin{:});
            AD = model.AutoDiffBackend;
            prf.useDiagonalReduction = isa(AD, 'DiagonalAutoDiffBackend') && AD.useMex;
        end
        
        function weights = evaluateOnDomain(prop, model, state) %#ok
            % Get state with the full primary variable set as AD
            [stateAD, ncell, ncomp] = getStateAD(model, state);
            % Timestep and mass at current and previous timestep
            dt    = state.reductionFactorProps.dt;
            mass  = model.getProps(stateAD, 'ComponentTotalMass');
            mass0 = state.reductionFactorProps.mass0;
            % Pressure at current and previous iteration
            p      = value(state.pressure);
            % Weights at current and previous iteration
            w      = getWeights(prop, mass, mass0, dt, ncell, ncomp);
            w_prev = state.reductionFactorProps.weights;
            % Get derivatives of weights wrt pressure
            next = struct('pressure', p, 'mass', {mass});
            dwdp = getWeightDerivatives(w, w_prev, next, state.reductionFactorProps);
            % Construct AD weights
            pAD = state.pressure;
            isAD = isa(pAD, 'ADI') && prop.includeDerivatives;
            weights = cell(ncomp, 1);
            for i = 1:ncomp
                wi = w(:, i);
                if isAD
                    dpi = dwdp(:, i);
                    if any(dpi)
                        Wp = model.AutoDiffBackend.convertToAD(wi, pAD);
                        Wp.jac{1} = sparse(1:ncell, 1:ncell, dpi, ncell, ncell);
                    else
                        Wp = wi;
                    end
                end
                weights{i} = Wp;
            end
        end
        
        function w = getWeights(prp, mass, mass0, dt, ncell, ncomp)
            % Compute mass difference
            acc = cellfun(@(m, m0) (m - m0)./dt, mass, mass0, 'UniformOutput', false);
            isD = cellfun(@isnumeric, acc);
            if all(isD)
                w = ones(ncell, ncomp);
                return;
            end
            ndof = ncell*ncomp;
            if prp.useDiagonalReduction
                % Special case: We can assume that all jacobians are diagonal.
                diags = cellfun(@(x) x.jac{1}.diagonal', acc, 'UniformOutput', false);
                M = vertcat(diags{:});
                ncell = size(M, 2);
                sz = repmat(ncomp, ncell, 1);
                Mi = invv(reshape(M, [], 1), sz);
                [i, j] = blockDiagIndex(sz, sz);
                b = zeros(ndof, 1);
                scale = 1/barsa;
                 % scale = sum(M(1:ncomp:end-ncomp+1, :), 1);
                b(1:ncomp:end-ncomp+1) = scale;
                M_inv = sparse(i, j, Mi, ndof, ndof);
                w = M_inv*b;
            else
                c = combineEquations(acc{:});
                % Get mass difference Jacobians
                A = c.jac{1};
                Safeguard against singular system (typically incomp/ weakly incomp)
                [~, jj, v] = find(A);
                ix  = jj <= ncell;
                tol = 1e-16;
                if ~any(ix) || norm(v(ix), inf) < tol
                    A = A + sparse(1:ncell, 1:ncell, 1, ndof, ndof);
                end
                b = zeros(1, ndof);
                b(1:ncell) = 1/barsa;
                % Compute weights
                w = b/A;
            end
            w = reshape(w', [], ncomp);
        end
    end
end

function [state, ncell, ncomp] = getStateAD(model, state)
    % Remove existing AD
    state = model.reduceState(state, true);
    % Get variables
    [vars, names, origin] = model.getPrimaryVariables(state);
    isP   = strcmp(names, 'pressure');
    origP = origin{isP};
    isW   = cellfun(@(x) ~strcmp(x, origP), origin);
    % Initialize AD
    [vars{~isW}] = model.AutoDiffBackend.initVariablesAD(vars{~isW});
    state = model.initStateAD(state, vars, names, origin);
    % Get various sizes needed later
    ncell = model.G.cells.num;
    ncomp = model.getNumberOfComponents();
end

function dwdp = getWeightDerivatives(w, w0, next, prev)
    % Compute weight derivatives by numerical differentiation
    if ~isempty(w0)
        p0 = prev.pressure;
        p  = next.pressure;
        dp = p - p0;
        dw = w - w0;
        dwdp = bsxfun(@rdivide, dw, dp);
        dwdp(~isfinite(dwdp)) = 0;
    else
        dwdp = 0*w;
    end
end

