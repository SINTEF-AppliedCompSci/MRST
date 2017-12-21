classdef OverallCompositionCompositionalModel < ThreePhaseCompositionalModel
    % Overall composition model for compositional problems
    %
    % SYNOPSIS:
    %   model = ThreePhaseCompositionalModel(G, rock, fluid, compFluid)
    %
    % DESCRIPTION:
    %   The overall composition model relies on primary variables
    %   pressure and overall compositions for all components. As a
    %   consequence, the the phase behavior and mixing is offloaded to the
    %   equation of state directly. A requirement for this model is that
    %   the equation of state is fullfilled at each Newton iteration.
    %
    % PARAMETERS:
    %   G         - Grid structure
    %   rock      - Rock structure for the reservoir
    %   fluid     - The flow fluid, containing relative permeabilities,
    %               surface densities and flow properties for the
    %               aqueous/water phase (if present)
    %   compFluid - CompositionalFluid instance describing the species
    %               present.
    %
    % RETURNS:
    %  model - Initialized class instance
    %
    % SEE ALSO:
    %   `ThreePhaseCompositionalModel`, `NaturalVariablesCompositionalModel`

    properties
        
    end
    
    methods
        function model = OverallCompositionCompositionalModel(varargin)
            model = model@ThreePhaseCompositionalModel(varargin{:});
        end
        
        
        function [xM, yM, sL, sV, rhoL, rhoV, muL, muV, report] = computeTwoPhaseFlowProps(model, state, p, temp, z)
            isADI = isa(p, 'ADI') || isa(temp, 'ADI') || any(cellfun(@(x) isa(x, 'ADI'), z));
            report = struct();
            if isADI
                t1 = tic();
                [x, y, L] = model.EOSModel.getPhaseFractionAsADI(state, p, temp, z);
                report.t_derivatives = toc(t1);
                t2 = tic();
                [Z_L, Z_V] = model.EOSModel.getCompressibility(state, p, temp, x, y, z);
                report.t_compressibility = toc(t2);
            else
                [x, y, L, Z_L, Z_V] = model.getProps(state, 'x', 'y', 'L', 'Z_L', 'Z_V');
                x = expandMatrixToCell(x);
                y = expandMatrixToCell(y);
                [report.t_derivatives, report.t_compressibility] = deal(0);
            end
            
            t1 = tic();
            xM = model.EOSModel.getMassFraction(x);
            yM = model.EOSModel.getMassFraction(y);
            report.t_massfraction = toc(t1);
            
            t2 = tic();
            rhoL = model.PropertyModel.computeDensity(p, x, Z_L, temp, true);
            rhoV = model.PropertyModel.computeDensity(p, y, Z_V, temp, false);
            report.t_density = toc(t2);
            
            [sL, sV] = model.EOSModel.computeSaturations(rhoL, rhoV, x, y, L, Z_L, Z_V);
            
            t3 = tic();
            if nargout > 6
                muL = model.PropertyModel.computeViscosity(p, x, Z_L, temp, true);
                muV = model.PropertyModel.computeViscosity(p, y, Z_V, temp, false);
            end
            report.t_viscosity = toc(t3);
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsCompositional(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            state0 = state;
            
            var0 = problem.primaryVariables;
            vars = var0;
            removed = false(size(vars));
            
            wix = strcmpi(vars, 'sW');
            if any(wix)
                state = model.updateStateFromIncrement(state, dx{wix}, problem, 'sW', inf, model.dsMaxAbs);
                removed(wix) = true;
                vars = vars(~wix);
            end
            
            % Components
            cnames = model.EOSModel.fluid.names;
            ncomp = numel(cnames);
            ok = false(ncomp, 1);

            z = state.components;
            rm = 0;

            for i = 1:ncomp
                name = lower(cnames{i});
                cix = strcmpi(var0, name);
                if any(cix)
                    z0 = z(:, i);

                    dz = dx{cix};
                    if isfinite(model.dzMaxAbs)
                        dz = sign(dz).*min(abs(dz), model.dzMaxAbs);
                    end
                    z(:, i) = min(max(z0 + dz, 0), 1);

                    ok(i) = true;
                    [vars, ix] = model.stripVars(vars, {name});
                    removed(~removed) = removed(~removed) | ix;

                    rm = rm - (z(:, i) - z0);
                end
            end
            
            if any(ok)
                % We had components as active variables somehow
                assert(nnz(~ok) == 1)
                z(:, ~ok) = min(max(z(:, ~ok) + rm, 0), 1);
                z = bsxfun(@rdivide, z, sum(z, 2));
                state.components = z;
                if model.water
                    v  = 1 - state.s(:, 1);
                    v0 = 1 - state0.s(:, 1);
                else
                    [v, v0] = deal(1);
                end

                state.dz = computeChange(z, state0.components, v, v0);
            end

            % Parent class handles almost everything for us
            problem.primaryVariables = vars;
            dx(removed) = [];
            [state, report] = updateState@ThreePhaseCompositionalModel(model, state, problem, dx, drivingForces);
            
            if problem.iterationNo == 1
                state.switched = false(model.G.cells.num, 1);
                state.switchCount = zeros(model.G.cells.num, 1);
            end
            twoPhase0 = state.L < 1 & state.L > 0;
            % Update saturations etc using flash routine
            state = model.computeFlash(state, problem.dt, problem.iterationNo);
            twoPhase = state.L < 1 & state.L > 0;
            switched = twoPhase0 ~= twoPhase;
            
            dispif(model.verbose > 1, '%d gas, %d oil, %d two-phase\n', nnz(state.L == 0), nnz(state.L == 1), nnz(twoPhase));
            state.switchCount = state.switchCount + double(switched);
            
            state.components = ensureMinimumFraction(state.components);
            state.x = ensureMinimumFraction(state.x);
            state.y = ensureMinimumFraction(state.y);
        end
    end
end

function dz = computeChange(z, z0, s, s0)
    z_prev = bsxfun(@times, s0, z0);
    z_curr = bsxfun(@times, s, z);
    dz = max(abs(z_curr - z_prev), [], 1);
end
