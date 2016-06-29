function [problem, state] = equationsOilWaterSurfactant(state0, state, model, ...
                                                      dt, drivingForces, ...
                                                      varargin)
%
%
% SYNOPSIS:
%   function [problem, state] = equationsOilWaterSurfactant(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: Assemble the linearized equations for an oil-water-surfactant
% system, computing both the residuals and the Jacobians. Returns the result as
% an instance of the class LinearizedProblem which can be solved using instances
% of LinearSolverAD.
%
% A description of the modeling equations can be found in the directory
% ad-eor/docs.
%
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD
%

    opt = struct('Verbose', mrstVerbose, ...
                 'reverseMode', false, ...
                 'resOnly', false, ...
                 'iteration', -1 );
    opt = merge_options(opt, varargin{:});

    W     = drivingForces.W;
    fluid = model.fluid;
    op    = model.operators;
    G     = model.G;

    % Properties at current timestep
    [p, sW, c, cmax, wellSol] = model.getProps(state, 'pressure', 'water', 'surfactant', ...
                                                      'surfactantmax', 'wellsol');

    % Properties at previous timestep
    [p0, sW0, c0, cmax0] = model.getProps(state0, 'pressure', 'water', 'surfactant', 'surfactantmax');

    pBH   = vertcat(wellSol.bhp);
    qWs   = vertcat(wellSol.qWs);
    qOs   = vertcat(wellSol.qOs);
    qWSft = vertcat(wellSol.qWSft);

    % Initialize independent variables.
    if ~opt.resOnly,
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode,
            [p, sW, c, qWs, qOs, qWSft, pBH] = initVariablesADI(p, sW, c, qWs, qOs, qWSft, pBH);
        else
            % Not implemented yet
            % zw = zeros(size(pBH));
            % [p0, sW0, c0, zw, zw, zw, zw] = ...
            %     initVariablesADI(p0, sW0, c0, zw, zw, zw, zw); %#ok
            % clear zw
        end
    end

    % We will solve for pressure, water saturation (oil saturation follows via
    % the definition of saturations), surfactant concentration and well rates +
    % bhp.
    primaryVars = {'pressure', 'sW', 'surfactant', 'qWs', 'qOs', 'qWSft', 'bhp'};



    % EQUATIONS ---------------------------------------------------------------

    % Compute fluxes and other properties for oil and water.
    [dpO, dpW, mobO, mobW, upco, upcw, bO, bW, pvMult, bO0, bW0, pvMult0, T] = ...
        computeFluxAndPropsOilWaterSurfactant(model, p0, p, sW, c, pBH, W);

    % Upstream weight b factors and multiply by interface fluxes to obtain the
    % fluxes at standard conditions.
    vO     = -op.faceUpstr(upco, mobO).*T.*dpO;
    vW     = -op.faceUpstr(upcw, mobW).*T.*dpW;
    bOvO   = op.faceUpstr(upco, bO).*vO;
    bWvW   = op.faceUpstr(upcw, bW).*vW;

    % Conservation of mass for water
    water = (op.pv/dt).*(pvMult.*bW.*sW - pvMult0.*bW0.*sW0) + op.Div(bWvW);

    % Conservation of mass for oil
    sO  = 1 - sW;
    sO0 = 1 - sW0;
    oil = (op.pv/dt).*(pvMult.*bO.*sO - pvMult0.*bO0.*sO0) + op.Div(bOvO);

    % Computation of adsoprtion term
    poro = model.rock.poro;
    ads  = computeEffAds(c, cmax, fluid);
    ads0 = computeEffAds(c0, cmax0, fluid);
    ads_term = fluid.rhoRSft.*((1-poro)./poro).*(ads - ads0);

    % Conservation of surfactant in water:
    mobSft = mobW.*c;
    vSft   = -op.faceUpstr(upcw, mobSft).*T.*dpW;
    bWvSft = op.faceUpstr(upcw, bW).*vSft;
    surfactant = (op.pv/dt).*((pvMult.*bW.*sW.*c - pvMult0.*bW0.*sW0.*c0) +  ads_term) + ...
        op.Div(bWvSft);

    if model.extraStateOutput
        sigma = fluid.ift(c);
    end

    eqs   = {water, oil, surfactant};
    names = {'water', 'oil', 'surfactant'};
    types = {'cell', 'cell', 'cell'};

    % Well conditions
    if ~isempty(W)
        wm = model.wellmodel;
        if ~opt.reverseMode
            wc   = vertcat(W.cells);
            pw   = p(wc);
            rhos = [fluid.rhoWS, fluid.rhoOS];
            bw   = {bW(wc), bO(wc)};
            mw   = {mobW(wc), mobO(wc)};
            s    = {sW(wc), sO(wc)};

            [cqs, weqs, ctrleqs, wc, state.wellSol] = wm.computeWellFlux(model, W, wellSol, pBH, {qWs, ...
                                qOs}, pw, rhos, bw, mw, s, {}, 'nonlinearIteration', opt.iteration);

            % Store the well equations (relate well bottom hole pressures to
            % influx).
            eqs(4:5) = weqs;
            % Store the control equations (trivial equations ensuring that each
            % well will have values corresponding to the prescribed value)
            eqs{7} = ctrleqs;
            % Add source terms to the equations. Negative sign may be
            % surprising if one is used to source terms on the right hand side,
            % but this is the equations on residual form.
            eqs{1}(wc) = eqs{1}(wc) - cqs{1};
            eqs{2}(wc) = eqs{2}(wc) - cqs{2};

            % surfactant well equations
            [wSft, wciSft, iInxW] = getWellSurfactant(W);
            cw        = c(wc);
            cw(iInxW) = wciSft;

            % Add surfactant
            bWqSft = cw.*cqs{1};
            eqs{3}(wc) = eqs{3}(wc) - bWqSft;

            % Well surfactant rate for each well is water rate in each perforation
            % multiplied with surfactant concentration in that perforated cell.
            perf2well = getPerforationToWellMapping(W);
            Rw = sparse(perf2well, (1:numel(perf2well))', 1, ...
                        numel(W), numel(perf2well));
            eqs{6} = qWSft - Rw*(cqs{1}.*cw);

            names(4:7) = {'waterWells', 'oilWells', 'surfactantWells', 'closureWells'};
            types(4:7) = {'perf', 'perf', 'perf', 'well'};
        else
            % not implemented yet
            % [eq, n, typ] = ...
            %     wm.createReverseModeWellEquations(model, state0.wellSol, p0);
            % % Add another equation for surfactant well rates
            % [eqs{4:7}] = deal(eq{1});
            % [names{4:7}] = deal(n{1});
            % [types{4:7}] = deal(typ{1});
        end
    else
        error('The surfactant model does not support senarios without wells now!');
    end

    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

function [wSft, wciSft, iInxW] = getWellSurfactant(W)
    if isempty(W)
        wSft = [];
        wciSft = [];
        iInxW = [];
        return
    end
    inj   = vertcat(W.sign)==1;
    surfactInj = cellfun(@(x)~isempty(x), {W(inj).surfact});
    wSft = zeros(nnz(inj), 1);
    wSft(surfactInj) = vertcat(W(inj(surfactInj)).surfact);
    wciSft = rldecode(wSft, cellfun(@numel, {W(inj).cells}));

    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw    = numel(W);
    perf2well = rldecode((1:nw)', nPerf);
    compi = vertcat(W.compi);
    iInx  = rldecode(inj, nPerf);
    iInx  = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end
