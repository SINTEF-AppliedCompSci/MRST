function [problem, state] = equationsPressureSaturationForOilWaterSurfactant(state0, state, model, dt, drivingForces, varargin)
%
%
% SYNOPSIS:
%   function [problem, state] = equationsPressureSaturationForOilWaterSurfactant(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: Assemble the linearized equations for an oil-water-surfactant
% system, computing both the residuals and the Jacobians. 
%    
% Only the mass conservation equations for the oil and water saturations are
% considered. The surfactant concentration is set as constant. 
%  
% This function is used to compute implicitly the pressure and saturation, when
% an implicit-explicit approach is taken
%
% The result is returned as an instance of the class LinearizedProblem which can
% be solved using instances of LinearSolverAD.
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
% SEE ALSO: LinearizedProblem, LinearSolverAD, equationsOilWaterSurfactant,
% ImplicitExplicitOilWaterSurfactantModel, PressureSaturationSurfactantModel,
% equationsSurfactantTransport

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
    [p, sW, cs, csmax, wellSol] = model.getProps(state, 'pressure', 'water', 'surfactant', ...
                                                      'surfactantmax', 'wellsol');

    % Properties at previous timestep
    [p0, sW0, cs0, csmax0] = model.getProps(state0, 'pressure', 'water', 'surfactant', 'surfactantmax');

    pBH    = vertcat(wellSol.bhp);
    qWs    = vertcat(wellSol.qWs);
    qOs    = vertcat(wellSol.qOs);

    % Initialize independent variables.
    if ~opt.resOnly,
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode,
            [p, sW, qWs, qOs, pBH] = initVariablesADI(p, sW, qWs, qOs, pBH);
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
    primaryVars = {'pressure', 'sW', 'qWs', 'qOs', 'bhp'};


    % EQUATIONS ---------------------------------------------------------------
    
    % Compute fluxes and other properties for oil and water.
    [dpO, dpW, mobO, mobW, upco, upcw, bO, bW, pvMult, bO0, bW0, pvMult0, T] = ...
        computeFluxAndPropsOilWaterSurfactant(model, p0, p, sW, c, pBH, W);

    % Upstream weight b factors and multiply by interface fluxes to obtain the
    % fluxes at standard conditions.
    vO   = -op.faceUpstr(upco, mobO).*T.*dpO;
    vW   = -op.faceUpstr(upcw, mobW).*T.*dpW;
    bOvO = op.faceUpstr(upco, bO).*vO;
    bWvW = op.faceUpstr(upcw, bW).*vW;

    % Conservation of mass for water
    water = (op.pv/dt).*( pvMult.*bW.*sW - pvMult0.*bW0.*sW0 ) + op.Div(bWvW);

    % Conservation of mass for oil
    sO  = 1 - sW;
    sO0 = 1 - sW0;
    oil = (op.pv/dt).*( pvMult.*bO.*sO - pvMult0.*bO0.*sO0 ) + op.Div(bOvO);

    eqs   = {water, oil};
    names = {'water', 'oil'};
    types = {'cell', 'cell'};

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
            eqs(3:4) = weqs;
            % Store the control equations (trivial equations ensuring that each
            % well will have values corresponding to the prescribed value)
            eqs{5} = ctrleqs;
            % Add source terms to the equations. Negative sign may be
            % surprising if one is used to source terms on the right hand side,
            % but this is the equations on residual form.
            eqs{1}(wc) = eqs{1}(wc) - cqs{1};
            eqs{2}(wc) = eqs{2}(wc) - cqs{2};

            names(3:5) = {'waterWells', 'oilWells', 'closureWells'};
            types(3:5) = {'perf', 'perf', 'well'};
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
