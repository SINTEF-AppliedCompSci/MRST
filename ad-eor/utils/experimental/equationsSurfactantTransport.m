function [problem, state] = equationsSurfactantTransport(state0, state, ...
                                                      model, dt, drivingForces, varargin)
%
%
% SYNOPSIS:
%   function [problem, state] = equationsSurfactantTransport(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: Assemble the linearized equations for an oil-water-surfactant
% system, computing both the residuals and the Jacobians. 
%    
% Only the mass conservation equation for the surfactant is considered. The
% pressure and saturations are set as constant.
%  
% This function is used to compute explicitly the concentration, when an
% implicit-explicit approach is taken.
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
% ImplicitExplicitOilWaterSurfactantModel, PressureSaturationSurfactantModel

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
                 'iteration', -1, ...
                 'solveExplicitly', true );
    opt = merge_options(opt, varargin{:});
    
    W = drivingForces.W;
    fluid = model.fluid;
    op    = model.operators;
    G = model.G;

    % Properties at current timestep
    [p, sW, cs, csmax, wellSol] = model.getProps(state, 'pressure', 'water', 'surfactant', ...
                                                      'surfactantmax', 'wellsol');

    % Properties at previous timestep
    [p0, sW0, cs0, csmax0] = model.getProps(state0, 'pressure', 'water', 'surfactant', 'surfactantmax');

    pBH   = vertcat(wellSol.bhp);
    qWs   = vertcat(wellSol.qWs);
    qOs   = vertcat(wellSol.qOs);
    qWSft = vertcat(wellSol.qWSft);

    % Initialize independent variables.
    if ~opt.resOnly,
        % ADI variables needed since we are not only computing residuals.
        if ~opt.reverseMode,
            [c_adi, qWSft] = initVariablesADI(c, qWSft);
        else
            % Not implemented
            % zw = zeros(size(pBH));
            % [p0, sW0, c0, zw, zw, zw, zw] = ...
            %     initVariablesADI(p0, sW0, c0, zw, zw, zw, zw); %#ok
            % clear zw
        end
    else
        c_adi = c;
    end

    % We will solve for pressure, water saturation (oil saturation follows via
    % the definition of saturations), surfactant concentration and well rates +
    % bhp.
    primaryVars = {'surfactant', 'qWSft'};

    c_impl = c_adi;
    if opt.solveExplicitly
        c_flux = c0;
    else
        c_flux = c_adi;
    end
    
    % EQUATIONS ---------------------------------------------------------------

    % Compute fluxes and other properties for oil and water.
    [dpO, dpW, mobO, mobW, upco, upcw, bO, bW, pvMult, bO0, bW0, pvMult0, T]  = ...
        computeFluxAndPropsOilWaterSurfactant(model, p0, p, sW, c_flux, pBH, W);

    mobSft = c_flux.*mobW;
    vSft   = -op.faceUpstr(upcw, mobSft).*T.*dpW;
    bWvSft = op.faceUpstr(upcw, bW).*vSft;

    poro = model.rock.poro;
    ads  = computeEffAds(cs_impl, csmax, fluid);
    ads0 = computeEffAds(cs0, csmax0, fluid);
    ads_term = fluid.rhoRSft.*((1-poro)./poro).*(ads - ads0);

    % Conservation of surfactant in water:
    surfactant = (op.pv/dt).*((pvMult.*bW.*sW.*c_impl - pvMult0.*bW0.*sW0.*c0) +  ads_term) + ...
        op.Div(bWvSft);

    if model.extraStateOutput
        sigma = fluid.ift(c_impl);
    end

    eqs   = {surfactant};
    names = {'surfactant'};
    types = {'cell'};

    % Well conditions
    if ~isempty(W)
        wm = model.wellmodel;
        if ~opt.reverseMode
            wc   = vertcat(W.cells);
            pw   = p(wc);
            rhos = [fluid.rhoWS, fluid.rhoOS];
            bw   = {bW(wc), bO(wc)};
            mw   = {mobW(wc), mobO(wc)};
            s    = {sW(wc), 1 - sW(wc)};

            [cqs, weqs, ctrleqs, wc, state.wellSol] = wm.computeWellFlux(model, W, wellSol, pBH, {qWs, ...
                                qOs}, pw, rhos, bw, mw, s, {}, 'nonlinearIteration', opt.iteration);

            % surfactant well equations
            [wSft, wciSft, iInxW] = getWellSurfactant(W);
            cw        = c_impl(wc);
            cw(iInxW) = wciSft;

            % Add surfactant
            bWqSft = cw.*cqs{1};
            eqs{1}(wc) = eqs{1}(wc) - bWqSft;

            % Well surfactant rate for each well is water rate in each perforation
            % multiplied with surfactant concentration in that perforated cell.
            perf2well = getPerforationToWellMapping(W);
            Rw = sparse(perf2well, (1:numel(perf2well))', 1, ...
                        numel(W), numel(perf2well));
            eqs{2} = qWSft - Rw*(cqs{1}.*cw);

            names(2) = {'surfactantWells'};
            types(2) = {'perf'};
            
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
