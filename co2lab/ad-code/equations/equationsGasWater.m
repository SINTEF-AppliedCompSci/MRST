function [problem, state] = equationsGasWater(model, state0, state, dt, drivingForces, varargin)
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

   
    opt = struct('Verbose'     , mrstVerbose , ...
                 'reverseMode' , false       , ...
                 'resOnly'     , false       , ...
                 'iteration'   , -1          , ...
                 'stepOptions' , []); % compatibility only
    opt = merge_options(opt, varargin{:});
    
    assert(isempty(drivingForces.src));  % unsupported
    W = drivingForces.W;
    s = model.operators;
    f = model.fluid;
    G = model.G;
    t = model.t;
    
    % Extract current and previous values of all variables to solve for
    [p, sG, wellSol] = model.getProps(state, 'pressure', 'sg', 'wellsol');
    [p0, sG0]        = model.getProps(state0, 'pressure', 'sg');
    
    bhp = vertcat(state.wellSol.bhp);
    qWs = vertcat(state.wellSol.qWs);
    qGs = vertcat(state.wellSol.qGs);

    % ------------------ Initialization of independent variables ------------------
    
    if ~opt.resOnly
        if ~opt.reverseMode
            [p, sG, qWs, qGs, bhp] = initVariablesADI(p, sG, qWs, qGs, bhp);
        else
            [p0, sG0, ~, ~, ~] = initVariablesADI(p0, sw0, 0*qWs, 0*qGs, 0*bhp);
        end
    end
        
    % ----------------------------------------------------------------------------
    
    % Check for p-dependent tran mult:
    trMult = 1;
    if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end;
    
    % Check for p-dependent porv mult:
    pvMult = 1; pvMult0 = 1;
    if isfield(f, 'pvMultR')
        pvMult  = f.pvMultR(p);
        pvMult0 = f.pvMultR(p0);
    end
    transMult = 1;
    if isfield(f, 'transMult')
        transMult = f.transMult(p);
    end

    trans = s.T .* transMult;
    
    % Check for capillary pressure
    pcWG = 0;
    if isfield(f, 'pcWG')
        pcWG = f.pcWG(sG);
    end
    % ----------------------------------------------------------------------------
    
    [krW, krG] = deal(f.krW(1-sG), f.krG(sG));
    
    % computing densities, mobilities and upstream indices
    [bW, mobW, fluxW] = compMFlux(p - pcWG, f.bW, f.muW, f.rhoWS, trMult, krW, s, trans, model);
    [bG, mobG, fluxG] = compMFlux(p,        f.bG, f.muG, f.rhoGS, trMult, krG, s, trans, model);

    bW0 = f.bW(p0 - pcWG, t);
    bG0 = f.bG(p0, t);
    
    % --------------------------- Continuity equations ---------------------------
    
    eqs{1} = (s.pv/dt) .* (pvMult .* bW .* (1-sG) - pvMult0 .* bW0 .* (1-sG0)) + s.Div(fluxW);
    eqs{2} = (s.pv/dt) .* (pvMult .* bG .* sG     - pvMult0 .* bG0 .* sG0    ) + s.Div(fluxG);
    
    % ---------------------------- Boundary conditions ----------------------------

    [bc_cells, b_fW, b_fG] = ...
        BCFluxes(G, s, p, t, f, drivingForces.bc, krW, krG, transMult, trMult);

    % @@ If a cell has more than one pressure boundary condition, only one
    % will be taken into account.
    eqs{1}(bc_cells) = eqs{1}(bc_cells) + b_fW;
    eqs{2}(bc_cells) = eqs{2}(bc_cells) + b_fG;
    
    % ------------------------------ Well equations ------------------------------

    primaryVars = {'pressure' , 'sG'   , 'qWs'      , 'qGs', 'bhp'};
    types = {'cell'           , 'cell' , 'perf'       , 'perf'     , 'well'};
    names = {'water'          , 'gas'  , 'waterWells' , 'gasWells' , 'closureWells'};

    sW = 1-sG;
    if ~isempty(W)
       wm = model.wellmodel;
       if ~opt.reverseMode
          wc = vertcat(W.cells);
          [cqs, weqs, ctrleqs, wc, state.wellSol] =                       ...
              wm.computeWellFlux(model, W, wellSol, bhp                 , ...
                                 {qWs, qGs}                             , ...
                                 p(wc)                                  , ...
                                 [model.fluid.rhoWS, model.fluid.rhoGS] , ...
                                 {bW(wc), bG(wc)}                       , ...
                                 {mobW(wc), mobG(wc)}                   , ...
                                 {sW(wc), sG(wc)}                       , ...
                                 {}                                     , ...
                                 'allowControlSwitching', false         , ...
                                 'nonlinearIteration', opt.iteration);
          % Store the separate well equations (relate well bottom hole
          % pressures to influx)
          eqs(3:4) = weqs;
          
          % Store the control equations (ensuring that each well has values
          % corresponding to the prescribed value)
          eqs{5} = ctrleqs;
          
          % Add source term to equations.  Negative sign may be surprising if
          % one is used to source terms on the right hand side, but this is
          % the equations on residual form
          eqs{1}(wc) = eqs{1}(wc) - cqs{1};
          eqs{2}(wc) = eqs{2}(wc) - cqs{2};
          
       else
          [eqs(3:5), names(3:5), types(3:5)] = ...
              wm.createReverseModeWellEquations(model, state0.wellSol, p0);%#ok
       end
    else
       %eqs(3:5) = {bhp, bhp, bhp}; % empty ADIs
    end

    % ----------------------------------------------------------------------------

    if isempty(W)
       % Remove names/types associated with wells, as no well exist
       types = types(1:2);
      names = names(1:2);
    end
    
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
    
end

% ============================= END MAIN FUNCTION =============================

function [cells, fluxW, fluxG] = BCFluxes(G, s, p, t, f, bc, krW, krG, transMult, trMult)
    
    if isempty(bc)
        % all boundary conditions are no-flow; nothing to do here.
        cells = []; fluxW = []; fluxG = []; return;
    end

    assert(all(strcmp(bc.type, 'pressure'))); % only supported type for now
    cells = sum(G.faces.neighbors(bc.face, :), 2);
    %assert(numel(unique(cells)) == numel(cells)); % multiple BC per cell not supported
    
    % prepare boundary-specific values
    bp   = p(cells); 
    bt   = t(cells);
    bdp  = bc.value - bp;
    bdz  = G.faces.centroids(bc.face, 3) - G.cells.centroids(cells,3);
    
    assert(isscalar(transMult)); % not implemented for face-wise values
    assert(isscalar(trMult));    % not implemented for face-wise values
    trans = s.T_all(bc.face) .* transMult; 
    krWf = trMult * krW(cells);
    krGf = trMult * krG(cells);
    
    g = norm(gravity); 
    
    [bw, rhoWf, mobW] = computeRhoMobBFace(bp, bt, f.bW, f.rhoWS, krWf, f.muW);
    [bg, rhoGf, mobG] = computeRhoMobBFace(bp, bt, f.bG, f.rhoGS, krGf, f.muG);
    
    dptermW = bdp - rhoWf .* g .* bdz;
    dptermG = bdp - rhoGf .* g .* bdz;
    
    % Adjust upstream-weighted mobilities to prevent gas from re-entering the domain
    ix = dptermG > 0;
    mobG(ix) = 0;
    mobW(ix) = trMult ./ f.muW(bp(ix), bt(ix));
 
    % compute fluxes
    fluxW = - bw .* mobW .* trans .* dptermW;
    fluxG = - bg .* mobG .* trans .* dptermG;
    
    %fprintf('%f, %f\n', max(abs(fluxW.val)), max(abs(fluxG.val)));
end

% ----------------------------------------------------------------------------

function [bb, brho, bmob] = computeRhoMobBFace(bp, bt, bfun, rhoS, krf, mufun)

    bb   = bfun(bp, bt);
    brho = bb .* rhoS;
    bmob = krf ./ mufun(bp, bt);

end

% ----------------------------------------------------------------------------
function [b, mob, flux] = compMFlux(p, bfun, mufun, rhoS, trMult, kr, s, trans, model)
    
    g   = norm(gravity);
    b   = bfun(p, model.t);
    mob = trMult .* kr ./ mufun(p, model.t);

    rhoFace = s.faceAvg(b*rhoS);
    
    dp   = s.Grad(p) - rhoFace .* model.getGravityGradient();
    upc  = (double(dp)<=0);
    flux = -s.faceUpstr(upc, b .* mob) .* trans .* dp;

end

% ----------------------------------------------------------------------------
function [wc, cqs] = checkForRepetitions(wc, cqs)
[c, ~, ic] = unique(wc, 'stable');
if numel(c) ~= numel(wc)
    A = sparse(ic, (1:numel(wc))', 1, numel(c), numel(wc));
    wc = c;
    for k=1:numel(cqs)
        cqs{k} = A*cqs{k};
    end
end
end
