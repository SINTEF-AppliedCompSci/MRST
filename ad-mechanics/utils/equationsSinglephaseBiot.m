function [eqs, state] = equationsSinglephaseBiot(p0, p, qWs, pBH, state, model, dt, mechTerm, ...
                                           drivingForces, varargin)

    % Note that state is given only for output
    opt = struct('iteration', -1, ...
                 'resOnly', false); % just to avoid warning
    opt = merge_options(opt, varargin{:});

    W = drivingForces.W;

    s = model.operators;
    G = model.G;
    f = model.fluid;
    rock = model.rock;

    %grav  = gravity;
    %gdz   = s.Grad(G.cells.centroids) * model.gravity';
    gdz   = s.Grad(G.cells.centroids) * model.getGravityVector()';
    %--------------------
    %check for p-dependent tran mult:
    trMult = 1;
    if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

    %check for p-dependent porv mult:
    pvMult = 1; pvMult0 = 1;
    if isfield(f, 'pvMultR')
        pvMult =  f.pvMultR(p);
        pvMult0 = f.pvMultR(p0);
    end
    transMult=1;
    if isfield(f, 'transMult')
        transMult=f.transMult(p);
    end

    trans=s.T.*transMult;
    % -------------------------------------------------------------------------
    % water props (calculated at oil pressure OK?)
    bW     = f.bW(p);
    rhoW   = bW.*f.rhoWS;
    % rhoW on face, avarge of neighboring cells (E100, not E300)
    rhoWf  = s.faceAvg(rhoW);
    mobW   = trMult./f.muW(p);
    dpW     = s.Grad(p) - rhoWf.*gdz;
    % water upstream-index
    upcw = (double(dpW)<=0);
    vW = - s.faceUpstr(upcw, mobW).*trans.*dpW;
    bWvW = s.faceUpstr(upcw, bW).*vW;

    if model.outputFluxes
        state = model.storeFluxes(state, vW, [], []);
    end

    if model.extraStateOutput
        state = model.storebfactors(state, bW, [], []);
        state = model.storeMobilities(state, mobW, [], []);
        state = model.storeUpstreamIndices(state, upcw, [], []);
    end

    % EQUATIONS ---------------------------------------------------------------

    eqs{1} = (1 ./ dt) .*                                                            ...
             ((rock.poro .* (G.cells.volumes .* pvMult)  + rock.alpha .* mechTerm.new) .* bW -      ...
              (rock.poro .* (G.cells.volumes .* pvMult0) + rock.alpha .* mechTerm.old) .* f.bW(p0)) ...
             + s.Div(bWvW);

    eqs = addFluxesFromSourcesAndBC(model, eqs, ...
                                    {p},...
                                    {rhoW},...
                                    {mobW}, ...
                                    {bW},  ...
                                    {ones(numel(p0),1)}, ...
                                    drivingForces);

    % well equations
    if ~isempty(W)
        wc    = vertcat(W.cells);
        pw   = p(wc);
        rhos = [f.rhoWS];
        bw   = {bW(wc)};
        mw   = {mobW(wc)};
        s = {1};

        wm = WellModel();
        [cqs, weqs, ctrleqs, wc, state.wellSol]  = wm.computeWellFlux(model, ...
                                                          W, state.wellSol, ...
                                                          pBH, {qWs}, pw, rhos, ...
                                                          bw, mw, s, {}, ...
                                                          'nonlinearIteration', ...
                                                          opt.iteration, ...
                                                          'referencePressureIndex', ...
                                                          1);
        eqs(2) = weqs;
        eqs{3} = ctrleqs;

        eqs{1}(wc) = eqs{1}(wc) - cqs{1};
    end
end
