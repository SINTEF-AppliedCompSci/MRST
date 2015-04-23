function [state, report] = incompMultiscale(state, CG, T, fluid, basis, varargin)
    opt = struct('getSmoother',    [], ...
                 'iterations',  0,...
                 'tolerance',   1e-6, ...
                 'useGMRES',    false, ...
                 'LinSolve',    @mldivide, ...
                 'Verbose',     true, ...
                 'reconstruct', true);
    [opt, incompOpt] = merge_options(opt, varargin{:});
    
    G = CG.parent;
    nc = G.cells.num;
    
    [A, q] = getSystemIncompTPFA(state, G, T, fluid, incompOpt{:});
    
    if numel(q) > nc
        [A, q, A_ww, A_wp, q_w] = eliminateWellEquations(A, q, nc);
        
        recover = @(p) recoverWellSolution(A_ww, A_wp, q_w, p);
    else
        recover = @(p) p;
    end
    
    [p_ms, report] = solveMultiscaleIteratively(A, q, basis, opt.getSmoother, ...
                                                             opt.tolerance,...
                                                             opt.iterations, ...
                                                             opt.LinSolve, ...
                                                             opt.useGMRES);

    
    state = setFluxes(state, CG, T, fluid, A, q, p_ms, recover, opt, incompOpt);
    
end

function state = setFluxes(state, CG, T, fluid, A, rhs, pressure, recover, opt, incompOpt)
    G = CG.parent;
    
    setFlux = @(p) ...
        incompTPFA(state, G, T, fluid, 'LinSolve', @(varargin) p, 'use_trans', numel(T) == G.faces.num, incompOpt{:});
    
    p_primal = recover(pressure);
    state = setFlux(p_primal);
    
    if opt.reconstruct
        sp = reconstructPressure(CG, pressure, A, rhs);
        sp = recover(sp);
        
        state_o = setFlux(sp);
        flux = state_o.flux;
        flux(CG.faces.fconn) = state.flux(CG.faces.fconn);
        state.flux = flux;
        state.reconstructedPressure = sp(1:G.cells.num);
        if isfield(state_o, 'wellSol')
            state.wellSol = state_o.wellSol;
        end
    end
end
