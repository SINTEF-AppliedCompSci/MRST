function cfl = estimateCFL(model, state, dt, W)
% Estimate CFL
    op = model.operators;
    N = op.N;
    u = abs(sum(state.flux(op.internalConn, :), 2));
    cfl_face = (u*dt)./op.pv(N);
    
    cfl = zeros(model.G.cells.num, 1);
    cfl(N(:, 1)) = cfl_face(:, 1);
    cfl(N(:, 2)) = max(cfl_face(:, 2), cfl(N(:, 2)));
    
    if nargin > 3
        for i = 1:numel(W)
            w = W(i);
            ws = state.wellSol(i);
            
            q = abs(sum(ws.flux, 2));
            
            cfl_w = q*dt./op.pv(w.cells);
            cfl(w.cells) = max(cfl(w.cells), cfl_w);
        end
    end
end