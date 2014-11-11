function basis = computeDefaultBasis(basis, G, state, system, W, fluid, pv, T, varargin)
    opt = struct('compIx',    [], ...
                 'linsolve',  @mldivide);
    opt = merge_options(opt, varargin{:});
    if isempty(opt.compIx)
        state = solveStationaryPressure(G, state, system, W, fluid, pv, T,  'computeBasisOnly', true, 'linsolve', opt.linsolve);
        basis = state.basis;
    else
        ix = opt.compIx;
        state = solveStationaryPressure(G, state, system, W(ix), fluid, pv, T,  'computeBasisOnly', true, 'linsolve', opt.linsolve);
        basis.B(:,ix) = state.basis.B;
        basis.R(ix,:) = state.basis.R;
    end
end


