function [flag_v, flag_g] = hybridUpwind(Gi, vT)
% Hybrid upwinding - two-phase only at the moment
    nPh = numel(Gi);
    assert(nPh == 2, ...
        'Hybrid upwinding currently only supported for two-phase problems.');
    
    % Total velocity for viscous flow, density sorting for buoyancy terms
    flag_v = repmat(vT > 0, 1, nPh);
    dp = Gi{2} - Gi{1};
    d_diff = dp <= 0;
    flag_g = [d_diff, ~d_diff];
end