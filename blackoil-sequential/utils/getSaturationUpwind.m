function [flag_v, flag_g] = getSaturationUpwind(upwtype, state, G, vT, T, K, upstr)
% Get upwind flags for viscous and gravity parts of flux function.
    switch lower(upwtype)
        case 'potential'
            % Standard PPU upwind
            flag = multiphaseUpwindIndices(G, vT, T, K, upstr);
            [flag_v, flag_g] = deal(flag);
        case 'hybrid'
            % Hybrid upwinding
            [flag_v, flag_g] = hybridUpwind(G, vT);
        case 'static'
            % Use whatever flags are kept in state, likely from the
            % pressure solver. Debug only.
            [flag_v, flag_g] = deal(state.upstreamFlag);
        otherwise
            error(['Unknown upwind type: ''', model.upwindType, '''']);
    end
end