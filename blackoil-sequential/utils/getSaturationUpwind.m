function [flag_v, flag_g] = getSaturationUpwind(upwtype, state, G, vT, T, K, upstr)
    switch lower(upwtype)
        case 'potential'
            flag = multiphaseUpwindIndices(G, vT, T, K, upstr);
            [flag_v, flag_g] = deal(flag);
        case 'hybrid'
            [flag_v, flag_g] = hybridUpwind(G, vT);
        case 'static'
            [flag_v, flag_g] = deal(state.upstreamFlag);
        otherwise
            error(['Unknown upwind type: ''', model.upwindType, '''']);
    end
end