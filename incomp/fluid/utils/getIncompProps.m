function varargout = getIncompProps(state, fluid)
    varargout = cell(1, nargout);
    if isfield(fluid, 'properties')
        [varargout{:}] = getPropsLegacy(state, fluid);
    else
        [varargout{:}] = getPropsAD(state, fluid);
    end
end

function [rho, kr, mu, dkr] = getPropsLegacy(state, fluid)
    [mu, rho] = fluid.properties(state);
    s         = fluid.saturation(state);
    if nargout < 4
        kr = fluid.relperm(s, state);
    else
        [kr, dkr] = fluid.relperm(s, state);
    end
end

function [rho, kr, mu, dkr] = getPropsAD(state, fluid)
    if isfield(state, 'pressure')
        p = state.pressure;
        p_avg = mean(p);
    else
        p_avg = 1*atm;
    end
    
    rho = [fluid.rhoWS*fluid.bW(p_avg), ...
           fluid.rhoOS*fluid.bO(p_avg)];
    mu = [fluid.muW(p), ...
           fluid.muO(p)];
    getDer = nargout > 3;
    s = state.s(:, 1);
    if getDer
        s = initVariablesAD_diagonal(s);
    end
    
    krW = fluid.krW(s);
    so = 1 - s;
    if isfield(fluid, 'krOW')
        krO = fluid.krOW(so);
    else
        krO = fluid.krO(so);
    end
    kr = [value(krW), value(krO)];
    if getDer
        dkr = [krW.jac{1}.diagonal, krO.jac{1}.diagonal]; 
    end
end