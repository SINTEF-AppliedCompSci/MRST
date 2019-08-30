function varargout = getIncompProps(state, fluid, varargin)
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

function [rho, kr, mu, dkr] = getPropsAD(state, fluid, varargin)
    nc = size(state.pressure, 1);
    opt = struct('pressure', repmat(1*atm, nc, 1));
    opt = merge_options(opt, varargin{:});
    
    rho = [fluid.rhoWS.*fluid.bW(opt.pressure), ...
           fluid.rhoOS.*fluid.bO(opt.pressure)];
    mu = [fluid.muW(opt.pressure), ...
           fluid.muO(opt.pressure)];
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