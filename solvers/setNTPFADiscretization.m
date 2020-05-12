function model = setNTPFADiscretization(model, varargin)

    % Set NTPFA discretization on a model
    isWrapper = isa(model, 'WrapperModel');
    if isWrapper
        m = model.parentModel;
    else
        m = model;
    end

    m = setNTPFA(m, varargin{:});

    if isWrapper
        model.parentModel = m;
    else
        model = m;
    end
end

function model = setNTPFA(model, varargin)
    require nfvm

    if isempty(model.FluxDiscretization)
        model = model.setupStateFunctionGroupings();
    end

    opt = struct('avgmpfa', false);
    [opt, extra] = merge_options(opt, varargin{:});
    
    if opt.avgmpfa
        ntpfa = AvgMPFA(model, extra{:});
    else
        ntpfa = NTPFA(model, extra{:});
    end
    
    % Discrete gradient
    fd = model.FluxDiscretization;
    dp = fd.getStateFunction('PressureGradient');
    dp.Grad = @(p) ntpfa.gradient(p);
    fd = fd.setStateFunction('PressureGradient', dp);

    % % Gravity potential difference
    % dg = fd.getStateFunction('GravityPotentialDifference');
    % dg.weight = Mg;
    % fd = fd.setStateFunction('GravityPotentialDifference', dg);

    model.FluxDiscretization = fd;

end
