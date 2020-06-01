function model = setAvgMPFADiscretization(model, varargin)

    % Set AvgMPFA discretization on a model
    isWrapper = isa(model, 'WrapperModel');
    if isWrapper
        m = model.parentModel;
    else
        m = model;
    end

    m = setAvgMPFA(m, varargin{:});

    if isWrapper
        model.parentModel = m;
    else
        model = m;
    end
end

function model = setAvgMPFA(model, varargin)

    require nfvm

    if isempty(model.FluxDiscretization)
        model = model.setupStateFunctionGroupings();
    end

    avgmpfa = AvgMPFA(model, varargin{:});

    % Discrete gradient
    fd = model.FluxDiscretization;
    dp = fd.getStateFunction('PressureGradient');
    dp.Grad = @(p) avgmpfa.gradient(p);
    fd = fd.setStateFunction('PressureGradient', dp);

    % % Gravity potential difference
    % dg = fd.getStateFunction('GravityPotentialDifference');
    % dg.weight = Mg;
    % fd = fd.setStateFunction('GravityPotentialDifference', dg);

    model.FluxDiscretization = fd;

end
