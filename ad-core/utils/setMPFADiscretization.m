function model = setMPFADiscretization(model)
    % Set MPFA discretization on a model
    isWrapper = isa(model, 'WrapperModel');
    if isWrapper
        m = model.parentModel;
    else
        m = model;
    end
    m = setMPFA(m);
    if isWrapper
        model.parentModel = m;
    else
        model = m;
    end
end

function model = setMPFA(model)
    require mpfa
    [~, M] = computeMultiPointTrans(model.G, model.rock); % From MPFA code
    Tv = M.rTrans; % Cells -> Inner faces
    Tg = M.rgTrans(model.operators.internalConn, :); % Inner -> Inner
    % Change sign and re-scale operators to fit with AD-OO
    % definition of gradient.
    scale = -1./(2*model.operators.T);
    MPFAGrad = Tv.*scale;
    Mg = -Tg.*scale/2;
    assert(all(M.N(:) == model.operators.N(:)), ...
        'Operator neighborship does not match MPFA neighborship. NNC?');
    if isempty(model.FluxDiscretization)
        model = model.setupStateFunctionGroupings();
    end
    % Discrete gradient
    fd = model.FluxDiscretization;
    dp = fd.getStateFunction('PressureGradient');
    dp.Grad = @(x) MPFAGrad*x;
    fd = fd.setStateFunction('PressureGradient', dp);
    % Gravity potential difference
    dg = fd.getStateFunction('GravityPotentialDifference');
    dg.weight = Mg;
    fd = fd.setStateFunction('GravityPotentialDifference', dg);
    
    model.FluxDiscretization = fd;
end
