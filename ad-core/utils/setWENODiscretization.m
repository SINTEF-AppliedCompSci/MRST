function model = setWENODiscretization(model, varargin)
    % Set WENO discretization on a model
    isWrapper = isa(model, 'WrapperModel');
    if isWrapper
        m = model.parentModel;
    else
        m = model;
    end
    m = setWENO(m, varargin{:});
    if isWrapper
        model.parentModel = m;
    else
        model = m;
    end
end

function model = setWENO(model, varargin)
    weno = WENOUpwindDiscretization(model, model.G.griddim, varargin{:}); % Create WENO
    if isempty(model.FluxDiscretization)
        model = model.setupStateFunctionGroupings();
    end
    fd = model.FluxDiscretization;
    fd = fd.setStateFunction('FaceMobility', FaceMobility(model, weno));
    fd = fd.setStateFunction('FaceComponentMobility', FaceComponentMobility(model, weno));
    model.FluxDiscretization = fd;
end
