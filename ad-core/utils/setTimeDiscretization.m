function model = setTimeDiscretization(model, type, varargin)
    % Set the discretization choice for a model
    switch lower(type)
        case {'fully-implicit', 'fim'}
            fb = FlowStateBuilder(varargin{:});
        case {'adaptive-implicit', 'aim'}
            fb = AdaptiveImplicitFlowStateBuilder(varargin{:});
        case {'explicit', 'impes'}
            fb = ExplicitFlowStateBuilder(varargin{:});
        otherwise
            error('%s not supported. Valid choices: FIM, AIM or IMPES');
    end

    if isa(model, 'SequentialPressureTransportModel')
        % Trigger validation - on wrapper
        model.transportModel = model.transportModel.validateModel();
        model.transportModel.parentModel = setFSB(model.transportModel.parentModel, fb);
    else
        model = setFSB(model, fb);
    end
end

function model = setFSB(model, fb)
    assert(isa(model, 'ExtendedReservoirModel'), ...
        'Only the "Generic" class of models support different temporal discretization');
    if isempty(model.FluxDiscretization)
        model = model.validateModel();
    end
    flux = model.FluxDiscretization;
    flux = flux.setFlowStateBuilder(fb);
    model.FluxDiscretization = flux;
end