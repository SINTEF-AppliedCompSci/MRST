function rmodel = getReservoirModel(model)
% Get the underlying reservoir model of a WrapperModel. If model is an
% instance of SequentialPressureTransportModel, the function returns the
% reservoir model of the pressure model.
    rmodel = model;
    if isa(rmodel, 'SequentialPressureTransportModel')
        % We have a sequential implicit model, use pressure model
        rmodel = model.pressureModel;
    end
    while isprop(rmodel, 'parentModel')
        % We have a WrapperModel - travers to the bottom of the hierarchy
        rmodel = rmodel.parentModel;
    end
end