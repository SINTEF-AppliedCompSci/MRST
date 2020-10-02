function model = setReservoirModel(model, rmodel)
    % Set the underlying reservoir model to a model.
    % This function probably belongs in ad-core
    if isa(model, 'SequentialPressureTransportModel')
        % Model is a sequential implicit model - set reservoir model for
        % pressure and transport by recursion
        model.pressureModel  = model.setReservoirModel(model.pressureModel , rmodel);
        model.transportModel = model.setReservoirModel(model.transportModel, rmodel);
        return
    elseif ~isprop(model, 'parentModel')
        % This model does not have a parent, simply return rmodel
        model = rmodel;
        return
    elseif ~isa(model.parentModel, 'ReservoirModel')
        % We have not reached the bottom, recursive call with parent
        model.parentModel = model.setReservoirModel(model.parentModel, rmodel);
        return
    end
    % We have reached the bottom - replace parent
    model.parentModel = rmodel;
end