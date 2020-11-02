function model = setReservoirModel(model, rmodel)
% Replace the underlying reservoir model of a WrapperModel with rmdoel. If
% model is an instance of SequentialPressureTransportModel, the function
% replaces the reservoir model of the pressure and transport model.
    if isa(model, 'SequentialPressureTransportModel')
        % Sequential implicit model - set reservoir model for both pressure
        % and transport
        model.pressureModel  = setReservoirModel(model.pressureModel , rmodel);
        model.transportModel = setReservoirModel(model.transportModel, rmodel);
        return
    elseif ~isprop(model, 'parentModel')
        % This model does not have a parent, simply return rmodel
        model = rmodel;
        return
    elseif ~isa(model.parentModel, 'ReservoirModel')
        % We have not reached the bottom, recursive call with parent
        model.parentModel = setReservoirModel(model.parentModel, rmodel);
        return
    end
    % We have reached the bottom - replace parent
    model.parentModel = rmodel;
end