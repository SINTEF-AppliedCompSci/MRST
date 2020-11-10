function rmodel = getReservoirModel(model)
% Get the underlying reservoir model from a model.
% This function probably belongs in ad-core

%{
#COPYRIGHT#
%}

    rmodel = model;
    if isa(rmodel, 'SequentialPressureTransportModel')
        % We have a sequential implicit model, use pressure parent model
        rmodel = model.pressureModel.parentModel;
    end
    while isprop(rmodel, 'parentModel')
        % We have a WrapperModel - travers to the bottom of the hierarchy
        rmodel = rmodel.parentModel;
    end
end