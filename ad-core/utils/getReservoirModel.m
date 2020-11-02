function rmodel = getReservoirModel(model)
%Get the underlying reservoir model of a WrapperModel
%
% SYNOPSIS:
%   rmodel = getReservoirModel(model)
%
% PARAMETERS:
%   model - Instance of WrapperModel class.
%
% RETURNS:
%   rmodel - ReservoirModel contained in wrapper model instance. If `model`
%            is an instance of `SequentialPressureTransportModel` the
%            function returns the `ReservoirModel` of the contained
%            `PressureModel`.
%
%            Otherwise, this is the top-most `ReservoirModel` - the one
%            that does not have a `parentModel`.
%
% EXAMPLE:
%   tmodel = TransportModel(model);     % Make transport model
%   rmodel = getReservoirModel(tmodel); % rmodel = model
%
% SEE ALSO:
%   `WrapperModel`, `ReservoirModel`, `setReservoirModel`
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