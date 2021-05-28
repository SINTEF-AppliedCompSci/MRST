function model = setReservoirModel(model, rmodel)
%Set the underlying reservoir model of a WrapperModel
%
% SYNOPSIS:
%   model = setReservoirModel(model, rmodel)
%
% PARAMETERS:
%   model  - Instance of WrapperModel class.
%   rmodel - Instance of ReservoirModel class to replace the current
%            underlying reservoir model of `model`
%
% RETURNS:
%   model - WrapperModel with updated Reservoir model. If `model`
%           is an instance of `SequentialPressureTransportModel` the
%           function replaces the `ReservoirModel` of the contained
%           `PressureModel` and `TransportModel`
%
% EXAMPLE:
%   tmodel = TransportModel(model);             % Make transport model
%   rmodel = getReservoirModel(tmodel);         % rmodel = model
%   rmodel.rock.perm = rmodel.rock.perm*1e3;    % Multiply permeability
%   rmodel = rmodel.setupOperators();           % Update model operators
%   tmodel = setReservoirModel(tmodel, rmodel); % Update tmodel
%
% SEE ALSO:
%   `WrapperModel`, `ReservoirModel`, `getReservoirModel`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
