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

    rmodel = model;
    if isa(rmodel, 'SequentialPressureTransportModel')
        % We have a sequential implicit model, use pressure model
        if isempty(model.parentModel)
            rmodel = model.pressureModel;
        end
    end
    while isprop(rmodel, 'parentModel')
        % We have a WrapperModel - travers to the bottom of the hierarchy
        rmodel = rmodel.parentModel;
    end
end
