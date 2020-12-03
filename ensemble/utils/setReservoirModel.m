function model = setReservoirModel(model, rmodel)
% Set the underlying reservoir model to a model.
% This function probably belongs in ad-core

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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