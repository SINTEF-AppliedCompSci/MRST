function model = setTimeDiscretization(model, type, varargin)
    % Set the discretization choice for a model

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

    switch lower(type)
        case {'fully-implicit', 'fim'}
            fb = ImplicitFlowStateBuilder(varargin{:});
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
    assert(isa(model, 'GenericReservoirModel'), ...
        'Only the "Generic" class of models support different temporal discretization');
    if isempty(model.FlowDiscretization)
        model = model.validateModel();
    end
    flux = model.FlowDiscretization;
    flux = flux.setFlowStateBuilder(fb);
    model.FlowDiscretization = flux;
end
