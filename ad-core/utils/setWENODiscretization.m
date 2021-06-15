function model = setWENODiscretization(model, varargin)
    % Set WENO discretization on a model

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

    isWrapper = isa(model, 'WrapperModel');
    if isWrapper
        m = model.parentModel;
    else
        m = model;
    end
    m = setWENO(m, varargin{:});
    if isWrapper
        model.parentModel = m;
    else
        model = m;
    end
end

function model = setWENO(model, varargin)
    weno = WENOUpwindDiscretization(model, model.G.griddim, varargin{:}); % Create WENO
    if isempty(model.FlowDiscretization)
        model = model.setupStateFunctionGroupings();
    end
    fd = model.FlowDiscretization;
    fd = fd.setStateFunction('FaceMobility', FaceMobility(model, weno));
    fd = fd.setStateFunction('FaceComponentMobility', FaceComponentMobility(model, weno));
    model.FlowDiscretization = fd;
end
