function model = setNTPFADiscretization(model, varargin)
% Set NTPFA discretization on a model

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    m = setNTPFA(m, varargin{:});

    if isWrapper
        model.parentModel = m;
    else
        model = m;
    end
end

function model = setNTPFA(model, varargin)

    require nfvm

    if isempty(model.FlowDiscretization)
        model = model.setupStateFunctionGroupings();
    end

    ntpfa = NTPFA(model, varargin{:});

    % Discrete gradient
    fd = model.FlowDiscretization;
    dp = fd.getStateFunction('PressureGradient');
    dp.Grad = @(p) ntpfa.gradient(p);
    fd = fd.setStateFunction('PressureGradient', dp);
    model.FlowDiscretization = fd;

    model.operators.Grad = dp.Grad;

    % Gravity
    if norm(model.getGravityVector()) > 0
        model.operators.gdz = setup_gdz(model, dp.Grad);
    else
        model.operators.gdz = zeros(size(model.operators.N, 1), 1);
    end

end
