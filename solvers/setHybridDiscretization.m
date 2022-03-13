function model = setHybridDiscretization(model, models, faceblocks, varargin)
% Set hybrid discretizations on a model

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    m2 = models;
    for k = 1:numel(models)
        if isa(models{k}, 'WrapperModel')
            m2{k} = models{k}.parentModel;
        end
    end

    m = setHybrid(m, m2, faceblocks, varargin{:});

    if isWrapper
        model.parentModel = m;
    else
        model = m;
    end

end

function model = setHybrid(model, models, faceblocks, varargin)

    if isempty(model.FlowDiscretization)
        model = model.setupStateFunctionGroupings();
    end

    fdmodels = cell(1, numel(models));

    for k = 1:numel(models)
        if isempty(models{k}.FlowDiscretization)
            models{k} = models{k}.setupStateFunctionGroupings();
        end
        fdmodels{k} = models{k}.FlowDiscretization;
    end

    hybrid = Hybrid(model, fdmodels, faceblocks);

    % Discrete gradient
    fd = model.FlowDiscretization;
    dp = fd.getStateFunction('PressureGradient');
    dp.Grad = @(p) hybrid.gradient(p);
    fd = fd.setStateFunction('PressureGradient', dp);

    % % Gravity potential difference
    % dg = fd.getStateFunction('GravityPotentialDifference');
    % dg.weight = Mg;
    % fd = fd.setStateFunction('GravityPotentialDifference', dg);

    model.FlowDiscretization = fd;

    model.operators.Grad = dp.Grad;

    % Gravity
    model.operators.gdz = setup_gdz(model, dp.Grad);

end
