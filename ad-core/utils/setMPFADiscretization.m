function model = setMPFADiscretization(model)
    % Set MPFA discretization on a model

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
    m = setMPFA(m);
    if isWrapper
        model.parentModel = m;
    else
        model = m;
    end
end

function model = setMPFA(model)
    require mpfa
    [~, M] = computeMultiPointTrans(model.G, model.rock); % From MPFA code
    Tv = M.rTrans; % Cells -> Inner faces
    Tg = M.rgTrans(model.operators.internalConn, :); % Inner -> Inner
    % Change sign and re-scale operators to fit with AD-OO
    % definition of gradient.
    T = getFaceTransmissibility(model.G, model.rock);
    scale = -1./(2*T(model.operators.internalConn));
    MPFAGrad = bsxfun(@times, Tv, scale);
    Mg = -bsxfun(@times, Tg, scale/2);
    assert(all(M.N(:) == model.operators.N(:)), ...
        'Operator neighborship does not match MPFA neighborship. NNC?');
    if isempty(model.FlowDiscretization)
        model = model.setupStateFunctionGroupings();
    end
    % Discrete gradient
    fd = model.FlowDiscretization;
    dp = fd.getStateFunction('PressureGradient');
    model.operators.mpfagrad = MPFAGrad;
    dp.Grad = @(x) MPFAGrad*x;
    fd = fd.setStateFunction('PressureGradient', dp);
    % Gravity potential difference
    dg = fd.getStateFunction('GravityPotentialDifference');
    dg.weight = Mg;
    fd = fd.setStateFunction('GravityPotentialDifference', dg);
    
    model.FlowDiscretization = fd;
end
