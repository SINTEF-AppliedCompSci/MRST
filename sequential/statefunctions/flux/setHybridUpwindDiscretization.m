function model = setHybridUpwindDiscretization(model)
    % Set HU discretization on a sequential model

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
        % Set it for transport model
        model = model.validateModel(); % Make sure that facility etc is set up.
        model.transportModel = setHybridUpwindDiscretization(model.transportModel);
    else
        isWrapper = isa(model, 'WrapperModel');
        if isWrapper
            assert(isa(model, 'TransportModel'));
            m = model.parentModel;
        else
            m = model;
        end
        m = setHybrid(m);
        if isWrapper
            model.parentModel = m;
        else
            model = m;
        end
    end
end

function model = setHybrid(model)
    if isempty(model.FlowDiscretization)
        model = model.setupStateFunctionGroupings();
    end
    % Discrete gradient
    fd = model.FlowDiscretization;
    ppu_upwind = PhasePotentialUpwindFlag(model);
    ppu_upwind.includeTotalVelocity = false;
    
    fd = fd.setStateFunction('PhaseUpwindFlag', PhaseUpwindFlagTotalVelocity(model));
    fd = fd.setStateFunction('PhaseUpwindFlagGravity', ppu_upwind);
    
    upwind = UpwindFunctionWrapperDiscretization(model);
    fd = fd.setStateFunction('FaceComponentMobilityGravity', FaceComponentMobility(model, upwind, 'PhaseUpwindFlagGravity'));
    fd = fd.setStateFunction('FaceMobilityGravity', FaceMobility(model, upwind, 'PhaseUpwindFlagGravity'));
    fd = fd.setStateFunction('FaceTotalMobilityGravity', FaceTotalMobility(model, 'FaceMobilityGravity'));

    fd = fd.setStateFunction('ComponentPhaseFlux', ComponentPhaseFluxFractionalFlowHybridUpwind(model));
    model.FlowDiscretization = fd;
end
