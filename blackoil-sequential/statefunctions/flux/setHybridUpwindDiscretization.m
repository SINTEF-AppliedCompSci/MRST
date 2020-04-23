function model = setHybridUpwindDiscretization(model)
    % Set HU discretization on a sequential model
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
    if isempty(model.FluxDiscretization)
        model = model.setupStateFunctionGroupings();
    end
    % Discrete gradient
    fd = model.FluxDiscretization;
    ppu_upwind = PhasePotentialUpwindFlag(model);
    ppu_upwind.includeTotalVelocity = false;
    
    fd = fd.setStateFunction('PhaseUpwindFlag', PhaseUpwindFlagTotalVelocity(model));
    fd = fd.setStateFunction('PhaseUpwindFlagGravity', ppu_upwind);
    
    upwind = UpwindFunctionWrapperDiscretization(model);
    fd = fd.setStateFunction('FaceComponentMobilityGravity', FaceComponentMobility(model, upwind, 'PhaseUpwindFlagGravity'));
    fd = fd.setStateFunction('FaceMobilityGravity', FaceMobility(model, upwind, 'PhaseUpwindFlagGravity'));
    fd = fd.setStateFunction('FaceTotalMobilityGravity', FaceTotalMobility(model, 'FaceMobilityGravity'));

    fd = fd.setStateFunction('ComponentPhaseFlux', ComponentPhaseFluxFractionalFlowHybridUpwind(model));
    model.FluxDiscretization = fd;
end
