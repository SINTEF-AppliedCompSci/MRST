function [model, nls] = setupBlockTestSolver(model, type, mode)
    nls = getNonLinearSolver(model);
    % We already have a short initial time-step
    nls.timeStepSelector.firstRampupStepRelative = 1;
    arg = {'tolerance', 1e-3, 'relaxation', 'ilu0', 'solver', 'bicgstab'};
    linsolve = AMGCL_CPRSolverAD(arg{:});
    linsolve.doApplyScalingCPR = false;
    linsolve.strategy = 'amgcl';
    nls.LinearSolver = linsolve;
    switch type
        case 'legacy'
            model = ThreePhaseBlackOilModel(model.G, model.rock, model.fluid, 'disgas', false, 'vapoil', false);
        case 'sparse'
            model.AutoDiffBackend = AutoDiffBackend();
        case 'diagonal'
            model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true);
        case 'diagonal-block'
            model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true, 'deferredAssembly', true);
            lsolve_block = AMGCL_CPRSolverBlockAD(arg{:});
            nls.LinearSolver = lsolve_block;
    end
    model = model.setStateFunctionEvaluationMode(mode);
    model = model.validateModel();
    if strcmpi(mode, 'full')
        % Optional: Ignore some stuff
        filter = {'FlowDiscretization.PhaseFlux', 'FlowDiscretization.FaceMobility'};
        model = model.setupStateFunctionGraph('filter', filter);
    end
    ctm = model.FlowPropertyFunctions.ComponentTotalMass;
    ctm = ctm.setMinimumDerivatives([1e-10, 1e-6, 1e-6]);
    model.FlowPropertyFunctions.ComponentTotalMass = ctm;
end