function [model, nls] = setupBlockTestSolver(model, type, mode)
%Undocumented Utility Function

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
